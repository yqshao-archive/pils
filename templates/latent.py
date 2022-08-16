#!/usr/bin/env python
#
def save_tmp_ckpt(model_dir):
    import tempfile
    import tensorflow as tf
    tmpdir = tempfile.mkdtemp()
    vars = {}
    reader = tf.train.load_checkpoint(model_dir)
    dtypes = reader.get_variable_to_dtype_map()
    for key in dtypes.keys():
        if ('global_step' in key) or ('TRAIN_OP' in key):
            continue
        vars[key+':0'] = tf.Variable(reader.get_tensor(key))
    tmp_ckpt = tf.train.Checkpoint(**vars).save(f'{tmpdir}/ckpt')
    return tmp_ckpt


def bpnn_extract(tensors, network, idx=0):
    """Extract features from a PiNet

    Args:
        tensors (dict of tensors): dictionary of tensors
        network (PiNet): a PiNet network instance
        idx (int): index of structure

    Return:
        dict of np.array: latent variables
    """
    arrays = {}
    tensors = network.preprocess(tensors)
    tensors = network.fingerprint(tensors)
    for k, v in network.feed_forward.ff_layers.items():
        prop = tensors['elem_fps'][k]
        for i, layer in enumerate(v[:-1]):
            prop = layer(prop)
            arrays[f'{k}_ff{i}'] = prop
        arrays[f'{k}_final'] = prop * v[-1].kernel[None,:,0]
    return arrays



def pinet_extract(tensors, network, idx=0):
    """Extract features from a PiNet

    Args:
        tensors (dict of tensors): dictionary of tensors
        network (PiNet): a PiNet network instance
        idx (int): index of structure

    Return:
        dict of np.array: latent variables
    """
    import numpy as np
    import tensorflow as tf

    tensors = network.preprocess(tensors)
    fc = network.cutoff(tensors["dist"])
    basis = network.basis_fn(tensors["dist"], fc=fc)
    output = 0.0
    arrays = {}
    final = []  # aggretes the

    ind_2 = tensors["ind_2"]
    ind_i = ind_2[:, 0]
    ind_j = ind_2[:, 1]
    for i in range(network.depth):
        arrays[f"gc{i}_0"] = tensors[
            "prop"
        ].numpy()  # gc{i}_0 initial property each GC block
        # manual propogation to get the latent variables
        ptmp = tensors["prop"]
        for j, layer in enumerate(network.gc_blocks[i].pp_layer.dense_layers):
            ptmp = layer(ptmp)
            arrays[f"gc{i}_pp{j}"] = ptmp.numpy()

        prop_i = tf.gather(ptmp, ind_i)
        prop_j = tf.gather(ptmp, ind_j)
        inter = tf.concat([prop_i, prop_j], axis=-1)
        inter = network.gc_blocks[i].pi_layer.ff_layer(inter)
        arrays[f"W_{i}"] = inter

        prop = network.gc_blocks[i]([tensors["ind_2"], tensors["prop"], basis])
        arrays[f"gc{i}_pi"] = prop.numpy()
        ptmp = prop
        for j, layer in enumerate(network.out_layers[i].ff_layer.dense_layers):
            ptmp = layer(ptmp)
            arrays[f"gc{i}_out{j}"] = ptmp.numpy()
        final.append(
            (ptmp * network.out_layers[i].out_units.kernel[None, :, 0]).numpy()
        )
        output = network.out_layers[i]([tensors["ind_1"], prop, output])
        tensors["prop"] = network.res_update[i]([tensors["prop"], prop])

    arrays["final"] = np.concatenate(final, axis=1)
    return arrays


def extract_feature(model_dir, ds_path, take, batch):
    """Extract Features from a Trained Model

    Args:
       model_dir (str): model directory
       ds_path (str): dataset path (in runner format)

    Returns:
       dictionary of numpy arrays with labels
    """
    import yaml
    import numpy as np
    import tensorflow as tf
    from pinn.io import load_runner, sparse_batch
    from pinn import get_network

    tf.keras.backend.clear_session()
    with open(f'{model_dir}/params.yml') as f:
        nn_spec = yaml.safe_load(f)['network']
    if nn_spec['name']=='PiNet':
        extract = pinet_extract
    elif nn_spec['name']=='BPNN':
        nn_spec['params']['use_jacobian'] = False
        extract = bpnn_extract
    else:
        raise "Unknow network architecture"

    # converts the checkpoint to eager-mode-compatible format
    tmp_ckpt = save_tmp_ckpt(model_dir)

    # setup the network
    dataset = load_runner(ds_path).take(take).apply(sparse_batch(batch))
    network = get_network(nn_spec)
    network(next(iter(dataset)))
    ckpt = tf.train.Checkpoint(**{v.name: v for v in network.variables})

    try:
        ckpt.restore(tmp_ckpt).assert_consumed()
    except:
        print('Inconsistent checkpoint, printing the first 10 variables')
        reader = tf.train.load_checkpoint(tmp_ckpt)
        shapes = reader.get_variable_to_shape_map()
        print(f"Checkpoint at '{tmp_ckpt}':")
        for key in sorted(shapes.keys())[:10]:
            print(f"  |- key='{key.replace('.S','/')}', shape={shapes[key]}")
        print(f"Variables in graph")
        for v in sorted(network.variables, key=lambda x: x.name)[:10]:
            print(f"  |- key='{v.name}', shape={v.shape}")
        print("Exiting, check the checkpoints!")
        return

    # the real extraction process
    iterator = iter(dataset)
    data = extract(next(iterator), network, 0)
    for i, next_data in enumerate(iterator):
        print(f"\\r{i+1}", end="")
        datum = extract(next_data, network, i + 1)
        for k in data.keys():
            data[k] = np.concatenate([data[k], datum[k]], axis=0)
    return data

import re
import numpy as np

setup = {
  'take': 100, # log interval in steps
  'batch': 10, # pressure in bar
  'key': 'final',
}

flags = {
  k: v for k,v in
    re.findall('--(.*?)[\\s,\\=]([^\\s]*)', "$flags")
}

setup.update(flags)
data = extract_feature("$model", "$ds", int(setup['take']), int(setup['batch']))
np.save('latent.npy', data[setup['key']])
