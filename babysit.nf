#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
This is a workflow for "babysitting" an active learning progress, several tools
are implemented as sub-workflows:

- ebmd: this use a biased ensemble md to extend a trajectory
- label: this filters a previously sampled traj and labels it with cp2k
- retrain: this generates a new generation

The naming pattern:

- models: models/$base-gen$g-seed$s
- resampled-traj: trajs/ebmd/$base-gen$g/$i
- labelled-ds: datasets/

where $base is basename of the model, $g is generation, $s is random seed for
data splitting, $i is filename of inital structure

*/

params.init = 'skel/init/*.xyz'
params.base = 'pils-v5-ekf-v3'
params.gen = 0
params.bias = 'nobias'
params.cp2k_aux = 'skel/cp2k-aux/*'

include { aseMD } from './tips/nextflow/ase.nf' addParams(publish: 'trajs/emd')

workflow ebmd {
  presets = [
    ['nobias', ''],
    ['bias1', ' --bias heaviside --kb 1'],
    ['bias001', ' --bias heaviside --kb 0.01'],
  ]
  md = '--ensemble nvt --T 340 --t 1 --dt 0.5 --log-every 10'

  channel.fromPath("models/$params.base-gen$params.gen-seed*/model", type:'dir') \
    | map {model -> [(model.parent.name =~ /(.*)-seed\d/)[0][1], model]}  \
    | groupTuple                                                          \
    | combine( channel.fromPath(params.init) ) \
    | combine( channel.fromList(presets)) \
    | map { base, model, init, preset, bias -> \
           ["$params.base-gen$params.gen-$preset/$init.baseName", model, init, md+bias]} \
    | aseMD
}


include { cp2kGenInp } from './tips/nextflow/cp2k.nf' addParams(publish: 'trajs/cp2k-sp')
include { cp2k } from './tips/nextflow/cp2k.nf' addParams(publish: 'trajs/cp2k-sp')
workflow label {
  strategies = [
    ['uniform', '-f asetraj --subsample --strategy uniform --nsample 50'],
    ['forces', '-f asetraj --subsample --strategy sorted --nsample 50'],
  ]

  skel = file('./skel/cp2k/singlepoint.inp')
  channel.fromPath("trajs/emd/$params.base-gen$params.gen-$params.bias/*/asemd.traj") \
    | map {traj -> ["${traj.parent.name}", traj]} \
    | combine( channel.fromList(strategies) ) \
    | map {name, traj, tag, flag ->  \
           ["$params.base-gen$params.gen-$params.bias/$name/$tag", skel, traj, flag]} \
    | cp2kGenInp \
    | flatMap {name, inps -> inps.collect {["$name/$it.baseName", it]}}
    | map {name, inp -> [name, inp, file(params.cp2k_aux)]}
    | cp2k
}


process checkConverge {
  label 'tips'

  input:
  tuple val(name), path(idx), path(logs), path(traj)

  output:
  tuple val(name), path('next.xyz'), stdout

  script:
  ftol = params.ftol
  etol = params.etol
  """
  #!/usr/bin/env python
  import numpy as np
  from ase.io import read
  from tips.io import load_ds

  eunit = 27.2114 # hartree to eV, hardcoded for CP2K for now
  idx = [int(i) for i in np.loadtxt("$idx")]
  logs = load_ds(read("$logs", index=':'), fmt='ase')
  traj = load_ds("$traj", fmt='asetraj')
  e_label = np.array([datum['energy']/len(datum['elem']) for datum in logs])*eunit
  f_label = np.array([datum['force'] for datum in logs])*eunit
  e_pred = np.array([traj[i]['energy']/len(traj[i]['elem']) for i in idx])
  f_pred = np.array([traj[i]['force'] for i in idx])

  ecnt = np.sum(np.abs(e_pred-e_label)>$etol)
  fcnt = np.sum(np.any(np.abs(f_pred-f_label)>$ftol,axis=(1,2)))
  converged = (ecnt==0) and (fcnt==0)

  ermse = np.mean((e_pred-e_label)**2)
  frmse = np.mean((f_pred-f_label)**2)

  if converged:
      msg = f'Converged; will restart from latest frame.'
      traj[-1:].convert('next.xyz', fmt='extxyz')
  else:
      msg = f'energy: {ecnt}/{len(idx)} failed, rmse={ermse:.2e}; force {fcnt}/{len(idx)} failed, rmse={frmse:.2e}.'
      traj[:1].convert('next.xyz', fmt='extxyz')
  print(msg)
  """
}

process mergeDS {
  input: tuple val(name), val(idx), path(logs, stageAs:'*.log')
  output: tuple val(name), path('merged.idx'), path('merged.xyz')

  script:
  """
  printf "${idx.join('\\n')}" > merged.idx
  tips convert ${logs.join(' ')} -f cp2klog -o merged
  """
}

process mkNewDS {
  input: tuple val(name), val(conv), path(new), path(newDS), val(newFlag), val(oldFlag)
  output: tuple val(name), path('mixDS.{yml,tfr}')

  script
  """
  tips subsample $newDS -f pinn -o newDs -of pinn $newFlag
  tips subsample $oldDS -f pinn -o oldDs -of pinn $oldFlag
  tips convert newDS.yml oldDS.yml -f pinn -of mixDS.yml
  rm {new,old}DS.*
  """
}

params.ftol = 0.200
params.etol = 0.005
params.initSize = 6
params.old_flag = '--nsample 24000'
params.new_flag = '--nsample 600'

workflow mkNextGenInp {
  base = params.base
  gen = params.gen
  bias = params.bias

  channel.fromPath("./trajs/cpsp/$base-gen$gen-$bias/*/uniform/*/cp2k.log") \
    | map {[it.parent.parent.parent.name, it.parent.name.toInteger(), it]} \
    | groupTuple(size:40, sort:true) \
    | mergeDS \
    | set {logs} // init, idx, log

  channel.fromPath("./trajs/emd/$base-gen$gen-$bias/*/asemd.traj") \
    | map {[it.parent.name, it]} \
    | set {trajs}

  // new init
  logs \
    | join (trajs) \
    | checkConverge \
    | view {name, init, msg -> "[gen$gen-$name] ${msg.trim()}"}\
    | map {name, init, msg -> name, init} \
    | groupTuple {size:params.initSize} \
    | set {newInit}

  // new converge flag
  checkConverge.out \
    | map {name, init, msg -> name, msg}
    | groupTuple {size:params.initSize} \
    | {name, msgs-> [name, ]} \
    | set {newConv}

  // new training set


}
