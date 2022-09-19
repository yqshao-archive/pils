params.pinn_models = './models/*/model'
params.remd_tag = 'nvt-340k-100ps'
params.remd_flags = '--ensemble nvt --T 340 --t 100 --dt 0.5 --log-every 20'

workflow remd {
  channel.fromPath(params.pinn_models, type:'dir') \
    | map {model -> [model.parent.name, model]}    \
    | set {models}

  // models
  //   | combine(channel.fromPath(params.asemd_init)) \
  //   | map { name, mode, init -> ["$name-$init.baseName", mode, init, params.ase_flags]} \
  //   | aseMD

  models
    | map {name, model -> [(name =~ /(.*)\d/)[0][1], model]} \
    | groupTuple                                              \
    | combine(channel.fromPath(params.asemd_init))            \
    | map { name, model, init -> ["$params.remd_tag/$init.baseName", model, init, params.remd_flags]} \
    | aseEMD
}


//========================================//
// Extract Feature and Perform Clustering //
//========================================//
params.models = 'exp/adam-singlet-run1/models/gen6/*/model'
params.latent_ds = 'datasets/pils-v6.{yml,tfr}'
params.latent_flags = '--take 100'

include { extract_latent } from './nextflow/latent.nf'

workflow latent {
  Channel.fromPath(params.models, type:'dir')     \
    | combine(Channel.fromFilePairs(params.latent_ds)) \
    | map {model, dsname, ds ->                           \
           ["${model.parent.name}-${dsname}", model, ds, params.latent_flags]} \
    | extract_latent
}
