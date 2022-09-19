params.pinn_models = './models/*/model'
params.remd_tag = 'nvt-340k-100ps'
params.remd_flags = '--ensemble nvt --T 340 --t 100 --dt 0.5 --log-every 20'

workflow {
  channel.fromPath(params.pinn_models, type:'dir') \
    | map {model -> [model.parent.name, model]}    \
    | set {models}

  models
    | map {name, model -> [(name =~ /(.*)\d/)[0][1], model]} \
    | groupTuple                                              \
    | combine(channel.fromPath(params.asemd_init))            \
    | map { name, model, init -> ["$params.remd_tag/$init.baseName", model, init, params.remd_flags]} \
    | aseEMD
}
