params.pinn_models = './models/*/model'
params.remd_tag = 'nvt-340k-100ps'
params.remd_flags = '--ensemble nvt --T 340 --t 100 --dt 0.5 --log-every 20'

workflow prod {
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

params.cp2k_inp = 'skel/cp2k/nvt-10ps.inp'
params.cp2k_aux = 'skel/cp2k-aux/*'
params.cp2k_geo = 'skel/init/trial-geo/*.xyz'
params.cp2k_from = '10'
include { cp2kMD; cp2k } from './tips/nextflow/cp2k.nf' addParams(publish: 'exp/trial-adam/cp2k-vali')

workflow vali_init {
  channel.fromPath(params.cp2k_geo) \
    | map {it -> \
           ["nvt-0-10ps/${it.baseName}", file(params.cp2k_inp), it, "--fmt asetraj --idx -1"]} \
    | cp2kMD
}

workflow vali_cont {
  int from = params.cp2k_from.toInteger()
  channel.fromPath("exp/trial-adam/cp2k-vali/nvt-${from-10}-${from}ps/*/cp2k-md-1_${2000*from}.restart") \
    | map {it -> \
           ["nvt-${from}-${from+10}ps/${it.parent.name}", it, file(params.cp2k_aux)] }\
    | cp2k
}
