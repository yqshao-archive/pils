params.exp = 'exp/prod-adam-run2'
params.gen = '30'
params.init = 'skel/init/lmp-geo/{a32b32i0,a16b16i16,a0b0i32}-r{1.08,1.16}.xyz'
params.md_from = 0
params.md_tag = 'nvt-340k-5ns'
params.md_flags = '--ensemble nvt --T 340 --t 5000 --dt 0.5 --log-every 200'

include { aseMD } from './tips/nextflow/ase.nf' addParams(publish: "$params.exp/prod/gen$params.gen")

workflow prod_init {
  channel.value(file("$params.exp/models/gen$params.gen/model*/model", type:'dir')) \
    | combine(channel.fromPath(params.init)) \
    | map { model, init -> ["${params.md_tag}-0/${init.baseName}", model, init, params.md_flags]} \
    | aseMD
}

workflow prod_cont {
  prev = params.md_from.toInteger()
  next = prev + 1
  restart = channel.fromPath("$params.exp/prod/gen$params.gen/$params.md_tag-$prev/*/asemd.traj")
  channel.value(file("$params.exp/models/gen$params.gen/model*/model", type:'dir')) \
    | combine(restart) \
    | map { model, init -> ["${params.md_tag}-${next}/${init.parent.name}", model, init, params.md_flags]} \
    | aseMD
}

params.cp2k_inp = 'skel/cp2k/nvt-10ps.inp'
params.cp2k_aux = 'skel/cp2k-aux/*'
params.cp2k_geo = 'skel/init/prod-geo/*.xyz'
params.cp2k_from = '10'
include { cp2kMD; cp2k } from './tips/nextflow/cp2k.nf' addParams(publish: 'exp/prod-adam-run2/cp2k-vali')

workflow vali_init {
  channel.fromPath(params.cp2k_geo) \
    | map {it -> \
           ["nvt-0-10ps/${it.baseName}", file(params.cp2k_inp), it, "--fmt asetraj --idx -1"]} \
    | cp2kMD
}

workflow vali_cont {
  int from = params.cp2k_from.toInteger()
  channel.fromPath("exp/prod-adam-run2/cp2k-vali/nvt-${from-10}-${from}ps/*/cp2k-md-1_${2000*from}.restart") \
    | map {it -> \
           ["nvt-${from}-${from+10}ps/${it.parent.name}", it, file(params.cp2k_aux)] }\
    | cp2k
}
