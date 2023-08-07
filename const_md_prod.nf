params.exp = 'exp/transfer'
params.gen = '35'
params.init = 'skel/restart/*.xyz'
params.md_from = 0
params.md_tag = 'nvt-340k-5ns'
params.md_flags = '--ensemble nvt --t 5000 --dt 0.5 --log-every 200'

include { constMD } from './nf/const_md.nf'

workflow {
  channel.value(file("$params.exp/models/gen$params.gen/model*/model", type:'dir')) \
    | combine(channel.fromPath(params.init)) \
    | map { model, init -> [model, init, (init =~ /(\d+)k/)[0][1]] } \
    | map { model, init, T -> ["nvt-${T}k-5ns-0/${init.baseName}", model, init, "${params.md_flags} --T ${T}", file('./py', type:'dir')]} \
    | constMD
}