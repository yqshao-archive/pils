#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.ds       = 'datasets/pils-50ps.{yml,tfr}'
params.take     = 100
params.cp2k_aux = 'skel/cp2k-aux/*'

include { cp2k; cp2kGenInp } from './tips/nextflow/cp2k.nf' addParams(publish: "ck-cp2k-frc")

process subsample {
  label 'tips'
  publishDir "ck-cp2k-frc"

  input:
    path(ds)

  output:
    path('ds.traj')

  script:
  """
  #!/usr/bin/env python
  from tips.io import load_ds

  ds = load_ds("${ds[0].baseName}.yml", fmt='pinn')
  ds[:$params.take].convert('ds.traj', fmt='asetraj')
  """
}

workflow {
  ch_inps = Channel.fromList(
    [
      ['sp-prod', file('./skel/cp2k/singlepoint.inp')],
      ['sp-tight', file('./skel/cp2k/sp-tight.inp')],
      ['sp-hicut', file('./skel/cp2k/sp-hicut.inp')]
    ]
  )
  channel.of(file(params.ds)) \
    | subsample \
    | combine(ch_inps) \
    | map {ds, name, inp -> [name, inp, ds, '-f asetraj --subsample --psample 100']} \
    | cp2kGenInp \
    | flatMap {name, inps -> inps.collect{ ["$name/$it.baseName", it] }} \
    | map {name, inp -> [name, inp, file(params.cp2k_aux)]} \
    | cp2k
}
