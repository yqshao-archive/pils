nextflow.enable.dsl=2

params.publish = 'analyses'

process compute_msd {
  publishDir "$params.publish/$name"
  input:
  tuple val(name), path(ds, stageAs:'ds??/*'), path(lib), val(flags)

  output:
  path '*.{npy,dat}'

  script:
  """
  python -m ${lib}.tcorr msd $ds $flags
  """
}

process compute_rdf {
  publishDir "$params.publish/$name"
  input:
  tuple val(name), path(ds, stageAs:'ds??/*'), path(lib), val(flags)

  output:
  path '*.{npy,dat}'

  script:
  """
  python -m ${lib}.rdf $ds $flags
  """
}

process compute_hbnet {
  publishDir "$params.publish/$name"
  input:
  tuple val(name), path(ds, stageAs:'ds??/*'), path(lib), val(flags)

  output:
  path '*.{npy,dat}'

  script:
  """
  python -m ${lib}.tcorr hbnet $ds $flags
  """
}
