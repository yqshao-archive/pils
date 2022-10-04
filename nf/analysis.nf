nextflow.enable.dsl=2

params.publish = 'analyses'

process compute_diff {
  publishDir "$params.publish/$name"
  input:
  tuple val(name), path(ds, stageAs:'ds??/*'), path(lib), val(flags)

  output:
  path '*.{npy,dat}'

  script:
  """
  python -m ${lib}.transport diff $ds $flags
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
  python -m ${lib}.2drdf $ds $flags
  """
}
