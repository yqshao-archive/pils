nextflow.enable.dsl=2

process compute_diff {
  publishDir "analyses/diff/$name"
  input:
  tuple val(name), path(ds, stageAs:'ds??/*'), path(lib), val(flags)

  output:
  path '*.{npy,dat}'

  script:
  """
  python -m ${lib}.transport diff $ds $flags
  """
}
