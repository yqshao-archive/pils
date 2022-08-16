// This file holds the custom analyis modules used in this project
nextflow.enable.dsl=2

params.publish = 'latent'

process extract_latent {
  tag "$name"
  label "latent"
  publishDir "$params.publish/$name"

  input:
  tuple val(name), path(model), path(ds), val(flags)

  output:
  tuple val(name), path('latent.npy')

  script:
  template 'latent.py'
}
