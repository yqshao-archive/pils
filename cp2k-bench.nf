#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { cp2kGenInp } from './tips/nextflow/cp2k.nf' addParams(publish: './')

workflow {
  ch_opts = Channel.fromList(
    [
      [ 'N4-T32',  '--nodes 4 --ntasks-per-node 32 --cpus-per-task 1'],
      [ 'N8-T32',  '--nodes 8 --ntasks-per-node 32 --cpus-per-task 1'],
      ['N16-T32', '--nodes 16 --ntasks-per-node 32 --cpus-per-task 1'],
      ['N32-T32', '--nodes 32 --ntasks-per-node 32 --cpus-per-task 1'],
   ]
  )

  ch_tags = Channel.fromList(['nvt-scan', 'nvt-hicut'])
  ch_geos = Channel.fromPath('trajs/lmp/*a16b16*/equi.dump')
  ch_aux = Channel.value('skel/cp2k-aux/*')

  ch_geos\
    | combine(ch_tags) \
    | map {geo, tag -> \
           ["bench/$tag-${geo.parent.name}", file("skel/cp2k/$tag-100fs.inp"), geo, \
            "--fmt lammps-dump --idx -1 --emap ${file("trajs/build/${geo.parent.name}/system.data")}"]} \
    | cp2kGenInp \
    | combine(ch_aux) \
    | combine(ch_opts) \
    | map {name, inp, aux, opt, copt -> \
           [name, inp, file(aux), opt, copt]} \
    | cp2kBench
}

process cp2kBench {
  tag "$name"
  scratch '/scratch/local'
  module 'CP2K/9.1-psmp-PLUMED'
  publishDir "$name-$opt"
  clusterOptions "-A SNIC2022-1-27 $copt"

  input:
    tuple val(name), path(input), path(aux), val(opt),  val(copt)

  output:
    tuple val(name), path('cp2k.log'), emit:logs
    tuple val(name), path('*.{ener,xyz,stress,cell}'), emit:traj, optional: true
    tuple val(name), path('*.restart'), emit:restart, optional:true

  script:
    """
    #!/bin/bash
    $params.cp2k_cmd -i $input | tee cp2k.log
    """
}
