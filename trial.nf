#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
 This is the workflow for an trial run, it has the following entry points
 `-entry lmpinit` generates the initial geometry guess with moltemplate and lammps
 `-entry cp2k` generates the cp2k trajectories or extends the trajectory
 `-entry train` trains the inital MD models and perform some basic analyses
*/

//======================================//
// THE GENERATION OF INITIAL GEOMETRIES //
//======================================//
include { moltemplate; tag2inp } from './nextflow/build.nf' addParams(publish: 'trajs/build')
include { lammpsMD } from './tips/nextflow/lammps.nf' addParams(publish: 'trajs/lmp')

params.tags = 'a32b32i0,a0b0i32,a16b16i16'
params.rhos = '1.1551,1.0753' // g cm^-3
params.lmp_inp = 'skel/lmp/init.in'

workflow lmpinit {
  channel.fromList(params.tags.tokenize(','))               \
    | combine (channel.fromList(params.rhos.tokenize(','))) \
    | map { tag, rho -> tag2inp(tag, rho) }                 \
    | moltemplate                                           \
    | map {name, aux -> [name, file(params.lmp_inp), aux]}  \
    | lammpsMD
}


//==========================================//
// THE GENERATION OF THE INITAL TRANING SET //
//==========================================//
params.cp2k_inp = 'skel/cp2k/nvt-10ps.inp'
params.cp2k_aux = 'skel/cp2k-aux/*'
params.cp2k_geo = 'trajs/lmp/*/equi.dump'
params.cp2k_from = '10'
include { cp2kMD; cp2k } from './tips/nextflow/cp2k.nf' addParams(publish: 'trajs/cp2k')

workflow cp2kinit {
  channel.fromPath(params.cp2k_geo)                    \
    | map {it ->                                       \
           ["nvt-0-10ps/${it.parent.name}", file(params.cp2k_inp), it, \
            "--fmt lammps-dump --idx -1 --emap ${file("trajs/build/${it.parent.name}/system.data")}"]} \
    | cp2kMD
}

workflow cp2krestart {
  int from = params.cp2k_from.toInteger()
  channel.fromPath("trajs/cp2k/nvt-${from-10}-${from}ps/*/cp2k-md-1_${2000*from}.restart") \
    | map {it ->                                                                    \
           ["nvt-${from}-${from+10}ps/${it.parent.name}", it, file(params.cp2k_aux)] }\
    | cp2k
}

//===============================================//
// THE INITIAL TRAINING AND EVALUATION OF MODELS //
//===============================================//
include { pinnTrain } from './tips/nextflow/pinn.nf' addParams(publish: 'models')
include { aseMD; aseMD as aseEMD } from './tips/nextflow/ase.nf' addParams(publish: 'trajs/ase')
params.datasets = './datasets/*.{yml,tfr}'
params.pinn_inp = './skel/pinn/*.yml'
params.asemd_init = './skel/init/*.xyz'
params.repeats = 5
params.adam_flags = '--train-steps 5000000 --log-every 10000 --ckpt-every 100000 --batch 1 --max-ckpts 1 --shuffle-buffer 1000 --init'
params.ekf_flags = '--train-steps 500000 --log-every 1000 --ckpt-every 10000 --batch 1 --max-ckpts 1 --shuffle-buffer 1000 --init'
params.ase_flags = '--ensemble nvt --T 300 --t 1 --dt 0.5 --log-every 20'

workflow train {
  // combine
  channel.fromFilePairs(params.datasets)         \
    | combine(channel.fromPath(params.pinn_inp)) \
    | combine(channel.of(1..params.repeats))     \
    | map {name, ds, inp, seed ->                \
           ["$name-$inp.baseName-$seed", ds, inp,\
            inp.baseName.contains('ekf') ?       \
            "--seed $seed $params.ekf_flags":    \
            "--seed $seed $params.adam_flags"]}  \
    | pinnTrain // ^- change pinn flags according to input

  pinnTrain.out.model                              \
    | combine(channel.fromPath(params.asemd_init)) \
    | map { name, mode, init -> ["$name-$init.baseName", mode, init, params.ase_flags]} \
    | aseMD

  pinnTrain.out.model                                         \
    | map {name, model -> [(name =~ /(.*)-\d/)[0][1], model]} \
    | groupTuple                                              \
    | combine(channel.fromPath(params.asemd_init))            \
    | map { name, mode, init -> ["$name-$init.baseName", mode, init, params.ase_flags]} \
    | aseEMD
}


//========================================//
// Extract Feature and Perform Clustering //
//========================================//
params.models = 'models-old/*/model'
params.latent_ds = 'datasets/samples.data'
params.latent_flags = '--take 100'

include { extract_latent } from './nextflow/latent.nf'

workflow latent {
  Channel.fromPath(params.models, type:'dir')     \
    | combine(Channel.fromPath(params.latent_ds)) \
    | map {model, ds ->                           \
           ["${model.parent.name}-$ds.baseName", model, ds, params.latent_flags]} \
    | extract_latent
}
