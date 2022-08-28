#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
This is a workflow for "babysitting" an active learning progress, several tools
are implemented as sub-workflows:

- ebmd: this use a biased ensemble md to extend a trajectory
- label: this filters a previously sampled traj and labels it with cp2k
- retrain: this generates a new generation

The naming pattern:

- models: models/$base-gen$g-seed$s
- resampled-traj: trajs/ebmd/$base-gen$g/$i
- labelled-ds: datasets/

where $base is basename of the model, $g is generation, $s is random seed for
data splitting, $i is filename of inital structure

*/

params.init = 'skel/init/*.xyz'
params.base = 'pils-v5-ekf-v3'
params.gen = 0

include { aseMD } from './tips/nextflow/ase.nf' addParams(publish: 'trajs/emd')

workflow ebmd {
  presets = [
    ['nobias', ''],
    ['bias1', ' --bias heaviside --kb 1'],
    ['bias001', ' --bias heaviside --kb 0.01'],
    // ['bias01', ' --bias heaviside --kb 0.1'],
    // ['bias3', ' --bias heaviside --kb 3'],
    //['bias10', ' --bias heaviside --kb 10'],
    // ['bias30', ' --bias heaviside --kb 30'],
    // ['bias100', ' --bias heaviside --kb 100'],
    // ['bias500', ' --bias heaviside --kb 500'],
  ]
  md = '--ensemble nvt --T 340 --t 20 --dt 0.5 --log-every 20'

  channel.fromPath("models/$params.base-gen$params.gen-seed*/model", type:'dir') \
    | map {model -> [(model.parent.name =~ /(.*)-seed\d/)[0][1], model]}  \
    | groupTuple                                                          \
    | combine( channel.fromPath(params.init) ) \
    | combine( channel.fromList(presets)) \
    | map { base, model, init, preset, bias -> \
           ["$params.base-gen$params.gen-$preset/$init.baseName", model, init, md+bias]} \
    | aseMD
}

params.bias = 'nobias'
params.cp2k_aux = 'skel/cp2k-aux/*'

include { cp2kGenInp } from './tips/nextflow/cp2k.nf' addParams(publish: 'geninp')
include { cp2k } from './tips/nextflow/cp2k.nf' addParams(publish: 'trajs/cpsp')
workflow label {
  strategies = [
    ['uniform', '-f asetraj --subsample --strategy uniform --nsample 40'],
    ['forces', '-f asetraj --subsample --strategy sorted --nsample 40'],
  ]

  skel = file('./skel/cp2k/singlepoint.inp')
  // all trajs from the  preset
  channel.fromPath("trajs/emd/$params.base-gen$params.gen-$params.bias/*/asemd.traj") \
    | map {traj -> ["${traj.parent.name}", traj]} \
    | combine( channel.fromList(strategies) ) \
    | map {name, traj, tag, flag ->  \
           ["$params.base-gen$params.gen-$params.bias/$name/$tag", skel, traj, flag]} \
    | cp2kGenInp \
    | flatMap {name, inps -> inps.collect {["$name/$it.baseName", it]}}
    | map {name, inp -> [name, inp, file(params.cp2k_aux)]}
    | cp2k
}

