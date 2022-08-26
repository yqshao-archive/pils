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
    ['bias10', ' --bias heaviside --kb 10'],
    ['bias100', ' --bias heaviside --kb 100'],
    ['bias500', ' --bias heaviside --kb 500'],
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
