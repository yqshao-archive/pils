#!/usr/bin/env nextflow

//========================================//
// Extract Feature and Perform Clustering //
//========================================//
params.models = 'exp/adam-singlet-run1/models/gen6/*/model'
params.latent_ds = 'datasets/pils-v6.{yml,tfr}'
params.latent_flags = '--take 100'

include { extract_latent } from './nf/latent.nf'

workflow latent {
  Channel.fromPath(params.models, type:'dir')     \
    | combine(Channel.fromFilePairs(params.latent_ds)) \
    | map {model, dsname, ds ->                           \
           ["${model.parent.name}-${dsname}", model, ds, params.latent_flags]} \
    | extract_latent
}

//========================//
// Diffusion Coefficients //
//========================//
include { compute_diff; compute_rdf } from './nf/analysis.nf'

workflow diffusion {
  lib = file('py', type:'dir')
  msd_50ps = Channel.of(
    ['msd-10-60ps', '-w 40 -ts 10 -te 60 -s 0.01'],
    // ['msd-10-20ps', '-w 9 -ts 10 -te 20 -s 0.01'],
    // ['msd-20-30ps', '-w 9 -ts 20 -te 30 -s 0.01'],
    // ['msd-30-40ps', '-w 9 -ts 30 -te 40 -s 0.01'],
    // ['msd-40-50ps', '-w 9 -ts 40 -te 50 -s 0.01'],
    // ['msd-50-60ps', '-w 9 -ts 50 -te 60 -s 0.01'],
  )
  msd_5ns = Channel.of(['msd-0-5ns', '-w 4000 -te 5000 -s 1'])

  Channel.fromPath('trajs/cp2k/nvt*ps/*{1.0753,1.1551}/', type:'dir') \
    | map {traj -> [traj.name, traj]} \
    | groupTuple \
    | map {name, paths -> [(name=~/(a\d+b\d+i\d+-rho.*)/)[0][1], paths.sort()]} \
    | combine(msd_50ps) \
    | map {name, ds, msd, flag -> ["cp2k/$name-$msd", ds, lib, "-dt 0.0005 $flag"]} \
    | set {ds_cp2k}

  Channel.fromPath('exp/trial-adam/prod/gen31/nvt-340k-5ns-0/*/asemd.traj') \
    | map {path -> [path.parent.name, path]} \
    | combine(msd_50ps.concat(msd_5ns)) \
    | map {name, ds, msd, flag -> ["trial/$name-$msd", ds, lib, "-dt 0.01 $flag"]} \
    | set {ds_trial}

  ds_cp2k | concat(ds_trial) | compute_diff
}

workflow rdf {
  lib = file('py', type:'dir')
  rdf_50ps = Channel.of(
    ['rdf-10-60ps', '-ts 10 -te 60 -s 0.01'],
  )
  rdf_5ns = Channel.of(['rdf-0-5ns', '-te 5000 -s 1'])

  Channel.fromPath('trajs/cp2k/nvt*ps/*{1.0753,1.1551}/', type:'dir') \
    | map {traj -> [traj.name, traj]} \
    | groupTuple \
    | map {name, paths -> [(name=~/(a\d+b\d+i\d+-rho.*)/)[0][1], paths.sort()]} \
    | combine(rdf_50ps) \
    | map {name, ds, rdf, flag -> ["cp2k/$name-$rdf", ds, lib, "-dt 0.0005 $flag"]} \
    | set {ds_cp2k}

  Channel.fromPath('exp/trial-adam/prod/gen31/nvt-340k-5ns-0/*/asemd.traj') \
    | map {path -> [path.parent.name, path]} \
    | combine(rdf_50ps.concat(rdf_5ns)) \
    | map {name, ds, rdf, flag -> ["trial/$name-$rdf", ds, lib, "-dt 0.01 $flag"]} \
    | set {ds_trial}

  ds_cp2k | concat(ds_trial) | compute_rdf
}
