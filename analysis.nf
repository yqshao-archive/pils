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

Channel.fromPath('trajs/cp2k/nvt*ps/*{1.0753,1.1551}/', type:'dir') \
  | map {traj -> [traj.name, traj]} \
  | groupTuple \
  | map {name, paths -> [(name=~/(a\d+b\d+i\d+-rho.*)/)[0][1], paths.sort()]} \
  | set {cp2k_traj}

Channel.fromPath('trajs/cp2k/nvt*ps/*{1.0753,1.1551}/', type:'dir') \
  | map {traj -> [traj.name, traj]} \
  | groupTuple \
  | map {name, paths -> [(name=~/(a\d+b\d+i\d+-rho.*)/)[0][1], paths.sort()]} \
  | set {cp2k_vali}

Channel.fromPath('exp/trial-adam/prod/gen31/nvt-340k-5ns-0/*/asemd.traj') \
  | map {path -> [path.parent.name, path]} \
  | set {pinn_prod}

workflow msd {
  lib = file('py', type:'dir')
  msd_50ps = Channel.of(['msd-10-60ps', '-w 40 -ts 10 -te 60 -s 0.01'])
  msd_20ps = Channel.of(['msd-0-20ps', '-w 10 -ts 0 -te 20 -s 0.01'])
  msd_5ns = Channel.of(['msd-0-5ns', '-w 4000 -te 5000 -s 1'])

  cp2k_traj
    | combine(msd_50ps) \
    | map {name, ds, msd, flag -> ["cp2k/$name-$msd", ds, lib, "-dt 0.0005 $flag"]} \
    | set {msd_cp2k}

  cp2k_vali
    | combine(msd_50ps.concat(msd_20ps)) \
    | map {name, ds, msd, flag -> ["vali/$name-$msd", ds, lib, "-dt 0.0005 $flag"]} \
    | set {msd_vali}

  pinn_prod
    | combine(msd_50ps.concat(msd_5ns)) \
    | map {name, ds, msd, flag -> ["trial/$name-$msd", ds, lib, "-dt 0.01 $flag"]} \
    | set {msd_prod}

  msd_cp2k | concat(msd_vali) | concat (msd_prod) | compute_diff
}

workflow rdf {
  lib = file('py', type:'dir')
  rdf_50ps = Channel.of(['rdf-10-60ps', '-ts 10 -te 60 -s 0.01'])
  rdf_20ps = Channel.of(['rdf-0-20ps', '-ts 0 -te 20 -s 0.01'])
  rdf_5ns = Channel.of(['rdf-0-5ns', '-te 5000 -s 1'])

  cp2k_traj
    | combine(rdf_50ps) \
    | map {name, ds, rdf, flag -> ["cp2k/$name-$rdf", ds, lib, "-dt 0.0005 $flag"]} \
    | set {rdf_cp2k}

  cp2k_vali
    | combine(rdf_50ps.concat(rdf_20ps)) \
    | map {name, ds, rdf, flag -> ["vali/$name-$rdf", ds, lib, "-dt 0.0005 $flag"]} \
    | set {rdf_vali}

  pinn_prod
    | combine(rdf_50ps.concat(rdf_5ns)) \
    | map {name, ds, rdf, flag -> ["trial/$name-$rdf", ds, lib, "-dt 0.01 $flag"]} \
    | set {rdf_prod}

  rdf_cp2k | concat(rdf_vali) | concat (rdf_prod) | compute_rdf
}
