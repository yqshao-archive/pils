#!/usr/bin/env nextflow

params.proj = 'exp/prod-adam-run2'
params.gens = '0,14,30' // all gen to latent analysis
params.gen = '30' //final production gen
params.latent_ds = 'datasets/pils-50ps.{yml,tfr}'
params.latent_flags = '--take 100'

//========================================//
// Extract Feature and Perform Clustering //
//========================================//
include { extract_latent } from './nf/latent.nf' addParams(publish: "$params.proj/analyses/latent")

workflow latent {
  Channel.fromList(params.gens.tokenize(',')
             .collect{g->[g, file("$params.proj/models/gen$g/model1/model", type:'dir')]}) \
    | combine(Channel.fromFilePairs(params.latent_ds)) \
    | map {g, model, dsname, ds ->                           \
           ["gen$g", model, ds, params.latent_flags]} \
    | extract_latent
}

//========================//
// Diffusion Coefficients //
//========================//
include { compute_msd; compute_rdf; compute_hbnet } from './nf/analysis.nf' addParams(publish: "$params.proj/analyses")

Channel.fromPath('trajs/cp2k/nvt*ps/*{1.0753,1.1551}/', type:'dir') \
  | map {traj -> [traj.name, traj]} \
  | groupTuple \
  | map {name, paths -> [(name=~/(a\d+b\d+i\d+-rho.*)/)[0][1], \
                         paths.sort{p-> (p=~/-(\d+)ps/)[0][1].toInteger()}]} \
  | set {cp2k_traj}

Channel.fromPath('exp/prod-adam-run2/cp2k-vali/nvt*ps/*/', type:'dir') \
  | map {traj -> [traj.name, traj]} \
  | groupTuple \
  | map {name, paths -> [(name=~/(a\d+b\d+i\d+-r.*)/)[0][1], \
                         paths.sort{p-> (p=~/-(\d+)ps/)[0][1].toInteger()}]} \
  | set {cp2k_vali}

Channel.fromPath("$params.proj/prod/gen$params.gen/nvt-*k-5ns-0/*/asemd.traj") \
  | map {path -> ["${path.parent.parent.name}/${path.parent.name}", path]} \
  | set {pinn_prod}

workflow msd {
  lib = file('py', type:'dir')
  msd0 = Channel.of(['msd-10-110ps', '-w 40 -ts 10 -te 110 -s 0.01'])
  msd1 = Channel.of(['msd-10-50ps', '-w 30 -ts 10 -te 50 -s 0.01'])
  msd2 = Channel.of(['msd-5-10ns', '-w 4000 -te 5000 -s 1'])
  msd3 = Channel.of(['msd-10-110ps', '-w 40 -ts 10 -te 110 -s 0.1'])

  cp2k_traj \
    | combine(msd0) \
    | map {name, ds, msd, flag -> ["cp2k/$name/$msd", ds, lib, "-dt 0.0005 $flag"]} \
    | set {msd_cp2k}

  cp2k_vali \
    | combine(msd1) \
    | map {name, ds, msd, flag -> ["vali/$name/$msd", ds, lib, "-dt 0.0005 $flag"]} \
    | set {msd_vali}

  pinn_prod \
    | combine(msd2.concat(msd3)) \
    | map {name, ds, msd, flag -> ["prod/$name/$msd", ds, lib, "-dt 0.1 $flag"]} \
    | set {msd_prod}

  msd_cp2k | concat(msd_vali) | concat (msd_prod) | compute_msd
}

workflow rdf {
  lib = file('py', type:'dir')
  rdf0 = Channel.of(['rdf-10-110ps', '-ts 10 -te 110 -s 0.01'])
  rdf1 = Channel.of(['rdf-10-50ps', '-ts 10 -te 50 -s 0.01'])
  rdf2 = Channel.of(['rdf-5-10ns', '-te 5000 -s 1'])
  rdf3 = Channel.of(['rdf-10-110ns', '-ts 10 -te 110 -s 0.1'])

  cp2k_traj \
    | combine(rdf0) \
    | map {name, ds, rdf, flag -> ["cp2k/$name/$rdf", ds, lib, "-dt 0.0005 $flag"]} \
    | set {rdf_cp2k}

  cp2k_vali \
    | combine(rdf1) \
    | map {name, ds, rdf, flag -> ["vali/$name/$rdf", ds, lib, "-dt 0.0005 $flag"]} \
    | set {rdf_vali}

  pinn_prod \
    | combine(rdf2.concat(rdf3)) \
    | map {name, ds, rdf, flag -> ["prod/$name/$rdf", ds, lib, "-dt 0.1 $flag"]} \
    | set {rdf_prod}

  rdf_cp2k | concat(rdf_vali) | concat (rdf_prod) | compute_rdf
}

workflow hbnet {
  lib = file('py', type:'dir')
  flag = '-w 30 -s 0.1'

  cp2k_traj \
    | map {name, ds -> ["cp2k/$name/", ds, lib, "-dt 0.0005 $flag"]} \
    | set {hbnet_cp2k}

  cp2k_vali \
    | map {name, ds -> ["vali/$name/", ds, lib, "-dt 0.0005 $flag"]} \
    | set {hbnet_vali}

  pinn_prod \
    | map {name, ds -> ["prod/$name/", ds, lib, "-dt 0.1 $flag"]} \
    | set {hbnet_prod}

  hbnet_cp2k | concat(hbnet_vali) | concat (hbnet_prod) | compute_hbnet
}

workflow {
  latent()
  msd()
  rdf()
  hbnet()
}
