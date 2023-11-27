#!/usr/bin/env nextflow

params.proj = 'exp/scan'
params.cp2k = 'nvt-scan'
params.gens = '0,25,49' // all gen to latent analysis
params.gen = '49' //final production gen
params.latent_ds = 'datasets/scan.{yml,tfr}'
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

Channel.fromPath("trajs/cp2k/${params.cp2k}*ps/*{1.0753,1.1551}/", type:'dir') \
  | map {traj -> [traj.name, traj]} \
  | groupTuple \
  | map {name, paths -> [(name=~/(a\d+b\d+i\d+-rho.*)/)[0][1], \
                         paths.sort{p-> (p=~/-(\d+)ps/)[0][1].toInteger()}]} \
  | set {cp2k_traj}

Channel.fromPath("$params.proj/prod/gen$params.gen/nvt-*k*-0/*/asemd.traj") \
  | map {path -> ["${path.parent.parent.name}/${path.parent.name}", path]} \
  | set {pinn_prod}

workflow msd {
  lib = file('py', type:'dir')
  msd1 = Channel.of(['msd-10-110ps', '-w 40 -ts 10 -te 110 -s 0.01'])
  msd2 = Channel.of(['msd-0-5ns', '-w 4000 -te 5000 -s 1'])

  cp2k_traj \
    | combine(msd1) \
    | map {name, ds, msd, flag -> ["cp2k/${params.cp2k}/$name/$msd", ds, lib, "-dt 0.0005 $flag"]} \
    | set {msd_cp2k}

  pinn_prod \
    | combine(msd2) \
    | map {name, ds, msd, flag -> ["prod/$name/$msd", ds, lib, "-dt 0.1 $flag"]} \
    | set {msd_prod}

  msd_prod | concat(msd_cp2k) | compute_msd
}

workflow rdf {
  lib = file('py', type:'dir')
  rdf1 = Channel.of(['rdf-10-110ps', '-ts 10 -te 110 -s 0.01'])
  rdf2 = Channel.of(['rdf-0-5ns', '-te 5000 -s 1'])

  cp2k_traj \
    | combine(rdf1) \
    | map {name, ds, rdf, flag -> ["cp2k/${params.cp2k}/$name/$rdf", ds, lib, "-dt 0.0005 $flag"]} \
    | set {rdf_cp2k}
  pinn_prod \
    | combine(rdf2) \
    | map {name, ds, rdf, flag -> ["prod/$name/$rdf", ds, lib, "-dt 0.1 $flag"]} \
    | set {rdf_prod}

  rdf_prod | concat(rdf_cp2k) | compute_rdf
}

workflow hbnet {
  lib = file('py', type:'dir')
  flag = '-w 30 -s 0.1'

  cp2k_traj \
    | map {name, ds -> ["cp2k/$name/", ds, lib, "-dt 0.0005 $flag"]} \
    | set {hbnet_cp2k}

  pinn_prod \
    | map {name, ds -> ["prod/$name/", ds, lib, "-dt 0.1 $flag"]} \
    | set {hbnet_prod}

  hbnet_prod | concat(hbnet_cp2k) | compute_hbnet
}

workflow {
  // latent()
  msd()
  rdf()
  hbnet()
}
