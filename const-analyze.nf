#!/usr/bin/env nextflow

params.exp = "exp/scan"

Channel.fromPath("${params.exp}/const/nvt-340k-5ns-0/*/asemd.traj") \
  | map {path -> ["${path.parent.name}", path]} \
  | set {const_prod}

include { compute_msd } from './nf/analysis.nf' addParams(publish: "$params.exp/const/analyses")

workflow msd {
  lib = file('py', type:'dir')
  msd2 = Channel.of(['msd-0-5ns', '-w 4000 -te 5000 -s 1'])

  const_prod \
    | combine(msd2) \
    | map {name, ds, msd, flag -> ["$name/$msd", ds, lib, "-dt 0.1 $flag"]} \
    | set {msd_prod}

  msd_prod | compute_msd
}

workflow {
  msd()
}

