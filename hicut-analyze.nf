// temporary analysis for hi-cut cases

params.proj = 'preview'
params.prefix = 'nvt-hicut'

include { compute_msd; compute_rdf; compute_hbnet } from './nf/analysis.nf' addParams(publish: "$params.proj/analyses/$params.prefix")

Channel.fromPath("trajs/cp2k/${params.prefix}*ps/*{1.0753,1.1551}/", type:'dir') \
  | map {traj -> [traj.name, traj]} \
  | groupTuple \
  | map {name, paths -> [(name=~/(a\d+b\d+i\d+-rho.*)/)[0][1], \
                         paths.sort{p-> (p=~/-(\d+)ps/)[0][1].toInteger()}]} \
  | set {cp2k_traj}

workflow msd {
  lib = file('py', type:'dir')
  msd0 = Channel.of(['msd-10-110ps', '-w 40 -ts 10 -te 110 -s 0.01'])

  cp2k_traj \
    | combine(msd0) \
    | map {name, ds, msd, flag -> ["cp2k/$name/$msd", ds, lib, "-dt 0.0005 $flag"]} \
    | set {msd_cp2k}

  msd_cp2k | compute_msd
}

workflow rdf {
  lib = file('py', type:'dir')
  rdf0 = Channel.of(['rdf-10-110ps', '-ts 10 -te 110 -s 0.01'])

  cp2k_traj \
    | combine(rdf0) \
    | map {name, ds, rdf, flag -> ["cp2k/$name/$rdf", ds, lib, "-dt 0.0005 $flag"]} \
    | set {rdf_cp2k}

  rdf_cp2k | compute_rdf
}

workflow hbnet {
  lib = file('py', type:'dir')
  flag = '-w 30 -s 0.1'

  cp2k_traj \
    | map {name, ds -> ["cp2k/$name/", ds, lib, "-dt 0.0005 $flag"]} \
    | set {hbnet_cp2k}

  hbnet_cp2k | compute_hbnet
}

workflow {
  msd()
  rdf()
  hbnet()
}
