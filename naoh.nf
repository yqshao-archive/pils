#!/usr/bin/env nextflow

process tame_onsager {
  publishDir "data/naoh_onsager/$traj.baseName"

  input:
  path traj

  output:
  path "corr*.dat"

  """
  tame onsager $traj/CSVR-5ns-0.01ps-E50/RUN{0..3}/prod.dump.xz --tags '1,1 2,2 3,3 1,2 1,3 2,3' -dt 0.5 --max-dt 150.0 --seg 4000
  """
}

process tame_cell {
  publishDir "data/naoh_onsager/$traj.baseName"

  input:
  path traj

  output:
  path "cell.dat"

  """
  #!/usr/bin/env python
  import numpy as np
  from tame.io import load_traj

  traj = load_traj("$traj/CSVR-5ns-0.01ps-E50/RUN0/prod.dump.xz")
  np.savetxt('cell.dat', traj.data['cell'].val)
  """
}

workflow {
  ch_trajs = Channel.fromPath( '/home/yunqi/proj_cond/2019_NaOH/traj/naoh_h2o-runner/*', type:'dir' )
  // ch_trajs | tame_onsager
  ch_trajs | tame_cell
}
