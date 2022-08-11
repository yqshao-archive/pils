// This file holds the custom analyis modules used in this project
nextflow.enable.dsl=2

// parameters
params.ase_trajs = './trajs/*/asemd.traj'

// processes
process ase_traj_fmax {
  // simply extracts the max force components in each frame
  publishDir "$publish"
  label 'local'
  cache false

  input:
    path traj
    val publish

  output:
    path 'fmax.dat'

  script:
    """
    #!/usr/bin/env python
    import numpy as np
    from ase.io import read

    traj = read("$traj", index=':')
    fmax = np.array([np.abs(a.get_forces()).max() for a in traj])
    np.savetxt('fmax.dat', fmax)
    """
}

// entry points
workflow fmax {
  ch = channel.fromPath(params.ase_trajs)
    .multiMap {
      traj: it
      pub: "analyses/fmax/$it.parent.name"
    }
  ase_traj_fmax(ch.traj, ch.pub)
}
