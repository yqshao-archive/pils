// This file holds the custom analyis modules used in this project
nextflow.enable.dsl=2

process hbonds {
  // count molecules, hydrogen bond (networks), and local geometry of X-OAc-X
  publishDir "analyses/hbonds/$name"
  label 'local'

  input:
  tuple val(name), path(dataset)

  output:
  path '*.{npy,dat}'

  script:
  """
  #!/usr/bin/env python
  import numpy as np
  from ase import Atoms, neighborlist
  from ase.io import write
  from scipy import sparse
  from tips.io import load_ds

  def process_datum(idx, datum):
      atoms = Atoms(datum["elem"], positions=datum["coord"], cell=datum["cell"], pbc=True)
      atoms.wrap()

      (heavy,) = np.where(atoms.numbers != 1)
      (hydro,) = np.where(atoms.numbers == 1)
      natoms = len(atoms)
      heavya2r = lambda x: heavy.tolist().index(x)

      cutoff = {
          ("H", "C"): 2,
          ("H", "N"): 2,
          ("H", "O"): 2,
          ("C", "C"): 2,
          ("C", "N"): 2,
          ("C", "O"): 2,
      }

      # build initial nl, see ase doc [ase/neighborlsit] ===================================
      nl_i, nl_j, nl_d = neighborlist.neighbor_list("ijd", atoms, cutoff, self_interaction=False)
      conMat = sparse.dok_matrix((natoms, natoms), dtype=np.int8)
      conMat[nl_i, nl_j] = 1  # we have several running indices here prefixed by (nl, mol, h)
      conMat[nl_j, nl_i] = 1  # v---- shamelessly taken from the ase documentation
      n_mol, mol_assign = sparse.csgraph.connected_components(conMat[heavy, :][:, heavy])
      symbs = [str(atoms[heavy][mol_assign == mol_i].symbols) for mol_i in range(n_mol)]
      cnt_mol = [symbs.count(k) for k in set(symbs)]
      assert set(symbs) == set(["C2O2", "CNC2NC"]), f"broken molecules: {set(symbs)} @ {idx}"
      type_mol = np.array(["N" in symb for symb in symbs], int)
      # ======================= designate the barebond with the notion 0->[H]OAc; 1: [H]C1Mim

      # attaching hydrogen to the heavy atoms, here we count the ions ======================
      h_n1a = np.array([nl_j[nl_i == h_ia][np.argmin(nl_d[nl_i == h_ia])] for h_ia in hydro])
      # ^----- naming scheme n1: first neighbor, a: abolute indices (instead of relative)
      (hacid_r,) = np.where(atoms.numbers[h_n1a] == 8)
      (hbase_r,) = np.where(atoms.numbers[h_n1a] == 7)
      hactive_r = np.concatenate([hacid_r, hbase_r])
      assert len(hactive_r) == 32, f"wrong number of active hydrogen ({len(hactive_r)}) @ {idx}"
      # =====================================================================================

      # build the HB network ================================================================
      hb_conn = sparse.dok_matrix((n_mol, n_mol), dtype=np.int8)
      hb_via = np.zeros((n_mol, n_mol))
      for h_ir, h_ia in zip(hactive_r, hydro[hactive_r]):
          if len(nl_j[nl_i == h_ia]) > 1: # nl is in abolute indices
              this_h = nl_i == h_ia
              h_n2a = nl_j[this_h][np.argsort(nl_d[this_h])[1]]
              if atoms.numbers[h_n2a] in [7, 8]:
                  mol1 = mol_assign[heavya2r(h_n1a[h_ir])]
                  mol2 = mol_assign[heavya2r(h_n2a)]
                  hb_conn[mol1, mol2] = 1 # <- connectitivy is symmetric
                  hb_conn[mol2, mol1] = 1 # v- via gives the distance to H, asymmetric
                  hb_via[mol1, mol2] = nl_d[(nl_i == h_ia) & (nl_j == h_n1a[h_ir])]
                  hb_via[mol2, mol1] = nl_d[(nl_i == h_ia) & (nl_j == h_n2a)]
      n_hb, hb_assign = sparse.csgraph.connected_components(hb_conn)
      # =====================================================================================

      # here one process the HB netwrok, first find the xax sites ===========================
      (center,) = np.where((hb_conn.toarray().sum(axis=0).reshape(-1) == 2) & (~type_mol))
      if center.shape[0] == 0:
          xax_info = np.zeros([0,5]) # edge cases
      else:
          left, right = np.array([np.where(row != 0)[0] for row in hb_conn.toarray()[center]]).T
          baa = type_mol[left] & (~type_mol[right])
          aab = (~type_mol[left]) & type_mol[right]
          left[baa], right[baa] = right[baa], left[baa] # swap the asymmetric configurations
          type_triplet = type_mol[left] + type_mol[right]
          xax_info = np.array([type_triplet, # type, d1, d2, d3, d4
                               hb_via[left, center],hb_via[center, left],
                               hb_via[center, right],hb_via[right, center]]).T
      # =====================================================================================

      # here one summarized the populations counts ==========================================
      pop_info = np.zeros([10]) # frame_idx, n_mol, n_ion, HB2, HB3, ....
      pop_info[:3] = [idx, len(hacid_r), len(hbase_r)]
      for hb_i in range(n_hb):
          size_hb = sum(hb_assign == hb_i)
          if size_hb > 1:
              pop_info[size_hb+1] += size_hb
      # =====================================================================================

      return xax_info, pop_info[None, :]

  dataset = "$dataset"
  if 'traj' in dataset:
      ds = load_ds(dataset, fmt="asetraj")
  else:
      ds = load_ds(dataset, fmt="cp2k", cp2k_frc="", cp2k_ener="")
  all_xax, all_pop = np.zeros([0,5]), np.zeros([0,10])

  faulty_frames = []
  for idx, datum in enumerate(ds[::]):
      try:
          xax_info, pop_info = process_datum(idx, datum)
          all_xax = np.append(all_xax, xax_info, axis=0)
          all_pop = np.append(all_pop, pop_info, axis=0)
      except:
          faulty_frames.append(idx)
          # if len(faulty_frames)>=1000:
          #     print(f"too many faulty frames: {faulty_frames}... aborting")
          #     break

  np.save('xax.npy', all_xax)
  np.save('pop.npy', all_pop)
  np.savetxt('faulty.dat', faulty_frames)
  """
}

params.cp2k_projs = './trajs/cp2k/nvt*ps/hoac-*/'
params.al_trajs = './trajs/al-adam1-sin-run1-gen25/nvt-340k-512ps/*/asemd.traj'
params.which = 'al'

workflow {
  if (params.which=='al') {
     channel.fromPath(params.al_trajs) \
       | map {traj -> ["${traj.parent.parent.parent.name}/${traj.parent.parent.name}/${traj.parent.name}", traj]} \
       | hbonds
  }
  if (params.which=='cp2k') {
     channel.fromPath(params.cp2k_projs, type:'dir') \
       | map {proj -> ["${proj.parent.name}/${proj.name}", proj]} \
       | hbonds
  }
}
