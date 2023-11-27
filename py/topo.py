#!/usr/bin/env python

import numpy as np

def mkconst(atoms):
    """
    returns 
    """
    
    import numpy as np
    from ase import Atoms, neighborlist
    from scipy import sparse
    from ase.constraints import Hookean

    (heavy,) = np.where(atoms.numbers != 1)
    (hydro,) = np.where(atoms.numbers == 1)
    (nitro,) = np.where(atoms.numbers == 7)
    (oxyge,) = np.where(atoms.numbers == 8)
    O_act = oxyge # all oxygen are active
    natoms = len(atoms)
    heavy2r = {k:i for i,k in enumerate(heavy)}

    rc = 5

    cutoff = {
        ("H", "C"): 2.0,
        ("H", "N"): rc,
        ("H", "O"): rc,
        ("C", "C"): 2.0,
        ("C", "N"): 2.0,
        ("C", "O"): 2.0,
    }

    # build initial nl, see ase doc [ase/neighborlsit] ===================================
    nl_i, nl_j, nl_d = neighborlist.neighbor_list("ijd", atoms, cutoff, self_interaction=False)
    conMat = sparse.dok_matrix((natoms, natoms), dtype=np.int8)
    conMat[nl_i, nl_j] = 1  # v ---- mostly taken from the ase documentation
    conMat[nl_j, nl_i] = 1  #
    n_mol, mol_assign = sparse.csgraph.connected_components(conMat[heavy, :][:, heavy])
    mol_sets = [heavy[mol_assign == mol_i] for mol_i in range(n_mol)]
    CN_N = np.squeeze(np.asarray(conMat[nitro, :][:, heavy].sum(axis=1)))
    N_act = nitro[CN_N==2]
    ALL_act = np.concatenate([O_act, N_act])
    mol_acts = [np.intersect1d(mol_set, ALL_act) for mol_set in mol_sets]
    type_mol = np.array([2-len(act) for act in mol_acts])

    #active protons:
    sel0 = [np.where(nl_i == h_ia)[0] for h_ia in hydro]
    h_n0a = np.array([nl_j[_sel][np.argmin(nl_d[_sel])] for _sel in sel0])
    H_act = hydro[np.in1d(h_n0a, ALL_act)]
    H_lig = h_n0a[np.in1d(h_n0a, ALL_act)]
    H_type = atoms.numbers[H_lig]
    
    # make constraints
    all_bonds = []
    cnt = 0 
    for i, j, j_type in zip(H_act, H_lig, H_type):
        # rb = 1.06 if j_type==7 else 1.01
        # rb = atoms.get_distance(i, j, mic=True)
        all_bonds.append(Hookean(int(i),int(j),5.,rt=1.5))
        if j_type==7:
            cnt += 1
    atoms.set_constraint(all_bonds)
    return atoms, cnt/len(H_act)

def mktopo(datum, level=0, hbcut=2.25):
    from ase import Atoms, neighborlist
    from scipy import sparse
    atoms = Atoms(datum["elem"], positions=datum["coord"], cell=datum["cell"], pbc=True)

    (heavy,) = np.where(atoms.numbers != 1)
    (hydro,) = np.where(atoms.numbers == 1)
    (nitro,) = np.where(atoms.numbers == 7)
    (oxyge,) = np.where(atoms.numbers == 8)
    O_act = oxyge # all oxygen are active
    natoms = len(atoms)
    heavy2r = {k:i for i,k in enumerate(heavy)}

    rc = 5 if level<=2 else 6.0

    cutoff = {
        ("H", "C"): 2.0,
        ("H", "N"): rc,
        ("H", "O"): rc,
        ("C", "C"): 2.0,
        ("C", "N"): 2.0,
        ("C", "O"): 2.0,
    }

    # build initial nl, see ase doc [ase/neighborlsit] ===================================
    nl_i, nl_j, nl_d = neighborlist.neighbor_list("ijd", atoms, cutoff, self_interaction=False)
    conMat = sparse.dok_matrix((natoms, natoms), dtype=np.int8)
    conMat[nl_i, nl_j] = 1  # v ---- mostly taken from the ase documentation
    conMat[nl_j, nl_i] = 1  #
    n_mol, mol_assign = sparse.csgraph.connected_components(conMat[heavy, :][:, heavy])
    mol_sets = [heavy[mol_assign == mol_i] for mol_i in range(n_mol)]
    CN_N = np.squeeze(np.asarray(conMat[nitro, :][:, heavy].sum(axis=1)))
    N_act = nitro[CN_N==2]
    ALL_act = np.concatenate([O_act, N_act])
    mol_acts = [np.intersect1d(mol_set, ALL_act) for mol_set in mol_sets]
    type_mol = np.array([2-len(act) for act in mol_acts])

    # -- zeros pass, tag active protons:
    sel0 = [np.where(nl_i == h_ia)[0] for h_ia in hydro]
    h_n0a = np.array([nl_j[_sel][np.argmin(nl_d[_sel])] for _sel in sel0])
    H_act = hydro[np.in1d(h_n0a, ALL_act)]

    level0 = (H_act, O_act, N_act)
    if level==0:
        return level0

    # -- first pass, tag D(oners)
    sel1 = [np.where(nl_i == h_ia)[0] for h_ia in H_act]
    h_n1a = np.array([nl_j[_sel][np.argmin(nl_d[_sel])] for _sel in sel1])
    (hacid_r,) = np.where(atoms.numbers[h_n1a] == 8)
    (hbase_r,) = np.where(atoms.numbers[h_n1a] == 7)
    h_n1e = atoms.numbers[h_n1a]
    h_n1d = np.array([nl_d[_sel][np.argmin(nl_d[_sel])] for _sel in sel1])
    D_mol = np.array([mol_assign[heavy2r[n1a]] for n1a in h_n1a])
    h_act1 = [mol_acts[_di] for _di in D_mol]
    D = np.concatenate([np.array(D_mol)[:,None], atoms.positions[h_n1a]], axis=1)
    A_mol = np.setdiff1d(np.arange(n_mol), D_mol)
    # -- second pass, make candidate neighbor for A
    sel2 = [_sel[~np.in1d(nl_j[_sel], _act)] for _sel, _act in zip(sel1, h_act1)]
    h_n2a = np.array([nl_j[_sel][np.argmin(nl_d[_sel])] for _sel in sel2])
    hb_conn = sparse.dok_matrix((n_mol, n_mol), dtype=np.int8)
    hb_via = np.zeros((n_mol, n_mol))
    for h_ir, h_ia in enumerate(H_act):
        this_h_to_act = (nl_i == h_ia) & np.in1d(nl_j, ALL_act)
        if np.sum(this_h_to_act) > 1:
            if np.sort(nl_d[this_h_to_act])[1]<hbcut:
                h_n2a = nl_j[this_h_to_act][np.argsort(nl_d[this_h_to_act])[1]]
                mol1 = mol_assign[heavy2r[h_n1a[h_ir]]]
                mol2 = mol_assign[heavy2r[h_n2a]]
                hb_conn[mol1, mol2] = 1 # <- connectitivy is symmetric
                hb_conn[mol2, mol1] = 1 # v- via gives the distance to H, asymmetric
                hb_via[mol1, mol2] = nl_d[(nl_i == h_ia) & (nl_j == h_n1a[h_ir])]
                hb_via[mol2, mol1] = nl_d[(nl_i == h_ia) & (nl_j == h_n2a)]
    n_hb, hb_assign = sparse.csgraph.connected_components(hb_conn)

    # here one summarized the populations counts ==========================================
    pop_info = np.zeros([10]) # frame_idx, n_mol, n_ion, HB2, HB3, ....
    pop_info[:2] = [len(hacid_r), len(hbase_r)]
    pair_info = []
    for hb_i in range(n_hb):
        size_hb = sum(hb_assign == hb_i)
        if size_hb > 1:
            pop_info[size_hb] += size_hb
        if size_hb == 2:
            pair_idx = np.where(hb_assign == hb_i)[0]
            pair_type = np.sum(type_mol[pair_idx])
            pair_info.append([pair_type, *pair_idx])
    level1 = (pop_info, np.array(pair_info),
              np.array([atoms.numbers[h_n1a]==7, h_n1a]).T)
    if level==1:
        return level0, level1


def cktopo(this_topo, prev_topo, keep=False):
    """ Checks the consistency between two topologies"""
    import warnings

    if prev_topo is None:
        return this_topo

    topo_keys = ['h_act', 'o_act', 'n_act']
    for (key, this_t,  prev_t) in zip(topo_keys, this_topo, prev_topo):
        if not keep:
            assert np.all(this_t == prev_t), f'{key} is not consistent: {this_t} vs {prev_t}'
        elif not np.all(this_t == prev_t):
            warnings.warn(f'{key} is not consistent: {this_t} vs {prev_t}')

    if keep:
        return prev_topo
    else:
        return this_topo

