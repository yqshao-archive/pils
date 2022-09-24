#!/usr/bin/env python

import numpy as np

def mktopo(datum, level=0):
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

    rc = 4.0 if level<=2 else 6.0

    cutoff = {
        ("H", "C"): 2.0,
        ("H", "N"): 4.0,
        ("H", "O"): 4.0,
        ("C", "C"): 2.0,
        ("C", "N"): 2.0,
        ("C", "O"): 2.0,
    }

    # build initial nl, see ase doc [ase/neighborlsit] ===================================
    nl_i, nl_j, nl_d = neighborlist.neighbor_list("ijd", atoms, cutoff, self_interaction=False)
    conMat = sparse.dok_matrix((natoms, natoms), dtype=np.int8)
    conMat[nl_i, nl_j] = 1  # we have several running indices here prefixed by (nl, mol, h)
    conMat[nl_j, nl_i] = 1  # v---- shamelessly taken from the ase documentation
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

    if level==0:
        return H_act, O_act, N_act


def cktopo(this_topo, prev_topo, idx=None):
    """ Checks the consistency between two topologies"""
    if prev_topo is None:
        return this_topo
    topo_keys = ['h_act', 'o_act', 'n_act']
    for (key, this_t,  prev_t) in zip(topo_keys, this_topo, prev_topo):
        assert np.all(this_t == prev_t), f'{key} is not consistent: {this_t} vs {prev_t}'
    return this_topo
