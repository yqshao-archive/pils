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
    (nitro,) = np.where(atoms.numbers == 7)
    (oxyge,) = np.where(atoms.numbers == 8)
    O_act = oxyge # all oxygen are active
    natoms = len(atoms)
    heavy2r = lambda x: heavy.tolist().index(x)

    cutoff = {
        ("H", "C"): 2,
        ("H", "N"): 3,
        ("H", "O"): 3,
        ("C", "C"): 2,
        ("C", "N"): 2,
        ("C", "O"): 2,
    }

    # build initial nl, see ase doc [ase/neighborlsit] ===================================
    nl_i, nl_j, nl_d = neighborlist.neighbor_list("ijd", atoms, cutoff, self_interaction=False)
    conMat = sparse.dok_matrix((natoms, natoms), dtype=np.int8)
    conMat[nl_i, nl_j] = 1  # we have several running indices here prefixed by (nl, mol, h)
    conMat[nl_j, nl_i] = 1  # v---- shamelessly taken from the ase documentation

    mol_sets = [heavy[mol_assign == mol_i] for mol_i in range(n_mol)]
    CN_N = np.squeeze(np.asarray(conMat[nitro, :][:, heavy].sum(axis=1)))
    N_act = nitro[CN_N==2]
    ALL_act = np.concatenate([O_act, N_act])
    mol_acts = [np.intersect1d(mol_set, ALL_act) for mol_set in mol_sets]
    if check: assert check_mols(mol_sets), str(mol_sets)

    # -- zeros pass, tag active protons:
    sel0 = [np.where(nl_i == h_ia)[0] for h_ia in hydro]
    h_n0a = np.array([nl_j[_sel][np.argmin(nl_d[_sel])] for _sel in sel0])
    H_act = hydro[np.in1d(h_n0a, ALL_act)]

    # -- first pass, tag D(oners)
    sel1 = [np.where(nl_i == h_ia)[0] for h_ia in H_act]
    h_n1a = np.array([nl_j[_sel][np.argmin(nl_d[_sel])] for _sel in sel1])
    (hacid_r,) = np.where(atoms.numbers[h_n1a] == 8)
    (hbase_r,) = np.where(atoms.numbers[h_n1a] == 7)
    hactive_r = np.concatenate([hacid_r, hbase_r])
    h_n1e = atoms.numbers[h_n1a]
    h_n1d = np.array([nl_d[_sel][np.argmin(nl_d[_sel])] for _sel in sel1])
    D_mol = np.array([mol_assign[heavy2r[n1a]] for n1a in h_n1a])
    h_act1 = [mol_acts[_di] for _di in D_mol]
    D = np.concatenate([np.array(D_mol)[:,None], atoms.positions[h_n1a]], axis=1)
    A_mol = np.setdiff1d(np.arange(n_mol), D_mol)

    # -- second pass, make candidate neighbor for A
    sel2 = [_sel[~np.in1d(nl_j[_sel], _act)] for _sel, _act in zip(sel1, h_act1)]
    h_n2a = np.array([nl_j[_sel][np.argmin(nl_d[_sel])] for _sel in sel2])

    # build the HB network ================================================================
    hb_conn = sparse.dok_matrix((n_mol, n_mol), dtype=np.int8)
    hb_via = np.zeros((n_mol, n_mol))
    for h_ir, h_ia in zip(hactive_r, hydro[hactive_r]):
        if len(nl_j[nl_i == h_ia]) > 1: # nl is in abolute indices
            this_h = nl_i == h_ia
            h_n2a = nl_j[this_h][np.argsort(nl_d[this_h])[1]]
            if atoms.numbers[h_n2a] in [7, 8]:
                mol1 = mol_assign[heavy2r(h_n1a[h_ir])]
                mol2 = mol_assign[heavy2r(h_n2a)]
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
    ds = load_ds(dataset, fmt="asetraj", index='::10')
else:
    ds = load_ds(f"{dataset}/cp2k-md", fmt="cp2k", cp2k_frc="", cp2k_ener="")
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
