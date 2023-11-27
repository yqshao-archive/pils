#!/usr/bin/env python

import click
import numpy as np
from tips.io import load_ds
from tips.cli.common import load_opts
from .topo import mktopo, cktopo

def mic_dist(pos1, pos2, cell):
    """Minimal image convention, assuming orthogonal cell"""
    diff = pos1[:, None, :] - pos2[None,:,:]
    diff -= np.rint(diff/cell[None,None,:])*cell[None,None,:]
    dist = np.sqrt(np.sum(diff**2, axis=2))
    return dist

@click.command(name="rdf", short_help="self diffusion coefficients")
@click.argument("dataset", nargs=-1)
@load_opts
@click.option('-s', '--stride', default=0.1, help="stride in trajectory reading")
@click.option('-ts', '--t-start', default=0.0, help="starting time")
@click.option('-te', '--t-end', default=None, type=float, help="ending time, default is None")
@click.option('-tc', '--t-check', default=1.0, help="interval to check the topology")
@click.option('-dt', default=0.0005, help="time step of traj, default is 0.5fs")
@click.option('-rmin', default=0.0, help="Histogram bin min")
@click.option('-rmax', default=5, help="Histogram bin max")
@click.option('-rbin', default=0.025, help="Histogram bin width")
def main(
        dataset, fmt, emap, # load_opts
        stride, t_start, t_end, t_check, dt,
        rmin, rmax, rbin, # histogram option
):
    """main script for generating the 2d-rdf

    Takes a TIPS recognizable trajectory and returns the
    2D-RDF: NHN_rdf.dat, OHN_rdf.dat, and OHO_rdf.dat;
    1D-RDF: HN_rdf.dat, HO_rdf.dat.
    """
    from itertools import combinations, product

    stride = round(stride/dt)
    start  = round(t_start/dt)
    end    = round(t_end/dt if t_end else -1)
    check  = round(t_check/dt)
    if len(dataset)==1 and dataset[0].endswith('.traj'):
        ds = load_ds(dataset[0], fmt='asetraj', index=f'{start}:{end}:{stride}')
    else:
        ds = load_ds(f'{dataset[0]}/cp2k-md', fmt='cp2k')
        for restart in dataset[1:]:
            ds = ds.join(load_ds(f'{restart}/cp2k-md', fmt='cp2k')[1:])
        ds = ds[start:end:stride]

    bins = np.arange(rmin, rmax+1e-6, rbin)
    hn_cnt, ho_cnt, nhn_cnt, nho_cnt, oho_cnt = [0]*5

    prev_topo = None
    for i, datum in enumerate(ds):
        if i%check == 0:
            this_topo = mktopo(datum, level=0)
            prev_topo = cktopo(this_topo, prev_topo, keep=True)
            h_act, o_act, n_act = prev_topo

        h_pos = datum['coord'][h_act]
        o_pos = datum['coord'][o_act]
        n_pos = datum['coord'][n_act]

        cell = np.diag(datum['cell'])
        hn_dist = mic_dist(h_pos, n_pos, cell)
        ho_dist = mic_dist(h_pos, o_pos, cell)

        hn_cnt += np.histogram(hn_dist, bins)[0]
        ho_cnt += np.histogram(ho_dist, bins)[0]
        n_idx = np.arange(len(n_act))
        o_idx = np.arange(len(o_act))
        ni, nj = np.array([*combinations(n_idx, 2)]).T
        oi, oj = np.array([*combinations(o_idx, 2)]).T
        nhn_cnt += np.histogram2d(hn_dist[:,ni].flatten(), hn_dist[:,nj].flatten(), bins)[0]
        oho_cnt += np.histogram2d(ho_dist[:,oi].flatten(), ho_dist[:,oj].flatten(), bins)[0]
        ni, oj = np.array([*product(n_idx, o_idx)]).T
        nho_cnt += np.histogram2d(hn_dist[:,ni].flatten(), ho_dist[:,oj].flatten(), bins)[0]

    vol = 4*np.pi/3*(bins[1:]**3 - bins[:-1]**3)
    rmid = (bins[1:]+bins[:-1])/2
    hn_cnt=hn_cnt/(i*len(n_act)*len(h_act)*vol/np.prod(cell))
    ho_cnt=ho_cnt/(i*len(n_act)*len(h_act)*vol/np.prod(cell))
    nhn_cnt=nhn_cnt/i
    nho_cnt=nho_cnt/i
    oho_cnt=oho_cnt/i
    np.savetxt('HN_rdf.dat', np.array([rmid, hn_cnt]).T)
    np.savetxt('HO_rdf.dat', np.array([rmid, ho_cnt]).T)
    np.savetxt('NHN_rdf.dat', nhn_cnt)
    np.savetxt('NHO_rdf.dat', nho_cnt)
    np.savetxt('OHO_rdf.dat', oho_cnt)

if __name__ == "__main__":
    main()
