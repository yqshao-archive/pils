#!/usr/bin/env python

"""
For computing the diffusion and onsager coefficients

The main script acts as an entry point, see available commands in:
python transport.py --help
"""

import click
import numpy as np
from tips.io import load_ds
from tips.cli.common import load_opts
from .hbnet import mktopo, cktopo

def rmcom(old_pos, new_pos, cell, masses):
    if old_pos is None:
        return new_pos
    else:
        new_pos -= np.rint((new_pos-old_pos)/cell)*cell
        new_pos -= np.average(new_pos-old_pos, weights=masses, axis=0)
        return new_pos

def mkmsd(cache, msd, cnt, pos, window):
    """ updates the msd information with new position data"""
    if cache is None:
        # create the empty ones on the first pass
        cache = np.full([window, pos.shape[0], 3], np.nan)
        cache[0, :, :] = pos
        msd = np.zeros(window)
        cnt = np.zeros(window)
    else:
        new_msd = np.sum((pos[None,:,:] - cache)**2, axis=2)
        cache = np.concatenate([[pos], cache[:-1]], axis=0)
        msd += np.nansum(new_msd, axis=1)
        cnt += np.sum(~np.isnan(new_msd), axis=1)
    return cache, msd, cnt

@click.command(name="diff", short_help="self diffusion coefficients")
@click.argument("dataset", nargs=-1)
@load_opts
@click.option('-w', '--window', default=10.0, help="correlation time window")
@click.option('-s', '--stride', default=0.1, help="stride in trajectory reading")
@click.option('-ts', '--t-start', default=0.0, help="starting time")
@click.option('-te', '--t-end', default=None, type=float, help="ending time, default is None")
@click.option('-tc', '--t-check', default=1.0, help="interval to check the topology")
@click.option('-dt', default=0.0005, help="time step of traj, default is 0.5fs")
def diff(
        dataset, fmt, emap, # load_opts
        window, stride, t_start, t_end, t_check, dt,
):
    from ase.data import atomic_masses as am
    stride = round(stride/dt)
    window = round(window/dt/stride)
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

    masses = np.array([am[i] for i in ds[0]['elem']])
    prev_topo, prev_pos = None, None
    h_cache, o_cache, n_cache = [None]*3
    h_msd, o_msd, n_msd = [None]*3
    h_cnt, o_cnt, n_cnt = [None]*3
    for i, datum in enumerate(ds):
        if i%check == 0:
            this_topo = mktopo(datum, level=0)
            prev_topo = cktopo(this_topo, prev_topo)
            h_act, o_act, n_act = this_topo
        cell = np.diag(datum['cell'])[None,:]
        prev_pos = rmcom(prev_pos, datum['coord'], cell, masses)
        h_pos = prev_pos[h_act]
        o_pos = prev_pos[o_act]
        n_pos = prev_pos[n_act]
        h_cache, h_msd, h_cnt = mkmsd(h_cache, h_msd, h_cnt, h_pos, window)
        o_cache, o_msd, o_cnt = mkmsd(o_cache, o_msd, o_cnt, o_pos, window)
        n_cache, n_msd, n_cnt = mkmsd(n_cache, n_msd, n_cnt, n_pos, window)
    t = np.arange(1, 1+window)*stride*dt
    np.save('msd.npy', np.array([t, h_msd/h_cnt, o_msd/o_cnt, n_msd/n_cnt]).T)


@click.group()
def cli():
    """Transport coefficients in PILS"""
    pass
cli.add_command(diff)


if __name__ == "__main__":
    cli()
