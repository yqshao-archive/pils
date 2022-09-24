#!/usr/bin/env python

import click
from tips.cli.common import load_opts

def datum2hists(bins):
    from hbnet import datum2hstar, build_nl

    return non


@load_opts
@click.option("-rmin", default=0)
@click.option("-rmax", default=5)
@click.option("-rbin", default=0.5)
def main(
        dataset, fmt, emap, # load_opts
        rmin, rmix, rbin, # histogram option
):
    """main script for generating the 2d-rdf

    Takes a TIPS recognizable trajectory and returns three files:
    `NHN_rdf.dat`, `OHN_rdf.dat`, and `OHO_rdf.dat`.
    """


    np.savetxt()
