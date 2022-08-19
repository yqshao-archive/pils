#!/usr/bin/env python

def piview(atoms,
           draw_bonds=False, bond_color='#000000', bond_cutoff=1.2,
           width=200, height=200, download=None,
           props=None, pmin=None, pmax=None, pmap='RdBu_r',
           inter=None, imin=None, imax=None, imap='RdBu_r', itol=0.2,
           border_color='black',
           **kwargs):
    """
    Args:
        atoms: ASE atoms object
        draw_bonds: draw covalent bonds
        width, height: viewport size, in px
        download: save the picture automatically (set a string as the filename)
        props, pmin, pmax, pmap:
           property values, when set, atoms will be colored with props
        inter, imin, imax, imap:
           interaction values, when set, half bonds will be drawn
        kwargs:
           additional attributes for PiView, they are set as variable for now
           - zoom (default: 1): scale the canvas, default fills the canvas;
           - margin (default: 0.5): controls the margin of canvas, in Ang;
           - quality (default: 24): subdivision of balls and bonds;
           - atomScale (default: 0.24): scale the radius of atoms;
           - radiusType (default: 'vdW'): atom radius type (covalent or vdW);
           - edge (default: 0.1): width of the black edge;
           - bondRadi (default: 0.15): bond width, in Ang.
    """
    import numpy as np
    import json, matplotlib.cm, random, string, os
    from matplotlib.colors import to_hex
    from IPython.core.display import HTML
    from ase.data.colors import jmol_colors
    from ase.neighborlist import neighbor_list, natural_cutoffs

    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    with open(f'{dir_path}/PiView.js') as f:
        jscript = f.read()

    # jsonify atoms
    natoms = len(atoms)
    atoms.positions = atoms.positions-atoms.positions.mean(axis=0, keepdims=True)
    atomsDict = {
        'elems': atoms.numbers.tolist(),
        'coord': atoms.positions.tolist()}

    if props is not None:
        if (pmin is None) or (pmax is None):
            pmin = -np.max(np.abs(props))
            pmax = np.max(np.abs(props))
        norm = matplotlib.colors.Normalize(vmin=pmin, vmax=pmax)
        cmap = matplotlib.cm.get_cmap(pmap)
        props = [to_hex(c) for c in cmap(norm(props))]

    if inter is not None:
        np.fill_diagonal(inter, 0)
        if (imin is None) or (pmax is None):
            imin = -np.max(np.abs(inter))
            imax =  np.max(np.abs(inter))
        norm = matplotlib.colors.Normalize(vmin=imin, vmax=imax)
        cmap = matplotlib.cm.get_cmap(imap)
        assert (inter.shape[0]==natoms) and (inter.shape[1]==natoms)
        inter = [(mi, mj, to_hex(cmap(norm(inter[mi, mj]))))
                 for mi in range(natoms) for mj in range(natoms)
                 if (mi!=mj) and np.abs(inter[mi, mj])>itol*np.abs(inter).max()]
    elif draw_bonds:
        # color scheme
        if bond_color == 'jmol':
            bcolor = lambda i: to_hex(jmol_colors[atoms.numbers[i]])
        else:
            bcolor = lambda i: bond_color
        nl =  neighbor_list('ijd', atoms, natural_cutoffs(atoms, mult=bond_cutoff))
        inter = [(int(i), int(j), bcolor(i)) for i,j in zip(nl[0],nl[1])]

    atomsDict = json.dumps(atomsDict)
    props = json.dumps(props)
    inter = json.dumps(inter)
    download = f'download={json.dumps(download)}' if download else ''

    alphabet = string.ascii_lowercase + string.digits
    uuid = ''.join(random.choices(alphabet, k=8))
    setvars = ''.join(f'var {k}={json.dumps(v)};' for k,v in kwargs.items())
    html = f"""
    <html>
    <head>
    </head>
    <body>
    <div>
       <piview {download} id='{uuid}' cstyle='border: 0.2em solid; border-color: {border_color};' height='{height}' width='{width}' atoms='{atomsDict}' props='{props}' inter='{inter}'></piview>
       <script type="module">
           var uuid='{uuid}';
           {setvars}
           {jscript}
       </script>
    </div>
    </body>
    </html>
    """
    return HTML(html)
