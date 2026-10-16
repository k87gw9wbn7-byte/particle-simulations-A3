A3 OVITO snapshot update

Copy these files into the root of your existing A3 project, replacing the same-named files.
This is an update package, not a complete submission: keep your existing code and data.

A3.ipynb and A3.html now embed four OVITO 3.14.1 / OSPRay renderings of the actual final PDB frames.
D1: N=5, aAB=25, chain bonds shown. D2: N=3,5,10, aAB=28, A blue / B orange.
All use t=2400. Task cells, numerical output and simulation files are unchanged.

The script scripts/render_d_ovito.py regenerates the images with ovito==3.14.1 and Pillow.
The normal numerical analysis only loads cached snapshots, so OVITO is optional unless rerendering.

These changes have not been pushed to GitHub.
