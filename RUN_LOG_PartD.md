# Part D runs and analysis

Added on 2026-10-15. The notebook and code are a local working copy of
`k87gw9wbn7-byte/particle-simulations-A3`, cloned at commit
`4a6d4ae0e0af2b77f254fe5299d266a4f66917c5`. No changes were sent to GitHub.

## Build and tests

Every production case is a separate compilation of `code/*.c` with GCC-compatible
arguments, with the settings supplied to `setparameters.c` as compile definitions.
Each run writes `data/<tag>_manifest.json`: exact command, compiler, machine,
source hashes, physical settings, seed, duration and return status. The
manifest and the run log are the record of what actually ran.

Commands from this package's root:

```sh
python3 scripts/run_d.py check --name d_chain_check
(cd code && sh tests/run_checks.sh)
```

The new fixture verifies end-to-end bond unwrapping across periodic boundaries,
including a chain whose unwrapped end-to-end distance exceeds half a box.
Existing analytical force, finite-difference, virial, noise, chain initialization,
neighbour-list, momentum, RDF and restart checks also passed after the change.

## Production

All cases use density 3, imposed kBT=1, mass=1, cutoff=1, C=2, gamma=4.5,
sigma=3, dt=0.04, both conservative and DPD thermostat forces on, and seed 23.
All nonbonded pair factors remain 1, including bonded beads. There are no angle
or torsion forces. All chains contain only one bead type.

D1 uses 3000 beads in a cubic box of side 10, N=5, aAB=25, and a homogeneous
start. It runs for 60,000 steps, discarding 10,000. A preliminary 20,000-step
D1 run was superseded by the longer run; its outputs were replaced.

```sh
python3 scripts/run_d.py d1 --steps 60000 --equil 10000
```

Every blend uses 6000 beads in a 20x10x10 box, with initially separated A and B
halves along x. The initial runs use 60,000 steps, discarding 20,000. Samples are taken
every 50 steps and blocks span 2500 steps. Production therefore has 16 blocks
of duration 100. PDB files contain the initial and final configurations.

D2:

```sh
python3 scripts/run_d.py blend --n 3 --aab 28
python3 scripts/run_d.py blend --n 5 --aab 28
python3 scripts/run_d.py blend --n 10 --aab 28
```

D3 (the N=3, aAB=28 case above is reused):

```sh
python3 scripts/run_d.py blend --n 3 --aab 30
python3 scripts/run_d.py blend --n 3 --aab 32
python3 scripts/run_d.py blend --n 5 --aab 27
python3 scripts/run_d.py blend --n 5 --aab 28.5
python3 scripts/run_d.py blend --n 5 --aab 30
python3 scripts/run_d.py blend --n 10 --aab 26
python3 scripts/run_d.py blend --n 10 --aab 26.5
python3 scripts/run_d.py blend --n 10 --aab 27
```

The blend commands were dispatched with at most two simulations running at once.
The D1 run overlapped with the first blends. Each simulation remains a single
independent process with its own outputs; there is no threaded force loop.

After initial diagnostics, the N=5 cases at aAB=27 and 30 were extended:

```sh
python3 scripts/run_d.py blend --n 5 --aab 27 --steps 120000 --equil 40000
python3 scripts/run_d.py blend --n 5 --aab 30 --steps 120000 --equil 40000
```

The former initially had a production-half difference of 3.02 combined SE;
the latter had positive minority concentration in only 7 of 16 base blocks.
The extended runs supersede those initial outputs. Their current manifests
record the longer durations. They have 32 production blocks. During these
extensions at most three simulations overlapped.

To reproduce the initial complete sequence serially:

```sh
python3 scripts/run_d.py all
```

On Windows use `py` instead of `python3`, and add `--cc` with your GCC path.
The scripts pass compiler arguments directly, including paths with spaces.
The run helpers require only Python's standard library.

## Analysis and snapshots

```sh
python3 scripts/post_d.py all
```

This requires numpy and matplotlib and reads only saved run data. The notebook
contains the same analysis functions so its D cells do not depend on `__file__`
or on launching production runs. Run its D1 cell before D2 and D3.

Mean squared chain size is calculated in C by summing consecutive minimum-image
bond vectors. Production block means are used for uncertainty. D3 extracts the
A bead fraction in the registered B-rich central bulk window, uses Eq. (18)
to obtain chi*N, then divides by N. A block jackknife is used for the nonlinear
quantity, with reblocking factors 1,2,4. The largest supported estimated error
is retained. Correlation, window, drift, symmetry and slab diagnostics are
saved in the D tables. A warning calls for investigation, not automatic acceptance.

The D snapshots are plotted from final simulated PDB coordinates using
matplotlib. They are not OVITO screenshots. C/O labels mean bead types A/B,
not chemical carbon/oxygen. The files can also be opened in OVITO; D1's black
bonds in the rendered image connect consecutive beads, omitting boundary-spanning
segments rather than drawing them through the periodic box.

The comparison uses the analytic published Eq. (27), not invented digitized
Fig. 8 data; the source and convention are in `data/reference/D_reference.md`.

## Remaining submission work

The earlier B prose tables disagree with the committed B data, and C2's
convergence statement needs care. These pre-existing issues were not silently
rewritten while adding D. Group identity and lessons learned need the group's
own input. Review AI/source declarations and understand the implementation for
the oral. The course's `check_my_submission.py` and a fresh released notebook
are not supplied here, so no claim is made that the course submission validator
has passed. The HTML is a local notebook export; use the course's export/check
workflow when preparing the final hand-in.
