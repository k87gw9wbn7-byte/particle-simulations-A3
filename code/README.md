# Assignment 3 — Part A working version

Part A (A1–A3) is implemented. Parts B–D are unfinished. This is not a complete
assignment submission. The group details and your own review must still be added.
The user has compiled and run this version on an M1 Mac (macOS 27.0.1), passed
the focused checks, and opened/played the smoke trajectory in OVITO Basic 3.16.1.
The notebook now displays those Mac outputs.

This code starts from the group's uploaded A2-11.zip and the PBS MD starter by
E.A.J.F. Peters, J.T. Padding and Jeroen Hofman. The original course-code usage
restrictions still apply. See the notebook Sources section for AI disclosure.

## First run on your Mac

Open Terminal in this code folder, then run each command separately:

```sh
gcc -O3 *.c -o dpd -lm
./dpd > ../data/a_smoke_run.log 2> ../data/a_smoke_stderr.log
tail -n 3 ../data/a_smoke_thermo.csv
```

The data directory must exist next to code. The default is 3000 monomers in a
10×10×10 box, rho=3, aAA=aAB=aBB=25, dt=0.04, 1000 steps. This is a short
implementation check. Running it replaces the supplied a_smoke outputs. Make a
copy of this project first if you want to retain the current Mac results.
The Mac's C rand() implementation can give different trajectories for the same
seed; agreement should be statistical, not byte-for-byte.

Open `../data/a_smoke.pdb` in OVITO. A beads are labelled C and B beads O, solely
for colouring: these are coarse-grained beads, not carbon and oxygen atoms.
The file has periodic box dimensions and five frames. The user opened it in
OVITO and confirmed playback. A full visual defect assessment and equilibration
assessment are not claimed by this check.

## Where to learn and edit

1. `setparameters.c`: all settings. Recompile after each change. DPD reduced units
   have rc=m=kBT=1. `conservative_on` controls C; `thermostat_on` controls D and R
   together. `noise_sigma` is derived from sqrt(2 gamma kBT), so it is 3.
2. `forces.c`: the pair loop, one random number and one equal/opposite force per
   unordered interacting pair. Only C contributes to the pair energy/virial.
   The attractive zero-rest-length spring has C=2. No angle, torsion or Berendsen
   force remains.
3. `main.c`: half kick → drift → wrap → neighbour update → force → half kick.
   D uses the velocity after the first half kick (Groot–Warren lambda=1/2).
   R is drawn initially and once at each new force evaluation, never redrawn
   between the two kicks that use the same stored force.
4. `initialise.c`: `chain_length` beads and N−1 bonds per chain. All beads of a
   chain have one type. A chains come first. The A-chain count is rounded from
   `fraction_A * number_of_chains`. Random-walk step length is 0.8; it is not
   the spring rest length. `separated_start=1` confines initial A/B chains to
   their left/right halves by rejecting crossing steps. This initial state
   must subsequently equilibrate. Non-bonded exclusion factors are all 1.
5. `analysis.c`: AA, AB and BB RDFs, dr=0.02, up to 3. A separate neighbour list
   extends to 3 and is rebuilt only when sampling. Small boxes with fewer than
   three RDF cells per dimension use an exact pair loop. All pairs, including
   intrachain pairs, are counted once. Normalisation uses the exact shell volume
   and possible pair counts NA(NA−1)/2, NA NB, NB(NB−1)/2. An absent pair species
   has NaN g(r), because its normalisation is undefined. A run with no samples
   likewise writes NaNs. Block and total histograms include counts and frames.

The original vector helpers, neighbour-list algorithm, topology derivation,
position update, half kick and wrapping are retained. Unused angle/dihedral
connectivity is still built, as specified by the assignment.

## Reproduce the Part A checks

Run from code:

```sh
sh tests/run_checks.sh > ../data/a_checks.log 2> ../data/a_checks_stderr.log
```

This checks pair formulas/switches, finite-difference forces and virial, noise
moments, chain topology and initial positions, momentum, force neighbours against
all pairs after dynamics, RDF neighbours against all pairs, ideal uniform RDF
normalisation, small-box RDF fallback and restart I/O. The check executable is
temporary and deleted. The uniform RDF check uses 100 independent uniform
configurations of 500 monomers, with 150 A and 350 B. It is a normalisation unit
check, not the thermostatted ideal-gas MD experiment required by task B2.

Optional memory/undefined-behaviour check (supported by Clang/GCC):

```sh
ASAN_OPTIONS=detect_leaks=0 CFLAGS='-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra' sh tests/run_checks.sh > ../data/a_sanitizer.log 2> ../data/a_sanitizer_stderr.log
```

Leak checking was disabled because the assistant's sandbox does not allow the
process inspection LeakSanitizer requires. Address/undefined-behaviour checks
were run; no leak-check result is claimed.

## Restart and later work

The A2 binary layout is retained: bead count, positions, velocities, forces.
The loader now rejects a mismatched count rather than reallocating arrays behind
an already allocated neighbour list. It does not store box, chain settings,
parameters or RNG state. Use matching settings and an A3-generated state;
restarting redraws the initial random forces and does not reproduce a bitwise
continuation. The binary format uses native size_t and doubles; generate your own
restart on the machine where it will be used. The supplied smoke restart is from
the user's M1 Mac. Do not use an A2 physical state as a DPD equilibrium state.

For each future production run, save the exact parameter configuration and use a
new output prefix/filenames. Do not interpret the smoke RDF as a converged paper
comparison. Velocity histograms, composition profiles, diagonal pressure tensor,
chain-size analysis and B–D production results still need implementation/work.

The included HTML was exported with nbconvert. The course's nb2html.py and
check_my_submission.py were not supplied here. Use both official helpers on the
finished assignment, as the notebook instructs. No official submission-checker
pass is claimed.
