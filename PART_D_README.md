# Your Part D working package

Open `A3.html` to read the report, or `A3.ipynb` to inspect and rerun the analysis.
The added D results use simulations run from this package's code. Original A–C
content is retained, with an updated status note and explicit source declaration.
`RUN_LOG_PartD.md` records the commands and remaining submission work.

What each part is checking:

- D1: Do the springs produce sensible chain sizes? For five ideal beads linked
  by four springs the expected mean squared end-to-end distance is 6. The
  repelling DPD beads need not give exactly that value.
- D2: At the same unlike-bead repulsion, does making chains longer strengthen
  separation? Compare all three profiles on the same axes and inspect the
  actual simulated snapshots.
- D3: Convert the small minority concentration in each bulk phase into chi*N,
  divide by N to obtain chi, and compare with the paper. Mixing becomes less
  favourable for longer chains because their beads cannot mix independently.

Important oral-exam points:

1. N counts beads per chain; total bead density remains 3 in every case.
2. The whole box is always 50:50 A/B. The chi formula needs the local minority
   bead fraction in a flat bulk region, not that overall fraction.
3. Eq. (18) returns chi*N for chains. Forgetting to divide by N mislabels chi.
4. Two interfaces occur because the box is periodic.
5. Starting from separated halves does not prove equilibrium. Check the
   production histories and the diagnostic flags as well as the pictures.
6. Error bars describe finite sampling. They do not account for every finite-box
   or bulk-window bias.
7. The snapshot colours are A/B bead labels, not chemical atom identities.

The current C files build from `code/` with `gcc -O3 *.c -o dpd -lm`. With no
run-specific definitions that command retains the old short Part A default;
use `scripts/run_d.py` to reproduce the D settings without editing the defaults.
Production runs can take several minutes each. Notebook cells only analyze data.

This is a working package, not a final submission: review the earlier B number
mismatches, C2's convergence qualification, group information, lessons learned,
and the course submission checker before handing in.
