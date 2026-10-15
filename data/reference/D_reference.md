# Part D reference check

Groot and Warren, J. Chem. Phys. 107, 4423–4435 (1997), Eq. (27) and Fig. 8.
Author-uploaded full text accessed 2026-10-15:
https://www.researchgate.net/publication/257993073_Dissipative_particle_dynamics_Bridging_the_gap_between_atomistic_and_mesoscopic_simulation

In the paper, chi*N*kBT/Delta_a = (0.306 +/- 0.003)*N for short polymers.
The simulation uses reduced imposed kBT=1. The comparison in D3 therefore
plots chi*N/Delta_a versus N against (0.306 +/- 0.003)*N.
These are the analytic published fit and its stated coefficient uncertainty,
not newly digitized Fig. 8 points. N=10 is at the edge of the stated fit range.

Eq. (18), for equal lengths N, is chi*N = ln[(1-phi_A)/phi_A]/(1-2*phi_A).
The local B-rich bulk bead fraction is required, not the global 50:50 composition
and not the fraction of whole chains whose centres lie in the window.
