#ifndef INITIALISE_H_
#define INITIALISE_H_
#include "structs.h"
/* Consecutive beads within each chain; chain_length=1 has no bonds. */
void initialise_bond_connectivity(struct Parameters *, struct Vectors *);
/* Original A2 topology builder. Angle/dihedral lists are retained but unused
 * by forces. DPD excludes no 1-2, 1-3 or 1-4 pairs. */
void initialise_structure(struct Parameters *, struct Vectors *, struct Nbrlist *);
void initialise_types(struct Parameters *, struct Vectors *);
void initialise(struct Parameters *, struct Vectors *, struct Nbrlist *, size_t *, double *);
/* Isotropic random walks; optional separation into left/right half-boxes. */
void initialise_positions(struct Parameters *, struct Vectors *);
void initialise_velocities(struct Parameters *, struct Vectors *);
#endif
