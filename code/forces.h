#ifndef FORCES_H_
#define FORCES_H_
#include "structs.h"
/* Kernels add forces/unnormalised virials and return conservative energy.
 * calculate_forces first clears accumulators, then divides virials by 3V. */
double calculate_forces(struct Parameters *, struct Nbrlist *, struct Vectors *);
double calculate_forces_nb(struct Parameters *, struct Nbrlist *, struct Vectors *);
double calculate_forces_bond(struct Parameters *, struct Vectors *);
double dpd_noise(void);
#endif
