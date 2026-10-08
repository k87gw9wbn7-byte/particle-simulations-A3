#ifndef DYNAMICS_H_
#define DYNAMICS_H_
#include "structs.h"
void update_positions(struct Parameters *, struct Nbrlist *, struct Vectors *);
double update_velocities_half_dt(struct Parameters *, struct Nbrlist *, struct Vectors *);
void boundary_conditions(struct Parameters *, struct Vectors *);
void temperature_per_type(struct Parameters *, struct Vectors *, double *);
#endif
