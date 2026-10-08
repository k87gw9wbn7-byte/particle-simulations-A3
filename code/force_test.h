#ifndef FORCE_TEST_H_
#define FORCE_TEST_H_

#include "structs.h"

/**
 * @brief The type of a force routine: a function that fills the force array
 * and returns the potential energy. The tests take the routine to check as a
 * function pointer of this type, so @ref calculate_forces as well as a wrapper
 * around a single force term can be tested.
 */
typedef double (*ForceKernel)(struct Parameters*, struct Nbrlist*, struct Vectors*);

/**
 * @brief Check the accumulated virial against a finite difference of the
 * potential energy under isotropic scaling of the system.
 *
 * All positions and the box are scaled by a factor lambda close to one. For
 * such a scaling the virial obeys W = sum(f.r) = -dU/dlambda at lambda = 1, so
 * the analytically accumulated virial must match the numerical derivative of
 * the potential energy. Prints the result and returns the relative error.
 *
 * @param kernel force routine to test, e.g. calculate_forces
 * @param p_parameters simulation parameters (L is temporarily scaled)
 * @param p_nbrlist neighbor list (rebuilt during and after the test)
 * @param p_vectors particle arrays (restored on return)
 * @return double relative error between analytical and finite-difference virial
 */
double virial_test(ForceKernel kernel,
                   struct Parameters *p_parameters,
                   struct Nbrlist *p_nbrlist,
                   struct Vectors *p_vectors);

/**
 * @brief Check the analytical force on particle i against a central finite
 * difference of the potential energy.
 *
 * The particle is displaced by +delta and -delta along each axis and the force
 * component is compared with -(U(+delta) - U(-delta))/(2 delta). A correct
 * force routine must reproduce minus the gradient of its potential energy;
 * this catches most mistakes in a force implementation. Use it on every force
 * term you write (set Parameters.force_test in setparameters.c).
 *
 * A force smaller than the rounding noise of the energy difference cannot be
 * tested this way and is reported as untested. This happens, for instance, on
 * a perfect lattice, where forces cancel by symmetry.
 *
 * @param kernel force routine to test, e.g. calculate_forces
 * @param i index of the particle to test
 * @param p_parameters simulation parameters
 * @param p_nbrlist neighbor list (updated as the particle is displaced)
 * @param p_vectors particle arrays (restored on return)
 * @return double relative error between analytical and finite-difference force
 */
double forces_test(ForceKernel kernel, int i,
                   struct Parameters *p_parameters,
                   struct Nbrlist *p_nbrlist,
                   struct Vectors *p_vectors);


/** @brief Bond-stretch term only, with the signature of calculate_forces (for forces_test). */
double forces_bond_only(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

double forces_conservative_only(struct Parameters *, struct Nbrlist *, struct Vectors *);
#endif
