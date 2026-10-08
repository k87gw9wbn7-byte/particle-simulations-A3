#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "force_test.h"
#include "nbrlist.h"
#include "forces.h"

double forces_test(ForceKernel kernel, int i,
                   struct Parameters *p_parameters,
                   struct Nbrlist *p_nbrlist,
                   struct Vectors *p_vectors)
{
    // Allocate temporary force buffer (zeroed)
    Vec3D *f_comp = calloc(p_parameters->num_part, sizeof(Vec3D));
    if (!f_comp)
    {
        fprintf(stderr, "forces_test: allocation failed\n");
        return -1.0;
    }
    Vec3D *f_saved = p_vectors->f;
    p_vectors->f = f_comp;

    Vec3D *r = p_vectors->r;
    Vec3D r_saved = r[i];
    Vec3D *dr_i = &(p_vectors->dr[i]);

    const double delta = 1e-6;
    Vec3D f_fd = (Vec3D){0.0, 0.0, 0.0}; // initialize all components
    double e_scale = 0.0; // magnitude of the energies that are differenced

    // Finite difference force components
    for (int dim = 0; dim < 3; dim++)
    {
        Vec3D dr = (Vec3D){0.0, 0.0, 0.0};
        if (dim == 0)
            dr.x = delta;
        else if (dim == 1)
            dr.y = delta;
        else
            dr.z = delta;

        // E(+delta)
        *dr_i = r[i];
        r[i] = v3_add(r_saved, dr);
        *dr_i = v3_sub(r[i], *dr_i);
        update_nbrlist(p_parameters, p_vectors, p_nbrlist);
        double e_plus = kernel(p_parameters, p_nbrlist, p_vectors);

        // E(-delta)
        *dr_i = r[i];
        r[i] = v3_sub(r_saved, dr);
        *dr_i = v3_sub(r[i], *dr_i);
        update_nbrlist(p_parameters, p_vectors, p_nbrlist);
        double e_minus = kernel(p_parameters, p_nbrlist, p_vectors);

        if (fabs(e_plus) > e_scale) e_scale = fabs(e_plus);
        if (fabs(e_minus) > e_scale) e_scale = fabs(e_minus);

        double component = -(e_plus - e_minus) / (2.0 * delta);
        if (dim == 0)
            f_fd.x = component;
        else if (dim == 1)
            f_fd.y = component;
        else
            f_fd.z = component;
    }

    // Recompute analytical forces at the unperturbed position
    *dr_i = r[i];
    r[i] = r_saved;
    *dr_i = v3_sub(r[i], *dr_i);
    update_nbrlist(p_parameters, p_vectors, p_nbrlist);
    *dr_i = (Vec3D){0.0, 0.0, 0.0}; // reset displacement: otherwise later update_nbrlist calls re-apply it to the stored pair vectors
    kernel(p_parameters, p_nbrlist, p_vectors); // fills f_comp with forces at original position

    Vec3D f_an = f_comp[i];
    Vec3D df = v3_sub(f_an, f_fd);

    // The finite difference of two nearly equal energies is only meaningful
    // above its own rounding noise, which is of the order of the machine
    // precision times the magnitude of the energies, divided by delta. Forces
    // below that resolution carry no information: report them as untested
    // instead of as a relative error of order one. This happens for instance on
    // a perfect lattice, where all forces cancel by symmetry.
    double noise = 16.0 * DBL_EPSILON * e_scale / delta;
    double norm_an = v3_norm(f_an);
    double norm_fd = v3_norm(f_fd);
    double scale = (norm_an > norm_fd) ? norm_an : norm_fd;

    if (scale == 0.0 || scale < noise)
    {
        printf("Particle %d: force %g is below the finite-difference resolution %g, not tested\n",
               i, scale, noise);
        p_vectors->f = f_saved;
        free(f_comp);
        return 0.0;
    }

    double result = v3_norm(df) / scale;

    printf("Particle %d: Force analytical: (%g, %g, %g), Force FD: (%g, %g, %g), rel.err=%g\n",
           i, f_an.x, f_an.y, f_an.z, f_fd.x, f_fd.y, f_fd.z, result);

    // Restore original force pointer
    p_vectors->f = f_saved;
    free(f_comp);

    return result;
}
double virial_test(ForceKernel kernel,
                   struct Parameters *p_parameters,
                   struct Nbrlist *p_nbrlist,
                   struct Vectors *p_vectors)
{
    // Validate the accumulated virial against a numerical derivative of the
    // potential energy under isotropic scaling of all positions and the box:
    // r -> lambda r, L -> lambda L. For such scaling W = sum(f.r) = -dU/dlambda
    // at lambda = 1, so 3V (press_vir_nb + press_vir_bnd) must match -dU/dlambda.
    size_t num_part = p_parameters->num_part;
    Vec3D *r = p_vectors->r;
    Vec3D *r_saved = malloc(num_part * sizeof(Vec3D));
    if (!r_saved)
    {
        fprintf(stderr, "virial_test: allocation failed\n");
        return -1.0;
    }
    memcpy(r_saved, r, num_part * sizeof(Vec3D));
    Vec3D L_saved = p_parameters->L;

    const double h = 1e-6;
    double U[2];
    for (int sgn = 0; sgn < 2; ++sgn)
    {
        double lambda = (sgn == 0) ? 1.0 + h : 1.0 - h;
        for (size_t i = 0; i < num_part; ++i)
            r[i] = v3_scl(lambda, r_saved[i]);
        p_parameters->L = v3_scl(lambda, L_saved);
        build_nbrlist(p_parameters, p_vectors, p_nbrlist); // full rebuild: L changed
        U[sgn] = kernel(p_parameters, p_nbrlist, p_vectors);
    }

    // Restore the original state and recompute forces and virials there
    memcpy(r, r_saved, num_part * sizeof(Vec3D));
    p_parameters->L = L_saved;
    build_nbrlist(p_parameters, p_vectors, p_nbrlist);
    kernel(p_parameters, p_nbrlist, p_vectors);
    free(r_saved);

    double W_fd = -(U[0] - U[1]) / (2.0 * h);
    double V = L_saved.x * L_saved.y * L_saved.z;
    double W_an = 3.0 * V * (p_vectors->press_vir_nb + p_vectors->press_vir_bnd);
    // Note: angle and dihedral potentials depend only on directions, which are
    // invariant under isotropic scaling, so their exact virial is zero. To avoid
    // a meaningless relative error when W is (near) zero, the difference is
    // normalized by the energy scale in that case.
    double scale = fmax(fabs(W_fd), fabs(U[0]));
    double result = fabs(W_an - W_fd) / ((scale > 0.0) ? scale : 1.0);

    printf("Virial test: analytical W = %g, finite-difference W = %g, rel.err=%g\n",
           W_an, W_fd, result);
    return result;
}

/* The original finite-difference and virial tests above are retained.
 * These wrappers isolate each conservative potential and exclude D/R. */
double forces_bond_only(struct Parameters *p, struct Nbrlist *n, struct Vectors *v)
{
    int c = p->conservative_on, t = p->thermostat_on;
    p->conservative_on = 0; p->thermostat_on = 0;
    double e = calculate_forces(p, n, v);
    p->conservative_on = c; p->thermostat_on = t;
    return e;
}
double forces_conservative_only(struct Parameters *p, struct Nbrlist *n, struct Vectors *v)
{
    int t = p->thermostat_on;
    size_t bonds = v->num_bonds;
    p->thermostat_on = 0; v->num_bonds = 0;
    double e = calculate_forces(p, n, v);
    p->thermostat_on = t; v->num_bonds = bonds;
    return e;
}
