#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "random.h"

/* One draw per unordered interacting pair per force evaluation. */
double dpd_noise(void)
{
    return sqrt(3.0)*(2.0*generate_uniform_random()-1.0);
}

double calculate_forces(struct Parameters *p, struct Nbrlist *n, struct Vectors *v)
{
    for (size_t i=0; i<p->num_part; ++i) v->f[i] = v3(0,0,0);
    v->press_vir_nb = v->press_vir_bnd = 0;
    double U = calculate_forces_bond(p,v) + calculate_forces_nb(p,n,v);
    double threeV = 3.0*p->L.x*p->L.y*p->L.z;
    v->press_vir_nb /= threeV;
    v->press_vir_bnd /= threeV;
    return U;
}

double calculate_forces_bond(struct Parameters *p, struct Vectors *v)
{
    double U = 0;
    for (size_t b=0; b<v->num_bonds; ++b) {
        size_t i=v->bonds[b].i, j=v->bonds[b].j;
        Vec3D rij = v3_min_image(v3_sub(v->r[i],v->r[j]),p->L);
        /* A2 harmonic bond, with r0=0: no division by r needed, including r=0. */
        Vec3D f = v3_scl(-p->k_bond,rij);
        v->f[i] = v3_add(v->f[i],f);
        v->f[j] = v3_sub(v->f[j],f);
        U += 0.5*p->k_bond*v3_dot(rij,rij);
        v->press_vir_bnd += v3_dot(rij,f);
    }
    return U;
}

double calculate_forces_nb(struct Parameters *p, struct Nbrlist *n, struct Vectors *v)
{
    double U = 0, rc2=p->r_cut*p->r_cut;
    if (!p->conservative_on && !p->thermostat_on) return 0;
    double random_scale = p->noise_sigma/sqrt(p->dt);
    for (size_t q=0; q<n->num_nbrs; ++q) {
        size_t i=n->nbr[q].i, j=n->nbr[q].j;
        Vec3D rij = v3_min_image(v3_sub(v->r[i],v->r[j]),p->L);
        double r2=v3_dot(rij,rij);
        if (r2 >= rc2) continue;
        if (r2 == 0.0) {
            fprintf(stderr,"DPD pair %zu,%zu exactly coincident: direction undefined.\n",i,j);
            exit(EXIT_FAILURE);
        }
        double r=sqrt(r2), w=1.0-r/p->r_cut;
        Vec3D e=v3_scl(1.0/r,rij);
        double fc=0, magnitude=0;
        if (p->conservative_on) {
            double a=p->a[pair_index(v->type[i],v->type[j])];
            fc = a*w;
            U += 0.5*a*p->r_cut*w*w;
            v->press_vir_nb += r*fc; /* C only, never D or R */
            magnitude += fc;
        }
        if (p->thermostat_on) {
            Vec3D vij=v3_sub(v->v[i],v->v[j]); /* predicted half-step velocities */
            magnitude -= p->gamma*w*w*v3_dot(e,vij);
            magnitude += random_scale*w*dpd_noise();
        }
        Vec3D f=v3_scl(magnitude,e);
        /* The neighbour list visits each unordered pair once; sharing this
         * force is exactly zeta_ij = zeta_ji and gives Newton's third law. */
        v->f[i]=v3_add(v->f[i],f);
        v->f[j]=v3_sub(v->f[j],f);
    }
    return U;
}
