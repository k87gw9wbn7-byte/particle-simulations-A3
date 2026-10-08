/* Converted from the supplied PBS A2 MD program: the velocity-Verlet
 * sequence and its dynamics routines are retained (GW lambda=1/2). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "setparameters.h"
#include "initialise.h"
#include "nbrlist.h"
#include "forces.h"
#include "dynamics.h"
#include "memory.h"
#include "fileoutput.h"
#include "force_test.h"
#include "analysis.h"

static void thermo(struct Parameters *p, struct Vectors *v, size_t step,
                   double time, double U, double K)
{
    double T=2*K/(3.0*p->num_part-3), tt[TYPES];
    if (!isfinite(U) || !isfinite(K)) {
        fprintf(stderr,"Non-finite energy at step %zu\n",step); exit(EXIT_FAILURE);
    }
    temperature_per_type(p,v,tt);
    printf("Step %zu time %.6g U %.9g K %.9g E %.9g T %.6g\n",step,time,U,K,U+K,T);
    record_to_csv(0,p,step,time,U,K,T,tt);
}

int main(void)
{
    struct Parameters p;
    struct Vectors v;
    struct Nbrlist n;
    size_t step=0;
    double time=0, K=0, U;
    set_parameters(&p);
    if (p.force_test) p.thermostat_on=0; /* no potential for D/R */
    srand(p.seed);
    alloc_memory(&p,&v,&n);
    if (p.load_restart) {
        load_restart(&p,&v);
        boundary_conditions(&p,&v);
        initialise_types(&p,&v);
        initialise_structure(&p,&v,&n);
    } else {
        initialise(&p,&v,&n,&step,&time);
        boundary_conditions(&p,&v);
    }
    build_nbrlist(&p,&v,&n);
    U=calculate_forces(&p,&n,&v); /* initial random force, then once per full step */
    printf("DPD: beads=%zu chain_length=%zu L=(%g,%g,%g) rho=%g dt=%g seed=%u\n",
           p.num_part,p.chain_length,p.L.x,p.L.y,p.L.z,
           p.num_part/(p.L.x*p.L.y*p.L.z),p.dt,p.seed);
    printf("C=%d D/R=%d a=(%g,%g,%g) gamma=%g sigma=%g bonds=%zu separated=%d\n",
           p.conservative_on,p.thermostat_on,p.a[0],p.a[1],p.a[2],p.gamma,p.noise_sigma,
           v.num_bonds,p.separated_start);
    if (p.force_test) {
        double maxerr=0;
        for (size_t i=0;i<p.num_part;i+=(size_t)p.force_test) {
            double e=forces_test(forces_conservative_only,(int)i,&p,&n,&v);
            if (e>maxerr || e<0) maxerr=e<0 ? INFINITY : e;
            e=forces_test(forces_bond_only,(int)i,&p,&n,&v);
            if (e>maxerr || e<0) maxerr=e<0 ? INFINITY : e;
        }
        double ve=virial_test(calculate_forces,&p,&n,&v);
        printf("Maximum force relative error %.9g; virial %.9g\n",maxerr,ve);
        free_memory(&v,&n);
        return (maxerr<1e-4 && ve>=0 && ve<1e-4) ? 0 : 1;
    }
    if (p.num_dt_pdb) record_trajectories_pdb(1,&p,&v,time);
    double tt[TYPES]={0};
    record_to_csv(1,&p,0,0,0,0,0,tt);
    for (size_t i=0;i<p.num_part;++i) K+=0.5*p.mass[v.type[i]]*v3_dot(v.v[i],v.v[i]);
    thermo(&p,&v,step,time,U,K);
    analysis_start(&p,&v);
    while (step<p.num_dt_steps) {
        ++step; time=step*p.dt;
        update_velocities_half_dt(&p,&n,&v);
        update_positions(&p,&n,&v);
        boundary_conditions(&p,&v);
        update_nbrlist(&p,&v,&n);
        U=calculate_forces(&p,&n,&v);
        K=update_velocities_half_dt(&p,&n,&v);
        analysis_update(&p,&v,step);
        if (p.num_dt_pdb && step%p.num_dt_pdb==0) record_trajectories_pdb(0,&p,&v,time);
        if (p.num_dt_restart && step%p.num_dt_restart==0) save_restart(&p,&v);
        if (step%p.num_dt_output==0) thermo(&p,&v,step,time,U,K);
    }
    analysis_finish(&p,&v);
    save_restart(&p,&v);
    free_memory(&v,&n);
    return 0;
}
