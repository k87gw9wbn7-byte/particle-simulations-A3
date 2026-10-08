#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "setparameters.h"

/* Edit here, then recompile. Default: a short Part A smoke test, NOT a
 * production result for B-D. All quantities are in DPD reduced units. */
void set_parameters(struct Parameters *p)
{
    memset(p, 0, sizeof(*p));
    p->L = v3(10, 10, 10);
    p->num_part = 3000;             /* rho = Nbeads / volume = 3 */
    p->chain_length = 1;            /* N=1 monomers; later 3, 5, 10 */
    p->fraction_A = 0.5;            /* rounded to nearest whole A chain */
    p->separated_start = 0;         /* 1: A left half, B right half */
    p->initial_bond_length = 0.8;    /* random-walk step, not rest length */
    p->seed = 13;

    p->conservative_on = 1;
    p->thermostat_on = 1;           /* dissipative AND random together */
    p->kT = 1.0;
    p->mass[0] = p->mass[1] = 1.0;
    p->r_cut = 1.0;
    p->r_shell = 0.4;
    p->gamma = 4.5;
    p->noise_sigma = sqrt(2.0*p->gamma*p->kT); /* sigma = 3 */
    p->a[pair_index(0,0)] = 25.0;
    p->a[pair_index(0,1)] = 25.0;
    p->a[pair_index(1,1)] = 25.0;
    p->k_bond = 2.0;
    p->r0_bond = 0.0;
    p->factor_12_nb = p->factor_13_nb = p->factor_14_nb = 1.0;

    p->dt = 0.04;
    p->num_dt_steps = 1000;         /* time 40; short implementation check */
    p->num_dt_output = 50;
    p->num_dt_pdb = 250;            /* 0 disables PDB */
    p->num_dt_restart = 1000;       /* 0 disables intermediate restarts */
    p->force_test = 0;              /* >0 tests every kth bead, then exits */
    p->load_restart = 0;
    p->rescale_output = 1.0;
    snprintf(p->filename_csv, sizeof(p->filename_csv), "../data/a_smoke_thermo");
    snprintf(p->filename_pdb, sizeof(p->filename_pdb), "../data/a_smoke");
    snprintf(p->filename_xyz, sizeof(p->filename_xyz), "../data/a_smoke");
    snprintf(p->restart_in_filename, sizeof(p->restart_in_filename), "../data/a_smoke.restart");
    snprintf(p->restart_out_filename, sizeof(p->restart_out_filename), "../data/a_smoke.restart");

    p->analysis_enabled = 1;
    p->equil_steps = 250;           /* discard steps <=250; not an equilibration claim */
    p->sample_steps = 25;
    p->block_steps = 250;
    p->rdf_max = 3.0;
    p->rdf_dr = 0.02;
    snprintf(p->analysis_prefix, sizeof(p->analysis_prefix), "../data/a_smoke");
    validate_parameters(p);
}

void validate_parameters(const struct Parameters *p)
{
    double minL = fmin(p->L.x, fmin(p->L.y, p->L.z));
    if (p->num_part < 2 || !p->chain_length || p->num_part % p->chain_length ||
        p->fraction_A < 0 || p->fraction_A > 1 || p->dt <= 0 || p->kT <= 0 ||
        p->r_cut <= 0 || p->r_shell < 0 || minL < 3*(p->r_cut+p->r_shell) ||
        p->mass[0] <= 0 || p->mass[1] <= 0 || p->gamma < 0 ||
        p->initial_bond_length <= 0 || p->initial_bond_length >= 0.5*minL ||
        p->k_bond < 0 || p->r0_bond != 0 || p->factor_12_nb != 1 ||
        p->factor_13_nb != 1 || p->factor_14_nb != 1 || !p->num_dt_output ||
        (p->analysis_enabled && (!p->sample_steps || !p->block_steps ||
          p->block_steps % p->sample_steps || p->rdf_dr <= 0 ||
          p->rdf_max <= 0 || p->rdf_max > 0.5*minL))) {
        fprintf(stderr, "Invalid DPD settings: check bead count, chains, units, intervals and cutoffs.\n");
        exit(EXIT_FAILURE);
    }
    if (fabs(p->noise_sigma*p->noise_sigma-2*p->gamma*p->kT) > 1e-10) {
        fprintf(stderr, "sigma^2 must equal 2*gamma*kT.\n"); exit(EXIT_FAILURE);
    }
}
