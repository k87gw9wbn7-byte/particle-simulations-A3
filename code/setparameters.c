#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "setparameters.h"

/* Each experiment is a separate compilation; physical settings still live here.
 * B_CASE=0: monomer fluid (B3 and restart for B1)
 * B_CASE=1: ideal gas (B2+B3)
 * B_CASE=2: conservative-only energy test (B1)
 * With no B_CASE macro, retain the supplied 1000-step Part A smoke settings. */
#ifndef B_CASE
#define B_CASE (-1)
#endif
#ifndef B_DT
#define B_DT 0.04
#endif
#ifndef B_END
#define B_END 200.0
#endif
#ifndef B_STEPS
#define B_STEPS 20000
#endif
#ifndef B_EQUIL
#define B_EQUIL 5000
#endif
#ifndef B_PREFIX
#define B_PREFIX "../data/b_fluid"
#endif
#ifndef B_RESTART_IN
#define B_RESTART_IN "../data/b_fluid.restart"
#endif

/* C_CASE: -1 keeps B behaviour, 1=C1, 2=C2, 3=C3/C4, 4=C4 control. */
#ifndef C_CASE
#define C_CASE (-1)
#endif
#ifndef C_AAB
#define C_AAB 37.0
#endif
#ifndef C_DT
#define C_DT 0.04
#endif
#ifndef C_STEPS
#define C_STEPS 120000
#endif
#ifndef C_EQUIL
#define C_EQUIL 20000
#endif
#ifndef C_BLOCK
#define C_BLOCK 5000
#endif
#ifndef C_SAMPLE
#define C_SAMPLE 50
#endif
#ifndef C_PREFIX
#define C_PREFIX "../data/c1"
#endif
#ifndef C_SEED
#define C_SEED 13
#endif
#ifndef C_PDB
#define C_PDB 5000
#endif

/* Part D: 1=homogeneous chain fluid, 2=separated blend. */
#ifndef D_CASE
#define D_CASE (-1)
#endif
#ifndef D_N
#define D_N 5
#endif
#ifndef D_AAB
#define D_AAB 25.0
#endif
#ifndef D_STEPS
#define D_STEPS 60000
#endif
#ifndef D_EQUIL
#define D_EQUIL 20000
#endif
#ifndef D_BLOCK
#define D_BLOCK 5000
#endif
#ifndef D_SEED
#define D_SEED 23
#endif
#ifndef D_PREFIX
#define D_PREFIX "../data/d1"
#endif

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
    p->velocity_enabled = 0;
    p->equil_steps = 250;           /* discard steps <=250; not an equilibration claim */
    p->sample_steps = 25;
    p->block_steps = 250;
    p->rdf_max = 3.0;
    p->rdf_dr = 0.02;
    snprintf(p->analysis_prefix, sizeof(p->analysis_prefix), "../data/a_smoke");
    if (B_CASE >= 0) {
        if (B_CASE > 2) { fprintf(stderr,"B_CASE must be 0, 1 or 2\n"); exit(EXIT_FAILURE); }
        p->chain_length = 1; p->separated_start = 0;
        p->conservative_on = (B_CASE != 1);
        p->thermostat_on = (B_CASE != 2);
        p->dt = B_DT;
        p->num_dt_steps = B_STEPS;
        p->equil_steps = B_EQUIL;
        p->sample_steps = 50; p->block_steps = 1000;
        p->num_dt_output = 50;
        p->num_dt_pdb = 0; p->num_dt_restart = 0;
        p->analysis_enabled = (B_CASE != 2);
        p->velocity_enabled = (B_CASE != 2);
        p->load_restart = (B_CASE == 2);
        if (B_CASE == 2) {
            double steps = B_END / p->dt;
            if (B_END <= 0 || fabs(steps-round(steps)) > 1e-8) {
                fprintf(stderr,"B_END must be a positive whole number of time steps\n"); exit(EXIT_FAILURE);
            }
            p->num_dt_steps = (size_t)llround(steps);
            p->num_dt_output = 1; /* B1: retain energy at every completed step. */
        }
        snprintf(p->filename_csv,sizeof p->filename_csv,"%s_thermo",B_PREFIX);
        snprintf(p->filename_pdb,sizeof p->filename_pdb,"%s",B_PREFIX);
        snprintf(p->filename_xyz,sizeof p->filename_xyz,"%s",B_PREFIX);
        snprintf(p->analysis_prefix,sizeof p->analysis_prefix,"%s",B_PREFIX);
        snprintf(p->restart_out_filename,sizeof p->restart_out_filename,"%s.restart",B_PREFIX);
        snprintf(p->restart_in_filename,sizeof p->restart_in_filename,"%s",B_RESTART_IN);
    }
    if (C_CASE >= 1) {
        if (C_CASE > 4 || B_CASE >= 0) {
            fprintf(stderr,"Select C_CASE 1..4 without B_CASE overrides.\n"); exit(EXIT_FAILURE);
        }
        p->L = C_CASE == 1 ? v3(10,10,10) : v3(20,10,10);
        p->num_part = C_CASE == 1 ? 3000 : 6000;
        p->chain_length = 1; p->fraction_A = 0.5;
        p->separated_start = (C_CASE != 1);
        p->seed = C_SEED;
        p->conservative_on = 1; p->thermostat_on = 1;
        p->a[pair_index(0,0)] = p->a[pair_index(1,1)] = 25;
        p->a[pair_index(0,1)] = C_CASE == 1 || C_CASE == 4 ? 25 : C_AAB;
        p->dt = C_DT; p->num_dt_steps = C_STEPS; p->equil_steps = C_EQUIL;
        p->sample_steps = C_SAMPLE; p->block_steps = C_BLOCK;
        p->num_dt_output = 250;
        p->num_dt_pdb = C_CASE == 1 ? 0 : C_PDB;
        p->num_dt_restart = 0; p->load_restart = 0;
        p->analysis_enabled = (C_CASE == 1);
        p->velocity_enabled = 0;
        p->slab_enabled = (C_CASE != 1); p->profile_dr = 0.25;
        snprintf(p->filename_csv,sizeof p->filename_csv,"%s_thermo",C_PREFIX);
        snprintf(p->filename_pdb,sizeof p->filename_pdb,"%s",C_PREFIX);
        snprintf(p->filename_xyz,sizeof p->filename_xyz,"%s",C_PREFIX);
        snprintf(p->analysis_prefix,sizeof p->analysis_prefix,"%s",C_PREFIX);
        snprintf(p->restart_out_filename,sizeof p->restart_out_filename,"%s.restart",C_PREFIX);
    }
    if (D_CASE >= 1) {
        if (D_CASE > 2 || B_CASE >= 0 || C_CASE >= 1 || D_N < 2) {
            fprintf(stderr,"Select D_CASE 1..2 and D_N>=2 without B/C overrides.\n");
            exit(EXIT_FAILURE);
        }
        p->L = D_CASE == 1 ? v3(10,10,10) : v3(20,10,10);
        p->num_part = D_CASE == 1 ? 3000 : 6000;
        p->chain_length = D_N; p->fraction_A = 0.5;
        p->separated_start = D_CASE == 2; p->seed = D_SEED;
        p->conservative_on = p->thermostat_on = 1;
        p->a[pair_index(0,0)] = p->a[pair_index(1,1)] = 25;
        p->a[pair_index(0,1)] = D_CASE == 1 ? 25 : D_AAB;
        p->dt = 0.04; p->num_dt_steps = D_STEPS; p->equil_steps = D_EQUIL;
        p->sample_steps = 50; p->block_steps = D_BLOCK;
        p->num_dt_output = 1000; p->num_dt_pdb = D_STEPS;
        p->num_dt_restart = 0; p->load_restart = 0;
        p->analysis_enabled = p->velocity_enabled = 0;
        p->slab_enabled = D_CASE == 2; p->profile_dr = 0.25;
        snprintf(p->filename_csv,sizeof p->filename_csv,"%s_thermo",D_PREFIX);
        snprintf(p->filename_pdb,sizeof p->filename_pdb,"%s",D_PREFIX);
        snprintf(p->filename_xyz,sizeof p->filename_xyz,"%s",D_PREFIX);
        snprintf(p->analysis_prefix,sizeof p->analysis_prefix,"%s",D_PREFIX);
        snprintf(p->restart_out_filename,sizeof p->restart_out_filename,"%s.restart",D_PREFIX);
    }
    validate_parameters(p);
}

void validate_parameters(const struct Parameters *p)
{
    if (p->slab_enabled && (p->profile_dr <= 0 || !p->sample_steps || !p->block_steps ||
        p->block_steps % p->sample_steps || p->equil_steps % p->block_steps ||
        p->num_dt_steps % p->block_steps || p->equil_steps >= p->num_dt_steps ||
        fabs(p->L.x/p->profile_dr-round(p->L.x/p->profile_dr)) > 1e-8)) {
        fprintf(stderr,"Slab settings require whole equal blocks and profile bins.\n"); exit(EXIT_FAILURE);
    }
    double minL = fmin(p->L.x, fmin(p->L.y, p->L.z));
    if (p->num_part < 2 || !p->chain_length || p->num_part % p->chain_length ||
        p->fraction_A < 0 || p->fraction_A > 1 || p->dt <= 0 || p->kT <= 0 ||
        p->r_cut <= 0 || p->r_shell < 0 || minL < 3*(p->r_cut+p->r_shell) ||
        p->mass[0] <= 0 || p->mass[1] <= 0 || p->gamma < 0 ||
        p->initial_bond_length <= 0 || p->initial_bond_length >= 0.5*minL ||
        p->k_bond < 0 || p->r0_bond != 0 || p->factor_12_nb != 1 ||
        p->factor_13_nb != 1 || p->factor_14_nb != 1 || !p->num_dt_output ||
        ((p->analysis_enabled || p->velocity_enabled) && (!p->sample_steps || !p->block_steps ||
          p->block_steps % p->sample_steps || p->rdf_dr <= 0 ||
          p->rdf_max <= 0 || p->rdf_max > 0.5*minL ||
          p->equil_steps >= p->num_dt_steps ||
          p->num_dt_steps-p->equil_steps < p->sample_steps))) {
        fprintf(stderr, "Invalid DPD settings: check bead count, chains, units, intervals and cutoffs.\n");
        exit(EXIT_FAILURE);
    }
    if (fabs(p->noise_sigma*p->noise_sigma-2*p->gamma*p->kT) > 1e-10) {
        fprintf(stderr, "sigma^2 must equal 2*gamma*kT.\n"); exit(EXIT_FAILURE);
    }
}
