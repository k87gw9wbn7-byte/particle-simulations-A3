#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"
#include "initialise.h"

// This function initializes the particle types. The type of a particle is the
// index used to look up its type-dependent mass and pair repulsion (see pair_index in structs.h for the pair parameters).
void initialise_types(struct Parameters *p, struct Vectors *v)
{
    size_t chains = p->num_part / p->chain_length;
    size_t nA = (size_t)llround(p->fraction_A * (double)chains);
    for (size_t i = 0; i < p->num_part; ++i)
        v->type[i] = (i / p->chain_length < nA) ? 0 : 1;
}

void initialise_bond_connectivity(struct Parameters *p, struct Vectors *v)
{
    size_t N = p->chain_length, chains = p->num_part / N;
    v->num_bonds = chains * (N - 1);
    v->bonds = v->num_bonds ? malloc(v->num_bonds * sizeof(*v->bonds)) : NULL;
    size_t b = 0;
    for (size_t c = 0; c < chains; ++c)
        for (size_t k = 0; k + 1 < N; ++k)
            v->bonds[b++] = (struct Bond){c*N+k, c*N+k+1};
}

// This function derives the complete molecular structure from the list of
// bonds set up in initialise_bond_connectivity: every angle triplet i-j-k
// (two bonds sharing atom j) and every dihedral quadruplet i-j-k-l (three
// consecutive bonds), plus the 1-2, 1-3 and 1-4 partner lists that the
// neighbor list uses to exclude or scale non-bonded interactions.
//
// The partner lists all use the same compact layout: the partners of particle
// i are stored in pairs12[head12[i]] up to (excluding) pairs12[head12[i+1]].
// Such a list is built in three passes: count the partners of each particle,
// turn the counts into starting positions by a running sum, then fill the
// slots. No searching or sorting is needed.
void initialise_structure(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
{
    initialise_bond_connectivity(p_parameters, p_vectors); // Initialize bonds

    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    size_t num_part = p_parameters->num_part;

    // Build the 1-2 partner list from the bonds. Each bond i-j contributes two
    // entries: j is a partner of i and i is a partner of j.
    size_t *cnt = (size_t *)calloc(num_part + 1, sizeof(size_t));
    for (size_t i = 0; i < num_bonds; ++i)
    {
        ++cnt[bonds[i].i + 1];
        ++cnt[bonds[i].j + 1];
    }
    size_t *head12 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
    head12[0] = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        head12[i] = cnt[i] + head12[i - 1]; // running sum: where the partners of particle i start
        cnt[i] = head12[i];                 // reused as the next free slot while filling
    }
    size_t *pairs12 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
    for (size_t i = 0; i < num_bonds; ++i)
    {
        pairs12[cnt[bonds[i].i]++] = bonds[i].j;
        pairs12[cnt[bonds[i].j]++] = bonds[i].i;
    }

    // An angle i-j-k is a pair of bonds sharing the central atom j, so each
    // pair of 1-2 partners of an atom defines one angle centred on it.
    size_t num_angles = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        size_t num_pairs = head12[i] - head12[i - 1];
        num_angles += (num_pairs * (num_pairs - 1)) / 2;
    }
    struct Angle *angles = (struct Angle *)malloc(num_angles * sizeof(struct Angle));
    size_t m = 0;
    for (size_t i = 0; i < num_part; ++i)
        for (size_t j = head12[i]; j < head12[i + 1]; ++j)
            for (size_t k = j + 1; k < head12[i + 1]; ++k)
            {
                angles[m].i = pairs12[j];
                angles[m].j = i;
                angles[m].k = pairs12[k];
                ++m;
            }

    // A dihedral i-j-k-l is built around a central bond j-k: i is any other
    // partner of j, and l any other partner of k.
    size_t num_dihdr = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        num_dihdr += (head12[bonds[i].i + 1] - head12[bonds[i].i] - 1) * (head12[bonds[i].j + 1] - head12[bonds[i].j] - 1);
    }
    struct Dihedral *dihedrals = (struct Dihedral *) malloc(num_dihdr * sizeof(struct Dihedral));
    size_t k = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        struct Dihedral dihdr;
        dihdr.j = bonds[i].i;
        dihdr.k = bonds[i].j;
        for (size_t p = head12[dihdr.j]; p < head12[dihdr.j + 1]; ++p)
        {
            if (pairs12[p] == dihdr.k)
                continue;
            dihdr.i = pairs12[p];
            for (size_t q = head12[dihdr.k]; q < head12[dihdr.k + 1]; ++q)
            {
                if (pairs12[q] != dihdr.j)
                {
                    dihdr.l = pairs12[q];
                    dihedrals[k++] = dihdr;
                }
            }
        }
    }
    // 1-3 partners are the two outer atoms of each angle. The list is only
    // needed when the 1-3 non-bonded interaction is scaled (factor != 1).
    if (p_parameters->factor_13_nb != 1.0)
    {
        for (size_t i = 0; i <= num_part; ++i)
            cnt[i] = 0;
        for (size_t i = 0; i < num_angles; ++i)
        {
            ++cnt[angles[i].i + 1];
            ++cnt[angles[i].k + 1];
        }
        size_t *head13 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
        head13[0] = 0;
        for (size_t i = 1; i <= num_part; ++i)
        {
            head13[i] = cnt[i] + head13[i - 1];
            cnt[i] = head13[i];
        }
        size_t *pairs13 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
        for (size_t i = 0; i < num_angles; ++i)
        {
            pairs13[cnt[angles[i].i]++] = angles[i].k;
            pairs13[cnt[angles[i].k]++] = angles[i].i;
        }
        p_nbrlist->head13 = head13;
        p_nbrlist->pairs13 = pairs13;
    }
    // 1-4 partners are the two outer atoms of each dihedral, built the same way.
    if (p_parameters->factor_14_nb != 1.0)
    {
        for (size_t i = 0; i <= num_part; ++i)
            cnt[i] = 0;
        for (size_t i = 0; i < num_dihdr; ++i)
        {
            ++cnt[dihedrals[i].i + 1];
            ++cnt[dihedrals[i].l + 1];
        }
        size_t *head14 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
        head14[0] = 0;
        for (size_t i = 1; i <= num_part; ++i)
        {
            head14[i] = cnt[i] + head14[i - 1];
            cnt[i] = head14[i];
        }
        size_t *pairs14 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
        for (size_t i = 0; i < num_dihdr; ++i)
        {
            pairs14[cnt[dihedrals[i].i]++] = dihedrals[i].l;
            pairs14[cnt[dihedrals[i].l]++] = dihedrals[i].i;
        }
        p_nbrlist->head14 = head14;
        p_nbrlist->pairs14 = pairs14;
    }

    p_vectors->num_angles = num_angles;
    p_vectors->angles = angles;
    p_vectors->num_dihedrals = num_dihdr;
    p_vectors->dihedrals = dihedrals;
    // The 1-2 partner list is only kept if the neighbor list needs it to scale
    // the non-bonded interaction of bonded pairs.
    if (p_parameters->factor_12_nb != 1.0)
    {
        p_nbrlist->head12 = head12;
        p_nbrlist->pairs12 = pairs12;
    }
    else
    {
        free(head12);
        free(pairs12);
    }
    free(cnt);
}


// This function initializes the simulation by calling subroutines to initialize
// particle types, positions, velocities, and bond connectivity.
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, size_t *p_step, double *p_time)
{
    initialise_types(p_parameters, p_vectors);  // Initialize particle types
    initialise_structure(p_parameters, p_vectors, p_nbrlist);  // Initialize structure (bonds, angles, dihedrals)
    srand(p_parameters->seed);  // Seed random number generator
    initialise_positions(p_parameters, p_vectors);  // Initialize particle positions
    initialise_velocities(p_parameters, p_vectors);  // Initialize particle velocities
    *p_step = 0;   // Initialize step to zero
    *p_time = 0.0; // Initialize time to zero
    return;
}


/* Independent random-walk directions; reject steps crossing the assigned
 * half-box for a separated start. Periodic wrapping applies in all axes. */
void initialise_positions(struct Parameters *p, struct Vectors *v)
{
    size_t N = p->chain_length;
    for (size_t c = 0; c < p->num_part/N; ++c) {
        size_t first = c*N;
        int t = v->type[first];
        Vec3D r = v3(generate_uniform_random()*p->L.x,
                     generate_uniform_random()*p->L.y,
                     generate_uniform_random()*p->L.z);
        if (p->separated_start) r.x = 0.5*r.x + t*0.5*p->L.x;
        v->r[first] = v3_in_box(r, p->L);
        for (size_t k = 1; k < N; ++k) {
            do {
                double z = 2.0*generate_uniform_random()-1.0;
                double phi = 2.0*PI*generate_uniform_random();
                double s = sqrt(fmax(0.0, 1.0-z*z));
                Vec3D dr = v3_scl(p->initial_bond_length,
                                  v3(s*cos(phi), s*sin(phi), z));
                r = v3_in_box(v3_add(v->r[first+k-1], dr), p->L);
            } while (p->separated_start && ((r.x < 0.5*p->L.x) != (t == 0)));
            v->r[first+k] = r;
        }
    }
}

// This function initializes the velocities of particles based on the Maxwell-Boltzmann distribution.
// The total momentum is also removed to ensure zero total momentum (important for stability).
void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    /// \todo Use the type-dependent mass, and remove the total momentum rather
    /// than the average velocity, once the particles have different masses
    //double sqrtktm = sqrt(p_parameters->kT / p_parameters->mass);
    
    struct Vec3D p_tot = {0.0, 0.0, 0.0};  // Total momentum
    Vec3D *v = p_vectors->v;  // Pointer to particle velocities
    double m_tot = 0.0;     // total mass
    // Assign random velocities to each particle
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        double m = p_parameters->mass[p_vectors->type[i]];
        double sqrtktm = sqrt(p_parameters->kT / m);
        v[i] = v3_scl(sqrtktm, v3(gauss(), gauss(), gauss()));
        p_tot = v3_add(p_tot, v3_scl(m, v[i]));
        m_tot += m;
    }

    // Remove the centre-of-mass velocity, so that the total momentum is zero
    struct Vec3D v_cm = v3_scl(1.0 / m_tot, p_tot);
    for (size_t i = 0; i < p_parameters->num_part; i++)
        v[i] = v3_sub(v[i], v_cm);
}

