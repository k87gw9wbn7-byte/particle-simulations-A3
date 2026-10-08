#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

// Allocate memory for the cell-linked list based on the number of particles and grid size.
// The grid size is determined by dividing the simulation box into cells, each larger than the interaction range.
void alloc_celllist(struct Parameters *p_parameters, struct Celllist *p_celllist)
{
    struct Index3D size_grid;
    size_t Mtot;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell;
    const size_t num_part = p_parameters->num_part;
    
    // Determine the size of the grid based on cutoff distance and box dimensions
    size_grid.i = floor(p_parameters->L.x / rlist);
    size_grid.j = floor(p_parameters->L.y / rlist);
    size_grid.k = floor(p_parameters->L.z / rlist);
    p_celllist->size_grid.i = size_grid.i;
    p_celllist->size_grid.j = size_grid.j;
    p_celllist->size_grid.k = size_grid.k;

    // Calculate total number of cells and allocate memory for the cell lists
    Mtot = size_grid.i * size_grid.j * size_grid.k;
    p_celllist->head = (size_t *)malloc(Mtot * sizeof(size_t));
    p_celllist->num_cells_max = Mtot;
    p_celllist->num_cells = Mtot;
    p_celllist->num_part_max = num_part;
    p_celllist->particle2cell = (size_t *)malloc(num_part * sizeof(size_t));
    p_celllist->list = (size_t *)malloc(num_part * sizeof(size_t));
}

// Free the memory allocated for the cell-linked list structures.
void free_celllist(struct Celllist *p_celllist)
{
    free(p_celllist->head);
    p_celllist->head = NULL;
    free(p_celllist->particle2cell);
    p_celllist->particle2cell = NULL;
    free(p_celllist->list);
    p_celllist->list = NULL;
}

// Build the cell-linked list by assigning particles to grid cells based on their positions.
// This is the first step before building the neighbor list.
void build_celllist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Celllist *p_celllist)
{
    struct Index3D size_grid, indx;
    size_t Mtot, icell;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell;
    struct Vec3D mL;
    const size_t num_part = p_parameters->num_part;
    size_t *particle2cell, *head, *celllist;
    struct Vec3D *r;

    size_grid.i = floor(p_parameters->L.x / rlist);
    size_grid.j = floor(p_parameters->L.y / rlist);
    size_grid.k = floor(p_parameters->L.z / rlist);
    p_celllist->size_grid.i = size_grid.i;
    p_celllist->size_grid.j = size_grid.j;
    p_celllist->size_grid.k = size_grid.k;
    Mtot = size_grid.i * size_grid.j * size_grid.k;
    if (Mtot > p_celllist->num_cells_max)
    {
        p_celllist->head = (size_t *)realloc(p_celllist->head, Mtot * sizeof(size_t));
        p_celllist->num_cells_max = Mtot;
    }
    if (num_part > p_celllist->num_part_max)
    {
        p_celllist->particle2cell = (size_t *)realloc(p_celllist->particle2cell, num_part * sizeof(size_t));
        p_celllist->list = (size_t *)realloc(p_celllist->list, num_part * sizeof(size_t));
        p_celllist->num_part_max = num_part;
    }
    mL.x = ((double)size_grid.i) / p_parameters->L.x;
    mL.y = ((double)size_grid.j) / p_parameters->L.y;
    mL.z = ((double)size_grid.k) / p_parameters->L.z;
    head = p_celllist->head;
    for (icell = 0; icell < Mtot; ++icell)
        head[icell] = SIZE_MAX;
    particle2cell = p_celllist->particle2cell;
    celllist = p_celllist->list;
    r = p_vectors->r;
    for (size_t i = 0; i < num_part; ++i)
    {
        // Wrap the position into the box before computing the cell index: an
        // out-of-box position would give an out-of-range index and corrupt memory.
        struct Vec3D rw = v3_in_box(r[i], p_parameters->L);
        indx.i = floor(rw.x * mL.x);
        indx.j = floor(rw.y * mL.y);
        indx.k = floor(rw.z * mL.z);
        // Guard against the floating-point edge case r == L after wrapping
        if (indx.i >= size_grid.i) indx.i = size_grid.i - 1;
        if (indx.j >= size_grid.j) indx.j = size_grid.j - 1;
        if (indx.k >= size_grid.k) indx.k = size_grid.k - 1;
        icell = indx.i + size_grid.i * (indx.j + indx.k * size_grid.j);
        particle2cell[i] = icell;
        celllist[i] = head[icell];
        head[icell] = i;
    }
}

// Allocate memory for the neighbor list. This includes arrays for storing neighbors and bonded interactions.
void alloc_nbrlist(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist)
{
    double Ndouble = (double)p_parameters->num_part;
    double rlist = p_parameters->r_cut + p_parameters->r_shell;
    // Estimate the number of pairs from the average density: each particle sees
    // on average rho * (4/3) pi rlist^3 others, and each pair is stored once
    // (hence the factor 0.5). The list grows automatically if this is too small.
    double Nnbr = 0.5 * Ndouble * (4.0 / 3.0 * PI * rlist * rlist * rlist) / (p_parameters->L.x * p_parameters->L.y * p_parameters->L.z);
    size_t num_nbrs_max = ceil(0.6 * Ndouble * (Nnbr + 2.0 * sqrt(Nnbr)));
    p_nbrlist->p_celllist = (struct Celllist *)malloc(sizeof(struct Celllist));
    alloc_celllist(p_parameters, p_nbrlist->p_celllist);
    p_nbrlist->num_nbrs_max = num_nbrs_max;
    p_nbrlist->num_nbrs = 0; // no pairs until build_nbrlist fills the list
    p_nbrlist->nbr = (struct Pair *)malloc(num_nbrs_max * sizeof(struct Pair));
    p_nbrlist->dr = (struct DeltaR *)malloc(p_parameters->num_part * sizeof(struct DeltaR));
    p_nbrlist->head12 = NULL;
    p_nbrlist->pairs12 = NULL;
    p_nbrlist->head13 = NULL;
    p_nbrlist->pairs13 = NULL;
    p_nbrlist->head14 = NULL;
    p_nbrlist->pairs14 = NULL;
}

// Free the memory allocated for the neighbor list and associated structures.
void free_nbrlist(struct Nbrlist *p_nbrlist)
{
    free_celllist(p_nbrlist->p_celllist);
    free(p_nbrlist->p_celllist);
    p_nbrlist->p_celllist = NULL;
    free(p_nbrlist->nbr);
    p_nbrlist->nbr = NULL;
    free(p_nbrlist->dr);
    p_nbrlist->dr = NULL;
    free(p_nbrlist->head12);
    p_nbrlist->head12 = NULL;
    free(p_nbrlist->pairs12);
    p_nbrlist->pairs12 = NULL;
    free(p_nbrlist->head13);
    p_nbrlist->head13 = NULL;
    free(p_nbrlist->pairs13);
    p_nbrlist->pairs13 = NULL;
    free(p_nbrlist->head14);
    p_nbrlist->head14 = NULL;
    free(p_nbrlist->pairs14);
    p_nbrlist->pairs14 = NULL;
}

// Check if particles i and j are 1-2 connected (i.e., bonded).
int is_connected_12(size_t i, size_t j, struct Nbrlist *p_nbrlist)
{
    size_t bgn = p_nbrlist->head12[i];
    size_t nd = p_nbrlist->head12[i+1];
    size_t * pairs12 = p_nbrlist->pairs12;
    for(size_t k=bgn; k< nd; ++k)
        if (pairs12[k]==j) return 1;
    return 0;
}

// Check if particles i and j are 1-3 connected (e.g., part of an angle interaction).
int is_connected_13(size_t i, size_t j, struct Nbrlist *p_nbrlist)
{
    size_t bgn = p_nbrlist->head13[i];
    size_t nd = p_nbrlist->head13[i+1];
    size_t * pairs13 = p_nbrlist->pairs13;
    for(size_t k=bgn; k< nd; ++k)
        if (pairs13[k]==j) return 1;
    return 0;
}

// Check if particles i and j are 1-4 connected (e.g., part of a dihedral interaction).
int is_connected_14(size_t i, size_t j, struct Nbrlist *p_nbrlist)
{
    size_t bgn = p_nbrlist->head14[i];
    size_t nd = p_nbrlist->head14[i+1];
    size_t * pairs14 = p_nbrlist->pairs14;
    for(size_t k=bgn; k< nd; ++k)
        if (pairs14[k]==j) return 1;
    return 0;
}

/**
 * @brief How strongly the non-bonded interaction of a pair is scaled, given how
 * many bonds apart the two atoms are.
 *
 * The closest relation decides: a pair that is both 1-3 and 1-4 connected, as
 * happens in a ring, counts as 1-3. Each test is only made when its factor
 * differs from one, so a force field that does not care about connectivity, as
 * in DPD, pays nothing for these lines.
 */
static inline double nb_factor(size_t i, size_t j, struct Nbrlist *p_nbrlist,
                               int check_12, double f12, int check_13, double f13,
                               int check_14, double f14)
{
    if (check_12 && is_connected_12(i, j, p_nbrlist)) return f12;
    if (check_13 && is_connected_13(i, j, p_nbrlist)) return f13;
    if (check_14 && is_connected_14(i, j, p_nbrlist)) return f14;
    return 1.0;
}

// Build the neighbor list: all pairs of particles closer than r_cut + r_shell,
// each stored once. The cell-linked-list restricts the candidate pairs to
// particles in the same or adjacent cells; without it every pair in the system
// would have to be tested, at a cost growing quadratically with num_part.
// Pairs whose non-bonded interaction is excluded by their connectivity
// (factor 0) are left out of the list; a scaled pair carries its factor.
void build_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
{
    struct Index3D indx, indx_nbr;
    size_t icell, inbr;
    size_t num_nbrs, num_nbrs_max;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell; /* the radius for inclusion in the list is r_cut + r_shell */
    const double rlist_sq = rlist * rlist;
    struct Vec3D ri;
    struct Vec3D *r = p_vectors->r;
    struct DeltaR rij;
    const struct DeltaR dr = {{0.0, 0.0, 0.0}, 0.0};
    size_t *head, *particle2cell, *celllist;
    // Half of the 26 neighboring cells. Visiting only 13 stores each pair once:
    // the pair (i,j) with j in the opposite cell is found when the loop reaches
    // particle j, whose half-stencil contains the cell of i.
    const int nbr_indcs[13][3] = {{0, 0, 1}, {0, 1, -1}, {0, 1, 0}, {0, 1, 1}, {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0}, {1, 0, 1}, {1, 1, -1}, {1, 1, 0}, {1, 1, 1}};
    size_t num_part = p_parameters->num_part;
    struct Pair *nbr;

    // First build a cell-linked-list
    build_celllist(p_parameters, p_vectors, p_nbrlist->p_celllist);

    /* Use the cell-linked-list to build the neighbor list: a pair is included
       when its distance is less than r_cut + r_shell. */
    struct Index3D size_grid = p_nbrlist->p_celllist->size_grid;
    num_nbrs_max = p_nbrlist->num_nbrs_max;
    num_nbrs = 0;
    nbr = p_nbrlist->nbr;
    particle2cell = p_nbrlist->p_celllist->particle2cell;
    head = p_nbrlist->p_celllist->head;
    celllist = p_nbrlist->p_celllist->list;
    const double f12 = p_parameters->factor_12_nb;
    const double f13 = p_parameters->factor_13_nb;
    const double f14 = p_parameters->factor_14_nb;
    // A pair whose factor is one is treated like any other pair, so there is
    // nothing to look up: the test is skipped entirely.
    const int check_12 = (f12 != 1.0);
    const int check_13 = (f13 != 1.0);
    const int check_14 = (f14 != 1.0);
    for (size_t i = 0; i < num_part; ++i)
    {
        // Find neighbors of particle i in its own cell. Starting the walk at
        // celllist[i], the particles inserted before i, stores each pair once.
        ri = r[i];
        for (size_t j = celllist[i]; j != SIZE_MAX; j = celllist[j])
        {
            rij.v = v3_sub(ri, r[j]);
            rij.sq = v3_dot(rij.v, rij.v);
            if (rij.sq < rlist_sq)
            {
                double factor = nb_factor(i, j, p_nbrlist, check_12, f12, check_13, f13, check_14, f14);
                if (factor == 0.0)
                    continue; // fully excluded, so it does not belong in the list
                if (num_nbrs >= num_nbrs_max)
                {
                    num_nbrs_max += 5 * p_parameters->num_part;
                    nbr = (struct Pair *)realloc(nbr, num_nbrs_max * sizeof(struct Pair));
                    p_nbrlist->nbr = nbr;
                    p_nbrlist->num_nbrs_max = num_nbrs_max;
                } 
                nbr[num_nbrs].i = i;
                nbr[num_nbrs].j = j;
                nbr[num_nbrs].factor = factor;
                num_nbrs++;
            }
        }

        // Next find neighbors of particle i in the 13 cells of its half-stencil
        icell = particle2cell[i];
        indx.i = icell % size_grid.i;
        icell = icell / size_grid.i;
        indx.j = icell % size_grid.j;
        indx.k = icell / size_grid.j;
        for (size_t k = 0; k < 13; ++k)
        {
            ri = p_vectors->r[i];
            indx_nbr.i = (indx.i + nbr_indcs[k][0]);
            indx_nbr.j = (indx.j + nbr_indcs[k][1]);
            indx_nbr.k = (indx.k + nbr_indcs[k][2]);
            // The if-statements below implement the periodic boundary
            // conditions: a neighbor cell beyond the edge of the grid wraps to
            // the other side, and the position of particle i is shifted by a
            // box length so plain subtraction gives the correct distance to
            // the particles in that cell. (Note that an unsigned index that
            // would go negative wraps to SIZE_MAX.)
            if (indx_nbr.i == SIZE_MAX)
            {
                ri.x += p_parameters->L.x;
                indx_nbr.i = size_grid.i - 1;
            }
            else if (indx_nbr.i >= size_grid.i)
            {
                ri.x -= p_parameters->L.x;
                indx_nbr.i = 0;
            }
            if (indx_nbr.j == SIZE_MAX)
            {
                ri.y += p_parameters->L.y;
                indx_nbr.j = size_grid.j - 1;
            }
            else if (indx_nbr.j >= size_grid.j)
            {
                ri.y -= p_parameters->L.y;
                indx_nbr.j = 0;
            }
            if (indx_nbr.k == SIZE_MAX)
            {
                ri.z += p_parameters->L.z;
                indx_nbr.k = size_grid.k - 1;
            }
            else if (indx_nbr.k >= size_grid.k)
            {
                ri.z -= p_parameters->L.z;
                indx_nbr.k = 0;
            }
            inbr = indx_nbr.i + size_grid.i * (indx_nbr.j + indx_nbr.k * size_grid.j);
            for (size_t j = head[inbr]; j != SIZE_MAX; j = celllist[j])
            {
                // Test the distance before asking whether the pair is excluded:
                // most candidates in a neighbouring cell are out of range, and
                // the exclusion lists are searched, so asking the other way
                // round searches them for pairs that are thrown away anyway.
                rij.v = v3_sub(ri, p_vectors->r[j]);
                rij.sq = v3_dot(rij.v, rij.v);
                if (rij.sq < rlist_sq)
                {
                    double factor = nb_factor(i, j, p_nbrlist, check_12, f12, check_13, f13, check_14, f14);
                    if (factor == 0.0)
                        continue; // fully excluded, so it does not belong in the list
                    if (num_nbrs >= num_nbrs_max)
                    {
                        num_nbrs_max += 5 * p_parameters->num_part;
                        nbr = (struct Pair *)realloc(nbr, num_nbrs_max * sizeof(struct Pair));
                        p_nbrlist->nbr = nbr;
                        p_nbrlist->num_nbrs_max = num_nbrs_max;
                    }
                    nbr[num_nbrs].i = i;
                    nbr[num_nbrs].j = j;
                    nbr[num_nbrs].factor = factor;
                    num_nbrs++;
                }
            }
        }
    }
    p_nbrlist->num_nbrs = num_nbrs;
    p_nbrlist->dr = (struct DeltaR *)realloc(p_nbrlist->dr, num_part * sizeof(struct DeltaR));
    for (size_t i = 0; i < num_part; ++i) /*initialize particle displacements (with respect to creation time) to zero */
        p_nbrlist->dr[i] = dr;
}

// This function decides, every time step, whether the neighbor list is still
// valid. The list contains all pairs that were within r_cut + r_shell at build
// time, so it misses a pair within r_cut only after two particles have
// together bridged the margin r_shell — which cannot happen before one of them
// has moved r_shell/2. Once a particle has, the list is rebuilt.
int update_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
{
    const double dr_sq_max = 0.25 * (p_parameters->r_shell * p_parameters->r_shell);
    int isRebuild = 0;
    // Check whether any particle has moved more than r_shell/2 since the build
    for (size_t i = 0; i < p_parameters->num_part; ++i)
        if ((p_nbrlist->dr[i].sq) > dr_sq_max)
        {
            isRebuild = 1;
            break; // one particle is enough, the whole list is rebuilt anyway
        }
    if (isRebuild) // rebuild neighbor list
        build_nbrlist(p_parameters, p_vectors, p_nbrlist);
    // Nothing to do when no rebuild is needed: the pair list holds only the
    // indices of the pairs, and the force routine takes the connecting vector
    // from the current positions.
    return isRebuild;
}
