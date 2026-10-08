#ifndef NBRLIST_H_
#define NBRLIST_H_

#include <stddef.h>

/**
 * @brief Allocate arrays needed to store the cell-linked-list data
 *
 * @param p_parameters used members: L, r_cut, r_shell, num_part
 * @param p_celllist the cell-linked-list to allocate
 */
void alloc_celllist(struct Parameters *p_parameters, struct Celllist *p_celllist);

/**
 * @brief Free arrays used for the cell-linked-list
 *
 * @param p_celllist the cell-linked-list to free
 */
void free_celllist(struct Celllist *p_celllist);

/**
 * @brief Build the cell-linked-list
 * 
 * @param p_parameters 
 * @param p_vectors 
 * @param p_celllist 
 */
void build_celllist(struct Parameters *p_parameters, struct Vectors * p_vectors, struct Celllist *p_celllist);

/**
 * @brief Allocate arrays needed to store the neighbor list 
 * 
 * @param p_parameters 
 * @param p_nbrlist 
 */
void alloc_nbrlist(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist);

/**
 * @brief Free arrays use to store the neighbor list
 * 
 * @param p_nbrlist 
 */
void free_nbrlist(struct Nbrlist *p_nbrlist);


/**
 * @brief is_connected_12 returns 1 if i and j are 1-2 (directly bonded) connected
 *
 * Used while building the neighbor list to scale the non-bonded interaction of
 * bonded pairs by @ref Parameters.factor_12_nb (0 excludes the pair entirely).
 *
 * @param i particle index
 * @param j particle index
 * @param p_nbrlist pointer to neighbor list, used members: head12, pairs12
 * @return int 1 if particles i and j are 1-2 connected, 0 otherwise
 */
int is_connected_12(size_t i, size_t j, struct Nbrlist *p_nbrlist);

/**
 * @brief is_connected_13 returns 1 if i and j are 1-3 connected (two bonds
 * apart, the outer atoms of an angle)
 *
 * @param i particle index
 * @param j particle index
 * @param p_nbrlist pointer to neighbor list, used members: head13, pairs13
 * @return int 1 if particles i and j are 1-3 connected, 0 otherwise
 */
int is_connected_13(size_t i, size_t j, struct Nbrlist *p_nbrlist);

/**
 * @brief is_connected_14 returns 1 if i and j are 1-4 connected (three bonds
 * apart, the outer atoms of a dihedral)
 *
 * @param i particle index
 * @param j particle index
 * @param p_nbrlist pointer to neighbor list, used members: head14, pairs14
 * @return int 1 if particles i and j are 1-4 connected, 0 otherwise
 */
int is_connected_14(size_t i, size_t j, struct Nbrlist *p_nbrlist);

/**
 * @brief Build the neighbor list: all pairs closer than r_cut + r_shell.
 *
 * A cell-linked-list is built first, so only particles in neighboring cells
 * have to be considered as candidates; this makes the cost of the build grow
 * linearly with the number of particles instead of quadratically.
 *
 * @param p_parameters used members: r_cut, r_shell, factor_12_nb, factor_13_nb, factor_14_nb
 * @param p_vectors used members: r
 * @param p_nbrlist the neighbor list to (re)build
 */
void build_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

/**
 * @brief Rebuild the neighbor list when it can no longer be trusted.
 *
 * The list contains all pairs that were closer than r_cut + r_shell when it
 * was built. It remains valid as long as no particle has moved further than
 * r_shell/2 since then (two particles approaching each other can close the
 * gap from both sides); once one has, the list is rebuilt.
 *
 * @param p_parameters used members: r_shell
 * @param p_vectors used members: r
 * @param p_nbrlist used members: dr (per-particle displacement since the last build)
 * @return int 1 if the list was rebuilt, 0 if it was still valid.
 */
int update_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

#endif /* NBRLIST */
