#ifndef TYPES_MD_H_
#define TYPES_MD_H_
#include <math.h>
#include "constants.h"

/* This header file contains definitions of struct types used in the molecular
   dynamics code, and the small vector helpers (v3_*) used throughout. */

/**
 * @brief A 3D vector.
 *
 * All positions, velocities and forces are stored as arrays of Vec3D. The
 * v3_* helper functions below implement the usual vector operations; they are
 * 'static inline' so the compiler expands them in place and they cost the same
 * as writing the component arithmetic out by hand.
 */
typedef struct Vec3D
{
    double x, y, z; //!< The three Cartesian components of the vector
} Vec3D;

/** @brief Construct the vector (x, y, z). */
static inline Vec3D v3(double x, double y, double z)
{
    return (Vec3D){x, y, z};
}

/** @brief Vector addition: a + b. */
static inline Vec3D v3_add(Vec3D a, Vec3D b)
{
    return (Vec3D){a.x + b.x, a.y + b.y, a.z + b.z};
}

/** @brief Vector subtraction: a - b. */
static inline Vec3D v3_sub(Vec3D a, Vec3D b)
{
    return v3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/** @brief Component-wise product: (a.x b.x, a.y b.y, a.z b.z). */
static inline Vec3D v3_mul(Vec3D a, Vec3D b)
{
    return v3(a.x * b.x, a.y * b.y, a.z * b.z);
}

/** @brief Component-wise division: (a.x / b.x, a.y / b.y, a.z / b.z). */
static inline Vec3D v3_div(Vec3D a, Vec3D b)
{
    return v3(a.x / b.x, a.y / b.y, a.z / b.z);
}

/** @brief Multiplication by a scalar: s a. */
static inline Vec3D v3_scl(double s, Vec3D a)
{
    return v3(s * a.x, s * a.y, s * a.z);
}

/** @brief Dot (inner) product a . b. */
static inline double v3_dot(Vec3D a, Vec3D b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/** @brief The length (Euclidean norm) |a| of a vector. */
static inline double v3_norm(Vec3D a)
{
    return sqrt(v3_dot(a, a));
}

/** @brief Cross product a x b. */
static inline Vec3D v3_cross(Vec3D a, Vec3D b)
{
    return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

/**
 * @brief The minimum image of a connecting vector: the shortest vector between
 * two particles, taking the periodic images of the box into account.
 *
 * All particles are inside the box, so a connecting vector is never longer than
 * one box length and shifting it by at most one box does the job. That makes a
 * comparison enough, where the general expression r - L*floor(r/L + 0.5) would
 * need a division. The division is what costs: the compiler turns floor into a
 * few instructions, but a double division is an order of magnitude slower still.
 *
 * This relies on the positions being inside the box, which @ref initialise and
 * @ref boundary_conditions guarantee. Use @ref v3_in_box, which makes no such
 * assumption, for a vector that may be many box lengths long.
 */
static inline Vec3D v3_min_image(Vec3D r, Vec3D L)
{
    if (r.x >= 0.5 * L.x) r.x -= L.x; else if (r.x < -0.5 * L.x) r.x += L.x;
    if (r.y >= 0.5 * L.y) r.y -= L.y; else if (r.y < -0.5 * L.y) r.y += L.y;
    if (r.z >= 0.5 * L.z) r.z -= L.z; else if (r.z < -0.5 * L.z) r.z += L.z;
    return r;
}

/**
 * @brief Wrap a position into the box [0,L) in each direction, however far
 * outside it is. More general, but also more expensive, than @ref v3_min_image.
 */
static inline Vec3D v3_in_box(Vec3D r, Vec3D L)
{
    r.x -= L.x * floor(r.x / L.x);
    r.y -= L.y * floor(r.y / L.y);
    r.z -= L.z * floor(r.z / L.z);
    return r;
}

/**
 * @brief Map an (unordered) pair of particle types to an index in an array of
 * per-pair parameters, such as @ref Parameters.a.
 *
 * The interaction parameters are symmetric, a(i,j) = a(j,i), so
 * only the lower triangle of the type-by-type matrix is stored, packed row by
 * row: NUM_TYPES*(NUM_TYPES+1)/2 entries instead of NUM_TYPES^2.
 */
static inline int pair_index(int i, int j)
{
    return (i>j? (i*(i+1))/2 + j: (j*(j+1))/2 + i);
}

/**
 * @brief All run settings in one struct, set by @ref set_parameters. There are
 * no input files: changing a run means editing setparameters.c and recompiling.
 */
struct Parameters
{
    size_t num_part, chain_length;
    double fraction_A, initial_bond_length;
    int separated_start, conservative_on, thermostat_on;
    unsigned seed;
    double dt, kT, r_cut, r_shell, gamma, noise_sigma;
    Vec3D L;
    double mass[TYPES], a[TYPES*(TYPES+1)/2]; /* AA, AB, BB */
    double k_bond, r0_bond;
    double factor_12_nb, factor_13_nb, factor_14_nb;
    size_t num_dt_steps, num_dt_output, num_dt_pdb, num_dt_restart;
    int force_test, load_restart;
    double rescale_output;
    char filename_pdb[1024], filename_xyz[1024], filename_csv[1024];
    char restart_in_filename[1024], restart_out_filename[1024];
    int analysis_enabled;
    size_t equil_steps, sample_steps, block_steps;
    double rdf_max, rdf_dr;
    char analysis_prefix[1024];
};

/**
 * @brief A 3D vector stored together with its square length. Used for
 * connecting vectors and accumulated displacements, where the square length is
 * needed anyway for comparisons against a cut-off distance.
 */
struct DeltaR
{
    Vec3D v;   //!< x, y and z components
    double sq; //!< square length v.v
};

/**
 * @brief Struct to store i, j, k indices of a 3D grid
 * 
 */
struct Index3D
{
    size_t i, j, k; //!< 3 indices: i,j, k
};

/**
 * @brief Struct to store indices of bonded particles i-j
 * 
 */
struct Bond
{
    size_t i,j;
};

/**
 * @brief Struct to store indices of particles in an angle i-j-k
 * 
 */
struct Angle
{
    size_t i,j,k;
};

/**
 * @brief Struct to store indices of particles in a dihedral i-j-k-l
 * 
 */
struct Dihedral
{
    size_t i,j,k,l;
};

/**
 * @brief One entry of the neighbor list: a pair of particles within interaction
 * range. Only the indices are stored; the connecting vector is computed from
 * the current positions in the force loop.
 */
struct Pair
{
    size_t i, j;   //!< indices of the two particles forming the pair
    double factor; //!< scaling of the non-bonded interaction of this pair (1 for an ordinary pair, e.g. @ref Parameters.factor_14_nb for a 1-4 connected one)
};

/**
 * @brief All per-particle data of the simulation, stored as parallel arrays:
 * element i of each array belongs to particle i. Also carries the molecular
 * topology (bonds, angles, dihedrals) and the accumulated pressure
 * contributions of the current time step.
 */
struct Analysis;
struct Vectors
{
    struct Analysis *analysis;
    size_t size;                //!< allocated length of the particle arrays (can be > num_part)
    size_t num_bonds;           //!< number of bonds
    size_t num_angles;          //!< number of angles 
    size_t num_dihedrals;       //!< number of dihedrals
    int    *type;               //!< particle types (index into the per-type parameter arrays)
    struct Vec3D *r;            //!< positions
    struct Vec3D *dr;           //!< displacements of the last time step
    struct Vec3D *v;            //!< velocities
    struct Vec3D *f;            //!< forces
    struct Bond *bonds;         //!< bonds
    struct Angle *angles;       //!< angles
    struct Dihedral *dihedrals; //!< dihedrals
    double press_kin;           //!< kinetic (ideal-gas) contribution to the pressure, set in @ref update_velocities_half_dt
    double press_vir_nb;        //!< non-bonded virial contribution to the pressure, set in @ref calculate_forces
    double press_vir_bnd;       //!< bonded (bond + angle + dihedral) virial contribution to the pressure, set in @ref calculate_forces
};

/**
 * @brief A cell-linked-list: the box is divided into cells at least as large as
 * the interaction range, so all neighbors of a particle are found in its own
 * cell and the 26 surrounding ones.
 *
 * The particles of one cell form a linked list: head[icell] is the first
 * particle of cell icell, list[i] the particle after particle i, and SIZE_MAX
 * marks the end. This stores an arbitrary number of particles per cell in two
 * fixed-size arrays, with no per-cell allocation.
 */
struct Celllist
{
    size_t *head; //!< head[icell] is the index of the first particle in cell icell, or SIZE_MAX for an empty cell
    size_t *list; //!< list[i] is the next particle in the same cell as particle i; SIZE_MAX marks the end of the list
    size_t *particle2cell; //!< particle2cell[i] is the cell index of particle i
    size_t num_cells, num_cells_max, num_part_max; //!< number of cells used, and the numbers of cells and particles allocated for
    struct Index3D size_grid; //!< number of cells in each direction
};

/**
 * @brief The (Verlet) neighbor list: all particle pairs closer than
 * r_cut + r_shell when the list was last built.
 *
 * The extra margin r_shell allows the same list to be reused for many steps.
 * The accumulated displacement of each particle since the last build is
 * tracked in dr; once any particle has moved more than r_shell/2 the list is
 * no longer guaranteed to contain all pairs within r_cut and is rebuilt
 * (see @ref update_nbrlist).
 *
 * The head12/pairs12 arrays (and 13, 14 alike) store which particle pairs are
 * 1-2 (directly bonded), 1-3 and 1-4 connected, in compact form: the bonded
 * partners of particle i are pairs12[head12[i]] up to pairs12[head12[i+1]].
 * They are only built when the corresponding scaling factor differs from one
 * (see @ref Parameters.factor_12_nb).
 */
struct Nbrlist
{
    struct Celllist *p_celllist;   //!< cell-linked-list used to build the neighbor list
    size_t num_nbrs, num_nbrs_max; //!< number of pairs in the list, and the number allocated for
    struct Pair *nbr;              //!< the pairs of the neighbor list
    struct DeltaR *dr;             //!< per-particle displacement since the list was last built; drives the rebuild decision
    size_t *head12, *pairs12;          //!< 1-2 (bonded) partners of each particle
    size_t *head13, *pairs13;          //!< 1-3 partners (two bonds apart) of each particle
    size_t *head14, *pairs14;          //!< 1-4 partners (three bonds apart) of each particle
};

#endif /* TYPES_MD_H_ */
