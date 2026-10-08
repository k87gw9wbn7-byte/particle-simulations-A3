#ifndef DPD_ANALYSIS_H_
#define DPD_ANALYSIS_H_
#include <stdio.h>
#include "structs.h"
/* Separate RDF neighbour list: never changes the force list or draws noise. */
struct Analysis {
    struct Parameters rdf_parameters;
    struct Nbrlist rdf_list;
    int use_list;
    size_t bins, frames, block_frames, block, count[TYPES];
    unsigned long long *hist, *block_hist; /* packed bin*3 + pair_index */
    FILE *block_file;
};
void analysis_alloc(const struct Parameters *, struct Vectors *);
void analysis_free(struct Vectors *);
void analysis_start(const struct Parameters *, struct Vectors *);
void analysis_update(const struct Parameters *, struct Vectors *, size_t);
void analysis_finish(const struct Parameters *, struct Vectors *);
/* Exposed for the independent uniform-position normalisation check. */
void rdf_sample(const struct Parameters *, struct Vectors *);
#endif
