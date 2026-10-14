#ifndef SLAB_H
#define SLAB_H
#include <stdio.h>
#include "structs.h"
struct SlabAnalysis {
    size_t bins, frames, block_frames, block, aligned_frames;
    double *hist, *lab, *block_hist, *block_lab;
    double sumP[3], sumOrder;
    FILE *profiles, *pressure, *pressure_blocks;
};
void slab_start(const struct Parameters *p, struct SlabAnalysis *a);
void slab_update(const struct Parameters *p, const struct Vectors *v,
                 struct SlabAnalysis *a, size_t step);
void slab_finish(const struct Parameters *p, struct SlabAnalysis *a);
#endif
