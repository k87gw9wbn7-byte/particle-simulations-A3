#ifndef CHAINSTATS_H
#define CHAINSTATS_H
#include <stdio.h>
#include "structs.h"
struct ChainStats { FILE *file; size_t frames, block; double re2, b2; };
double chain_re2(const struct Parameters *, const struct Vectors *, size_t);
void chainstats_start(const struct Parameters *, struct ChainStats *);
void chainstats_update(const struct Parameters *, const struct Vectors *, struct ChainStats *, size_t);
void chainstats_finish(struct ChainStats *);
#endif
