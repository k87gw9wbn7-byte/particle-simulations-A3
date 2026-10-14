#ifndef DPD_VELOCITY_H
#define DPD_VELOCITY_H
#include <stdio.h>
#include "structs.h"
/* Component histogram [-5,5), width .05; speed histogram [0,8), width .05. */
#define V_BINS 200
#define S_BINS 160
#define V_MIN (-5.0)
#define V_WIDTH 0.05
struct VelocityAnalysis {
    FILE *blocks;
    size_t frames,block_frames,block;
    unsigned long long hist[3][V_BINS],speed[S_BINS],under[3],over[3],speed_over;
    double sum[3],sum2[3],block_sum[3],block_sum2[3];
};
void velocity_start(const struct Parameters *,struct VelocityAnalysis *);
void velocity_update(const struct Parameters *,const struct Vectors *,struct VelocityAnalysis *,size_t);
void velocity_finish(const struct Parameters *,struct VelocityAnalysis *);
#endif
