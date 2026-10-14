#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "setparameters.h"
#include "velocity.h"
#define CHECK(x) do { if(!(x)){fprintf(stderr,"Velocity check failed at line %d\n",__LINE__);return 1;} }while(0)
int main(void){
    struct Parameters p;set_parameters(&p);p.num_part=2;p.velocity_enabled=1;
    p.equil_steps=0;p.sample_steps=1;p.block_steps=2;
    snprintf(p.analysis_prefix,sizeof p.analysis_prefix,"../data/b_velocity_fixture");
    Vec3D vdata[2]={ {1,2,3},{-1,-2,-3} };
    struct Vectors v={0};v.v=vdata;struct VelocityAnalysis a;
    velocity_start(&p,&a);velocity_update(&p,&v,&a,0);CHECK(a.frames==0);
    velocity_update(&p,&v,&a,1);CHECK(a.frames==1 && a.block_frames==1);
    CHECK(a.sum[0]==0 && a.sum[1]==0 && a.sum[2]==0);
    CHECK(a.sum2[0]==2 && a.sum2[1]==8 && a.sum2[2]==18);
    CHECK(a.hist[0][80]==1 && a.hist[0][120]==1 && a.speed[74]==2);
    velocity_update(&p,&v,&a,2);CHECK(a.frames==2 && a.block_frames==0 && a.block==1);
    velocity_finish(&p,&a);
    puts("PASS velocity bins, moments, sample timing and blocks; fixture expected raw kT=14/3, corrected=28/3");
    return 0;
}
