#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "forces.h"
#include "setparameters.h"
#include "slab.h"
#define CHECK(x) do {if(!(x)){fprintf(stderr,"Tensor test failed line %d\n",__LINE__);return 1;}} while(0)
static int close_to(double a,double b){return fabs(a-b)<1e-10;}
int main(void){
    struct Parameters p;set_parameters(&p);p.num_part=2;p.thermostat_on=0;
    Vec3D r[2]={{.3,.4,0},{0,0,0}},vel[2]={{1,2,3},{-1,-2,-3}},f[2];int types[2]={0,1};
    struct Pair pair={0,1,1};struct Nbrlist n={0};n.nbr=&pair;n.num_nbrs=1;
    struct Vectors v={0};v.r=r;v.v=vel;v.f=f;v.type=types;
    calculate_forces(&p,&n,&v);
    CHECK(close_to(v.virial_diag.x,2.25)&&close_to(v.virial_diag.y,4)&&close_to(v.virial_diag.z,0));
    CHECK(close_to((v.virial_diag.x+v.virial_diag.y+v.virial_diag.z)/3000,v.press_vir_nb));
    /* D/R affects the force array but cannot enter the stored C virial. */
    p.thermostat_on=1;srand(13);calculate_forces(&p,&n,&v);CHECK(close_to(v.virial_diag.x,2.25));
    p.thermostat_on=0;p.conservative_on=0;struct Bond bond={0,1};v.bonds=&bond;v.num_bonds=1;
    calculate_forces(&p,&n,&v);CHECK(close_to(v.virial_diag.x,-.18)&&close_to(v.virial_diag.y,-.32));
    v.num_bonds=0;calculate_forces(&p,&n,&v);CHECK(close_to(v.virial_diag.x,0)&&close_to(v.virial_diag.y,0));
    /* Known velocities give P=(.002,.008,.018), sigma=5*(.002-.013)=-.055. */
    p.slab_enabled=1;p.profile_dr=.25;p.sample_steps=1;p.block_steps=2;p.equil_steps=0;
    snprintf(p.analysis_prefix,sizeof p.analysis_prefix,"../data/tensor_fixture");
    struct SlabAnalysis a;slab_start(&p,&a);slab_update(&p,&v,&a,1);slab_update(&p,&v,&a,2);
    CHECK(a.frames==2&&a.block==1);slab_finish(&p,&a);
    puts("PASS oblique conservative virial, D/R exclusion, spring virial, reset, kinetic tensor and sampler fixture");return 0;
}
