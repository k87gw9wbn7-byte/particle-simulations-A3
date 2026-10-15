#include <math.h>
#include <stdio.h>
#include "chainstats.h"
int main(void) {
 struct Parameters p={0};struct Vectors v={0};
 Vec3D positions[4]={{9,1,1},{1,1,1},{3,1,1},{5,1,1}};
 p.L=v3(10,10,10);p.chain_length=4;v.r=positions;
 if(fabs(chain_re2(&p,&v,0)-36)>1e-12)return 1;
 positions[0]=v3(9.6,1,1);positions[1]=v3(.2,1,1);positions[2]=v3(.8,1,1);p.chain_length=3;
 if(fabs(chain_re2(&p,&v,0)-1.44)>1e-12)return 2;
 puts("PASS end-to-end bond unwrapping, including endpoint separation beyond half a box");return 0;
}
