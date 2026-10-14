#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "slab.h"
static FILE *output(const struct Parameters *p,const char *suffix) {
    char name[1100]; snprintf(name,sizeof name,"%s%s",p->analysis_prefix,suffix);
    FILE *f=fopen(name,"w"); if(!f){perror(name);exit(EXIT_FAILURE);} return f;
}
void slab_start(const struct Parameters *p, struct SlabAnalysis *a) {
    memset(a,0,sizeof *a); if(!p->slab_enabled)return;
    a->bins=(size_t)llround(p->L.x/p->profile_dr);
    a->hist=calloc(2*a->bins,sizeof(double)); a->lab=calloc(2*a->bins,sizeof(double));
    a->block_hist=calloc(2*a->bins,sizeof(double)); a->block_lab=calloc(2*a->bins,sizeof(double));
    if(!a->hist||!a->lab||!a->block_hist||!a->block_lab){perror("slab allocation");exit(EXIT_FAILURE);}
    a->profiles=output(p,"_profile_blocks.csv");
    fputs("block,start_step,end_step,t_start,t_end,production,frames,aligned_fraction,order,x,rho_A,rho_B,phi_A,rho_A_lab,rho_B_lab\n",a->profiles);
    a->pressure=output(p,"_pressure_samples.csv");
    fputs("step,time,production,Pxx,Pyy,Pzz,sigma_two_interfaces,order,shift,aligned\n",a->pressure);
    a->pressure_blocks=output(p,"_pressure_blocks.csv");
    fputs("block,start_step,end_step,t_start,t_end,production,frames,Pxx,Pyy,Pzz,sigma_two_interfaces,order\n",a->pressure_blocks);
}
static void finish_block(const struct Parameters *p,struct SlabAnalysis *a,size_t end) {
    if(!a->block_frames)return;
    size_t start=end-p->block_steps;
    int prod=start>=p->equil_steps;
    double normal=a->block_frames*p->profile_dr*p->L.y*p->L.z;
    for(size_t b=0;b<a->bins;++b) {
        double ra=a->block_hist[2*b]/normal,rb=a->block_hist[2*b+1]/normal;
        fprintf(a->profiles,"%zu,%zu,%zu,%.17g,%.17g,%d,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                a->block,start,end,start*p->dt,end*p->dt,prod,a->block_frames,
                (double)a->aligned_frames/a->block_frames,a->sumOrder/a->block_frames,
                (b+0.5)*p->profile_dr,ra,rb,ra+rb>0?ra/(ra+rb):NAN,
                a->block_lab[2*b]/normal,a->block_lab[2*b+1]/normal);
    }
    double px=a->sumP[0]/a->block_frames,py=a->sumP[1]/a->block_frames,pz=a->sumP[2]/a->block_frames;
    fprintf(a->pressure_blocks,"%zu,%zu,%zu,%.17g,%.17g,%d,%zu,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            a->block,start,end,start*p->dt,end*p->dt,prod,a->block_frames,
            px,py,pz,0.5*p->L.x*(px-0.5*(py+pz)),a->sumOrder/a->block_frames);
    fflush(a->profiles); fflush(a->pressure_blocks);
    memset(a->block_hist,0,2*a->bins*sizeof(double));memset(a->block_lab,0,2*a->bins*sizeof(double));
    memset(a->sumP,0,sizeof a->sumP);a->sumOrder=0;a->block_frames=0;a->aligned_frames=0;++a->block;
}
void slab_update(const struct Parameters *p,const struct Vectors *v,struct SlabAnalysis *a,size_t step) {
    if(!p->slab_enabled||step==0)return;
    if(step%p->sample_steps==0) {
        double re=0,im=0,k=2*PI/p->L.x;
        Vec3D kinetic=v3(0,0,0);
        for(size_t i=0;i<p->num_part;++i) {
            double sign=v->type[i]==0?1:-1;
            re+=sign*cos(k*v->r[i].x);im+=sign*sin(k*v->r[i].x);
            kinetic=v3_add(kinetic,v3_scl(p->mass[v->type[i]],v3_mul(v->v[i],v->v[i])));
        }
        double order=hypot(re,im)/p->num_part;
        int align=p->a[pair_index(0,1)]>25 && order>=0.15;
        /* Pure A-left/B-right: phase=pi/2 and A centre=Lx/4. Remove only
         * translation along x; never rotate the slab or change particle data. */
        double shift=align?atan2(im,re)/k-p->L.x/4:0;
        int prod=step>p->equil_steps;
        for(size_t i=0;i<p->num_part;++i) {
            double x=v->r[i].x-shift;x-=p->L.x*floor(x/p->L.x);
            size_t b=(size_t)(x/p->profile_dr),bl=(size_t)(v->r[i].x/p->profile_dr);
            if(b>=a->bins)b=a->bins-1;
            if(bl>=a->bins)bl=a->bins-1;
            size_t t=(size_t)v->type[i];
            ++a->block_hist[2*b+t];++a->block_lab[2*bl+t];
            if(prod){++a->hist[2*b+t];++a->lab[2*bl+t];}
        }
        double V=p->L.x*p->L.y*p->L.z;
        Vec3D P=v3_scl(1/V,v3_add(kinetic,v->virial_diag));
        double trace=(P.x+P.y+P.z)/3;
        double scalar=v3_dot(kinetic,v3(1,1,1))/(3*V)+v->press_vir_nb+v->press_vir_bnd;
        if(fabs(trace-scalar)>1e-9*fmax(1,fabs(scalar))) {
            fprintf(stderr,"Pressure trace disagrees with scalar virial\n");exit(EXIT_FAILURE);
        }
        fprintf(a->pressure,"%zu,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",
                step,step*p->dt,prod,P.x,P.y,P.z,0.5*p->L.x*(P.x-0.5*(P.y+P.z)),order,shift,align);
        ++a->block_frames;if(prod)++a->frames;
        a->aligned_frames+=(size_t)align;a->sumOrder+=order;
        a->sumP[0]+=P.x;a->sumP[1]+=P.y;a->sumP[2]+=P.z;
    }
    if(step%p->block_steps==0)finish_block(p,a,step);
}
void slab_finish(const struct Parameters *p,struct SlabAnalysis *a) {
    if(!p->slab_enabled)return;
    FILE *f=output(p,"_profile.csv");
    fputs("frames,x,rho_A,rho_B,phi_A,rho_A_lab,rho_B_lab\n",f);
    double normal=a->frames*p->profile_dr*p->L.y*p->L.z;
    for(size_t b=0;b<a->bins;++b) {
        double ra=a->hist[2*b]/normal,rb=a->hist[2*b+1]/normal;
        fprintf(f,"%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",a->frames,(b+0.5)*p->profile_dr,
                ra,rb,ra+rb>0?ra/(ra+rb):NAN,a->lab[2*b]/normal,a->lab[2*b+1]/normal);
    }
    fclose(f);fclose(a->profiles);fclose(a->pressure);fclose(a->pressure_blocks);
    free(a->hist);free(a->lab);free(a->block_hist);free(a->block_lab);
    printf("Slabs: %zu production frames, %zu total monitoring blocks\n",a->frames,a->block);
}
