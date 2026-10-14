#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "velocity.h"
static FILE *output(const struct Parameters *p,const char *suffix){
    char path[1100];snprintf(path,sizeof path,"%s%s",p->analysis_prefix,suffix);
    FILE *f=fopen(path,"w");if(!f){perror(path);exit(EXIT_FAILURE);}return f;
}
static void moment_header(FILE *f){
    fprintf(f,"block,frames,beads_per_frame,n_per_component,mean_vx,mean_vy,mean_vz,var_vx,var_vy,var_vz,kT_variance,kT_dof_corrected\n");
}
static void moments(FILE *f,const struct Parameters *p,size_t block,size_t frames,
                    const double sum[3],const double sum2[3]){
    double n=(double)p->num_part*frames,mean[3],var[3],T=0;
    for(int k=0;k<3;k++){
        mean[k]=sum[k]/n;var[k]=sum2[k]/n-mean[k]*mean[k];
        T+=p->mass[0]*var[k]/3;
    }
    fprintf(f,"%zu,%zu,%zu,%.0f,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
        block,frames,p->num_part,n,mean[0],mean[1],mean[2],var[0],var[1],var[2],T,T*p->num_part/(p->num_part-1.0));
}
void velocity_start(const struct Parameters *p,struct VelocityAnalysis *a){
    memset(a,0,sizeof *a);if(!p->velocity_enabled)return;
    if(p->mass[0]!=p->mass[1]){fprintf(stderr,"B3 pooled velocities require equal masses\n");exit(EXIT_FAILURE);}
    a->blocks=output(p,"_velocity_blocks.csv");moment_header(a->blocks);
}
static void finish_block(const struct Parameters *p,struct VelocityAnalysis *a){
    if(!a->block_frames)return;
    moments(a->blocks,p,a->block,a->block_frames,a->block_sum,a->block_sum2);
    fflush(a->blocks);a->block++;a->block_frames=0;
    memset(a->block_sum,0,sizeof a->block_sum);memset(a->block_sum2,0,sizeof a->block_sum2);
}
void velocity_update(const struct Parameters *p,const struct Vectors *v,struct VelocityAnalysis *a,size_t step){
    if(!p->velocity_enabled || step<=p->equil_steps)return;
    size_t elapsed=step-p->equil_steps;
    if(elapsed%p->sample_steps==0){
        for(size_t i=0;i<p->num_part;i++){
            double c[3]={v->v[i].x,v->v[i].y,v->v[i].z};
            for(int k=0;k<3;k++){
                if(!isfinite(c[k])){fprintf(stderr,"Nonfinite sampled velocity\n");exit(EXIT_FAILURE);}
                a->sum[k]+=c[k];a->sum2[k]+=c[k]*c[k];
                a->block_sum[k]+=c[k];a->block_sum2[k]+=c[k]*c[k];
                if(c[k]<V_MIN)a->under[k]++;
                else if(c[k]>=V_MIN+V_BINS*V_WIDTH)a->over[k]++;
                else {
                    size_t b=(size_t)((c[k]-V_MIN)/V_WIDTH);
                    if(b<V_BINS)a->hist[k][b]++;else a->over[k]++;
                }
            }
            double speed=v3_norm(v->v[i]);
            if(speed>=S_BINS*V_WIDTH)a->speed_over++;
            else {size_t b=(size_t)(speed/V_WIDTH);if(b<S_BINS)a->speed[b]++;else a->speed_over++;}
        }
        a->frames++;a->block_frames++;
    }
    if(elapsed%p->block_steps==0)finish_block(p,a);
}
void velocity_finish(const struct Parameters *p,struct VelocityAnalysis *a){
    if(!p->velocity_enabled)return;
    if(!a->frames){fprintf(stderr,"No velocity samples: reduce equil_steps or increase run length\n");exit(EXIT_FAILURE);}
    finish_block(p,a);fclose(a->blocks);a->blocks=NULL;
    FILE *f=output(p,"_velocity_moments.csv");moment_header(f);
    moments(f,p,0,a->frames,a->sum,a->sum2);fclose(f);
    double n=(double)p->num_part*a->frames,den=n*V_WIDTH;
    f=output(p,"_velocity_hist.csv");
    fprintf(f,"v_lo,v_hi,v,p_vx,p_vy,p_vz,p_pool,count_vx,count_vy,count_vz,n_per_component,frames\n");
    for(size_t b=0;b<V_BINS;b++){
        double lo=V_MIN+b*V_WIDTH;unsigned long long *h0=a->hist[0],*h1=a->hist[1],*h2=a->hist[2];
        fprintf(f,"%.10g,%.10g,%.10g,%.12g,%.12g,%.12g,%.12g,%llu,%llu,%llu,%.0f,%zu\n",
            lo,lo+V_WIDTH,lo+V_WIDTH/2,h0[b]/den,h1[b]/den,h2[b]/den,
            (h0[b]+h1[b]+h2[b])/(3*den),h0[b],h1[b],h2[b],n,a->frames);
    }
    fclose(f);f=output(p,"_speed_hist.csv");
    fprintf(f,"s_lo,s_hi,s,p_speed,count,n_samples,frames\n");
    for(size_t b=0;b<S_BINS;b++){
        double lo=b*V_WIDTH;
        fprintf(f,"%.10g,%.10g,%.10g,%.12g,%llu,%.0f,%zu\n",lo,lo+V_WIDTH,lo+V_WIDTH/2,a->speed[b]/den,a->speed[b],n,a->frames);
    }
    fclose(f);f=output(p,"_velocity_tails.csv");
    fprintf(f,"channel,underflow,overflow,n_samples\n");
    for(int k=0;k<3;k++)fprintf(f,"v%d,%llu,%llu,%.0f\n",k,a->under[k],a->over[k],n);
    fprintf(f,"speed,0,%llu,%.0f\n",a->speed_over,n);fclose(f);
    printf("Velocity statistics: %zu frames, %zu blocks, component bin width=%g\n",a->frames,a->block,V_WIDTH);
}
