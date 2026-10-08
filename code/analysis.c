#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "analysis.h"
#include "nbrlist.h"

static FILE *open_output(const char *name)
{
    FILE *f=fopen(name,"w");
    if (!f) { perror(name); exit(EXIT_FAILURE); }
    return f;
}

void analysis_alloc(const struct Parameters *p, struct Vectors *v)
{
    v->analysis=NULL;
    if (!p->analysis_enabled) return;
    struct Analysis *a=calloc(1,sizeof(*a));
    if (!a) { perror("analysis"); exit(EXIT_FAILURE); }
    v->analysis=a;
    a->bins=(size_t)ceil(p->rdf_max/p->rdf_dr);
    a->hist=calloc(3*a->bins,sizeof(*a->hist));
    a->block_hist=calloc(3*a->bins,sizeof(*a->block_hist));
    if (!a->hist || !a->block_hist) { perror("RDF histogram"); exit(EXIT_FAILURE); }
    a->rdf_parameters=*p;
    a->rdf_parameters.r_cut=p->rdf_max;
    a->rdf_parameters.r_shell=0;
    /* The inherited half-stencil requires >=3 cells in every direction.
     * Small boxes use an exact i<j loop, still restricted to r<=Lmin/2. */
    double minL=fmin(p->L.x,fmin(p->L.y,p->L.z));
    a->use_list=(minL >= 3*p->rdf_max);
    if (a->use_list) alloc_nbrlist(&a->rdf_parameters,&a->rdf_list);
}

void analysis_free(struct Vectors *v)
{
    struct Analysis *a=v->analysis;
    if (!a) return;
    if (a->block_file) fclose(a->block_file);
    if (a->use_list) free_nbrlist(&a->rdf_list);
    free(a->hist); free(a->block_hist); free(a); v->analysis=NULL;
}

void analysis_start(const struct Parameters *p, struct Vectors *v)
{
    struct Analysis *a=v->analysis;
    if (!a) return;
    for (size_t i=0;i<p->num_part;++i) ++a->count[v->type[i]];
    char path[1100];
    snprintf(path,sizeof(path),"%s_rdf_blocks.csv",p->analysis_prefix);
    a->block_file=open_output(path);
    fprintf(a->block_file,"block,frames,r_lo,r_hi,r,g_AA,g_AB,g_BB,count_AA,count_AB,count_BB\n");
}

static void count_pair(const struct Parameters *p, struct Vectors *v, size_t i, size_t j)
{
    struct Analysis *a=v->analysis;
    Vec3D d=v3_min_image(v3_sub(v->r[i],v->r[j]),p->L);
    double r2=v3_dot(d,d);
    if (r2 >= p->rdf_max*p->rdf_max) return;
    size_t b=(size_t)(sqrt(r2)/p->rdf_dr);
    if (b>=a->bins) return;
    size_t k=3*b+(size_t)pair_index(v->type[i],v->type[j]);
    ++a->hist[k]; ++a->block_hist[k];
}

void rdf_sample(const struct Parameters *p, struct Vectors *v)
{
    struct Analysis *a=v->analysis;
    if (!a) return;
    if (a->use_list) {
        build_nbrlist(&a->rdf_parameters,v,&a->rdf_list);
        for (size_t q=0;q<a->rdf_list.num_nbrs;++q)
            count_pair(p,v,a->rdf_list.nbr[q].i,a->rdf_list.nbr[q].j);
    } else {
        for (size_t i=0;i<p->num_part;++i)
            for (size_t j=i+1;j<p->num_part;++j) count_pair(p,v,i,j);
    }
    ++a->frames; ++a->block_frames;
}

static void write_hist(FILE *f, const struct Parameters *p, const struct Analysis *a,
                       const unsigned long long *hist, size_t frames, int block)
{
    double na=(double)a->count[0], nb=(double)a->count[1];
    double pairs[3]={na*(na-1)/2,na*nb,nb*(nb-1)/2};
    double V=p->L.x*p->L.y*p->L.z;
    for (size_t b=0;b<a->bins;++b) {
        double lo=b*p->rdf_dr, hi=fmin((b+1)*p->rdf_dr,p->rdf_max);
        double shell=(4.0*PI/3.0)*(hi*hi*hi-lo*lo*lo), g[3];
        for (int t=0;t<3;++t) {
            double expected=frames*pairs[t]*shell/V;
            g[t]=expected>0 ? hist[3*b+t]/expected : NAN;
        }
        if (block) fprintf(f,"%zu,",a->block);
        fprintf(f,"%zu,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%llu,%llu,%llu\n",
                frames,lo,hi,0.5*(lo+hi),g[0],g[1],g[2],
                hist[3*b],hist[3*b+1],hist[3*b+2]);
    }
}

static void finish_block(const struct Parameters *p, struct Vectors *v)
{
    struct Analysis *a=v->analysis;
    if (!a->block_frames) return;
    write_hist(a->block_file,p,a,a->block_hist,a->block_frames,1);
    fflush(a->block_file);
    memset(a->block_hist,0,3*a->bins*sizeof(*a->block_hist));
    a->block_frames=0; ++a->block;
}

void analysis_update(const struct Parameters *p, struct Vectors *v, size_t step)
{
    if (!v->analysis || step<=p->equil_steps) return;
    size_t elapsed=step-p->equil_steps;
    if (elapsed%p->sample_steps==0) rdf_sample(p,v);
    if (elapsed%p->block_steps==0) finish_block(p,v);
}

void analysis_finish(const struct Parameters *p, struct Vectors *v)
{
    struct Analysis *a=v->analysis;
    if (!a) return;
    finish_block(p,v);
    char path[1100];
    snprintf(path,sizeof(path),"%s_rdf.csv",p->analysis_prefix);
    FILE *f=open_output(path);
    fprintf(f,"frames,r_lo,r_hi,r,g_AA,g_AB,g_BB,count_AA,count_AB,count_BB\n");
    write_hist(f,p,a,a->hist,a->frames,0);
    fclose(f); fclose(a->block_file); a->block_file=NULL;
    printf("RDF: %zu frames, %zu blocks, dr=%g, rmax=%g, NA=%zu, NB=%zu\n",
           a->frames,a->block,p->rdf_dr,p->rdf_max,a->count[0],a->count[1]);
}
