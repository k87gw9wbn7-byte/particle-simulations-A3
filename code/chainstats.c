#include <stdlib.h>
#include <string.h>
#include "chainstats.h"
/* Sum consecutive minimum-image bonds. Applying a single minimum image to
 * the two endpoints would fold extended chains back into the periodic box. */
double chain_re2(const struct Parameters *p, const struct Vectors *v, size_t c) {
    Vec3D R=v3(0,0,0); size_t first=c*p->chain_length;
    for(size_t j=1;j<p->chain_length;++j)
        R=v3_add(R,v3_min_image(v3_sub(v->r[first+j],v->r[first+j-1]),p->L));
    return v3_dot(R,R);
}
void chainstats_start(const struct Parameters *p, struct ChainStats *a) {
    memset(a,0,sizeof *a); if(p->chain_length<2)return;
    char name[1100];snprintf(name,sizeof name,"%s_chain_blocks.csv",p->analysis_prefix);
    a->file=fopen(name,"w");if(!a->file){perror(name);exit(EXIT_FAILURE);}
    fputs("block,t_start,t_end,production,frames,chains,Re2,bond2\n",a->file);
}
void chainstats_update(const struct Parameters *p,const struct Vectors *v,struct ChainStats *a,size_t step) {
    if(!a->file)return;
    if(step%p->sample_steps==0) {
        size_t chains=p->num_part/p->chain_length; double R=0,b=0;
        for(size_t c=0;c<chains;++c)R+=chain_re2(p,v,c);
        for(size_t k=0;k<v->num_bonds;++k) {
            Vec3D d=v3_min_image(v3_sub(v->r[v->bonds[k].i],v->r[v->bonds[k].j]),p->L);
            b+=v3_dot(d,d);
        }
        a->re2+=R/chains;a->b2+=b/v->num_bonds;++a->frames;
    }
    if(step%p->block_steps==0 && a->frames) {
        size_t start=step-p->block_steps;
        fprintf(a->file,"%zu,%.17g,%.17g,%d,%zu,%zu,%.17g,%.17g\n",
            a->block,start*p->dt,step*p->dt,start>=p->equil_steps,a->frames,
            p->num_part/p->chain_length,a->re2/a->frames,a->b2/a->frames);
        fflush(a->file);++a->block;a->frames=0;a->re2=a->b2=0;
    }
}
void chainstats_finish(struct ChainStats *a) { if(a->file)fclose(a->file); }
