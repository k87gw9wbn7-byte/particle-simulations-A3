/* Focused Part A verification, separate executable (not part of *.c build). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "setparameters.h"
#include "forces.h"
#include "initialise.h"
#include "memory.h"
#include "nbrlist.h"
#include "force_test.h"
#include "dynamics.h"
#include "analysis.h"
#include "random.h"
#include "fileoutput.h"
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#x); exit(1); } } while(0)
static void near(double a,double b,double tol) { CHECK(isfinite(a)); CHECK(fabs(a-b)<=tol); }
static void setup(struct Parameters *p,struct Vectors *v,struct Nbrlist *n,size_t N,size_t length)
{
    set_parameters(p); p->num_part=N; p->chain_length=length; p->analysis_enabled=0;
    alloc_memory(p,v,n);
    size_t step; double time;
    initialise(p,v,n,&step,&time);
}
static void pair_tests(void)
{
    struct Parameters p; struct Vectors v; struct Nbrlist n;
    setup(&p,&v,&n,2,1);
    v.r[0]=v3(5.5,5,5); v.r[1]=v3(5,5,5);
    v.v[0]=v3(1,0,0); v.v[1]=v3(0,0,0);
    p.a[0]=25; p.a[1]=37; p.a[2]=43;
    build_nbrlist(&p,&v,&n); CHECK(n.num_nbrs==1);
    for(int a=0;a<2;++a) for(int b=0;b<2;++b) {
        v.type[0]=a;v.type[1]=b;
        double rep=p.a[pair_index(a,b)];
        for(int c=0;c<2;++c) for(int t=0;t<2;++t) {
            p.conservative_on=c;p.thermostat_on=t;
            srand(97);double z=dpd_noise();srand(97);
            double e=calculate_forces(&p,&n,&v);
            double expected=c*rep*0.5+t*(-p.gamma*0.25+3*0.5*z/sqrt(p.dt));
            near(v.f[0].x,expected,1e-12);near(v.f[1].x,-expected,1e-12);
            near(v3_norm(v3_add(v.f[0],v.f[1])),0,1e-12);
            near(e,c*0.125*rep,1e-12);
            near(v.press_vir_nb,c*0.25*rep/3000,1e-12);
        }
    }
    p.conservative_on=1;p.thermostat_on=0;
    /* A near-boundary pair is separated by 0.5 using minimum image. */
    v.r[0]=v3(0.2,5,5);v.r[1]=v3(9.7,5,5);build_nbrlist(&p,&v,&n);
    near(calculate_forces(&p,&n,&v),0.125*43,1e-10);
    near(v.f[0].x,0.5*43,1e-10);
    v.r[0]=v3(1,5,5);v.r[1]=v3(0,5,5);build_nbrlist(&p,&v,&n);
    near(calculate_forces(&p,&n,&v),0,0);near(v3_norm(v.f[0]),0,0);
    /* Dissipation extracted by averaging two equal-noise runs with +/- v. */
    v.r[0]=v3(5.5,5,5);v.r[1]=v3(5,5,5);build_nbrlist(&p,&v,&n);
    p.conservative_on=0;p.thermostat_on=1;
    v.v[0]=v3(1,0,0); srand(8);calculate_forces(&p,&n,&v);double fp=v.f[0].x;
    v.v[0]=v3(-1,0,0);srand(8);calculate_forces(&p,&n,&v);double fm=v.f[0].x;
    near((fp-fm)/2,-4.5*0.25,1e-12);
    free_memory(&v,&n);puts("PASS pair forces, switches, types, cutoff, periodicity, momentum, conservative-only energy/virial");
}
static void finite_difference_tests(void)
{
    struct Parameters p;struct Vectors v;struct Nbrlist n;
    setup(&p,&v,&n,2,2);p.thermostat_on=0;
    v.r[0]=v3(5.2,5.3,5.1);v.r[1]=v3(5,5,5);
    build_nbrlist(&p,&v,&n);
    CHECK(forces_test(forces_conservative_only,0,&p,&n,&v)<1e-7);
    CHECK(forces_test(forces_bond_only,0,&p,&n,&v)<1e-7);
    CHECK(virial_test(calculate_forces,&p,&n,&v)<1e-7);
    /* Bonded DPD pair is retained: U = C*r^2/2 + U_C. */
    double r=v3_norm(v3_sub(v.r[0],v.r[1]));
    near(calculate_forces(&p,&n,&v),r*r+12.5*(1-r)*(1-r),1e-12);
    p.conservative_on=0;v.r[0]=v.r[1];
    near(calculate_forces(&p,&n,&v),0,0);near(v3_norm(v.f[0]),0,0);
    free_memory(&v,&n);puts("PASS finite-difference conservative/spring forces, virial, bonded repulsion, zero-length spring");
}
static void noise_test(void)
{
    const size_t n=200000;double s=0,ss=0;
    srand(123);
    for(size_t i=0;i<n;++i) {double z=dpd_noise();CHECK(fabs(z)<=sqrt(3.0));s+=z;ss+=z*z;}
    double mean=s/n,var=ss/n-mean*mean;
    printf("Noise samples=%zu mean=%.9g variance=%.9g\n",n,mean,var);
    CHECK(fabs(mean)<0.015);CHECK(fabs(var-1)<0.015);puts("PASS uniform noise moments");
}
static void chain_tests(void)
{
    size_t lengths[]={1,3,5,10};
    for(size_t l=0;l<4;++l) for(int separated=0;separated<2;++separated) {
        struct Parameters p;struct Vectors v;struct Nbrlist n;
        setup(&p,&v,&n,120,lengths[l]);p.separated_start=separated;
        initialise_positions(&p,&v);
        size_t N=p.chain_length;
        CHECK(v.num_bonds==(120/N)*(N-1));
        CHECK(v.num_angles==(120/N)*(N>2?N-2:0));
        CHECK(v.num_dihedrals==(120/N)*(N>3?N-3:0));
        for(size_t i=0;i<120;++i) {
            CHECK(v.type[i]==v.type[(i/N)*N]);
            if(separated) CHECK((v.r[i].x<5)==(v.type[i]==0));
        }
        for(size_t b=0;b<v.num_bonds;++b) {
            size_t i=v.bonds[b].i,j=v.bonds[b].j;CHECK(j==i+1);CHECK(i/N==j/N);
            near(v3_norm(v3_min_image(v3_sub(v.r[i],v.r[j]),p.L)),0.8,1e-12);
        }
        build_nbrlist(&p,&v,&n);calculate_forces(&p,&n,&v);
        Vec3D momentum=v3(0,0,0);
        for(int step=0;step<50;++step) {
            update_velocities_half_dt(&p,&n,&v);update_positions(&p,&n,&v);
            boundary_conditions(&p,&v);update_nbrlist(&p,&v,&n);
            CHECK(isfinite(calculate_forces(&p,&n,&v)));
            CHECK(isfinite(update_velocities_half_dt(&p,&n,&v)));
        }
        for(size_t i=0;i<120;++i) momentum=v3_add(momentum,v.v[i]);
        CHECK(v3_norm(momentum)<1e-9);
        /* Compare the reused/rebuilt force list with a fresh all-pairs sum. */
        p.thermostat_on=0;
        size_t saved_bonds=v.num_bonds;v.num_bonds=0;
        double actual=calculate_forces(&p,&n,&v), expected=0;
        Vec3D brute[120];for(size_t i=0;i<120;++i) brute[i]=v3(0,0,0);
        for(size_t i=0;i<120;++i) for(size_t j=i+1;j<120;++j) {
            Vec3D d=v3_min_image(v3_sub(v.r[i],v.r[j]),p.L);
            double r=v3_norm(d);
            if(r<p.r_cut) {
                double a=p.a[pair_index(v.type[i],v.type[j])],w=1-r/p.r_cut;
                expected+=0.5*a*p.r_cut*w*w;
                Vec3D f=v3_scl(a*w/r,d);brute[i]=v3_add(brute[i],f);brute[j]=v3_sub(brute[j],f);
            }
        }
        near(actual,expected,1e-9);
        for(size_t i=0;i<120;++i) near(v3_norm(v3_sub(v.f[i],brute[i])),0,1e-9);
        v.num_bonds=saved_bonds;free_memory(&v,&n);
    }
    puts("PASS chains N=1,3,5,10; homogeneous/separated starts; lengths; 50-step dynamics; momentum; neighbour forces vs all pairs");
}
static void rdf_test(void)
{
    struct Parameters p;struct Vectors v;struct Nbrlist n;
    set_parameters(&p);p.num_part=500;p.fraction_A=0.3;p.equil_steps=0;
    p.sample_steps=1;p.block_steps=10;p.conservative_on=0;
    snprintf(p.analysis_prefix,sizeof(p.analysis_prefix),"../data/a_uniform");
    alloc_memory(&p,&v,&n);size_t step;double t;initialise(&p,&v,&n,&step,&t);
    analysis_start(&p,&v);
    struct Analysis *a=v.analysis;
    unsigned long long *brute=calloc(3*a->bins,sizeof(*brute));CHECK(brute);
    for(size_t frame=1;frame<=100;++frame) {
        initialise_positions(&p,&v); /* independent uniform positions, NOT an MD run */
        if(frame==1) {
            build_nbrlist(&p,&v,&n);
            size_t force_pairs=n.num_nbrs;
            for(size_t i=0;i<p.num_part;++i) for(size_t j=i+1;j<p.num_part;++j) {
                double r=v3_norm(v3_min_image(v3_sub(v.r[i],v.r[j]),p.L));
                if(r<p.rdf_max) ++brute[3*(size_t)(r/p.rdf_dr)+pair_index(v.type[i],v.type[j])];
            }
            analysis_update(&p,&v,frame);
            CHECK(n.num_nbrs==force_pairs);
            for(size_t k=0;k<3*a->bins;++k) CHECK(brute[k]==a->hist[k]);
        } else analysis_update(&p,&v,frame);
    }
    double na=a->count[0],nb=a->count[1],pairs[3]={na*(na-1)/2,na*nb,nb*(nb-1)/2};
    for(int type=0;type<3;++type) {
        unsigned long long count=0;
        for(size_t b=25;b<a->bins;++b) count+=a->hist[3*b+type];
        double expected=100*pairs[type]*(4*PI/3)*(27-0.125)/1000;
        double g=count/expected;printf("Uniform integrated g[%d] for 0.5<=r<3: %.9g\n",type,g);
        CHECK(fabs(g-1)<0.02);
    }
    analysis_finish(&p,&v);free(brute);free_memory(&v,&n);
    /* Small-box fallback: same exact histogram and rmax=3 in a box of 6. */
    p.L=v3(6,6,6);p.num_part=60;
    snprintf(p.analysis_prefix,sizeof(p.analysis_prefix),"../data/a_smallbox");
    alloc_memory(&p,&v,&n);initialise(&p,&v,&n,&step,&t);analysis_start(&p,&v);
    CHECK(!v.analysis->use_list);rdf_sample(&p,&v);
    unsigned long long total=0,expected=0;
    for(size_t k=0;k<3*v.analysis->bins;++k) total+=v.analysis->hist[k];
    for(size_t i=0;i<p.num_part;++i) for(size_t j=i+1;j<p.num_part;++j)
        if(v3_norm(v3_min_image(v3_sub(v.r[i],v.r[j]),p.L))<3) ++expected;
    CHECK(total==expected);analysis_finish(&p,&v);free_memory(&v,&n);
    puts("PASS RDF all-pairs comparison, force-list independence, ideal uniform normalisation, small-box fallback");
}
static void restart_test(void)
{
    struct Parameters p;struct Vectors v;struct Nbrlist n;setup(&p,&v,&n,20,5);
    snprintf(p.restart_out_filename,sizeof(p.restart_out_filename),"../data/a_test.restart");
    snprintf(p.restart_in_filename,sizeof(p.restart_in_filename),"../data/a_test.restart");
    build_nbrlist(&p,&v,&n);calculate_forces(&p,&n,&v);
    Vec3D r=v.r[7],vel=v.v[7],f=v.f[7];save_restart(&p,&v);
    v.r[7]=v.v[7]=v.f[7]=v3(0,0,0);load_restart(&p,&v);
    near(v3_norm(v3_sub(r,v.r[7])),0,0);near(v3_norm(v3_sub(vel,v.v[7])),0,0);
    near(v3_norm(v3_sub(f,v.f[7])),0,0);remove(p.restart_out_filename);
    free_memory(&v,&n);puts("PASS restart round trip");
}
int main(void)
{
    pair_tests();finite_difference_tests();noise_test();chain_tests();rdf_test();restart_test();
    puts("ALL PART A CHECKS PASSED");return 0;
}
