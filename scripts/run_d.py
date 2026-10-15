#!/usr/bin/env python3
"""Run the assignment's D1/D2/D3 cases; record exact settings and source hashes."""
import argparse,datetime,hashlib,json,os,pathlib,platform,subprocess,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
SERIES={3:[3,5,7],5:[2,3.5,5],10:[1,1.5,2]}
def tag(n,a):return f'd_n{n}_a{a:g}'.replace('.','p')
def run(case,n=5,aab=25,steps=None,equil=None,block=2500,seed=23,cc='gcc',name=None):
    steps=steps or 60000
    equil=equil if equil is not None else (10000 if case=='d1' else 20000)
    if steps<=equil or equil<0 or equil%block or steps%block or block%50:raise ValueError('Use whole blocks, steps>equil>=0; block divisible by 50')
    label=name or ('d1' if case=='d1' else tag(n,aab))
    if pathlib.Path(label).name!=label:raise ValueError('Simple filename prefix required')
    build=ROOT/'build';build.mkdir(exist_ok=True)
    sources=sorted((ROOT/'code').glob('*.c'))
    args=[cc,'-std=c11','-O3','-Wall','-Wextra','-Wpedantic','-Werror','-I',str(ROOT/'code')]
    if case=='check':
        sources=[ROOT/'code/chainstats.c',ROOT/'code/tests/test_chainstats.c']
    else:
        args += [f'-DD_CASE={1 if case=="d1" else 2}',f'-DD_N={n}',f'-DD_AAB={aab}',f'-DD_STEPS={steps}',f'-DD_EQUIL={equil}',f'-DD_BLOCK={block}',f'-DD_SEED={seed}',f'-DD_PREFIX="../data/{label}"']
    exe=build/(label+('.exe' if os.name=='nt' else ''))
    args += [str(s) for s in sources]+['-lm','-o',str(exe)]
    print('Compile',label,flush=True);subprocess.run(args,check=True)
    print('Run',label,flush=True);t=time.perf_counter()
    with (ROOT/'data'/f'{label}_run.log').open('w') as f:
        result=subprocess.run([str(exe)],cwd=ROOT/'code',stdout=f,stderr=subprocess.STDOUT)
    m=dict(case=case,tag=label,chain_length=n,aab=25 if case=='d1' else aab,dt=.04,steps=steps,equil=equil,block=block,sample=50,seed=seed,kT=1,C=2,
           L=[10,10,10] if case=='d1' else [20,10,10],num_part=3000 if case=='d1' else 6000,
           compile_command=args,run_command=[str(exe)],working_directory=str(ROOT/'code'),
           compiler=subprocess.check_output([cc,'--version'],text=True).splitlines()[0],machine=platform.platform(),
           utc_finished=datetime.datetime.now(datetime.timezone.utc).isoformat(),run_seconds=time.perf_counter()-t,returncode=result.returncode,
           source_sha256={str(s.relative_to(ROOT)):hashlib.sha256(s.read_bytes()).hexdigest() for s in sources})
    (ROOT/'data'/f'{label}_manifest.json').write_text(json.dumps(m,indent=2))
    if result.returncode:raise RuntimeError(f'{label} failed; inspect its run log')
    print(f'Finished {label}: {m["run_seconds"]:.1f}s',flush=True)
    return m
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('case',choices=['check','d1','blend','all'])
    p.add_argument('--n',type=int,choices=[3,5,10],default=5);p.add_argument('--aab',type=float,default=28)
    p.add_argument('--steps',type=int);p.add_argument('--equil',type=int);p.add_argument('--block',type=int,default=2500)
    p.add_argument('--seed',type=int,default=23);p.add_argument('--cc',default=os.environ.get('CC','gcc'));p.add_argument('--name')
    a=p.parse_args();kw=vars(a).copy();case=kw.pop('case');kw['n']=kw.pop('n')
    if case=='all':
        if a.name:p.error('--name requires an individual case')
        run('check',cc=a.cc,name='d_chain_check');run('d1',cc=a.cc)
        cases=[(n,28) for n in [3,5,10]]+[(n,25+d) for n,ds in SERIES.items() for d in ds]
        for n,ab in dict.fromkeys(cases):run('blend',**(kw|dict(n=n,aab=ab)))
    else:run(case,**kw)
if __name__=='__main__':main()
