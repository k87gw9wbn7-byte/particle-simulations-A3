#!/usr/bin/env python3
"""Compile/run C1-C4 with GCC. Run from the project root on any platform."""
import argparse,datetime,hashlib,json,pathlib,platform,subprocess,time,os
ROOT=pathlib.Path(__file__).resolve().parents[1]
def run(case,cc='gcc',aab=37,steps=120000,equil=20000,block=5000,sample=50,dt=.04,seed=13,pdb=5000,name=None):
    mapping={'c1':1,'c2':2,'c3':3,'control':4,'checks':-1,'tensor-check':-1}
    tag=name or ('c3_a'+f'{aab:g}'.replace('.','p') if case=='c3' else 'c4_control' if case=='control' else case)
    if pathlib.Path(tag).name!=tag:raise ValueError('name must be a simple output prefix')
    if not steps>equil>=0 or dt<=0 or block<1 or sample<1 or block%sample or equil%block or steps%block:
        raise ValueError('Require steps>equil>=0, dt>0, equal whole blocks and block divisible by sample')
    (ROOT/'data').mkdir(exist_ok=True);(ROOT/'build').mkdir(exist_ok=True)
    exe=ROOT/'build'/(tag+'.exe');sources=sorted((ROOT/'code').glob('*.c'))
    args=[cc,'-std=c11','-O3','-Wall','-Wextra','-Wpedantic','-Werror','-I',str(ROOT/'code')]
    if case in ['checks','tensor-check']:
        sources=[s for s in sources if s.name!='main.c']+[ROOT/'code/tests'/('test_dpd.c' if case=='checks' else 'test_tensor.c')]
    else:
        args += [f'-DC_CASE={mapping[case]}',f'-DC_AAB={aab}',f'-DC_DT={dt}',f'-DC_STEPS={steps}',f'-DC_EQUIL={equil}',f'-DC_BLOCK={block}',f'-DC_SAMPLE={sample}',f'-DC_SEED={seed}',f'-DC_PDB={pdb}',f'-DC_PREFIX="../data/{tag}"']
    args += [str(s) for s in sources]+['-lm','-o',str(exe)]
    cwd=ROOT/'code'
    if case in ['checks','tensor-check']:
        cwd=ROOT/'data/validation/code';cwd.mkdir(parents=True,exist_ok=True);(cwd.parent/'data').mkdir(exist_ok=True)
    print('Compile',tag,flush=True);t=time.perf_counter();subprocess.run(args,check=True);compile_seconds=time.perf_counter()-t
    print('Run',tag,'(progress in data/'+tag+'_run.log)',flush=True);t=time.perf_counter()
    with (ROOT/'data'/f'{tag}_run.log').open('w') as f:r=subprocess.run([str(exe)],cwd=cwd,stdout=f,stderr=subprocess.STDOUT)
    metadata={'case':case,'tag':tag,'aab':25 if case in ['c1','control'] else aab,'seed':seed,'dt':dt,'steps':steps,'equil':equil,'block':block,'sample':sample,
              'L':[10,10,10] if case=='c1' else [20,10,10],'num_part':3000 if case=='c1' else 6000,
              'compile_command':args,'run_command':[str(exe)],'working_directory':str(cwd),'compiler':subprocess.check_output([cc,'--version'],text=True).splitlines()[0],
              'machine':platform.platform(),'utc_finished':datetime.datetime.now(datetime.timezone.utc).isoformat(),'compile_seconds':compile_seconds,'run_seconds':time.perf_counter()-t,'returncode':r.returncode,
              'source_sha256':{str(s.relative_to(ROOT)):hashlib.sha256(s.read_bytes()).hexdigest() for s in sources}}
    (ROOT/'data'/f'{tag}_manifest.json').write_text(json.dumps(metadata,indent=2))
    if r.returncode:raise RuntimeError('Run failed: inspect data/'+tag+'_run.log')
    print(f'Finished {tag}, {metadata["run_seconds"]:.2f} s',flush=True)
def main():
    p=argparse.ArgumentParser();p.add_argument('case',choices=['checks','tensor-check','c1','c2','c3','control','all'])
    p.add_argument('--cc',default=os.environ.get('CC','gcc'));p.add_argument('--aab',type=float,default=37)
    p.add_argument('--steps',type=int,default=120000);p.add_argument('--equil',type=int,default=20000)
    p.add_argument('--block',type=int,default=5000);p.add_argument('--sample',type=int,default=50)
    p.add_argument('--dt',type=float,default=.04);p.add_argument('--seed',type=int,default=13);p.add_argument('--pdb',type=int,default=5000);p.add_argument('--name')
    a=p.parse_args();kw={k:v for k,v in vars(a).items() if k!='case'}
    if a.case=='all':
        if a.name:p.error('--name is only for individual runs')
        run('checks',cc=a.cc);run('tensor-check',cc=a.cc)
        for case,val in [('c1',25),('c2',37),('c3',35),('c3',39),('c3',43),('c3',47),('control',25)]:
            run(case,**(kw|{'aab':val}))
    else:run(a.case,**kw)
if __name__=='__main__':main()
