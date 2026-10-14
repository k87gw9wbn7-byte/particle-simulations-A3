#!/usr/bin/env python3
"""Compile each experiment separately, run from code/, save exact commands/timings."""
import argparse,datetime,hashlib,json,os,pathlib,platform,subprocess,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def run(case,cc='gcc',dt=.04,steps=20000,equil=5000,end=200.,name=None,restart=None):
    mapping={'fluid':0,'ideal':1,'energy':2}
    tag=name or (f'b_energy_dt{dt:g}'.replace('.','p') if case=='energy' else f'b_{case}')
    if pathlib.Path(tag).name!=tag:raise ValueError('name must be a filename prefix without directories')
    build=ROOT/'build';build.mkdir(exist_ok=True);(ROOT/'data').mkdir(exist_ok=True)
    exe=build/(tag+'.exe');prefix='../data/'+tag
    command=[cc,'-std=c11','-O3','-Wall','-Wextra','-Wpedantic','-Werror','-I',str(ROOT/'code')]
    sources=sorted((ROOT/'code').glob('*.c'))
    if case in ['checks','velocity-check']:
        test='test_dpd.c' if case=='checks' else 'test_velocity.c'
        sources=[p for p in sources if p.name!='main.c']+[ROOT/'code'/'tests'/test]
    else:
        command+=[f'-DB_CASE={mapping[case]}',f'-DB_DT={dt}',f'-DB_STEPS={steps}',f'-DB_EQUIL={equil}',f'-DB_END={end}',f'-DB_PREFIX="{prefix}"']
        if case=='energy':
            state=pathlib.Path(restart).resolve() if restart else ROOT/'data'/'b_fluid.restart'
            if not state.is_file():raise FileNotFoundError(f'Run fluid first; restart missing: {state}')
            command+=[f'-DB_RESTART_IN="{state.as_posix()}"']
    command+=[str(p) for p in sources]+['-lm','-o',str(exe)]
    run_cwd=ROOT/'code'
    if case in ['checks','velocity-check']:
        run_cwd=ROOT/'data'/'validation'/'code';run_cwd.mkdir(parents=True,exist_ok=True)
        (run_cwd.parent/'data').mkdir(exist_ok=True)
    record={'case':case,'name':tag,'dt':dt,'steps':round(end/dt) if case=='energy' else steps,
            'discard_steps':0 if case=='energy' else equil,'end_time':end if case=='energy' else steps*dt,
            'machine':platform.platform(),'compiler':subprocess.check_output([cc,'--version'],text=True).splitlines()[0],
            'utc_started':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'compile_command':command,'run_command':[str(exe)],'working_directory':str(run_cwd)}
    if case=='energy':record['input_restart']=str(state);record['input_restart_sha256']=digest(state)
    print('Compile',tag,flush=True);t=time.perf_counter();subprocess.run(command,check=True);record['compile_seconds']=time.perf_counter()-t
    print('Run',tag,'-> data/'+tag+'_run.log',flush=True);t=time.perf_counter()
    with (ROOT/'data'/f'{tag}_run.log').open('w') as f:
        result=subprocess.run([str(exe)],cwd=run_cwd,stdout=f,stderr=subprocess.STDOUT)
    record['wall_seconds']=time.perf_counter()-t;record['returncode']=result.returncode
    record['source_sha256']={p.name:digest(p) for p in sources}
    (ROOT/'data'/f'{tag}_manifest.json').write_text(json.dumps(record,indent=2))
    if result.returncode:raise RuntimeError(f'{tag} failed: inspect data/{tag}_run.log')
    print(f'Finished {tag}, {record["wall_seconds"]:.2f} s',flush=True)
    return record
def main():
    p=argparse.ArgumentParser();p.add_argument('case',choices=['checks','velocity-check','fluid','ideal','energy','all'])
    p.add_argument('--cc',default=os.environ.get('CC','gcc'));p.add_argument('--dt',type=float,default=.04)
    p.add_argument('--steps',type=int,default=20000);p.add_argument('--equil',type=int,default=5000)
    p.add_argument('--end',type=float,default=200);p.add_argument('--name');p.add_argument('--restart')
    a=p.parse_args()
    if a.dt<=0 or a.end<=0 or a.steps<1 or a.equil<0:p.error('positive dt/end/steps and nonnegative equil required')
    if a.case=='all':
        run('checks',a.cc,name='b_checks');run('velocity-check',a.cc,name='b_velocity_check');run('fluid',a.cc);run('ideal',a.cc)
        for dt in [.04,.02,.01,.005]:run('energy',a.cc,dt=dt,end=a.end)
    else:run(a.case,a.cc,a.dt,a.steps,a.equil,a.end,a.name,a.restart)
if __name__=='__main__':main()
