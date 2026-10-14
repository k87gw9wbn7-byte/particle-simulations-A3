#!/usr/bin/env python3
"""B1: energy figure, relative-drift table and measured timestep scaling."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data'
DT=[.04,.02,.01,.005]
def tag(dt):return f'b_energy_dt{dt:g}'.replace('.','p')
def main():
    rows=[];runs={};initial=[];hashes=[]
    for dt in DT:
        name=tag(dt);d=np.genfromtxt(DATA/f'{name}_thermo.csv',delimiter=',',names=True);runs[dt]=d
        m=json.loads((DATA/f'{name}_manifest.json').read_text());hashes.append(m['input_restart_sha256'])
        E0=d['Etot'][0];initial.append(E0);err=(d['Etot']-E0)/abs(E0)
        assert np.isfinite(err).all() and np.isclose(d['time'][-1],200), 'Incomplete/wrong-duration run'
        slope=np.polyfit(d['time'],err,1)[0]
        rows.append((dt,len(d)-1,d['time'][-1],E0,err[-1],np.max(np.abs(err)),np.sqrt(np.mean(err**2)),slope))
    assert len(set(hashes))==1 and np.allclose(initial,initial[0],rtol=0,atol=1e-9), 'Different initial states'
    a=np.array(rows)
    np.savetxt(DATA/'B1_drift_table.csv',a,delimiter=',',comments='',fmt='%.12g',
        header='dt,steps,t_final,E0,final_signed_relative_drift,max_abs_relative_error,rms_relative_error,relative_slope_per_time')
    fig,axes=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained');d=runs[.01]
    for ax,k,label in zip(axes,['Etot','Ekin','Epot'],['Total','Kinetic','Potential']):
        ax.plot(d['time'],d[k],lw=.7);ax.set(ylabel=f'{label} energy');ax.grid(alpha=.2)
    axes[-1].set_xlabel('Time (DPD units)');fig.suptitle('B1: conservative-only monomers, dt = 0.01')
    fig.savefig(DATA/'B1_energy.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for dt,d in runs.items():ax[0].plot(d['time'],(d['Etot']-d['Etot'][0])/abs(d['Etot'][0]),lw=.8,label=f'dt={dt:g}')
    ax[0].set(xlabel='Time',ylabel='[E(t) - E(0)] / |E(0)|');ax[0].legend();ax[0].grid(alpha=.2)
    ax[1].loglog(a[:,0],a[:,5],'o-',label='Max |relative error|')
    ax[1].loglog(a[:,0],a[:,6],'s-',label='RMS relative error')
    ax[1].loglog(a[:,0],a[-1,5]*(a[:,0]/a[-1,0])**2,'--',label='dt² reference')
    ax[1].set(xlabel='dt',ylabel='Error over t = 200');ax[1].legend();ax[1].grid(alpha=.2)
    fig.savefig(DATA/'B1_drift.png',dpi=180);plt.close(fig)
    scaling={'fit_range_dt':a[1:,0].tolist(),
        'max_error_exponent_smallest_three':float(np.polyfit(np.log(a[1:,0]),np.log(a[1:,5]),1)[0]),
        'rms_error_exponent_smallest_three':float(np.polyfit(np.log(a[1:,0]),np.log(a[1:,6]),1)[0]),
        'pairwise_max_error_exponents':(np.log(a[:-1,5]/a[1:,5])/np.log(a[:-1,0]/a[1:,0])).tolist()}
    (DATA/'B1_scaling.json').write_text(json.dumps(scaling,indent=2))
    print('dt | final signed drift | max absolute relative error | RMS relative error')
    for r in rows:print(f'{r[0]:.3g} | {r[4]:+.6g} | {r[5]:.6g} | {r[6]:.6g}')
    print('Observed scaling:',json.dumps(scaling))
if __name__=='__main__':main()
