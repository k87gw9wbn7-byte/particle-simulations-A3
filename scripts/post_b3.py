#!/usr/bin/env python3
"""B3: component and speed PDFs; measured temperature from unbinned moments."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
DATA=Path(__file__).resolve().parents[1]/'data'
def gaussian(v,T):return np.exp(-v*v/(2*T))/np.sqrt(2*np.pi*T) # m=1
def maxwell(s,T):return np.sqrt(2/np.pi)*s*s*np.exp(-s*s/(2*T))/T**1.5
def main():
    fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained');rows=[]
    for row,(name,label) in enumerate([('b_fluid','Conservative on'),('b_ideal','Conservative off')]):
        h=np.genfromtxt(DATA/f'{name}_velocity_hist.csv',delimiter=',',names=True)
        s=np.genfromtxt(DATA/f'{name}_speed_hist.csv',delimiter=',',names=True)
        m=np.genfromtxt(DATA/f'{name}_velocity_moments.csv',delimiter=',',names=True)
        b=np.atleast_1d(np.genfromtxt(DATA/f'{name}_velocity_blocks.csv',delimiter=',',names=True))
        T=float(m['kT_variance']);Tc=float(m['kT_dof_corrected']);SE=b['kT_variance'].std(ddof=1)/np.sqrt(len(b))
        assert np.isfinite(T) and (b['frames']==b['frames'][0]).all()
        assert np.isclose(b['kT_variance'].mean(),T,atol=1e-10)
        rows.append([label,1,T,SE,Tc,100*(T-1),float(m['var_vx']),float(m['var_vy']),float(m['var_vz']),int(m['frames'])])
        for k in ['vx','vy','vz']:ax[row,0].plot(h['v'],h['p_'+k],lw=.8,alpha=.65,label=k)
        ax[row,0].plot(h['v'],gaussian(h['v'],1),'k--',label='Imposed kT=1')
        ax[row,0].plot(h['v'],gaussian(h['v'],T),color='#bb4b28',lw=1.5,label=f'Gaussian, measured kT={T:.4f}')
        ax[row,0].set(title=label,xlabel='Velocity component',ylabel='Probability density',xlim=(-4,4));ax[row,0].legend(fontsize=8)
        ax[row,1].plot(s['s'],s['p_speed'],label='Measured speed PDF')
        ax[row,1].plot(s['s'],maxwell(s['s'],1),'k--',label='Maxwell, imposed kT=1')
        ax[row,1].plot(s['s'],maxwell(s['s'],T),color='#bb4b28',label='Maxwell, measured kT')
        ax[row,1].set(title=label,xlabel='Speed |v|',ylabel='Probability density',xlim=(0,5));ax[row,1].legend(fontsize=8)
        # Full normalisation includes explicitly recorded out-of-range samples.
        tails=np.genfromtxt(DATA/f'{name}_velocity_tails.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
        n=float(h['n_per_component'][0])
        for i,k in enumerate(['vx','vy','vz']):
            assert h['count_'+k].sum()+int(tails['underflow'][i])+int(tails['overflow'][i])==n
        assert s['count'].sum()+int(tails['overflow'][-1])==n
    fig.suptitle('B3: Maxwell–Boltzmann statistics; dt=0.04, gamma=4.5, sigma=3')
    fig.savefig(DATA/'B3_velocity.png',dpi=180);plt.close(fig)
    with (DATA/'B3_temperature_table.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['case','imposed_kT','kT_from_variance','block_SE','kT_COM_corrected','deviation_percent','var_vx','var_vy','var_vz','frames']);w.writerows(rows)
    for row in rows:print(row)
if __name__=='__main__':main()
