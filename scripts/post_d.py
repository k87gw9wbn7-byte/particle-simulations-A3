#!/usr/bin/env python3
"""Reproduce Part D tables, figures and convergence checks from saved data."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data'
SERIES={3:[3,5,7],5:[2,3.5,5],10:[1,1.5,2]}
def tag(n,a):return f'd_n{n}_a{a:g}'.replace('.','p')
def load(name):return np.atleast_1d(np.genfromtxt(DATA/name,delimiter=',',names=True))
def manifest(t):return json.loads((DATA/f'{t}_manifest.json').read_text())
def reblock(v):
    v=np.asarray(v,float);rows=[]
    for q in [1,2,4]:
        nb=len(v)//q
        if nb>=4:
            z=v[:nb*q].reshape(nb,q).mean(axis=1)
            rows.append((q,nb,float(z.std(ddof=1)/np.sqrt(nb))))
    if not rows:raise ValueError('Need at least four complete production blocks')
    return max(row[2] for row in rows),rows

def savefig(fig,name):
    fig.tight_layout();fig.savefig(DATA/name,dpi=180)
    try:
        from IPython import get_ipython
        from IPython.display import display
        if get_ipython() is not None:display(fig)
    except ImportError:pass
    plt.close(fig)
def write_table(name,rows):
    with (DATA/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def print_table(rows,keys):
    print(' | '.join(keys))
    for row in rows:print(' | '.join(f'{row[k]:.6g}' if isinstance(row[k],(float,np.floating)) else str(row[k]) for k in keys))

def snapshot(t,name,chains=False):
    # The C writer stores the initial and final configurations; use the final one.
    xyz=[];types=[];time=None
    for line in (DATA/f'{t}.pdb').read_text().splitlines():
        if line.startswith('MODEL'):xyz=[];types=[]
        elif line.startswith('REMARK TIME'):time=float(line.split('=')[1])
        elif line.startswith('HETATM'):
            xyz.append([float(line[30:38]),float(line[38:46]),float(line[46:54])]);types.append(line[12:16].strip())
    m=manifest(t);r=np.array(xyz);types=np.array(types);L=np.array(m['L'])
    if len(r)!=m['num_part'] or not np.isclose(time,m['steps']*m['dt']):raise ValueError('Missing final snapshot')
    fig=plt.figure(figsize=(9,5));ax=fig.add_subplot(111,projection='3d')
    for lab,color in [('C','#2577bb'),('O','#df6c31')]:
        z=r[types==lab];ax.scatter(*z.T,s=5 if chains else 8,c=color,alpha=.55 if chains else .8,label='A' if lab=='C' else 'B',depthshade=True,rasterized=True)
    if chains:
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        segments=[];n=m['chain_length']
        for c in range(len(r)//n):
            for j in range(n-1):
                p,q=r[c*n+j:c*n+j+2]
                if np.all(np.abs(q-p)<L/2):segments.append([p,q])
        ax.add_collection3d(Line3DCollection(segments,colors='black',linewidths=.45,alpha=.55))
    ax.set(xlim=(0,L[0]),ylim=(0,L[1]),zlim=(0,L[2]),xlabel='x / rc',ylabel='y / rc',zlabel='z / rc')
    ax.set_box_aspect(L);ax.view_init(elev=21,azim=-65);ax.legend(loc='upper right')
    ax.set_title(f'N={m["chain_length"]}, aAB={m["aab"]:g}, t={time:g}; final simulated coordinates')
    savefig(fig,name)

def d1():
    m=manifest('d1');b=load('d1_chain_blocks.csv');v=b['Re2'][b['production']==1]
    e,rb=reblock(v);ideal=3*(m['chain_length']-1)*m['kT']/m['C']
    half=len(v)//2;first,last=v[:half].mean(),v[half:].mean()
    row=dict(N=m['chain_length'],steps=m['steps'],discard_steps=m['equil'],blocks=len(v),Re2=float(v.mean()),Re2_SE=e,ideal_Re2=ideal,
             mean_bond2=float(b['bond2'][b['production']==1].mean()),first_half=float(first),second_half=float(last))
    write_table('D1_chain_size.csv',[row]);print_table([row],list(row))
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    axes[0].plot(b['t_end'],b['Re2'],'o-');axes[0].axhline(ideal,color='black',ls='--',label='Ideal springs, imposed kBT')
    axes[0].axvline(m['equil']*m['dt'],color='grey',ls=':');axes[0].set(xlabel='Time',ylabel='Block mean end-to-end distance squared');axes[0].legend()
    axes[1].plot([q*m['block']*m['dt'] for q,nb,se in rb],[se for q,nb,se in rb],'o-');axes[1].set(xlabel='Merged block duration',ylabel='SE of mean Re²')
    savefig(fig,'D1_chain_convergence.png');snapshot('d1','D1_snapshot.png',True)
    (DATA/'D1_summary.json').write_text(json.dumps(row,indent=2));return row

def profile(t,cut=None):
    m=manifest(t);raw=load(f'{t}_profile_blocks.csv');groups=[raw[raw['block']==b] for b in np.unique(raw['block'])]
    x=groups[0]['x'];a=np.array([g['rho_A'] for g in groups]);b=np.array([g['rho_B'] for g in groups]);end=np.array([g['t_end'][0] for g in groups]);start=np.array([g['t_start'][0] for g in groups])
    assert all(len(g)==len(x) and np.allclose(g['x'],x) for g in groups)
    assert len({g['frames'][0] for g in groups})==1
    L=m['L'][0];cut=m['equil']*m['dt'] if cut is None else cut;keep=start>=cut-1e-9
    scale=(x[1]-x[0])*m['L'][1]*m['L'][2]
    assert np.allclose(a.sum(axis=1)*scale,m['num_part']/2) and np.allclose(b.sum(axis=1)*scale,m['num_part']/2)
    def fraction(lo,hi,component):
        w=(x>=lo*L)&(x<hi*L)
        return (a if component=='A' else b)[:,w].sum(axis=1)/(a+b)[:,w].sum(axis=1)
    minority=fraction(.7,.8,'A');opposite=fraction(.2,.3,'B');wide=fraction(.65,.85,'A')
    return dict(m=m,x=x,a=a,b=b,phi=a/(a+b),end=end,keep=keep,cut=cut,v=minority,opposite=opposite,wide=wide,
                aligned=np.array([g['aligned_fraction'][0] for g in groups]))

def diagnostics(d):
    k=d['keep'];v=d['v'][k];se,rb=reblock(v);half=len(v)//2
    s1,_=reblock(v[:half]);s2,_=reblock(v[half:]);den=np.hypot(s1,s2)
    drift=float(abs(v[:half].mean()-v[half:].mean()));drift_z=drift/den if den>0 else float('inf')
    delta=v-d['wide'][k];delta_se,_=reblock(delta);sym=v-d['opposite'][k];symse,_=reblock(sym)
    x=d['x'];a=d['a'][k];b=d['b'][k];L=d['m']['L'][0]
    def subwindow(lo,hi):
        w=(x>=lo*L)&(x<hi*L)
        return a[:,w].sum(axis=1)/(a+b)[:,w].sum(axis=1)
    central=subwindow(.7,.75)-subwindow(.75,.8);central_se,_=reblock(central)
    flags=[]
    if drift_z>3:flags.append('production-half drift')
    if d['aligned'][k].min()<.95:flags.append('weak/absent x slab')
    if abs(delta.mean())>max(3*delta_se,.1*v.mean()):flags.append('bulk-window sensitivity')
    if abs(central.mean())>max(3*central_se,.1*v.mean()):flags.append('nonflat central bulk')
    if abs(sym.mean())>max(3*symse,.1*v.mean()):flags.append('A/B asymmetry')
    if not 0<v.mean()<.5:flags.append('no valid minority plateau')
    return dict(phi=float(v.mean()),phi_SE=se,drift=drift,drift_z=drift_z,window_difference=float(delta.mean()),window_difference_SE=delta_se,
                opposite_phi=float(d['opposite'][k].mean()),symmetry_difference_SE=symse,central_difference=float(central.mean()),central_difference_SE=central_se,zero_minority_blocks=int(np.count_nonzero(v<=0)),flags='; '.join(flags) or 'none'),rb

def d2():
    fig,axes=plt.subplots(1,2,figsize=(11,4));rows=[]
    for n in [3,5,10]:
        d=profile(tag(n,28));k=d['keep'];A=d['a'][k].mean(axis=0);B=d['b'][k].mean(axis=0);z,rb=diagnostics(d)
        axes[0].plot(d['x'],A/(A+B),label=f'N={n}')
        axes[1].plot(d['end'],d['v'],'o-',ms=3,label=f'N={n}')
        rows.append(dict(N=n,**z));snapshot(tag(n,28),f'D2_snapshot_N{n}.png')
    axes[0].set(xlabel='x / rc',ylabel='A bead fraction',ylim=(0,1));axes[0].legend()
    axes[1].set(xlabel='Time',ylabel='A fraction in B-rich central window');axes[1].legend();axes[1].axvline(d['cut'],color='grey',ls=':')
    savefig(fig,'D2_profiles.png');write_table('D2_diagnostics.csv',rows);print_table(rows,['N','phi','phi_SE','drift_z','flags']);return rows

def chiN(p):return np.log((1-p)/p)/(1-2*p)
def d3():
    rows=[];rbrows=[];fig,axes=plt.subplots(3,3,figsize=(12,10),sharex=True,sharey=True)
    fig2,axes2=plt.subplots(3,1,figsize=(9,9))
    for row,(n,deltas) in enumerate(SERIES.items()):
        for col,delta in enumerate(deltas):
            t=tag(n,25+delta);d=profile(t);k=d['keep'];v=d['v'][k];z,rbs=diagnostics(d)
            if not 0<z['phi']<.5:raise ValueError(f'{t}: cannot extract chi from absent coexistence')
            est=float(chiN(z['phi']));errors=[]
            for q,nb,se in rbs:
                values=v[:q*nb].reshape(nb,q).mean(axis=1);loo=(values.sum()-values)/(nb-1)
                if np.any(loo<=0):raise ValueError(f'{t}: insufficient minority counts')
                estimates=chiN(loo);err=float(np.sqrt((nb-1)/nb*np.sum((estimates-estimates.mean())**2)));errors.append(err)
                rbrows.append(dict(N=n,delta_a=delta,merge=q,blocks=nb,block_time=q*d['m']['block']*.04,phi_SE=se,chiN_SE=err))
            err=max(errors)
            rows.append(dict(N=n,delta_a=delta,aAB=25+delta,**z,chi=est/n,chi_SE=err/n,chiN=est,chiN_SE=err,ratio=est/delta,ratio_SE=err/delta,
                             paper_ratio=.306*n,discard_time=d['cut'],steps=d['m']['steps'],blocks=int(k.sum())))
            ax=axes[row,col];meanA=d['a'][k].mean(axis=0);meanB=d['b'][k].mean(axis=0)
            ax.plot(d['x'],meanA/(meanA+meanB));ax.axvspan(14,16,color='orange',alpha=.2);ax.axvspan(4,6,color='blue',alpha=.12)
            ax.set(title=f'N={n}, Δa={delta:g}',xlabel='x / rc',ylabel='A bead fraction',ylim=(0,1))
            axes2[row].plot(d['end'],d['v'],'o-',ms=3,label=f'Δa={delta:g}')
        axes2[row].axvline(d['cut'],color='grey',ls=':');axes2[row].set(xlabel='Time',ylabel='Minority bead fraction',title=f'N={n}');axes2[row].legend()
    savefig(fig,'D3_profiles.png');savefig(fig2,'D3_bulk_convergence.png')
    write_table('D3_table.csv',rows);write_table('D3_reblocking.csv',rbrows)
    fig,ax=plt.subplots(figsize=(8,5));xx=np.linspace(2,11,100)
    ax.plot(xx,.306*xx,'k--',label='Groot & Warren Eq. (27): 0.306 N')
    ax.fill_between(xx,.303*xx,.309*xx,color='grey',alpha=.2,label='Published coefficient ±0.003')
    for n in SERIES:
        group=[r for r in rows if r['N']==n]
        offsets=np.linspace(-.18,.18,len(group))
        for off,r in zip(offsets,group):
            ax.errorbar(n+off,r['ratio'],yerr=r['ratio_SE'],fmt='o',capsize=3,label=f'N={n}, Δa={r["delta_a"]:g}'+(' (check)' if r['flags']!='none' else ''))
    ax.set(xlabel='Chain length N (points offset slightly for visibility)',ylabel='χN / Δa',xticks=[3,5,10]);ax.legend(fontsize=8,ncol=2)
    savefig(fig,'D3_chi_comparison.png');print_table(rows,['N','delta_a','phi','phi_SE','chi','chi_SE','chiN','chiN_SE','drift_z','flags'])
    print('Flags require inspection or more sampling; agreement with a line does not establish equilibration.')
    return rows
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('part',choices=['d1','d2','d3','all']);a=p.parse_args()
    for name in ['d1','d2','d3'] if a.part=='all' else [a.part]:globals()[name]()
