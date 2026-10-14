#!/usr/bin/env python3
"""B2: ideal-gas partial RDF, block uncertainty, and integrated checks."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
DATA=Path(__file__).resolve().parents[1]/'data'
def main():
    g=np.genfromtxt(DATA/'b_ideal_rdf.csv',delimiter=',',names=True)
    b=np.genfromtxt(DATA/'b_ideal_rdf_blocks.csv',delimiter=',',names=True)
    ids=np.unique(b['block']);assert len(ids)>1
    frames=np.array([b['frames'][b['block']==i][0] for i in ids])
    assert (frames==frames[0]).all(), 'Use equal-size blocks for this simple SE calculation'
    fig,ax=plt.subplots(1,3,figsize=(12,3.6),sharey=True,layout='constrained');rows=[]
    for a,pair in zip(ax,['AA','AB','BB']):
        blocks=np.array([b['g_'+pair][b['block']==i] for i in ids]);sem=blocks.std(axis=0,ddof=1)/np.sqrt(len(ids))
        assert np.allclose(blocks.mean(axis=0),g['g_'+pair])
        a.plot(g['r'],g['g_'+pair],lw=1,label='Measured')
        a.fill_between(g['r'],np.maximum(0,g['g_'+pair]-2*sem),g['g_'+pair]+2*sem,alpha=.2,label='±2 block SE')
        a.axhline(1,color='black',ls='--',lw=1);a.set(title=pair,xlabel='r / rc',xlim=(0,3));a.grid(alpha=.2)
        keep=g['r_lo']>=.3;w=g['r_hi'][keep]**3-g['r_lo'][keep]**3
        block_mean=np.average(blocks[:,keep],axis=1,weights=w)
        rows.append([pair,np.average(g['g_'+pair][keep],weights=w),block_mean.std(ddof=1)/np.sqrt(len(ids)),len(ids)])
    ax[0].set_ylabel('g(r)');ax[-1].legend(fontsize=8);fig.suptitle('B2: DPD ideal gas, conservative force off')
    fig.savefig(DATA/'B2_ideal_rdf.png',dpi=180);plt.close(fig)
    with (DATA/'B2_rdf_summary.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['pair','shell_weighted_mean_0p3_to_3','block_SE','blocks']);w.writerows(rows)
    print('Ideal-gas RDF summary (r >= .3; all raw bins retained in CSV):')
    for row in rows:print(row)
if __name__=='__main__':main()
