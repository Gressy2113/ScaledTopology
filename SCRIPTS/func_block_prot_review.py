import plumed
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib

matplotlib.rcParams['svg.fonttype'] = 'none'

from SCRIPTS.func_block import reweight_3d, read_fes_3d, dG_calc
from SCRIPTS.func_plot import plot_prof
import warnings
warnings.filterwarnings("ignore")
import scipy.ndimage
from tqdm import tqdm

BIASF = 5
kBT=310*8.314462618*0.001 # kJ/mol


def plot_2d_d_cn(axis, k, COLVAR, dist1, dist2, cn, fes, Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN, FOLDERS):
    Ubias = -fes * (1-1/BIASF)

    WEIGHTS, ED = np.histogramdd(COLVAR[['cn', 'd_580', 'd_489']].to_numpy(), 
                                bins = [Nbins_CN, Nbins_D2, Nbins_D1],
                                range = (Bonds_CN, Bonds_D2, Bonds_D1),
                                density=True, 
                                )
    weighted_avg = np.sum(np.exp(1/kBT * Ubias) * WEIGHTS, axis=axis)
    weighted_avg[weighted_avg==0]=np.nan
    norm = np.sum(np.exp(1/kBT * (Ubias-np.max(Ubias))))
    fes_dens = -kBT * np.log(weighted_avg/norm)
    #fes_dens[fes_dens==np.nan] = np.nanmax(fes_dens)
    
    fes_dens = scipy.ndimage.gaussian_filter(fes_dens, 0.1)

    vmax = np.nanmax(fes_dens)//10 * 10

    if axis == 1: 
        d = dist1[0, 0, :]
    else: 
        d = dist2[0, :, 0]
    h = plt.contourf(d, cn[:, 0, 0], fes_dens, 
                cmap = 'jet', 
                origin = 'lower',
               levels = np.arange(0, vmax, 1), 
                )

    plt.contour(d, cn[:, 0, 0], fes_dens, 
                colors = 'k', 
                linewidths=0.5,
               levels = np.arange(0, vmax, 5)
                )
    plt.colorbar(h, label="Free energy, kJ/mol")
    plt.xlabel(f'L{axis}, nm')
    plt.ylabel("CN")
    plt.title(f'{FOLDERS[k]} L{axis}')
    
    plt.xlim(0.2, 1.5)
    plt.show()


def read_fes(k, FOLDERS):
    fes_2d = pd.read_csv(f'{FOLDERS[k]}/fes_dens_2D.csv', sep = ' ', header = None)
    fes_2d = fes_2d.fillna(np.nanmax(fes_2d))

    d1 = pd.read_csv(f'{FOLDERS[k]}/d1.csv', sep = ' ', header = None).T.rename(columns={0:'d', 1:'G'})
    d2 = pd.read_csv(f'{FOLDERS[k]}/d2.csv', sep = ' ', header = None).T.rename(columns={0:'d', 1:'G'})

    with open(f'{FOLDERS[k]}/reweight/data.txt', 'w') as f: 
        for i in range(fes_2d.shape[1]): 
            for j in range(fes_2d.shape[0]): 
                #print(i, j, d1.shape, d2.shape, fes_2d.shape)
                f.write(f"{d1['d'].iloc[i]} {d2['d'].iloc[j]} {fes_2d[i][j]}\n")
                
    return (d1, d2, fes_2d)

def plot_1d(d1, d2, k):    
    fig, ax = plt.subplots(figsize = (5.5, 3))
    dG1 = dG_calc(d1['d'], d1['G'], BOND_MAX=0.4, R_RES=0.9)
    dG2 = dG_calc(d2['d'], d2['G'], BOND_MAX=0.4, R_RES=0.9)
    print(dG1, dG2)
    plt.plot(d1['d'], d1['G'], '.-', label = 'D489', color = 'blue', 
             markeredgecolor = 'k', 
             markeredgewidth = 0.5
             )
    plt.plot(d2['d'], d2['G'], '.-', label = 'D580', color = 'tab:orange', 
             markeredgecolor = 'k', 
             markeredgewidth = 0.5
             )

    plt.ylim(-25, 20)
    ax.grid(True, which='major', linestyle=  '-')
    ax.grid(True, which='minor', linestyle=  '-', lw=0.2)
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(-100, 45, 10)))
    ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(-100, 45, 10/4)))

    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(0, 3, 0.5)))
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(0, 3, 0.1)))

    plt.xlim(0, 2)#, 1.5)
    plt.title(k)
    plt.legend()
    plt.xlabel('L, nm')
    plt.ylabel('Free Energy, kJ/mol')
    plt.savefig(f'IMAGES/{k}_reweighting_1d.svg', dpi=300, bbox_inches = 'tight')
    plt.show()

def plot_2d(d1, d2, fes_2d, k, prof_name, FOLDERS, colors):    
    vmax = np.max(fes_2d)//10 * 10
    h = plt.contourf(d1['d'], d2['d'], fes_2d, 
                cmap = 'jet', 
                origin = 'lower',
                aspect = 'equal',
                # vmin=0,
                # vmax=100
                levels = np.arange(0, vmax, 1), 
                # algorithm='threaded'
                #extend='max'
                )
    # # h.cmap.set_over('red')
    # # h.cmap.set_under('blue')

    plt.contour(d1['d'], d2['d'], fes_2d, 
                colors = 'k', 
                linewidths=0.5,
                levels = np.arange(0, vmax, 5)
                )
    if prof_name: 
        path = np.loadtxt(f'{FOLDERS[k]}/{prof_name['pref']}_pts_npts{prof_name['npts']}_stepmax{prof_name['stepmax']}.csv')
        plt.plot(path[:, 0], path[:, 1],  '-', color = 'k')
    plt.colorbar(h, label = 'FES, kJ/mol')
    plt.xlim(0.2, 1.5)
    plt.ylim(0.2, 1.5)
    plt.xlabel('L1, nm')
    plt.ylabel('L2, nm')
    plt.title(k)
    plt.savefig(f'IMAGES/{k}_reweighting_2d.svg', dpi=300, bbox_inches = 'tight')
    plt.show()
    
    if prof_name: 
        prof = np.loadtxt(f'{FOLDERS[k]}/{prof_name['pref']}_prof_path_npts{prof_name['npts']}_stepmax{prof_name['stepmax']}.csv')
        plot_prof(prof, color = colors[k], ylim = (-35, 25))
        plt.savefig(f'IMAGES/{k}_reweighting_path_1d.svg', dpi=300, bbox_inches = 'tight')
        plt.show()

def plot_stride_dG(k, ax, FOLDERS, colors):

    stride = pd.read_csv(f'{FOLDERS[k]}/stride_reweight/dG_bond.csv', sep = '\t', header = None, 
                            names = ['t', 'dG1', 'dG2'])
    stride['t'] /= 40
    # if FOLDERS[i].split('/')[0].split('.')[-1] == "multiwalkers": 

    plt.plot(stride['t'], stride['dG1'], '.-', color = colors[k], label = FOLDERS[k],
             markeredgecolor = 'k', markeredgewidth = 0.5)
    plt.plot(stride['t'], stride['dG2'], '.-',  color = colors[k], markeredgecolor = 'k', 
             markeredgewidth = 0.5)


    plt.xlabel('T, μs')
    plt.ylabel('ΔG, kJ/mol')
    plt.legend()

    plt.xlim(0, 2.5)
    plt.ylim(-25, 20)

    plt.grid(True, which='major', linestyle=  '-')
    plt.grid(True, which='minor', linestyle=  '-', lw=0.2)

    plt.axhline(0, color = 'k', lw = 0.5)
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(0, 3.1, 0.5)))
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(0, 3.1, 0.1)))

    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(-30, 30, 5)))
    ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(-30, 30, 2.5)))
    
def calc_dGmeanstd(COLVAR, 
                   dist1,dist2, cn, fes, 
                   Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN, FOLDER):
    #CVLR=COLVAR #.iloc[3* 1000**2+1:8* 1000**2+1:]
    N_blocks = np.arange(3, 25, 1)

    dG1_mean = np.zeros_like(N_blocks).astype(float)
    dG1_std = np.zeros_like(N_blocks).astype(float)

    dG2_mean = np.zeros_like(N_blocks).astype(float)
    dG2_std = np.zeros_like(N_blocks).astype(float)

    for i in tqdm(range(len(N_blocks))): 
        dG1_cur, dG2_cur = [], []
        dt = len(COLVAR)//N_blocks[i]
        for j in range(0, len(COLVAR), dt):
            dG1_, dG2_ = reweight_3d(COLVAR[j:j+dt], 
                                        dist1, dist2, cn, fes, 
                                        Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN, 
                                        plot=False)

            if np.isnan(dG1_) == False and dG1_ < np.inf:
                dG1_cur.append(dG1_)
            if np.isnan(dG2_) == False and dG2_ < np.inf:
                dG2_cur.append(dG2_)

        dG1_mean[i] = np.mean(dG1_cur)
        dG1_std[i] = np.std(dG1_cur)/np.sqrt(len(dG1_cur)) 
        
        dG2_mean[i] = np.mean(dG2_cur)
        dG2_std[i] = np.std(dG2_cur)/np.sqrt(len(dG2_cur)) 

    np.savetxt(f'{FOLDER}/block_analysis.csv', 
           np.vstack([N_blocks, dG1_mean, dG1_std, dG2_mean, dG2_std]).T, 
           header = 'N_blocks dG1_mean dG1_std dG2_mean dG2_std')
    
    return (N_blocks, dG1_mean, dG1_std, dG2_mean, dG2_std)
