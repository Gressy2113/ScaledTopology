import shutil
import sys
import pandas as pd
import plumed
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
matplotlib.rcParams['svg.fonttype'] = 'none'

import datetime

import warnings
warnings.filterwarnings("ignore")


FOLDER = sys.argv[1]
DATA_FOLDER=f'{FOLDER}/stride_reweight'
if not os.path.exists(DATA_FOLDER):
    os.mkdir(DATA_FOLDER)

BIASF = float(sys.argv[2])
Nmin, Nmax = int(sys.argv[3]), int(sys.argv[4])
STRIDE = np.arange(Nmin, Nmax)
T0=1000001

CV = ['dp', 'cn']
kBT=310*8.314462618*0.001 # kJ/mol BINS=[50, 50, 25]
BULK_MIN=2
BULK_MAX=2.5
BOND_MAX=0.4
K_RES = 100000 # kJ/mol/nm2
R_RES = 0.7 #nm
Vmol = 1.66 #nm3

def dG_calc(d, fesd, BOND_MAX=0.4):
    is_bulk=np.int_(np.logical_and(d>BULK_MIN,d<BULK_MAX))    
    is_bond=np.int_(d<BOND_MAX)

    hstep = d[2]-d[1]
    dW = kBT*np.log(np.sum(np.exp(-fesd/kBT)*is_bulk)*hstep/(np.sum(is_bulk)*hstep))
    L_bond = np.sum(np.exp(-fesd/kBT)*is_bond*hstep)#*BOND_MAX
    #print (L_bond)
    L_bulk = BULK_MAX-BULK_MIN
    #S_cyl = 2*np.pi*kBT/K_RES + np.pi*R_RES**2
    S_cyl = np.pi*R_RES**2 + 2*np.pi*R_RES*np.sqrt(np.pi*kBT/(2*K_RES)) + 2*np.pi*kBT/K_RES

    #print (2*np.pi*kBT/K_RES,np.pi*R_RES**2)
    dG0_cyl = dW - kBT*np.log(S_cyl*L_bond/Vmol)
    dW_bond = -kBT*np.log(L_bond/L_bulk)
    dGv = - kBT*np.log(S_cyl*L_bulk/Vmol)
    print ('dW_bulk= ' + str(dW) + ' kJ/mol')
    print ('dW_bond= ' + str(dW_bond) + ' kJ/mol')
    print ('dGv= ' + str
           (dGv) + ' kJ/mol')
    print ('dG0= ' + str(dG0_cyl) + ' kJ/mol')
    return (dW, dW_bond, dGv, dG0_cyl)




FilePath = f"{DATA_FOLDER}/dG_bond.csv" # replace the temp with your file path/name
if os.path.isfile(FilePath):
    modifiedTime = os.path.getmtime(FilePath) 

    timeStamp =  datetime.datetime.fromtimestamp(modifiedTime).strftime("%b-%d-%y-%H:%M:%S")
    os.rename(FilePath,FilePath+"_"+timeStamp)
    shutil.copyfile(FilePath+"_"+timeStamp, FilePath)
    
COLVAR = plumed.read_as_pandas(f"{FOLDER}/COLVAR") 

for n in STRIDE:
    fes = plumed.read_as_pandas(f"{FOLDER}/FES/{n}.dat")
    fes_2d = fes.pivot(index='dp', columns='cn', values='file.free')
    Nbins_DP, Nbins_CN=0,0
    Bonds_DP, Bonds_CN=[-1, -1], [-1, -1]
    with open(f'{FOLDER}/FES/{n}.dat', 'r') as f: 
        ss = f.readlines()
        for s in [ss[_] for _ in np.where([s[:2] == '#!' for s in ss])[0]]: 
            s1 = s.split()
            if s1[2] == f'nbins_{CV[0]}': 
                Nbins_DP = int(s1[3])
            elif s1[2] == f'nbins_{CV[1]}': 
                Nbins_CN = int(s1[3])
            elif s1[2] == f'min_{CV[0]}': 
                Bonds_DP[0] = float(s1[3])
            elif s1[2] == f'max_{CV[0]}':
                Bonds_DP[1] = float(s1[3])
            elif s1[2] == f'min_{CV[1]}': 
                Bonds_CN[0] = float(s1[3])
            elif s1[2] == f'max_{CV[1]}':
                Bonds_CN[1] = float(s1[3])
    print(Nbins_DP, Nbins_CN)
    print(Bonds_DP, Bonds_CN)
    
    Ubias = -fes_2d.to_numpy()*(1-1/BIASF) 
    DP = fes_2d.index.to_numpy()
    WEIGHTS, ED = np.histogramdd(COLVAR.iloc[T0:(n+1)*1000**2 // 10][['dp', 'cn']].to_numpy(), 
                                bins = [Nbins_DP, Nbins_CN],
                                range = (Bonds_DP, Bonds_CN),
                                density=True, 
                                )

    weighted_avg = np.sum(np.exp(1/kBT * Ubias) * WEIGHTS, axis=1)
    weighted_avg[weighted_avg==0]=np.nan
    norm = np.sum(np.exp(1/kBT * (Ubias-np.max(Ubias))))

    #print(weighted_avg, norm)
    fes_dens = -kBT * np.log(weighted_avg/norm)

    is_bulk=np.int_((BULK_MIN < DP) & (DP < BULK_MAX))
    shift = np.nansum(is_bulk*fes_dens)/np.nansum(is_bulk)
    fes_dens-=shift

    fig, ax = plt.subplots()


    plt.plot(DP, fes_dens, 'o-', label = 'reweight', markersize = 4, markeredgecolor = 'k')
    np.savetxt(f'{DATA_FOLDER}/prof_1D_{n}_reweight.dat', np.concatenate(([DP], [fes_dens])).T)

    # plt.plot(df[0], df[1], label = 'rew_my')

    plt.xlabel('D, nm')
    plt.ylabel('FES, kJ/mol')

    plt.ylim(-25, 20)
    ax.grid(True, which='major', linestyle=  '-')
    ax.grid(True, which='minor', linestyle=  '-', lw=0.2)
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(-100, 45, 10)))
    ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(-100, 45, 10/4)))

    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(0, 3, 0.5)))
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(0, 3, 0.1)))

    plt.xlim(0, 1.5)
    plt.title(FOLDER)

    dist_dens = np.sum(np.exp(-fes_2d/kBT),axis=1)/np.sum(np.exp(-fes_2d/kBT)).sum()
    DP = fes_2d.index.to_numpy()

    fes_1d_sum = -kBT*np.log(dist_dens)
    is_bulk=np.int_((BULK_MIN < DP) & (DP < BULK_MAX))
    shift = np.sum(is_bulk*fes_1d_sum)/np.sum(is_bulk) #np.mean(dist_fes[is_bulk==1]) #/np.sum(is_bulk)
    fes_1d_sum -= shift

    plt.plot(DP, fes_1d_sum, label = 'sum_hills')

    plt.legend()
    plt.savefig(f'{DATA_FOLDER}/{n}.png', bbox_inches = 'tight')
    plt.close()

    prof = pd.read_csv(f'{DATA_FOLDER}/prof_1D_{n}_reweight.dat', sep = ' ', header = None)
    
    dG = dG_calc(prof[0], prof[1], BOND_MAX=BOND_MAX)
    with open(FilePath, "a") as f:
        f.write(f"{n}\t{dG[3]}\n")


