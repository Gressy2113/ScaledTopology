import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import matplotlib
import plumed
import os

matplotlib.rcParams['svg.fonttype'] = 'none'

kBT=310*8.314462618*0.001 # kJ/mol BINS=[50, 50, 25]
BULK_MIN=2
BULK_MAX=2.5

K_RES = 100000 # kJ/mol/nm2
Vmol = 1.66 #nm3
BIASF=5

def dG_calc(d, fesd, BOND_MAX, R_RES):
    is_bulk=np.int_(np.logical_and(d>BULK_MIN,d<BULK_MAX))    
    is_bond=np.int_(d<BOND_MAX)

    hstep = d[2]-d[1]
    dW = kBT*np.log(np.nansum(np.exp(-fesd/kBT)*is_bulk)*hstep/(np.nansum(is_bulk)*hstep))
    L_bond = np.nansum(np.exp(-fesd/kBT)*is_bond*hstep)#*BOND_MAX
    L_bulk = BULK_MAX-BULK_MIN
    S_cyl = np.pi*R_RES**2 + 2*np.pi*R_RES*np.sqrt(np.pi*kBT/(2*K_RES)) + 2*np.pi*kBT/K_RES

    dG0_cyl = dW - kBT*np.log(S_cyl*L_bond/Vmol)
    dW_bond = -kBT*np.log(L_bond/L_bulk)
    dGv = - kBT*np.log(S_cyl*L_bulk/Vmol)
    return dG0_cyl

def reweight_2d(tstart, tfinal, COLVAR, dist, cn, fes, Nbins_D, Nbins_CN, Bonds_D, Bonds_CN, FOLDER, plot=False): 
    Ubias = -fes.T * (1-1/BIASF) #4/5 #-(1-1/5) * fes_2d.to_numpy()
    
    WEIGHTS, ED = np.histogramdd(COLVAR[['dp', 'cn']].to_numpy()[tstart:tfinal], 
                                bins = [Nbins_D, Nbins_CN],
                                range = (Bonds_D, Bonds_CN),
                                density=True, 
                                )
    if plot: 
        fig, ax = plt.subplots()
        sns.heatmap(WEIGHTS.T, cmap = cm.jet)
        ax.invert_yaxis()
        plt.show()
    
    ####1D####
    if plot:
        fig, ax = plt.subplots(figsize = (5.5, 3))

    weighted_avg = np.sum(np.exp(1/kBT * (Ubias-np.max(Ubias))) * WEIGHTS, axis=1)
    weighted_avg[weighted_avg==0]=np.nan
    norm = np.sum(np.exp(1/kBT * (Ubias-np.max(Ubias))))

    fes_dens = -kBT * np.log(weighted_avg/norm)

    is_bulk=np.int_((BULK_MIN < dist[0]) & (dist[0] < BULK_MAX))
    shift = np.nansum(is_bulk*fes_dens)/np.nansum(is_bulk)
    fes_dens-=shift
    dG_489 = dG_calc(dist[0], fes_dens, BOND_MAX=0.4, R_RES=0.7)

    if plot:
        plt.plot(dist[0], fes_dens, '.-', color = 'blue', label = f'dG={round(dG_489, 3)} kJ/mol')

    ###
    if plot:
        #plt.ylim(-25, 25)
        ax.grid(True, which='major', linestyle=  '-')
        ax.grid(True, which='minor', linestyle=  '-', lw=0.2)

        plt.xlim(0, 2)
        plt.legend()
        plt.xlabel('L, nm')
        plt.ylabel('Free Energy, kJ/mol')
        plt.savefig(f'{FOLDER}/reweighting_1d.svg', dpi=300, bbox_inches = 'tight')
        plt.show()
    
    if plot: 
        np.savetxt(f'{FOLDER}/prof_1D_reweight.dat', np.concatenate(([dist[0]], [fes_dens])).T)

    return (dG_489)

def reweight_3d(COLVAR,
                dist1, dist2, cn, fes, 
                Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN, 
                FOLDER = None, plot=False, FIG_FOLDER = 'FIG_REWEIGHT',FIG_NAME = 'reweight_2d'): 
    
    if not os.path.exists(FIG_FOLDER):
        os.mkdir(FIG_FOLDER)

    Ubias = -fes * (1-1/BIASF)

    WEIGHTS, ED = np.histogramdd(COLVAR[['cn', 'd_580', 'd_489']].to_numpy(), #[tstart:tfinal], 
                                bins = [Nbins_CN, Nbins_D2, Nbins_D1],
                                range = (Bonds_CN, Bonds_D2, Bonds_D1),
                                density=True, 
                                )
    weighted_avg = np.sum(np.exp(1/kBT * Ubias) * WEIGHTS, axis=0)
    weighted_avg[weighted_avg==0]=np.nan
    norm = np.sum(np.exp(1/kBT * (Ubias-np.max(Ubias))))
    fes_dens = -kBT * np.log(weighted_avg/norm)
    if plot: 
        fig, ax = plt.subplots()
        sns.heatmap(fes_dens, cmap = cm.jet)
        ax.invert_yaxis()
        plt.xlim(0, 100)
        plt.ylim(0, 100)
        plt.title(FOLDER)
        plt.savefig(f'{FIG_FOLDER}/{FIG_NAME}_2d.png', dpi=300, bbox_inches = 'tight')
        plt.show()
        np.savetxt(f'{FOLDER}/fes_dens_2D.csv', fes_dens)

    ####1D####
    if plot:
        fig, ax = plt.subplots(figsize = (5.5, 3))


    # Ubias_d = Ubias #-fes_dens * (1-1/BIASF)
    # WEIGHTS_d = np.nansum(WEIGHTS, axis=0)

    ###d1###
    dist = dist1[0, 0, :]
    dist_dens = np.sum(np.exp(1/kBT * Ubias) * WEIGHTS, axis=(0, 1))
    dist_fes = -kBT * np.log(dist_dens/norm)
    # dist_dens = np.nansum(np.exp(-fes_dens/kBT),axis=0) #/np.nansum(np.exp(-fes_dens/kBT))
    # dist_fes = -kBT*np.log(dist_dens)
    
    is_bulk=np.int_((BULK_MIN < dist) & (dist < BULK_MAX))
    shift = np.nansum(is_bulk*dist_fes)/np.nansum(is_bulk) #np.mean(dist_fes[is_bulk==1]) #/np.sum(is_bulk)
    dist_fes -= shift
    dG_489 = dG_calc(dist, dist_fes, BOND_MAX=0.4, R_RES=0.9)
    if plot:
        plt.plot(dist, dist_fes, '.-', label = f'D489\ndG={round(dG_489, 3)}', color = 'blue')
        
        np.savetxt(f'{FOLDER}/d1.csv', np.vstack((dist, dist_fes)))
    
    ###d2###
    dist = dist2[0, :, 0]
    # dist_dens = np.nansum(np.exp(1/kBT * Ubias_d) * WEIGHTS_d, axis=1)
    # norm = np.nansum(np.exp(1/kBT * (Ubias_d)))
    dist_dens = np.sum(np.exp(1/kBT * Ubias) * WEIGHTS, axis=(0, 2))

    dist_fes = -kBT * np.log(dist_dens/norm)
    # dist_dens = np.nansum(np.exp(-fes_dens/kBT),axis=1)/np.nansum(np.exp(-fes_dens/kBT))
    # dist_fes = -kBT*np.log(dist_dens)
    
    is_bulk=np.int_((BULK_MIN < dist) & (dist < BULK_MAX))
    shift = np.nansum(is_bulk*dist_fes)/np.nansum(is_bulk) #np.mean(dist_fes[is_bulk==1]) #/np.sum(is_bulk)
    dist_fes -= shift
    dG_580 = dG_calc(dist, dist_fes, BOND_MAX=0.4, R_RES=0.9)
    if plot:
        plt.plot(dist, dist_fes, '.-', label = f'D580\ndG={round(dG_580, 3)}', color = 'tab:orange')
        np.savetxt(f'{FOLDER}/d2.csv', np.vstack((dist, dist_fes)))

    ###
    if plot:
        #plt.ylim(-20, 20)
        ax.grid(True, which='major', linestyle=  '-')
        ax.grid(True, which='minor', linestyle=  '-', lw=0.2)

        plt.xlim(0, 2)
        plt.legend()
        plt.xlabel('L, nm')
        plt.ylabel('Free Energy, kJ/mol')
        plt.title(FOLDER)
        plt.savefig(f'{FIG_FOLDER}/{FIG_NAME}_1d.png', dpi=300, bbox_inches = 'tight')
        plt.show()
    return (dG_489, dG_580)

def read_fes_2d(FOLDER):    
    data = plumed.read_as_pandas(f"{FOLDER}/FES/fes.dat")

    CV = ['dp', 'cn']
    
    Nbins_DP, Nbins_CN=0,0
    Bonds_DP, Bonds_CN=[-1, -1], [-1, -1]

    with open(f'{FOLDER}/FES/fes.dat', 'r') as f: 
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
    print(Nbins_DP, Nbins_CN, Bonds_DP, Bonds_CN)

    dist = np.array(data[CV[0]]).reshape((Nbins_CN, Nbins_DP))
    cn = np.array(data[CV[1]]).reshape((Nbins_CN, Nbins_DP))

    fes    = np.array(data["file.free"]).reshape((Nbins_CN, Nbins_DP)) #* m

    return (dist, cn, fes, Nbins_DP, Nbins_CN, Bonds_DP, Bonds_CN)

def read_fes_3d(FOLDER, NAME='fes'):    
    data = plumed.read_as_pandas(f"{FOLDER}/FES/{NAME}.dat")

    CV = ['d_489', 'd_580', 'cn']
    
    Nbins_D1, Nbins_D2, Nbins_CN=0,0,0
    Bonds_D1, Bonds_D2, Bonds_CN=[-1, -1], [-1, -1], [-1, -1]
    with open(f'{FOLDER}/FES/{NAME}.dat', 'r') as f: 
        for s in f.readlines():
        #for s in [ss[_] for _ in np.where([s[:2] == '#!' for s in ss])[0]]: 
            if s[0] != '#': 
                break
            s1 = s.split()
            if s1[2] == f'nbins_{CV[0]}': 
                Nbins_D1 = int(s1[3])
            elif s1[2] == f'nbins_{CV[1]}': 
                Nbins_D2 = int(s1[3])
            elif s1[2] == f'nbins_{CV[2]}': 
                Nbins_CN = int(s1[3])
                
            elif s1[2] == f'min_{CV[0]}': 
                Bonds_D1[0] = float(s1[3])
            elif s1[2] == f'max_{CV[0]}':
                Bonds_D1[1] = float(s1[3])
            elif s1[2] == f'min_{CV[1]}': 
                Bonds_D2[0] = float(s1[3])
            elif s1[2] == f'max_{CV[1]}':
                Bonds_D2[1] = float(s1[3])
            elif s1[2] == f'min_{CV[2]}': 
                Bonds_CN[0] = float(s1[3])
            elif s1[2] == f'max_{CV[2]}':
                Bonds_CN[1] = float(s1[3])
    print(Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN)

    dist1 = np.array(data[CV[0]]).reshape((Nbins_CN, Nbins_D2, Nbins_D1))
    dist2 = np.array(data[CV[1]]).reshape((Nbins_CN, Nbins_D2, Nbins_D1))
    cn = np.array(data[CV[2]]).reshape((Nbins_CN, Nbins_D2, Nbins_D1))
    fes    = np.array(data["file.free"]).reshape((Nbins_CN, Nbins_D2, Nbins_D1)) 

    return (dist1, dist2, cn, fes,  Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN) 