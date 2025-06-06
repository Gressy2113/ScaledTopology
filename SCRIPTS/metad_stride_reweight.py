import shutil
import sys
import plumed
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
matplotlib.rcParams['svg.fonttype'] = 'none'

from func_block import reweight_3d, read_fes_3d
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

FilePath = f"{DATA_FOLDER}/dG_bond.csv" # replace the temp with your file path/name
if os.path.isfile(FilePath):
    modifiedTime = os.path.getmtime(FilePath) 

    timeStamp =  datetime.datetime.fromtimestamp(modifiedTime).strftime("%b-%d-%y-%H:%M:%S")
    os.rename(FilePath,FilePath+"_"+timeStamp)
    shutil.copyfile(FilePath+"_"+timeStamp, FilePath)
    
COLVAR = plumed.read_as_pandas(f"{FOLDER}/COLVAR") 

for n in STRIDE:
    # COLVAR.iloc[n*1000**2 // 10: (n+1)*1000**2 // 10]
    dist1, dist2, cn, fes, Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN = read_fes_3d(FOLDER, n)
    print(Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN)
    dG_489, dG_580 = reweight_3d(COLVAR.iloc[:(n+1)*1000**2 // 10], 
                dist1, dist2, cn, fes, 
                Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN, 
                DATA_FOLDER, plot=True, FIG_FOLDER=DATA_FOLDER, FIG_NAME=n
                )
    plt.close()
    
    with open(FilePath, "a") as f:
        f.write(f"{n}\t{dG_489}\t{dG_580}\n")


