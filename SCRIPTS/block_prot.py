import plumed
import sys

from func_block_prot_review import *

import warnings
warnings.filterwarnings("ignore")

'''
PARAMS
'''
FOLDER = sys.argv[1]
Tmin = int(sys.argv[2])

'''
RUN BLOCK ANALYSIS
'''
cvlr = plumed.read_as_pandas(f'{FOLDER}/COLVAR')
COLVAR = cvlr[(Tmin*1000**2 <= cvlr['time'])]

dist1, dist2, cn, fes, Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN = read_fes_3d(FOLDER)

calc_dGmeanstd(
    COLVAR, 
    dist1, dist2, cn, fes, 
    Nbins_D1, Nbins_D2, Nbins_CN, Bonds_D1, Bonds_D2, Bonds_CN,
    FOLDER
)