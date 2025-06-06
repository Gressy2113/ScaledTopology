#!/bin/bash

N_WALKERS=4
GPU_ID=1
NT=22
PINOFFSET=24
NT1=$((NT/N_WALKERS))
echo $NT1

for ((N = 0 ; N < $N_WALKERS ; N++ ))
do

rm -r walker_$N
mkdir -p walker_$N
cd walker_$N 
cp ../system.top . 
cp ../eq.$N.gro eq.gro
cp ../posres.gro .
cp ../init_eq.pdb .
cp ../md.mdp .

gmx_mpi grompp -f md -c eq.gro -r posres.gro -p system.top -n ../indata/grps.ndx -o md

echo "MOLINFO STRUCTURE=init_eq.pdb

# group with water oxygens
OW:  GROUP ATOMS=@water REMOVE=@hydrogens

# group with Mg ion
MG:  GROUP ATOMS=17511

# site: CG D489+580 + C=O 484+576
CG_489:  GROUP ATOMS=215
CG_580:  GROUP ATOMS=638
O_484:   GROUP ATOMS=144
O_576:   GROUP ATOMS=589

# distance between Mg and com
d_489:    DISTANCE ATOMS=MG,CG_489
d_580:    DISTANCE ATOMS=MG,CG_580
d_484:    DISTANCE ATOMS=MG,O_484
d_576:    DISTANCE ATOMS=MG,O_576


# wall to CD E588
CD_588: GROUP ATOMS=775
d_588:    DISTANCE ATOMS=MG,CD_588

cn:    COORDINATION GROUPA=MG GROUPB=OW R_0=0.275 NLIST NL_CUTOFF=0.6 NL_STRIDE=20


# this is the metadynamics
metadP: METAD ...
  ARG=d_489,d_580,cn
  SIGMA=0.05,0.05,0.1
  HEIGHT=0.3
  PACE=500
  TEMP=310
  BIASFACTOR=5
  FILE=../HILLS
# very permissive boundaries to avoid unexpected stops:
  GRID_MIN=0,0,0 GRID_MAX=8,8,9

  # WALKERS_DIR=..
  # WALKERS_N=$N_WALKERS 
  # WALKERS_ID=$N 

  WALKERS_MPI

...

lwall_588: LOWER_WALLS ARG=d_588 AT=0.500 EPS=0.1 EXP=1 KAPPA=500*(-1)


PRINT ARG=d_489,d_580,cn,d_484,d_576 FILE=COLVAR STRIDE=500
" > plumed.dat
cd ../

done

nohup mpirun -np $N_WALKERS gmx_mpi mdrun -plumed plumed.dat -v -deffnm md -multidir walker_* \
			-gpu_id $GPU_ID -pinoffset $PINOFFSET -cpi -cpo -cpt 15 > md.job & 
