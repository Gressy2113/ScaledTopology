#!/bin/bash

N_WALKERS=4
# GPU_ID=0,1
# NT=16
# PINOFFSET=0
# NT1=$((NT/N_WALKERS))
echo $NT1

MPIRUN=/usr/lib64/openmpi/bin/mpirun

for ((N = 0 ; N < $N_WALKERS ; N++ ))
do

rm -r walker_$N
mkdir -p walker_$N
cd walker_$N 
cp ../eq.$N.gro eq.gro
cp ../plmd.dat .
cp ../init.pdb .

gmx_mpi grompp -f ../md.mdp -c eq.gro -r ../posres.gro -p ../system.top -n ../indata/grps.ndx -o md

cd ../


done

# nohup taskset --cpu-list 0-$NT mpirun -np $N_WALKERS gmx_mpi mdrun -plumed plumed.dat -v -deffnm md -multidir walker_* \
# 			-gpu_id $GPU_ID -pinoffset $PINOFFSET -cpi -cpo -cpt 15 > md.job & 


nohup $MPIRUN -np $N_WALKERS gmx_mpi mdrun -plumed plmd.dat -v -deffnm md \
      -multidir walker_* \
      -nb gpu -bonded gpu -pme gpu -pmefft gpu \
      -gpu_id 0,1 -pin on -pinstride 1 \
      -cpi -cpo -cpt 15 > md.job & 


#nohup mpirun -np 4 gmx_mpi mdrun -plumed plumed.dat -v -deffnm md -multidir walker_* -gputasks 0011 -ntomp 11 -cpi -cpo -cpt 15 > md1.job & 
