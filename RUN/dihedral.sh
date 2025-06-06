#!bin/bash

# for SYS in S5-S6_adapted-CA/MD #S5-S6_adapted-MG/MD
# do
#     cd /home/veretenenko/TRPV6-Mg/ScaledTopology-final/$SYS
#     mkdir -p dihedral
#     nohup /nfs/belka2/soft/impulse/dev/inst/runtask.py -t ../../IMPULSE/calc_gp.mk \
#                 -f ../../IMPULSE/args_gp/args.ch1_489 > nohup_ch1_489.out & 

#     nohup /nfs/belka2/soft/impulse/dev/inst/runtask.py -t ../../IMPULSE/calc_gp.mk \
#                 -f ../../IMPULSE/args_gp/args.ch2_489 > nohup_ch2_489.out & 

#     nohup /nfs/belka2/soft/impulse/dev/inst/runtask.py -t ../../IMPULSE/calc_gp.mk \
#                 -f ../../IMPULSE/args_gp/args.ch1_580 > nohup_ch1_580.out & 

#     nohup /nfs/belka2/soft/impulse/dev/inst/runtask.py -t ../../IMPULSE/calc_gp.mk \
#                 -f ../../IMPULSE/args_gp/args.ch2_580 > nohup_ch2_580.out & 

# done

for SYS in S5-S6_adapted-CA/mw_0.8 #S5-S6_adapted-MG/mw_0.8
do
    cd /home/veretenenko/TRPV6-Mg/ScaledTopology-final/$SYS
    for i in 0 1 2 3
    do
        cd walker_$i
        mkdir -p dihedral
        nohup /nfs/belka2/soft/impulse/dev/inst/runtask.py -t ../../../IMPULSE/calc_gp.mk \
                    -f ../../../IMPULSE/args_gp/args.ch1_489 > nohup_ch1_489.out & 

        nohup /nfs/belka2/soft/impulse/dev/inst/runtask.py -t ../../../IMPULSE/calc_gp.mk \
                    -f ../../../IMPULSE/args_gp/args.ch2_489 > nohup_ch2_489.out & 

        nohup /nfs/belka2/soft/impulse/dev/inst/runtask.py -t ../../../IMPULSE/calc_gp.mk \
                    -f ../../../IMPULSE/args_gp/args.ch1_580 > nohup_ch1_580.out & 

        nohup /nfs/belka2/soft/impulse/dev/inst/runtask.py -t ../../../IMPULSE/calc_gp.mk \
                    -f ../../../IMPULSE/args_gp/args.ch2_580 > nohup_ch2_580.out & 
        cd ../
    done

done