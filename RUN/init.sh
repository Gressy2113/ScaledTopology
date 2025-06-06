#!/bin/bash

# SYSTEM=S5-S6_adapted-CA/mw_0.8
# FOLDER=/nfs/corona1/veretenenko/TRPV6-Mg/ScaledTopology-S5-S6/System.0.CA_bc.multiwalkers/METAD.3D.b5h0.3_new

# mkdir -p $SYSTEM
# cd $SYSTEM

# ln -s $FOLDER/HILLS
# ln -s $FOLDER/FES

# for w in 0 1 2 3
# do
#     mkdir -p walker_$w
#     ln -s $FOLDER/walker_$w/md.tpr walker_$w/md.tpr
#     ln -s $FOLDER/walker_$w/md.xtc walker_$w/md.xtc
#     ln -s $FOLDER/walker_$w/COLVAR.$w walker_$w/COLVAR.$w

# done


for folder in S5-S6_adapted-CA/mw_0.8 S5-S6_adapted-CA/mw_0.85 \
            S5-S6_apo-MG/mw_0.8 S5-S6_apo-CA/mw_0.85
do
    cd /home/veretenenko/TRPV6-Mg/ScaledTopology-final/$folder
    paste -d "\n" walker_0/COLVAR.0 \
                    walker_1/COLVAR.1 \
                    walker_2/COLVAR.2 \
                    walker_3/COLVAR.3 > COLVAR
    head -n -3 COLVAR > tmp && mv tmp COLVAR
    rm -f tmp
done