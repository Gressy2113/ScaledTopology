#!/bin/bash 

nohup python SCRIPTS/block_prot.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-MG/mw_0.8 0 > np1.out & 
nohup python SCRIPTS/block_prot.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-CA/mw_0.85 0 > np2.out & 
nohup python SCRIPTS/block_prot.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-CA/mw_0.8 0 > np3.out & 
