#!/bin/bash 


nohup python SCRIPTS/plotstring.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-MG/mw_0.8 200 200 att6 > np1.out & 
nohup python SCRIPTS/plotstring.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-CA/mw_0.85 200 200 att0 > np1.out & 
nohup python SCRIPTS/plotstring.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-CA/mw_0.8 200 200 att0 > np1.out & 
