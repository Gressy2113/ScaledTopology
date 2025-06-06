#!/bin/bash 


nohup python SCRIPTS/plotstring.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-MG/mw_0.8 200 200 fin > np1.out & 
nohup python SCRIPTS/plotstring.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-CA/mw_0.85 200 200 fin > np2.out & 
nohup python SCRIPTS/plotstring.py /home/veretenenko/TRPV6-Mg/ScaledTopology-final/S5-S6_adapted-CA/mw_0.8 200 200 fin > np3.out & 
