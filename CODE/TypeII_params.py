# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 10:00:02 2023

Optimum parameters for RC-LH1-LH2 in a 5800K spectrum

@author: btw430
"""

'''***************************************************************************
* Spectral parameters                                                      '''

sig=1.81E-18 #optical cross-section of the antenna (m^2)
B12=0.70 #ratio of Chlb/Chb peak amplitude 
lp1=799.0#B850 absorption peak (nm)
w1=10.5 #B850 Gaussian width (nm)
lp2=848.0 #B800 absorption peak
w2=12.0 #B850 Gaussian width (nm)
lpRC=870.0 #reaction centre
wRC=w2 #width is the same as Chla

G_params=[sig,B12,lp1,w1,lp2,w2,lpRC,wRC] #tuple of these parameters to send to the antenna model
'''*************************************************************************'''

'''***************************************************************************
* Dimensions parameters                                                     '''


N1=18.0 #number of accessory pigments in a single LHC
N2=9.0 #number of primary pigmets in a single LHCII
N_LHC=5.0 #The number of antenna complexes per reaction centre (does not have 
          #to be an integer)
a0 = 45.0

Size_params=[N1,N2,N_LHC,a0] #list of these parameters to send to the antenna model
#N_LHC appears twice as we want to use the firstinstance as a variable and the 
#fourth value as a reference constant

'''*************************************************************************'''

'''***************************************************************************
* Rate Parameters                                                          '''

k_diss=1.0/1.0E-9 #Chl excited state decay rate 
#k_trap=1.0/5.0E-12 #PSII trapping rate
k_trap=1/5.0E-12
k_con=1.0/10.0E-3 #PSII RC turnover rate
K_12=1.0/1.0E-12 #fundamental B800 - B850 hopping rate
K_2RC=1.0/50.0E-12 #fundanental B850-RC hopping rate


k_params=[k_diss,k_trap,k_con,K_12,K_2RC]

'''*************************************************************************'''

T=300.0