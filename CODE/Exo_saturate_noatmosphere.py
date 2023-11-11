# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 15:12:14 2023

@author: Ifigeneia
"""
import numpy as np
import Thermo_antenna as thermo
import TypeII_params
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as kB



Ts=['2300','2600','2800','3300','3800','4300','4800','5800']
filename='Scaled_Spectrum_PHOENIX_'

out = open('newParameters_no_atmosphere.txt','w')

    #loop over temperatures
l_temp=[] #x_axis for all spectra (nm)
flux_temp=[] #spectral irradiances (nm)
    
e_rate_Ts=[]
phi_F_Ts=[]
rate_Ts=[]
for temp in Ts:
    file_temp=filename+temp+'K.txt'
    fin=open(file_temp,'r')
    
    l=[] #x_axis for all spectra (nm)
    flux=[] #spectral irradiances (nm)
    for line in fin:
        line=line.rstrip()
        elements=line.split()
        l.append(float(elements[0]))
        flux.append(float(elements[1]))
        

    l=np.array(l)
    flux=np.array(flux)

        #import default params for PSII at T=300K ambient temerature
    G_params=TypeII_params.G_params    
    Size_params=TypeII_params.Size_params
    k_params=TypeII_params.k_params
    T=TypeII_params.T

    n=1 #start with N_LHC = 5
    Size_params[2]=n #re-assign number of antenna subunits
    e_rate_old=0.00001
    model=thermo.antenna_TypeII(l,flux,G_params,Size_params,k_params,T)
    e_rate_new=model['nu_e']

    phi_F=0.0
    while (e_rate_new-e_rate_old)/e_rate_old>=0.001:
        n=n+1
        Size_params[2]=n #re-assign number of antenna subunits
        e_rate_old=e_rate_new
        model=thermo.antenna_TypeII(l,flux,G_params,Size_params,k_params,T)
            
        if(model['nu_e']<100.0):
            e_rate_new=model['nu_e']
        else:
            e_rate_new=100.0
                
        phi_F_new=model['phi_F']

        #The current values of e_rate_new and phi_F_new will be the ones we want
    rate = e_rate_new/n
    e_rate_Ts.append(e_rate_new)
    phi_F_Ts.append(phi_F_new) 
    rate_Ts.append(e_rate_new/n)  
    out.write(str(temp)+' '+str(e_rate_new)+' '+str(rate)+' '+str(phi_F_new)+'\n') 
