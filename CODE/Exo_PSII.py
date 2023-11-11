# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:54:18 2023

USes Thermo_antenna to simulate PSII quantum efficiency in the habitable zones
of low mass stars. 

@author: Samir Chitnavis and C Duffy
"""

import numpy as np
import Thermo_antenna as thermo
import PSII_params
import matplotlib.pyplot as plt


Ts=['2300','2600','2800','3300','3800','4300','4800','5800']
filename='Inter_Spectrum_PHOENIX_'

e_rate=[] #photosynthetic performance variables
phi_F=[]
gamma1,gamma2=[],[]
DeltaH_12,DeltaS_12,DeltaG_12=[],[],[]
DeltaH_2RC,DeltaS_2RC,DeltaG_2RC=[],[],[]
k_12, k_2RC=[],[]
 
for Temp in Ts: #loop over stellar temperatures
    #import incident spectral flux
    l=[] #x_axis for all spectra (nm)
    irad_y=[] #spectral irradiances (nm)
    fin=open(filename+Temp+'K.txt','r')

    for line in fin:
        line=line.rstrip()
        elements=line.split()
        l.append(float(elements[0]))
        irad_y.append(float(elements[1]))
    
    l=np.array(l) #convert to numpy array to pass to functions in Thermo_antenna
    irad_y=np.array(irad_y)

    #import default params for PSII in 5800K irradiance at T=300K ambient temerature
    G_params=PSII_params.G_params    
    Size_params=PSII_params.Size_params
    k_params=PSII_params.k_params
    T=PSII_params.T
    
    #(1) Performance of PSII under different stars
    model=thermo.antenna(l,irad_y,G_params,Size_params,k_params,T)
    
    #save key parameters
    if model['nu_e']<=100.0:
        e_rate.append(model['nu_e'])
    else:
        e_rate.append(100.0)
    phi_F.append(model['phi_F'])
    gamma1.append(model['gamma'][0])
    gamma2.append(model['gamma'][1])
    DeltaH_12.append(model['DeltaH_12'][0])
    DeltaS_12.append(T*model['DeltaS_12'][0])
    DeltaG_12.append(model['DeltaG_12'][0])
    
    DeltaH_2RC.append(model['DeltaH_2RC'][0])
    DeltaS_2RC.append(T*model['DeltaS_2RC'][0])
    DeltaG_2RC.append(model['DeltaG_2RC'][0])
    
    k_12.append(model['k_12'])
    k_2RC.append(model['k_2RC'])

print('ne_e:', e_rate)
print('phi_e:', phi_F)

print('\nDeltaH_12: ',DeltaH_12)
print('\nDeltaS_12: ',DeltaS_12)
print('\nDeltaG_12: ',DeltaG_12)

print('\nDeltaH_2RC: ',DeltaH_2RC)
print('\nDeltaS_2RC: ',DeltaS_2RC)
print('\nDeltaG_2RC:',DeltaG_2RC)

print('\nk_12: ',k_12)
print('\nk_2RC: ',k_2RC)

#plot the rate and effiiency of PSII for different PS
fig = plt.figure()
plt.scatter(Ts,e_rate,color='maroon',marker='.')
plt.plot(Ts,e_rate,color='maroon',label='$\\nu_{e}$',linestyle='--')
plt.xlabel('$T_{s}$ (K)',fontsize=12)
plt.ylabel('$\\nu_{e}$ ($s^{-1}$)',fontsize=14)
plt.tick_params('x',labelsize=12)
plt.tick_params('y',labelsize=12)
plt.ylim(-5,105)
plt2=plt.twinx()
plt2.set_ylabel('$\phi_{F}$',fontsize=14)
plt2.tick_params('y',labelsize=12)
plt2.scatter(Ts,phi_F,color='darkblue',marker='.')
plt2.plot(Ts,phi_F,color='darkblue',label='$\phi_{F}$',linestyle='--')
plt2.set_ylim(-0.05,1.05)
plt.tight_layout()
fig.legend(frameon=True,fontsize=14,loc=(0.7,0.19))
plt.savefig('Exo_PSII_rate.pdf')
plt.show()

gamma_tot=[a+b for a,b in zip(gamma1,gamma2)]
plt.plot(Ts,gamma1,color='maroon',linestyle='--',label='$\gamma_{1}$')
plt.scatter(Ts,gamma1,color='maroon',marker='.')
plt.plot(Ts,gamma2,color='darkblue',linestyle='--',label='$\gamma_{2}$')
plt.scatter(Ts,gamma2,color='darkblue',marker='.')
plt.plot(Ts,gamma_tot,color='darkgreen',linestyle='--',label='$\gamma_{1}+\gamma_{2}$')
plt.scatter(Ts,gamma_tot,color='darkgreen',marker='.')
plt.xlabel('$T_{s}$ (K)',fontsize=14)
plt.tick_params('x',labelsize=12)
plt.ylabel('$\gamma_{i}$ (photon $s^{-1}$)',fontsize=14)
plt.tick_params('y',labelsize=12)
plt.ylim(0.0,132.0)
plt.legend(loc='upper left',fontsize=14)
plt.savefig('Exo_PSII_gamma.pdf')
plt.show()