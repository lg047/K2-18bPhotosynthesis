# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 12:17:29 2023

Calculates the efficiency of PSII-like photosystem with antennae of increasing 
sizes

@author: Samir Chitnavis and Chris Duffy
"""

import numpy as np
import Thermo_antenna as thermo
import PSII_params
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as kB


Ts=['2300','2600','2800','3300','3800','4300','4800','5800']
colours=['maroon','red','darkorange','gold','darkgreen','darkcyan','darkblue','fuchsia']
filename='Inter_Spectrum_PHOENIX_'
N_LHC=np.arange(1.0,56.0,5.0) #range of antenna subunits

e_rate=[]
phi_F=[]
DeltaS_GRC=[]
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
    
    #loop over N_LHCII
    e_rate_T=[] #save electron output rate
    phi_F_T=[] #save PSII yield
    DeltaS_GRC_T=[] #save free energy change
    for n in N_LHC:
        Size_params[2]=n #re-assign number of antenna subunits
        model=thermo.antenna(l,irad_y,G_params,Size_params,k_params,T)
        if model['nu_e']<=100.0:            
            e_rate_T.append(model['nu_e'])
        else: 
            e_rate_T.append(100.0)
            
        phi_F_T.append(model['phi_F'])
        DeltaS_GRC_T.append(model['DeltaG_2RC'][0])
        
    e_rate.append(e_rate_T)
    phi_F.append(phi_F_T)
    DeltaS_GRC.append(DeltaS_GRC_T)

#Plotting
#(1) electron production rate
plt.fill_between([-0.5,52.0],50.0,100.0,alpha=0.2,color='green',linewidth=0.0)
for i, temp in enumerate(Ts):    
    temp_label=str(temp)+'K'
    plt.plot(N_LHC,e_rate[i],color=colours[i],linestyle='-')
    plt.scatter(N_LHC,e_rate[i],label=temp_label,color=colours[i],marker='.')

#plt.hlines(125.0,-0.5,51.0,color='k',linestyle='--')
#plt.hlines(75.0,-0.5,51.0,color='k',linestyle='--')

plt.xlim(-0.5,72)
plt.ylim(-5.0,105.0)
plt.xlabel('$N_{LHC}$',fontsize=14)
plt.tick_params('x',labelsize=12)
plt.ylabel('$\\nu_{e}$ ($s^{-1}$)',fontsize=14)
plt.tick_params('y',labelsize=12)
plt.legend(fontsize=12,loc='upper right')
plt.savefig('Exo_N_LHC.pdf')
plt.show()

#(2) Quantum effiiency (independent of spectrum) and entropy
fig = plt.figure()
plt.vlines(5.0,-0.1,7.0,color='k',linestyle='--',alpha=0.2)
plt.hlines(0.85,-0.5,5.0,color='k',linestyle='--',alpha=0.2)
plt.plot(N_LHC,phi_F[-1],color='maroon',linestyle='-')
plt.scatter(N_LHC,phi_F[-1],color='maroon',marker='.',label='$\phi_{F}$')
plt.ylim(0.18,1.05)
plt.xlim(-0.5,53.0)
plt.xlabel('$N_{LHC}$',fontsize=14)
plt.ylabel('$\phi_{F}$',fontsize=14)
plt.tick_params('x',labelsize=12)
plt.tick_params('y',labelsize=12)

DG_kT=[]
for DG in DeltaS_GRC[0]:
    DG_kT.append(DG/(kB*300.0))    

plt2=plt.twinx()
plt2.hlines(4.15,5.0,60.0,color='k',linestyle='--',alpha=0.2)
plt2.set_ylabel('$\Delta F_{2,RC}/k_{B}T$',fontsize=14)
plt2.tick_params('y',labelsize=12)
plt2.set_ylim(1.9,10.2)
plt2.plot(N_LHC,DG_kT,color='darkblue',linestyle='-')
plt2.scatter(N_LHC,DG_kT,color='darkblue',marker='.',label='$\Delta F_{2,RC}$')

fig.legend(frameon=True,fontsize=14,loc=(0.66,0.69))
plt.savefig('Exo_phi_F.pdf')
plt.show()


