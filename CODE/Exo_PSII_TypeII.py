# -*- coding: utf-8 -*-
"""
Created on Friday 28th of July, 2023

USes Thermo_antenna to simulate turnover rate and quantum efficiency 
of a the PSII and Type II RC in the habitable zones
of low mass stars. 

@author: Samir Chitnavis and C Duffy
"""

import numpy as np
import Thermo_antenna as thermo
import PSII_params
import TypeII_params
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as kB


Ts=['2300','2600','2800','3300','3800','4300','4800','5800']
filename='Scaled_Spectrum_PHOENIX_'
T=300.0 #ambient temperature

nu_e_PSII_T, nu_e_TypeII_T=[], [] #photosynthetic performance variables
phi_e_PSII_T, phi_e_TypeII_T=[], []
 
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

    #(1) Plant PSII

    #import default params for PSII in 5800K irradiance at T=300K ambient temerature
    G_params_PSII=PSII_params.G_params    
    Size_params_PSII=PSII_params.Size_params
    k_params_PSII=PSII_params.k_params
    
    model_PSII=thermo.antenna(l,irad_y,G_params_PSII,Size_params_PSII,k_params_PSII,T)

    nu_e_PSII=model_PSII['nu_e']
    phi_e_PSII=model_PSII['phi_F']
    DG_2RC=[model_PSII['DeltaG_2RC'][0]/(kB*T),model_PSII['DeltaG_2RC'][1]/(kB*T)]
    tau_2RC=[1.0E12/model_PSII['k_2RC'][0],1.0E12/model_PSII['k_2RC'][1]]

    nu_e_PSII_T.append(nu_e_PSII)
    phi_e_PSII_T.append(phi_e_PSII)

    if(Temp=='5800'):
        print('PSII RC')
        print('Ts = ', str(Temp)+' K')
        print('nu_e = '+str(nu_e_PSII))
        print('phiF = '+str(phi_e_PSII))
        print('DG_2RC = ('+str(DG_2RC[0])+', '+str(DG_2RC[1])+')')
        print('tau_2RC = '+str(tau_2RC[0])+', '+str(tau_2RC[1]))
        print('\n')

    #(2) Purple bacterial Type-II

    G_params_TypeII=TypeII_params.G_params    
    Size_params_TypeII=TypeII_params.Size_params
    k_params_TypeII=TypeII_params.k_params

    model_TypeII=thermo.antenna_TypeII(l,irad_y,G_params_TypeII,Size_params_TypeII,k_params_TypeII,T)
    
    nu_e_TypeII=model_TypeII['nu_e']
    phi_e_TypeII=model_TypeII['phi_F']
    DG_2RC=[model_TypeII['DeltaG_2RC'][0]/(kB*T),model_TypeII['DeltaG_2RC'][1]/(kB*T)]
    tau_2RC=[1.0E12/model_TypeII['k_2RC'][0],1.0E12/model_TypeII['k_2RC'][1]]

    nu_e_TypeII_T.append(nu_e_TypeII)
    phi_e_TypeII_T.append(phi_e_TypeII)

    if(Temp=='5800'):
        print('Type-II RC')
        print('Ts = ', str(Temp)+' K')
        print('nu_e = '+str(nu_e_TypeII))
        print('phiF = '+str(phi_e_TypeII))
        print('DG_2RC = ('+str(DG_2RC[0])+', '+str(DG_2RC[1])+')')
        print('tau_2RC = '+str(tau_2RC[0])+', '+str(tau_2RC[1]))
        print('\n')
    

#plot the rate and effiiency of PSII and Type II for different Ts
#numerical Ts
Ts_num=[]
for Temp in Ts:
    Ts_num.append(float(Temp))

fig = plt.figure()
plt.scatter(Ts_num,nu_e_PSII_T,color='maroon',marker='.')
plt.plot(Ts_num,nu_e_PSII_T,color='maroon',label='PSII ($\phi_{e}=$'+str('%.*f' % (2, phi_e_PSII_T[0]))+')',linestyle='--')
plt.scatter(Ts_num,nu_e_TypeII_T,color='darkblue',marker='.')
plt.plot(Ts_num,nu_e_TypeII_T,color='darkblue',label='Type II ($\phi_{e}=$'+str('%.*f' % (2, phi_e_TypeII_T[0]))+')',linestyle='--')
plt.hlines(y=100.0,xmin=2000.0,xmax=6000.0,linestyle='--',color='grey')
plt.xticks(Ts_num,rotation=90)
plt.xlabel('$T_{s}$ (K)',fontsize=14)
plt.ylabel('$\\nu_{e}$ ($s^{-1}$)',fontsize=14)
plt.tick_params('x',labelsize=12)
plt.tick_params('y',labelsize=12)
plt.xlim(2000.0,6000.0)
plt.ylim(-5,110)

plt.legend(loc='lower right',fontsize=14)
plt.tight_layout()
plt.savefig('Exo_PSII_TypeII_rate.pdf')
plt.show()
'''

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

'''