# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 16:20:05 2023

@author: Ifigeneia
"""


import matplotlib.pyplot as plt #import plt function from matplotlib library
import numpy as np 


files=['earth','archean','archean_lower','archean_upper','archean_upperch4','archean_upperco2']
legends =['Modern Earth','Archean Earth','Archean lower limits','Archean upper limits','Archean high methane, low $CO_{2}$','Archean high $CO_{2}$, low methane']
#files = ['earth','water_0','water_2','water_4']
#files = ['earth','archean']
#files = ['earth','methane_05','methane_1','methane_5','methane_10']
#legends=['Modern Earth',r'$H_{2}O$ $0$%',r'$H_{2}O$ $2$%',r'$H_{2}O$ $4$%']
#legends=['Modern Earth', 'Archean Earth']
#legends=['Modern Earth',r'$CH_{4}$ $0.5$%',r'$CH_{4}$ $1$%',r'$CH_{4}$ $5$%',r'$CH_{4}$ $10$%']
lambda_lim=[400.0,750.0]
colours=['royalblue','gold','darkorange','red','maroon','purple']

for i, file in enumerate(files):
    
    transmissionData_wavelength = []
    transmissionData_transmittance = []

    fil = open(file+'.txt','r')
    lines_after_7 = fil.readlines()[7:]

    for line in lines_after_7: 

        line = line.rstrip()
        elements = line.split(' ')
        wavelength = float(elements[0])*1000 
        transmittance = float(elements[1])
        transmissionData_wavelength.append(wavelength)
        transmissionData_transmittance.append(transmittance)
  
    
    l=[] #x_axis for all spectra (nm)
    flux=[] #spectral irradiances (nm)
    fin=open('Scaled_Spectrum_PHOENIX_5800K.txt')

    for line in fin:
        line=line.rstrip()
        elements=line.split()
        l.append(float(elements[0]))
        flux.append(float(elements[1]))
    
    transmissionInterpolated = np.interp(l, transmissionData_wavelength, transmissionData_transmittance)
    print(transmissionInterpolated)


    difference = []
    new_flux = []
    count = 0
    for f in flux:
        new = f*transmissionInterpolated[count] #flux with the introduction of atmosphere
        new_flux.append(new)
        diff = (new-f)*100.0/f
        difference.append(diff)
        count= count + 1
    if i == 0:
        plt.plot(l,flux, label = 'No atmosphere',linestyle='-',linewidth=4,alpha=0.7,color='darkgrey')
    plt.plot(l,new_flux,label = legends[i], color=colours[i])
    plt.xlabel('Wavelength (nm)',fontsize=14)
    plt.ylabel('Incident spectral flux (Wm$^{-2}$nm$^{-1}$)',fontsize=14)


plt.legend(loc = 1)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.xlim(200,1400)
#plt.xlim(lambda_lim[0],lambda_lim[1])
plt.xlim(800.0,1100.0)
plt.title(r'$T_{s}=$5800 K',fontsize=14)
plt.savefig('archean_all_spectral_5800K.pdf')
plt.show()