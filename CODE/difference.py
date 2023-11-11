# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:01:50 2023

@author: Ifigeneia
"""


import matplotlib.pyplot as plt #import plt function from matplotlib library
import numpy as np 

#files = ['earth','water_0','water_2','water_4']
files = ['earth','archean']
#files = ['earth','methane_05','methane_1','methane_5','methane_10']
colours=['royalblue','gold','darkorange','red','maroon']
#colours=['royalblue','gold','darkorange','maroon']
#legends=['Modern Earth',r'$H_{2}O$ $0$%',r'$H_{2}O$ $2$%',r'$H_{2}O$ $4$%']
#legends=['Modern Earth',r'$CH_{4}$ $0.5$%',r'$CH_{4}$ $1$%',r'$CH_{4}$ $5$%',r'$CH_{4}$ $10$%']
legends=['Modern Earth', 'Archean Earth']
lambda_lim=[400.0,750.0]
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
    fin=open('Scaled_Spectrum_PHOENIX_2800K.txt')

    for line in fin:
        line=line.rstrip()
        elements=line.split()
        l.append(float(elements[0]))
        flux.append(float(elements[1]))
    
    
    transmissionInterpolated = np.interp(l, transmissionData_wavelength, transmissionData_transmittance)

    f_max=max(flux)
    difference = []
    for n,f in enumerate(flux):
        new = f*transmissionInterpolated[n]
        difference.append((new-f)/f_max*100.0)

    
    plt.plot(l,difference,label = legends[i], color=colours[i])
    plt.xlabel('Wavelength (nm)',fontsize=13)
    plt.ylabel('% flux reduction',fontsize=13)
plt.legend(loc=3)
#plt.xlim(200,1400)
plt.xlim(lambda_lim[0],lambda_lim[1])
#plt.xlim(800.0,1100.0)
plt.ylim(-100.0,0.0)
plt.title(r'$T_{s}=$2800 K')
plt.savefig('archean_diff_oxy_2800K.pdf')
plt.show()