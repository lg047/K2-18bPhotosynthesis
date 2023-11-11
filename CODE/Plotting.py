# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:51:23 2023

@author: Ifigeneia
"""

import numpy as np
import matplotlib.pyplot as plt


files = ['no_atmosphere','earth','water_0','water_2','water_4']
legends=['No atmosphere','Modern Earth',r'$H_{2}O$ $0$%',r'$H_{2}O$ $2$%',r'$H_{2}O$ $4$%']
#files=['earth','archean','archean_lower','archean_upper','archean_upperch4','archean_upperco2']
#legends =['Modern Earth','Archean Earth','Archean lower limits','Archean upper limits','Archean high methane, low $CO_{2}$','Archean high $CO_{2}$, low methane']
#files = ['no_atmosphere','earth','archean']
#legends = ['No atmosphere','Modern Earth', 'Archean Earth']
#files = ['no_atmosphere','earth','methane_05','methane_1','methane_5','methane_10']
#legends=['No atmosphere','Modern Earth',r'$CH_{4}$ $0.5$%',r'$CH_{4}$ $1$%',r'$CH_{4}$ $5$%',r'$CH_{4}$ $10$%']
#colours=['darkgrey','royalblue','gold','darkorange','red','maroon']
colours=['darkgrey','royalblue','gold','darkorange','maroon','purple']


for i,file in enumerate(files):
    
    temp = []
    e_rate = []
    rate = []
    phi_rate = []
    file_name='newParameters_'+file+'.txt'
    fil = open(file_name,'r')
    for n,line in enumerate(fil):
        temperature=["2300","2600","2800","3300","3800","4300","4800","5800"]
        line=line.rstrip()
        elements=line.split()
        temp.append(float(elements[0]))
        e_rate.append(float(elements[1]))
        rate.append(float(elements[2]))
        phi_rate.append(float(elements[3]))
        
    plt.scatter(temperature,phi_rate,label=legends[i],color=colours[i])
    plt.xlabel(r'$T_{s}$ $(K)$',fontsize=14)
    plt.ylabel(r'$\phi_{F}$',fontsize=14)
    #print(temp)
    #plt.xticks([2300,2600,2800,3300,3800,4300,4800,5800])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Type II')
plt.legend()
plt.savefig('water_typeII_phi.pdf')
plt.show()

for i,file in enumerate(files):
    temperature=["2300","2600","2800","3300","3800","4300","4800","5800"]
    temp = []
    e_rate = []
    rate = []
    phi_rate = []
    file_name='newParameters_'+file+'.txt'
    fil = open(file_name,'r')
    for line in fil:
        line=line.rstrip()
        elements=line.split()
        temp.append(float(elements[0]))
        e_rate.append(float(elements[1]))
        rate.append(float(elements[2]))
        phi_rate.append(float(elements[3]))
        
    plt.scatter(temperature,rate,label=legends[i],color=colours[i])
    plt.xlabel(r'$T_{s}$ $(K)$',fontsize=14)
    plt.ylabel(r'$\nu_{e}^{max}/N_{LHC}$ (s$^{-1}$)',fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Type II')
plt.legend()
plt.savefig('water_typeII__ratio.pdf')
plt.show()
    
for i,file in enumerate(files):
    temperature=["2300","2600","2800","3300","3800","4300","4800","5800"]
    temp = []
    e_rate = []
    rate = []
    phi_rate = []
    file_name='newParameters_'+file+'.txt'
    fil = open(file_name,'r')
    for line in fil:
        line=line.rstrip()
        elements=line.split()
        temp.append(float(elements[0]))
        e_rate.append(float(elements[1]))
        rate.append(float(elements[2]))
        phi_rate.append(float(elements[3]))

        
    plt.scatter(temperature,e_rate,label=legends[i],color=colours[i])
    plt.xlabel(r'$T_{s}$ $(K)$',fontsize=14)
    plt.ylabel(r'$\nu_{e}^{max}$ (s$^{-1}$)',fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.title('Type II')
plt.legend()
plt.savefig('water_typeII__nue.pdf')
plt.show()
    

