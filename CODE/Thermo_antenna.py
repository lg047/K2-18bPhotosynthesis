# -*- coding: utf-8 -*-
"""
Code to simulate the performance of a two-compartment themrodynamic antenna 
coupled to an irreversible exciton trap (reaction centre)

@author: Samir Chitnavis, Chris Duffy
"""

import numpy as np
from scipy.constants import h as h
from scipy.constants import c as c
from scipy.constants import Boltzmann as kB

#Absorption lines
#Returns a normalized Gaussian absorption line for a wavelength range, l (nm), 
#peak wavelength, l_p (nm), standard deviation, w (nm) 
def gauss_abs(l,lp,w):
    exponent=-((l-lp)**2)/(2.0*w**2)
    gauss_y=np.exp(exponent)
    #normalize
    N=np.trapz(gauss_y,l)
    return(gauss_y/N)    
    
#overlap integral calculator
def Olap_Int(l,f1,f2):
    olap=[a*b for a,b in zip(f1,f2)]
    return(np.trapz(olap,l))
    
#Enthalpy change on exciton transfer
#lp1 and lp2 are the peak wavelengths of the two antenna domains
#N12=n1/n2 where ni is the number of thermodynamically equivalent states in domain i 
#T is the temperature in K
def deltaG(lp1,lp2,N12,T):
    
    lp1=lp1*1.0E-9 #convert to SI units
    lp2=lp2*1.0E-9
    
    H_12=h*c*(lp1-lp2)/(lp1*lp2) #Enthalpy change on moving from domain 1 to domain 2
    H_21=h*c*(lp2-lp1)/(lp1*lp2) #Enthalpy change for the reverse transfer
    
    S_12=kB*np.log(1.0/N12) #Entropy change for transfer from domain 1 to domain 2
    S_21=kB*np.log(N12) #Entropy change for the reverse transfer
    
    G_12=H_12-(T*S_12)
    G_21=H_21-(T*S_21)
    
    return([H_12,H_21],[S_12,S_21],[G_12,G_21])

# Main antenna model
#Inputs:
    #l = a 1d np.array of wavelengths (nm) covering the irradiance window
    #Ip_y = the spectral irradiance for the points in l (in W per m^2 per nm).
    #G_params=(sig,B12,lp1,w1,lp2,w2)
        #sig = integrated cross-section (in m^2)
        #B12 = ration of the amplitudes of g1/g2
        #lpi = peak wavelength (nm) of peak i
        #wi = width (standard dev, nm) of peak i
    #k_params=(k_diss,k_trap,k_con), empirical rate constants for decay channels
        #k_diss = rate constant for dissipation of Chl excited states, generally ~1/4ns
        #k_trap = photochemical trapping rate in the RC, ~1/4ps
        #k_con = Qy reduction rate constant, generally ~1/10ms
        #K_12 = fundamental hopping rate from G1 to G2
        #K_23 = fundamental hopping rate from G1 to GRC (G3)
        #Size ratio of domain 1 and domain 2 
        #Size ration of domain 2 and the RC 
def antenna(l,Ip_y,G_params,Size_params,k_params,T):
    #spectral parameters
    sig,B12=G_params[0],G_params[1]
    lp1,w1=G_params[2],G_params[3]
    lp2,w2=G_params[4],G_params[5]
    lpRC,wRC=G_params[6],G_params[7] #for later developments involing spectral overlap
    
    #size parameters
    N1, N2, N_LHC=Size_params[0], Size_params[1], Size_params[2]
    
    #Decay rates
    k_diss, k_trap, k_con=k_params[0], k_params[1], k_params[2]
    K_12,K_2RC=k_params[3], k_params[4]

    #(1) Convert the spectral iradiance into the spectral photon flux
    fp_y=np.zeros(len(Ip_y)) #photon flux
    for i, item in enumerate(Ip_y):
        fp_y[i]=item*((l[i]*1.0E-9)/(h*c)) #factor of 1E-9 since l is in nm.
    
    #(2) Generate antenna absorption profile
    g1_l=gauss_abs(l,lp1,w1)
    g2_l=gauss_abs(l,lp2,w2)
    gRC_l=gauss_abs(l,lpRC,wRC) #for later developments involing spectral overlap
    
    #(3) Calculate gamma1 and gamma2 (photon input rate)
    #IT is dependent on the number of LHCs
    gamma1=N_LHC*sig*B12*Olap_Int(l,fp_y,g1_l)
    gamma2=N_LHC*sig*Olap_Int(l,fp_y,g2_l)
        
    #(4)Transfer rates
    #This first requires computing the free energy change associated with the 
    #transfer step 
    
    n1n2=N1/N2 #ratio of thermodynamically equivalent states in G1 and G2 
    thermo12=deltaG(lp1,lp2,n1n2,T) #thermodynamic parameters for G1 -> G2
    G_12,G_21=thermo12[2][0],thermo12[2][1]
    
    n2nRC=N_LHC*N2
    thermo2RC=deltaG(lp2,lpRC,n2nRC,T) #thermodynamic parameters for G2 -> GRC
    G_2RC,G_RC2=thermo2RC[2][0],thermo2RC[2][1]
    
    #(5) a: transfer between 1 -> 2
    if G_12==0.0:
        k_12, k_21=K_12, K_12
    elif G_12<0.0:
        k_12, k_21=K_12, K_12*np.exp(-G_21/(kB*T))
    elif G_12>0.0:
        k_12, k_21=K_12*np.exp(-G_12/(kB*T)), K_12
        
    #(5) b: transfer between 2 -> RC
    if G_2RC==0.0:
        k_2RC, k_RC2=K_2RC, K_2RC
    elif G_2RC<0.0:
        k_2RC, k_RC2=K_2RC, K_2RC*np.exp(-G_RC2/(kB*T))
    elif G_2RC>0.0:
        k_2RC, k_RC2=K_2RC*np.exp(-G_2RC/(kB*T)), K_2RC
        
    #(6) Build Transfer Matrix and Vector
    Kmat=np.zeros((4,4))
    Kmat[0][0], Kmat[0][1]=-(k_12+k_diss), k_21
    Kmat[1][0], Kmat[1][1], Kmat[1][2]=k_12, -(k_21+k_2RC+k_diss), k_RC2
    Kmat[2][1], Kmat[2][2]=k_2RC, -(k_RC2+k_trap+k_diss)
    Kmat[3][2], Kmat[3][3]= k_trap, -k_con
        
    gamma_vec=np.array([-gamma1, -gamma2, 0.0, 0.0])
    
    #(7) Solve the equation
    Kmat_inv=np.linalg.inv(Kmat) #take inverse
    Neq=np.zeros(4)
    for i in range(4):
        for j in range(4):
            Neq[i]=Neq[i]+Kmat_inv[i][j]*gamma_vec[j]
    
    #(8) Calculate the observables and output
    e_rate=k_con*Neq[3] #electron output rate
    phi_F=e_rate/(e_rate+k_diss*(Neq[0]+Neq[1])) #fluoresence quantum yield
    
    output_dict={
        'nu_e': e_rate,
        'phi_F': phi_F,
        'gamma': [gamma1,gamma2],
        'N_eq': Neq, 
        'K_mat': Kmat,
        'DeltaH_12': [thermo12[0][0],thermo12[0][1]],
        'DeltaS_12': [thermo12[1][0],thermo12[1][1]],
        'DeltaG_12': [thermo12[2][0],thermo12[2][1]],
        'DeltaH_2RC': [thermo2RC[0][0],thermo2RC[0][1]],
        'DeltaS_2RC': [thermo2RC[1][0],thermo2RC[1][1]],
        'DeltaG_2RC': [thermo2RC[2][0],thermo2RC[2][1]],
        'k_12': [k_12, k_21],
        'k_2RC': [k_2RC,k_RC2]
        }
    
    return(output_dict)
    
'''############################################################################
# Purple bacterial model                                                      #
############################################################################'''


# Main antenna model but for a basic purple bacterial (Typ-II RC) antenna
#Inputs:
    #l = a 1d np.array of wavelengths (nm) covering the irradiance window
    #Ip_y = the spectral irradiance for the points in l (in W per m^2 per nm).
    #G_params=(sig,B12,lp1,w1,lp2,w2)
        #sig = integrated cross-section (in m^2)
        #B12 = ration of the amplitudes of g1/g2
        #lpi = peak wavelength (nm) of peak i
        #wi = width (standard dev, nm) of peak i
    #Size_params
        #N1 = Number of B850 BCHl a's per LH2
        #N2 = Number of B800        ''
        #N_LHC = antenna size
        #N_LHC0 = reference antenna size
    #k_params=(k_diss,k_trap,k_con), empirical rate constants for decay channels
        #k_diss = rate constant for dissipation of Chl excited states, generally ~1/4ns
        #k_trap = photochemical trapping rate in the RC, ~1/4ps
        #k_con = Qy reduction rate constant, generally ~1/10ms
        #K_12 = fundamental hopping rate from G1 to G2
        #K_23 = fundamental hopping rate from G1 to GRC (G3)
def antenna_TypeII(l,Ip_y,G_params,Size_params,k_params,T):
    #spectral parameters
    sig,B12=G_params[0],G_params[1]
    lp1,w1=G_params[2],G_params[3]
    lp2,w2=G_params[4],G_params[5]
    lpRC,wRC=G_params[6],G_params[7] #for later developments involing spectral overlap
    
    #size parameters
    N1, N2, = Size_params[0], Size_params[1]
    N_LHC, a0 = Size_params[2], Size_params[3]
    
    #Decay rates
    k_diss, k_trap, k_con=k_params[0], k_params[1], k_params[2]
    K_12,K_2RC=k_params[3], k_params[4]

    #(1) Convert the spectral iradiance into the spectral photon flux
    fp_y=np.zeros(len(Ip_y)) #photon flux
    for i, item in enumerate(Ip_y):
        fp_y[i]=item*((l[i]*1.0E-9)/(h*c)) #factor of 1E-9 since l is in nm.
    
    #(2) Generate antenna absorption profile
    g1_l=gauss_abs(l,lp1,w1)
    g2_l=gauss_abs(l,lp2,w2)
    gRC_l=gauss_abs(l,lpRC,wRC) #for later developments involing spectral overlap
    
    #(3) Calculate gamma1 and gamma2 (photon input rate)
    #IT is dependent on the number of LHCs
    gamma1=N_LHC*sig*B12*Olap_Int(l,fp_y,g1_l)
    gamma2=N_LHC*sig*Olap_Int(l,fp_y,g2_l)
        
    #(4)Transfer rates
    #This first requires computing the free energy change associated with the 
    #transfer step 
    
    n1n2=N1/N2 #ratio of thermodynamically equivalent states in G1 and G2 
    thermo12=deltaG(lp1,lp2,n1n2,T) #thermodynamic parameters for G1 -> G2
    G_12,G_21=thermo12[2][0],thermo12[2][1]
    
    '''########################################################################
    # Phenomenological entropy function for the LH2-L1-RC transfer step.
    # (Improve later)
    ########################################################################'''

    thermo2RC=[[0.0,0.0],[0.0,0.0],[0.0,0.0]]

    #enthalpy
    thermo2RC[0][0]=h*c*(lp2-lpRC)/(lp2*lpRC)
    thermo2RC[0][1]=h*c*(lpRC-lp2)/(lp2*lpRC)


    #entropy
    thermo2RC[1][0]=kB*np.log(a0/(N2*N_LHC)) 
    thermo2RC[1][1]=kB*np.log((N2*N_LHC)/a0)
    
    #free energy
    thermo2RC[2][0]=thermo2RC[0][0]-(T*thermo2RC[1][0])        
    thermo2RC[2][1]=thermo2RC[0][1]-(T*thermo2RC[1][1])

    G_2RC,G_RC2=thermo2RC[2][0],thermo2RC[2][1]
    
    #(5) a: transfer between 1 -> 2
    if G_12==0.0:
        k_12, k_21=K_12, K_12
    elif G_12<0.0:
        k_12, k_21=K_12, K_12*np.exp(-G_21/(kB*T))
    elif G_12>0.0:
        k_12, k_21=K_12*np.exp(-G_12/(kB*T)), K_12
        
    #(5) b: transfer between 2 -> RC
    if G_2RC==0.0:
        k_2RC, k_RC2=K_2RC, K_2RC
    elif G_2RC<0.0:
        k_2RC, k_RC2=K_2RC, K_2RC*np.exp(-G_RC2/(kB*T))
    elif G_2RC>0.0:
        k_2RC, k_RC2=K_2RC*np.exp(-G_2RC/(kB*T)), K_2RC
        
    #(6) Build Transfer Matrix and Vector
    Kmat=np.zeros((4,4))
    Kmat[0][0], Kmat[0][1]=-(k_12+k_diss), k_21
    Kmat[1][0], Kmat[1][1], Kmat[1][2]=k_12, -(k_21+k_2RC+k_diss), k_RC2
    Kmat[2][1], Kmat[2][2]=k_2RC, -(k_RC2+k_trap+k_diss)
    Kmat[3][2], Kmat[3][3]= k_trap, -k_con
        
    gamma_vec=np.array([-gamma1, -gamma2, 0.0, 0.0])
    
    #(7) Solve the equation
    Kmat_inv=np.linalg.inv(Kmat) #take inverse
    Neq=np.zeros(4)
    for i in range(4):
        for j in range(4):
            Neq[i]=Neq[i]+Kmat_inv[i][j]*gamma_vec[j]
    
    #(8) Calculate the observables and output
    e_rate=k_con*Neq[3] #electron output rate
    phi_F=e_rate/(e_rate+k_diss*(Neq[0]+Neq[1])) #fluoresence quantum yield
    
    output_dict={
        'nu_e': e_rate,
        'phi_F': phi_F,
        'gamma': [gamma1,gamma2],
        'N_eq': Neq, 
        'K_mat': Kmat,
        'DeltaH_12': [thermo12[0][0],thermo12[0][1]],
        'DeltaS_12': [thermo12[1][0],thermo12[1][1]],
        'DeltaG_12': [thermo12[2][0],thermo12[2][1]],
        'DeltaH_2RC': [thermo2RC[0][0],thermo2RC[0][1]],
        'DeltaS_2RC': [thermo2RC[1][0],thermo2RC[1][1]],
        'DeltaG_2RC': [thermo2RC[2][0],thermo2RC[2][1]],
        'k_12': [k_12, k_21],
        'k_2RC': [k_2RC,k_RC2]
        }
    
    return(output_dict)
    




