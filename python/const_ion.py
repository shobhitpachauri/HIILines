from __future__ import division
import numpy as np
from scipy import interpolate
from const import *

#Recombination coefficients
alphaB_HI     = lambda T4:2.59*10**-13*T4**(-0.833-0.034*np.log(T4))                       #cm^3 s^-1
alphaB_HeI    = lambda T4:2.72*10**-13*T4**-0.789                                          #cm^3 s^-1 #Bruce Draine Eq 14.15
alphaB_Halpha = lambda T4:1.17*10**-13*T4**(-0.942-0.031*np.log(T4))                       #cm^3 s^-1
alphaB_Hbeta  = lambda T4:3.03*10**-14*T4**(-0.874-0.058*np.log(T4))                       #cm^3 s^-1
alpha1_HeI    = lambda T4:1.54*10**-13*T4**-0.486                                          #cm^3 s^-1 #Bruce Draine Eq 14.14
alphaA_HeI    = lambda T4:alphaB_HeI(T4)+alpha1_HeI(T4)
#Osterbrock Table A5.1
#alphaB_OIII   = 3.66*10**-12                                                               #cm^3 s^-1
#alphaB_OII    = 3.99*10**-13                                                               #cm^3 s^-1
#RR: Badnel2006 https://iopscience.iop.org/article/10.1086/508465
#DR: Badnell,N.R.1986,J.Phys.B,19,3827. 2006a
#RR and DR summarized in Aaron Smith COLT bitbucket
#https://bitbucket.org/aaron_smith/colt/src/fb0cd32aeadaedce637a2df46780b1a71a1d3864/src/rates.h
#########################Oxygen#########################################
def alpha_RR_OIII(T4):
    T   = T4*10000
    ST0 = np.sqrt(T / 0.1602)
    ST1 = np.sqrt(T / 4.377e6)
    Bp  = 0.7668 + 0.107 * np.exp(-139200./T)
    return 2.096e-9 / (ST0 * (1 + ST0)** (1-Bp) * (1 + ST1)**(1.+Bp))

def alpha_RR_OII(T4):
    T   = T4*10000
    ST0 = np.sqrt(T / 4.136)
    ST1 = np.sqrt(T / 4.214e6)
    Bp  = 0.6109 + 0.4093 * np.exp(-87700./T)
    return 6.622e-11 / (ST0 * (1 + ST0)**(1-Bp) * (1 + ST1)**(1+Bp))

def alpha_DR_OIII(T4):
    T   = T4*10000
    return T**-1.5 * (1.627e-7 * np.exp(-45.35/T) + 1.262e-7 * np.exp(-284.7/T) + 6.663e-7 * np.exp(-4166./T) + 3.925e-6 * np.exp(-28770./T) + 0.002406 * np.exp(-195300./T) + 0.001146 * np.exp(-364600./T))

def alpha_DR_OII(T4):
    T   = T4*10000
    return T**-1.5 * (5.629e-8 * np.exp(-5395./T) + 2.55e-7 * np.exp(-17700./T) + 0.0006173 * np.exp(-167100./T) + 0.0001627 * np.exp(-268700./T))

alphaB_OIII = lambda T4: alpha_RR_OIII(T4) + alpha_DR_OIII(T4)
alphaB_OII  = lambda T4: alpha_RR_OII(T4)  + alpha_DR_OII(T4)

#Charge transfer rate coefficient, Osterbrock Table A5.4
delta_OII     = 1.05*10**-9                                                                #cm^3 s^-1
delta_OI      = 1.04*10**-9                                                                #cm^3 s^-1
#Charge transfer OII + HI <-> OI + HII, Bruce Draine 14.7.1
k0_OI_ct      = lambda T4: 1.14*10**-9*T4**(0.4+0.018*np.log(T4))                          #cm^3 s^-1
k1_OI_ct      = lambda T4: 3.44*10**-10*T4**(0.451+0.036*np.log(T4))                       #cm^3 s^-1
k2_OI_ct      = lambda T4: 5.33*10**-10*T4**(0.384+0.024*np.log(T4))*np.exp(-97/T4/10000)  #cm^3 s^-1
k0r_OI_ct     = lambda T4: 8/5*k0_OI_ct(T4)*np.exp(-229/T4/10000)                          #cm^3 s^-1

###HI###
#Photoionization cross section
def sigma_HI(nu):
    sigma  = np.zeros(len(nu))   

    E               = h*nu*6.242*10**18
    idx_null        = np.where((E<13.6)|(E>5*10**4))[0]
    sigma[idx_null] = 0
    idx_valid       = np.where((E>=13.6)&(E<=5*10**4))[0]

    E0     = 4.298*10**-1 #eV
    sigma0 = 5.475*10**4 #Mb
    ya     = 3.288*10**1
    P      = 2.963
    yw     = 0
    y0     = 0
    y1     = 0

    x      = E[idx_valid]/E0-y0
    y      = np.sqrt(x**2+y1**2)
    
    sigma[idx_valid] = sigma0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**-P*10**-18 #cm^2
    return sigma

nu_Halpha    = 1.89*eV2J/h #Hz
nu_Hbeta     = 2.55*eV2J/h #Hz

###HeI###
# states: 0: 1^1S_0; 1: 2^1S_0; 2: 2^1P_1; 3: 2^3S_1                 #Osterbrock Table 2.5
#State degeneracy
g0_HeI       = 1
g1_HeI       = 1
g2_HeI       = 3
g3_HeI       = 3
#Spontaneous decay rate
A30_HeI      = 1.26*10**-4 #s^-1
#Collisional (de-)excitation coefficients                            #Osterbrock Table 2.5
T4_grid      = np.array([6000,8000,10000,15000,20000,25000])/10000
k31_grid     = np.array([1.95,2.45,2.60,3.05,2.55,2.68])*10**-8      #cm^3 s^-1
k31_HeI      = interpolate.interp1d(T4_grid,k31_grid,bounds_error=False, fill_value=(k31_grid[0],k31_grid[-1]))
k32_grid     = np.array([2.34,3.64,5.92,7.83,9.23,9.81])*10**-9      #cm^3 s^-1  
k32_HeI      = interpolate.interp1d(T4_grid,k32_grid,bounds_error=False, fill_value=(k32_grid[0],k32_grid[-1]))
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://aas.aanda.org/articles/aas/pdf/2000/18/ds9707.pdf
E30_HeI      = 19.8*eV2J
T4_grid      = 10**np.array([3.75,4.00,4.25,4.50,4.75,5.00,5.25,5.50,5.75])/10000
Omega03_grid = np.array([6.198,6.458,6.387,6.157,5.832,5.320,4.787,4.018,3.167])*10**-2
k30_grid     = 8.629*10**-8/np.sqrt(T4_grid)*Omega03_grid/g3_HeI     #cm^3 s^-1
k30_HeI      = interpolate.interp1d(T4_grid,k30_grid,bounds_error=False, fill_value=(k30_grid[0],k30_grid[-1]))
#Fraction of recombination radiation resulting in hydrogen ionization
p            = lambda ne,T4: 0.75*A30_HeI/(A30_HeI+ne*(k30_HeI(T4)+k31_HeI(T4)+k32_HeI(T4)))+0.25*2/3+0.75*ne*k32_HeI(T4)/(A30_HeI+ne*(k30_HeI(T4)+k31_HeI(T4)+k32_HeI(T4))) \
                           +(0.75*ne*k31_HeI(T4)/(A30_HeI+ne*(k30_HeI(T4)+k31_HeI(T4)+k32_HeI(T4)))+0.25*1/3)*0.56
#Photoionization cross section
def sigma_HeI(nu):
    sigma  = np.zeros(len(nu))

    E               = h*nu*6.242*10**18
    idx_null        = np.where((E<24.59)|(E>5*10**4))[0]
    sigma[idx_null] = 0
    idx_valid       = np.where((E>=24.59)&(E<=5*10**4))[0]


    E0     = 13.61 #eV
    sigma0 = 949.2 #Mb
    ya     = 1.469
    P      = 3.188
    yw     = 2.039
    y0     = 0.4434
    y1     = 2.136

    x      = E[idx_valid]/E0-y0
    y      = np.sqrt(x**2+y1**2)
    
    sigma[idx_valid] = sigma0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**-P*10**-18 #cm^2
    return sigma

def sigma_HeII(nu):
    sigma  = np.zeros(len(nu))

    E               = h*nu*6.242*10**18
    idx_null        = np.where((E<54.42)|(E>5*10**4))[0]
    sigma[idx_null] = 0
    idx_valid       = np.where((E>=54.42)&(E<=5*10**4))[0]


    E0     = 1.72    #eV
    sigma0 = 1.369e4 #Mb
    ya     = 32.88
    P      = 2.963
    yw     = 0
    y0     = 0
    y1     = 0

    x      = E[idx_valid]/E0-y0
    y      = np.sqrt(x**2+y1**2)

    sigma[idx_valid] = sigma0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**-P*10**-18 #cm^2
    return sigma

###OII### 
def sigma_OI(nu):
    sigma  = np.zeros(len(nu))

    E               = h*nu*6.242*10**18
    idx_null        = np.where((E<13.62)|(E>538))[0]
    sigma[idx_null] = 0
    idx_valid       = np.where((E>=13.62)&(E<=538))[0]

    E0     = 1.240 #eV
    sigma0 = 1.745*10**3 #Mb
    ya     = 3.784
    P      = 17.64
    yw     = 7.589*10**-2
    y0     = 8.698
    y1     = 1.271*10**-1

    x      = E[idx_valid]/E0-y0
    y      = np.sqrt(x**2+y1**2)
    sigma[idx_valid] = sigma0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**-P*10**-18 #cm^2
    return sigma

#State degeneracy                                                    
g0_OII  = 4
g1_OII  = 6
g2_OII  = 4
g3_OII  = 4
g4_OII  = 2
#Spontaneous decay rate                                             #Cloudy data source
A10_OII = 7.416e-06 + 3.382e-05
A20_OII = 1.414e-04 + 2.209e-05
A21_OII = 1.30e-07  + 1.49e-20
A30_OII = 5.22e-02  + 2.43e-07
A31_OII = 8.37e-03  + 9.07e-02
A32_OII = 1.49e-02  + 3.85e-02
A40_OII = 2.12e-02  + 3.72e-07
A41_OII = 8.34e-03  + 5.19e-02
A42_OII = 9.32e-03  + 7.74e-02
A43_OII = 1.41e-10  + 4.24e-24
#Level energy                                                       #Cloudy data source
E10_OII = 38575*kb                                                  #J
E20_OII = 38604*kb
E30_OII = 58225*kb
E40_OII = 58228*kb
E21_OII = E20_OII-E10_OII
E31_OII = E30_OII-E10_OII
E32_OII = E30_OII-E20_OII
E41_OII = E40_OII-E10_OII
E42_OII = E40_OII-E20_OII
E43_OII = E40_OII-E30_OII
#Level energy frequency
nu10_OII = E10_OII/h                                                #Hz
nu20_OII = E20_OII/h
nu21_OII = E21_OII/h
nu30_OII = E30_OII/h
nu31_OII = E31_OII/h
nu32_OII = E32_OII/h
nu40_OII = E40_OII/h
nu41_OII = E41_OII/h
nu42_OII = E42_OII/h
nu43_OII = E43_OII/h
#Collisional (de-)excitation coefficients                           #Bruce Draine Table F4
Omega10_OII = lambda T4:0.803*T4**(0.023-0.008*np.log(T4))
k10_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega10_OII(T4)/g1_OII #cm^3 s^-1
k01_OII     = lambda T4:g1_OII/g0_OII*k10_OII(T4)*np.exp(-E10_OII/kb/T4/10000)
Omega20_OII = lambda T4:0.550*T4**(0.054-0.004*np.log(T4))
k20_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega20_OII(T4)/g2_OII #cm^3 s^-1
k02_OII     = lambda T4:g2_OII/g0_OII*k20_OII(T4)*np.exp(-E20_OII/kb/T4/10000)
Omega21_OII = lambda T4:1.434*T4**(-0.176+0.004*np.log(T4))
k21_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega21_OII(T4)/g2_OII #cm^3 s^-1
k12_OII     = lambda T4:g2_OII/g1_OII*k21_OII(T4)*np.exp(-E21_OII/kb/T4/10000)
Omega30_OII = lambda T4:0.140*T4**(0.025-0.006*np.log(T4))
k30_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega30_OII(T4)/g3_OII #cm^3 s^-1
k03_OII     = lambda T4:g3_OII/g0_OII*k30_OII(T4)*np.exp(-E30_OII/kb/T4/10000)
Omega31_OII = lambda T4:0.349*T4**(0.060+0.052*np.log(T4))
k31_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega31_OII(T4)/g3_OII #cm^3 s^-1
k13_OII     = lambda T4:g3_OII/g1_OII*k31_OII(T4)*np.exp(-E31_OII/kb/T4/10000)
Omega32_OII = lambda T4:0.326*T4**(0.063+0.052*np.log(T4))
k32_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega32_OII(T4)/g3_OII #cm^3 s^-1
k23_OII     = lambda T4:g3_OII/g2_OII*k32_OII(T4)*np.exp(-E32_OII/kb/T4/10000)
Omega40_OII = lambda T4:0.283*T4**(0.023-0.004*np.log(T4))
k40_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega40_OII(T4)/g4_OII #cm^3 s^-1
k04_OII     = lambda T4:g4_OII/g0_OII*k40_OII(T4)*np.exp(-E40_OII/kb/T4/10000)
Omega41_OII = lambda T4:0.832*T4**(0.076+0.055*np.log(T4))
k41_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega41_OII(T4)/g4_OII #cm^3 s^-1
k14_OII     = lambda T4:g4_OII/g1_OII*k41_OII(T4)*np.exp(-E41_OII/kb/T4/10000)
Omega42_OII = lambda T4:0.485*T4**(0.059+0.052*np.log(T4))
k42_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega42_OII(T4)/g4_OII #cm^3 s^-1
k24_OII     = lambda T4:g4_OII/g2_OII*k42_OII(T4)*np.exp(-E42_OII/kb/T4/10000)
Omega43_OII = lambda T4:0.322*T4**(0.019+0.037*np.log(T4))
k43_OII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega43_OII(T4)/g4_OII #cm^3 s^-1
k34_OII     = lambda T4:g4_OII/g3_OII*k43_OII(T4)*np.exp(-E43_OII/kb/T4/10000)

#Five level rates
R01_OII = lambda ne,T4: ne*k01_OII(T4)
R02_OII = lambda ne,T4: ne*k02_OII(T4)
R03_OII = lambda ne,T4: ne*k03_OII(T4)
R04_OII = lambda ne,T4: ne*k04_OII(T4)
R10_OII = lambda ne,T4: ne*k10_OII(T4)+A10_OII
R12_OII = lambda ne,T4: ne*k12_OII(T4)
R13_OII = lambda ne,T4: ne*k13_OII(T4)
R14_OII = lambda ne,T4: ne*k14_OII(T4)
R20_OII = lambda ne,T4: ne*k20_OII(T4)+A20_OII
R21_OII = lambda ne,T4: ne*k21_OII(T4)+A21_OII
R23_OII = lambda ne,T4: ne*k23_OII(T4)
R24_OII = lambda ne,T4: ne*k24_OII(T4)
R30_OII = lambda ne,T4: ne*k30_OII(T4)+A30_OII
R31_OII = lambda ne,T4: ne*k31_OII(T4)+A31_OII
R32_OII = lambda ne,T4: ne*k32_OII(T4)+A32_OII
R34_OII = lambda ne,T4: ne*k34_OII(T4)
R40_OII = lambda ne,T4: ne*k40_OII(T4)+A40_OII
R41_OII = lambda ne,T4: ne*k41_OII(T4)+A41_OII
R42_OII = lambda ne,T4: ne*k42_OII(T4)+A42_OII
R43_OII = lambda ne,T4: ne*k43_OII(T4)+A43_OII

#Photoionization cross section
def sigma_OII(nu):
    #https://arxiv.org/abs/astro-ph/9601009v2
    sigma  = np.zeros(len(nu))

    E               = h*nu*6.242*10**18
    idx_null        = np.where((E<35.12)|(E>558.1))[0]
    sigma[idx_null] = 0
    idx_valid       = np.where((E>=35.12)&(E<=558.1))[0]


    E      = h*nu*6.242*10**18 #eV
    E0     = 1.386 #eV
    sigma0 = 5.967*10 #Mb
    ya     = 3.175*10
    P      = 8.943
    yw     = 1.934*10**-2
    y0     = 2.131*10
    y1     = 1.503*10**-2

    x      = E[idx_valid]/E0-y0
    y      = np.sqrt(x**2+y1**2)
    
    sigma[idx_valid] = sigma0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**-P*10**-18 #cm^2
    return sigma

###OIII###
#State degeneracy
g0_OIII  = 1
g1_OIII  = 3
g2_OIII  = 5
g3_OIII  = 5
g4_OIII  = 1
#Spontaneous decay rate                     #http://adsabs.harvard.edu/full/1968IAUS...34..143G
A10_OIII = 2.6*10**-5  #2.7*10**-5  #s^-1
A20_OIII = 3.5*10**-11 #3.1*10**-11 #s^-1
A21_OIII = 9.8*10**-5  #9.7*10**-5  #s^-1
A30_OIII = 1.9*10**-6               #s^-1
A31_OIII = 0.0071      #6.8*10**-3  #s^-1
A32_OIII = 0.021       #2.0*10**-2  #s^-1
A40_OIII = 0                        #s^-1
A41_OIII = 0.23                     #s^-1
A42_OIII = 7.1*10**-4               #s^-1
A43_OIII = 1.6                      #s^-1
#Level energy and frequency
E10_OIII  = 163*kb                  #J
E20_OIII  = 441*kb
E30_OIII  = 29169*kb
E40_OIII  = 61207*kb
E21_OIII  = E20_OIII-E10_OIII
E31_OIII  = E30_OIII-E10_OIII
E32_OIII  = E30_OIII-E20_OIII
E41_OIII  = E40_OIII-E10_OIII
E42_OIII  = E40_OIII-E20_OIII
E43_OIII  = E40_OIII-E30_OIII
nu10_OIII = E10_OIII/h              #s^-1 
nu20_OIII = E20_OIII/h
nu30_OIII = E30_OIII/h
nu40_OIII = E40_OIII/h
nu21_OIII = E21_OIII/h
nu31_OIII = E31_OIII/h
nu32_OIII = E32_OIII/h
nu41_OIII = E41_OIII/h
nu42_OIII = E42_OIII/h
nu43_OIII = E43_OIII/h
#Collisional (de-)excitation coefficients                           #Bruce Draine Table F2 
Omega10_OIII = lambda T4:0.522*T4**(0.033-0.009*np.log(T4))
k10_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega10_OIII(T4)/g1_OIII #cm^3 s^-1
k01_OIII     = lambda T4:g1_OIII/g0_OIII*k10_OIII(T4)*np.exp(-E10_OIII/kb/T4/10000)
Omega20_OIII = lambda T4:0.257*T4**(0.081+0.017*np.log(T4))
k20_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega20_OIII(T4)/g2_OIII #cm^3 s^-1
k02_OIII     = lambda T4:g2_OIII/g0_OIII*k20_OIII(T4)*np.exp(-E20_OIII/kb/T4/10000)
Omega21_OIII = lambda T4:1.23*T4**(0.053+0.007*np.log(T4))
k21_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega21_OIII(T4)/g2_OIII #cm^3 s^-1
k12_OIII     = lambda T4:g2_OIII/g1_OIII*k21_OIII(T4)*np.exp(-E21_OIII/kb/T4/10000)
Omega30_OIII = lambda T4:0.243*T4**(0.12+0.031*np.log(T4))
k30_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega30_OIII(T4)/g3_OIII #cm^3 s^-1
k03_OIII     = lambda T4:g3_OIII/g0_OIII*k30_OIII(T4)*np.exp(-E30_OIII/kb/T4/10000)
Omega31_OIII = lambda T4:0.243*T4**(0.12+0.031*np.log(T4))*3
k31_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega31_OIII(T4)/g3_OIII #cm^3 s^-1
k13_OIII     = lambda T4:g3_OIII/g1_OIII*k31_OIII(T4)*np.exp(-E31_OIII/kb/T4/10000)
Omega32_OIII = lambda T4:0.243*T4**(0.12+0.031*np.log(T4))*5
k32_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega32_OIII(T4)/g3_OIII #cm^3 s^-1
k23_OIII     = lambda T4:g3_OIII/g2_OIII*k32_OIII(T4)*np.exp(-E32_OIII/kb/T4/10000)
Omega40_OIII = lambda T4:0.0321*T4**(0.118+0.057*np.log(T4))
k40_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega40_OIII(T4)/g4_OIII #cm^3 s^-1
k04_OIII     = lambda T4:g4_OIII/g0_OIII*k40_OIII(T4)*np.exp(-E40_OIII/kb/T4/10000)
Omega41_OIII = lambda T4:0.0321*T4**(0.118+0.057*np.log(T4))*3
k41_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega41_OIII(T4)/g4_OIII #cm^3 s^-1
k14_OIII     = lambda T4:g4_OIII/g1_OIII*k41_OIII(T4)*np.exp(-E41_OIII/kb/T4/10000)
Omega42_OIII = lambda T4:0.0321*T4**(0.118+0.057*np.log(T4))*5
k42_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega42_OIII(T4)/g4_OIII #cm^3 s^-1
k24_OIII     = lambda T4:g4_OIII/g2_OIII*k42_OIII(T4)*np.exp(-E42_OIII/kb/T4/10000)
Omega43_OIII = lambda T4:0.523*T4**(0.210-0.099*np.log(T4))
k43_OIII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega43_OIII(T4)/g4_OIII #cm^3 s^-1
k34_OIII     = lambda T4:g4_OIII/g3_OIII*k43_OIII(T4)*np.exp(-E43_OIII/kb/T4/10000)
#Five level rates
R01_OIII = lambda ne,T4: ne*k01_OIII(T4)
R02_OIII = lambda ne,T4: ne*k02_OIII(T4)
R03_OIII = lambda ne,T4: ne*k03_OIII(T4)
R04_OIII = lambda ne,T4: ne*k04_OIII(T4)
R10_OIII = lambda ne,T4: ne*k10_OIII(T4)+A10_OIII
R12_OIII = lambda ne,T4: ne*k12_OIII(T4)
R13_OIII = lambda ne,T4: ne*k13_OIII(T4)
R14_OIII = lambda ne,T4: ne*k14_OIII(T4)
R20_OIII = lambda ne,T4: ne*k20_OIII(T4)+A20_OIII
R21_OIII = lambda ne,T4: ne*k21_OIII(T4)+A21_OIII
R23_OIII = lambda ne,T4: ne*k23_OIII(T4)
R24_OIII = lambda ne,T4: ne*k24_OIII(T4)
R30_OIII = lambda ne,T4: ne*k30_OIII(T4)+A30_OIII
R31_OIII = lambda ne,T4: ne*k31_OIII(T4)+A31_OIII
R32_OIII = lambda ne,T4: ne*k32_OIII(T4)+A32_OIII
R34_OIII = lambda ne,T4: ne*k34_OIII(T4)
R40_OIII = lambda ne,T4: ne*k40_OIII(T4)+A40_OIII
R41_OIII = lambda ne,T4: ne*k41_OIII(T4)+A41_OIII
R42_OIII = lambda ne,T4: ne*k42_OIII(T4)+A42_OIII
R43_OIII = lambda ne,T4: ne*k43_OIII(T4)+A43_OIII

#NII
#State degeneracy
g0_NII  = 1
g1_NII  = 3
g2_NII  = 5
g3_NII  = 5
g4_NII  = 1
#Spontaneous decay rate            #Cloudy
A10_NII = 2.08*10**-6               #s^-1
A20_NII = 1.12*10**-12             #s^-1
A21_NII = 7.46*10**-6              #s^-1
A30_NII = 5.25*10**-7              #s^-1
A31_NII = 9.22*10**-7+9.84*10**-4  #s^-1
A32_NII = 8.65*10**-6+2.91*10**-3  #s^-1
A40_NII = 0                        #s^-1
A41_NII = 3.18*10**-2              #s^-1
A42_NII = 1.55*10**-4              #s^-1
A43_NII = 1.14                     #s^-1
#Level energy and frequency
E10_NII  = 70*kb                   #J
E20_NII  = 188*kb
E30_NII  = 22037*kb
E40_NII  = 47033*kb
E21_NII  = E20_NII-E10_NII
E31_NII  = E30_NII-E10_NII
E32_NII  = E30_NII-E20_NII
E41_NII  = E40_NII-E10_NII
E42_NII  = E40_NII-E20_NII
E43_NII  = E40_NII-E30_NII
nu10_NII = E10_NII/h               #s^-1
nu20_NII = E20_NII/h
nu30_NII = E30_NII/h
nu40_NII = E40_NII/h
nu21_NII = E21_NII/h
nu31_NII = E31_NII/h
nu32_NII = E32_NII/h
nu41_NII = E41_NII/h
nu42_NII = E42_NII/h
nu43_NII = E43_NII/h
#Collisional (de-)excitation coefficients                           #Bruce Draine Table F2
Omega10_NII = lambda T4:0.431*T4**(0.099+0.014*np.log(T4))
k10_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega10_NII(T4)/g1_NII #cm^3 s^-1
k01_NII     = lambda T4:g1_NII/g0_NII*k10_NII(T4)*np.exp(-E10_NII/kb/T4/10000)
Omega20_NII = lambda T4:0.273*T4**(0.166+0.030*np.log(T4))
k20_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega20_NII(T4)/g2_NII #cm^3 s^-1
k02_NII     = lambda T4:g2_NII/g0_NII*k20_NII(T4)*np.exp(-E20_NII/kb/T4/10000)
Omega21_NII = lambda T4:1.15*T4**(0.137+0.024*np.log(T4))
k21_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega21_NII(T4)/g2_NII #cm^3 s^-1
k12_NII     = lambda T4:g2_NII/g1_NII*k21_NII(T4)*np.exp(-E21_NII/kb/T4/10000)
Omega30_NII = lambda T4:0.303*T4**(0.053+0.009*np.log(T4))
k30_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega30_NII(T4)/g3_NII #cm^3 s^-1
k03_NII     = lambda T4:g3_NII/g0_NII*k30_NII(T4)*np.exp(-E30_NII/kb/T4/10000)
Omega31_NII = lambda T4:0.909*T4**(0.053+0.010*np.log(T4))
k31_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega31_NII(T4)/g3_NII #cm^3 s^-1
k13_NII     = lambda T4:g3_NII/g1_NII*k31_NII(T4)*np.exp(-E31_NII/kb/T4/10000)
Omega32_NII = lambda T4:1.51*T4**(0.054+0.011*np.log(T4))
k32_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega32_NII(T4)/g3_NII #cm^3 s^-1
k23_NII     = lambda T4:g3_NII/g2_NII*k32_NII(T4)*np.exp(-E32_NII/kb/T4/10000)
Omega40_NII = lambda T4:0.0352*T4**(0.066+0.018*np.log(T4))
k40_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega40_NII(T4)/g4_NII #cm^3 s^-1
k04_NII     = lambda T4:g4_NII/g0_NII*k40_NII(T4)*np.exp(-E40_NII/kb/T4/10000)
Omega41_NII = lambda T4:0.105*T4**(0.070+0.021*np.log(T4))
k41_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega41_NII(T4)/g4_NII #cm^3 s^-1
k14_NII     = lambda T4:g4_NII/g1_NII*k41_NII(T4)*np.exp(-E41_NII/kb/T4/10000)
Omega42_NII = lambda T4:0.176*T4**(0.065+0.017*np.log(T4))
k42_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega42_NII(T4)/g4_NII #cm^3 s^-1
k24_NII     = lambda T4:g4_NII/g2_NII*k42_NII(T4)*np.exp(-E42_NII/kb/T4/10000)
Omega43_NII = lambda T4:0.806*T4**(-0.175-0.014*np.log(T4))
k43_NII     = lambda T4:8.629*10**-8/np.sqrt(T4)*Omega43_NII(T4)/g4_NII #cm^3 s^-1
k34_NII     = lambda T4:g4_NII/g3_NII*k43_NII(T4)*np.exp(-E43_NII/kb/T4/10000)
#Five level rates
R01_NII = lambda ne,T4: ne*k01_NII(T4)
R02_NII = lambda ne,T4: ne*k02_NII(T4)
R03_NII = lambda ne,T4: ne*k03_NII(T4)
R04_NII = lambda ne,T4: ne*k04_NII(T4)
R10_NII = lambda ne,T4: ne*k10_NII(T4)+A10_NII
R12_NII = lambda ne,T4: ne*k12_NII(T4)
R13_NII = lambda ne,T4: ne*k13_NII(T4)
R14_NII = lambda ne,T4: ne*k14_NII(T4)
R20_NII = lambda ne,T4: ne*k20_NII(T4)+A20_NII
R21_NII = lambda ne,T4: ne*k21_NII(T4)+A21_NII
R23_NII = lambda ne,T4: ne*k23_NII(T4)
R24_NII = lambda ne,T4: ne*k24_NII(T4)
R30_NII = lambda ne,T4: ne*k30_NII(T4)+A30_NII
R31_NII = lambda ne,T4: ne*k31_NII(T4)+A31_NII
R32_NII = lambda ne,T4: ne*k32_NII(T4)+A32_NII
R34_NII = lambda ne,T4: ne*k34_NII(T4)
R40_NII = lambda ne,T4: ne*k40_NII(T4)+A40_NII
R41_NII = lambda ne,T4: ne*k41_NII(T4)+A41_NII
R42_NII = lambda ne,T4: ne*k42_NII(T4)+A42_NII
R43_NII = lambda ne,T4: ne*k43_NII(T4)+A43_NII
#Cross section
def sigma_NI(nu):
    #https://arxiv.org/abs/astro-ph/9601009v2
    sigma           = np.zeros(len(nu))

    E               = h*nu*6.242*10**18
    idx_null        = np.where((E<14.53)|(E>404.8))[0]
    sigma[idx_null] = 0
    idx_valid       = np.where((E>=14.53)&(E<=404.8))[0]

    E      = h*nu*6.242*10**18 #eV
    E0     = 4.034             #eV
    sigma0 = 8.235*100         #Mb
    ya     = 8.033*10
    P      = 3.928
    yw     = 9.097*10**-2
    y0     = 8.598*10**-1
    y1     = 2.325

    x      = E[idx_valid]/E0-y0
    y      = np.sqrt(x**2+y1**2)

    sigma[idx_valid] = sigma0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**-P*10**-18 #cm^2
    return sigma

def sigma_NII(nu):
    sigma           = np.zeros(len(nu))

    E               = h*nu*6.242*10**18
    idx_null        = np.where((E<29.6)|(E>423.6))[0]
    sigma[idx_null] = 0
    idx_valid       = np.where((E>=29.6)&(E<=423.6))[0]

    E      = h*nu*6.242*10**18 #eV
    E0     = 6.128*10**-2      #eV
    sigma0 = 1.944             #Mb
    ya     = 8.163*10**2
    P      = 8.773
    yw     = 1.043*10
    y0     = 4.280*10**2
    y1     = 2.030*10

    x      = E[idx_valid]/E0-y0
    y      = np.sqrt(x**2+y1**2)

    sigma[idx_valid] = sigma0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**-P*10**-18 #cm^2
    return sigma

def sigma_NIII(nu):
    sigma           = np.zeros(len(nu))

    E               = h*nu*6.242*10**18
    idx_null        = np.where((E<47.45)|(E>447.3))[0]
    sigma[idx_null] = 0
    idx_valid       = np.where((E>=47.45)&(E<=447.3))[0]

    E      = h*nu*6.242*10**18 #eV
    E0     = 0.2420            #eV
    sigma0 = 0.9375            #Mb
    ya     = 278.8
    P      = 9.156
    yw     = 1.850
    y0     = 187.7
    y1     = 3.999

    x      = E[idx_valid]/E0-y0
    y      = np.sqrt(x**2+y1**2)

    sigma[idx_valid] = sigma0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**-P*10**-18 #cm^2
    return sigma

