#!/usr/bin/env python
# coding: utf-8

# ###### T.A. Maryam Azhin,
# 
# Department of Chemical and Materials Engineering, University of Alberta

# # Seminar 10. Isothermal effective diffusivity. Internal effectiveness factor. 

# **Q1. Calculation of effective diffusivity**. 
# 
# From “Introduction to Chemical Reactor Analysis”, R.E. Hayes, J.P. Mmbaga, 2nd Ed. CRC Press.
# 
# Calculate effective diffusivity for methane in air at $500 K$ and $120 kPa$ in a catalyst with a porosity of $0.4$ (uniform pore distribution), pellet density of $1,400 \frac{kg}{m^3}$ and a BET surface area of $180,000 \frac{m^2}{kg}$. Assume tortuosity factor of $4$.
# 
# Values for diffusion volumesome values, $νi$, for some molecules: 
# 
# |     |     |    |    |
# |:---:|----:|---:|---:|
# |H2   |7.07 | CO |18.9|
# |N2   |17.9 | CO2|26.9|
# |O2   |16.6 | H2O|12.7|
# |air  |20.1 | NH3|14.9|
# |CH4  |24.42| He |2.88|
# |Ar   |16.1 | N2O|35.9|

# **Answer to Q1.**
# 
# $T=500K$
# 
# $M_A(CH_4)=16\;\frac{g}{mol}$
# 
# $M_B(air)=29\;\frac{g}{mol}$
# 
# $P=120000 pa$
# 
# $\nu_A(CH_4)=24.42$
# 
# $\nu_B(air)=20.1$
# 
# Binary molecular diffusion in gases $[m2/s]$ can be calculated using the Fuller formula 
# 
# ($P$–pressure in $[Pa]$, $T$ – temperature in $[K]$, $M$ – molar mass in $[g/mol]$, $\nu_i$ is a diffusion volume.):
# 
# \begin{eqnarray}
# D_{AB}=\frac{1.013*10^{-2}*T^{1.75}(\frac{1}{MA}+\frac{1}{MB})^{0.5}}{P*((\sum \nu_i)_A^{1/3}+(\sum \nu_i)_B^{1/3})^2}=\frac{1.013*10^{-2}*500^{1.75}(\frac{1}{16}+\frac{1}{29})^{0.5}}{120000*(24.42^{1/3}+20.1^{1/3})^2}=4.4*10^{-5}\frac{m^2}{s}
# \end{eqnarray}
# 
# Which is a reasonable number for a gas phase.
# 
# Knudsen diffusivity $D_K$ for gas species varies with the pore diameter $d_p$ in $[m]$, $T$ the temperature in $[K]$, and the molar mass of the diffusion species, $M$, in $[g/mol]$. 
# 
# We can calculate $d_p$ from the following equation:
# 
# ($\phi_P$ – pellet porosity, $S$ – catalyst surface area in $[\frac{m^2}{kg}]$ (for example, from BET measurements) $\rho$ - pellet density in $[\frac{kg}{m^3}]$.)
# 
# \begin{eqnarray}
# d_p=\frac{4*\phi_p}{S*\rho}=\frac{4*0.4}{180000 [\frac{m^2}{kg}]*1400[\frac{kg}{m^3}]}=6.35*10^{-9} [m]\\
# \end{eqnarray}
# 
# The equivalent pore diameter is taken as the average diameter of the pores in the pellet, which gives a reasonable approximation provided that there is a relatively narrow pore size distribution. 
# 
# The Knudsen diffusion coefficient, $D_k$ in $[m^2/s]$ in a straight round pore can be calculated from the following formula: ($M=16 \frac{g}{mol}$ molecular mass of the diffusing species.)
# 
# \begin{eqnarray}
# D_k=48.5*d_p*(\frac{T}{M})^{0.5}=48.5*6.35*10^{-9}*(\frac{500}{16})^{0.5}=1.72*10^{-6}\frac{m^2}{s}
# \end{eqnarray}
# 
# $D_{pore}$ is diffusivity in a pore, which is a combination of molecular and Knudsen diffusivity. In general, the concentration dependence of the pore diffusion coefficient is not large, and in many cases, the following equation (Bosanquet formula) is used to compute the pore diffusion coefficient.
# 
# \begin{eqnarray}
# D_{pore}=(\frac{1}{D_{AB}}+\frac{1}{D_K})^{-1}=1.66*10^{-6}\frac{m^2}{s}
# \end{eqnarray}
# 
# Therefore, it is governed mostly by $D_k$ because of small pore size (width).
# 
# Effective diffusivity in a catalyst for unimodal pore size distribution, $D_{eff}$:
# 
# ($\phi_p$ – pellet porosity, $\tau$ – tortuosity factor, typically $3-4$, for calculation example refer to Fogler, Elements of CRE, example 12-1 page 817)
# 
# \begin{eqnarray}
# D_{eff}=\frac{\phi_pD_{pore}}{\tau}=\frac{0.4*1.66*10^{-6}}{4}=1.66*10^{-7}\frac{m^2}{s}
# \end{eqnarray}
# 
# Note that $D_{eff}$ is $265$ times lower than bulk diffusivity $D_{AB}$ and $D_{eff}$ does not depend on catalyst particle size.
# 
# FYI: When a catalyst pellet has two distinct regions of pore sizes (for example, micro or meso pores inside an individual catalyst grain and macro pores between catalyst  grains pressed into a catalyst pellet), then effective diffusivity is:
# 
# \begin{eqnarray}
# D_{eff}=D_{pore,M}\phi^2_{p,M}+[\frac{\phi^2_{p,Micro}(1+3\phi_{p,M})}{(1-\phi_{p,M})}]D_{pore,micro}
# \end{eqnarray}
# 
# Where $M$ refers to a large - pore region and "micro" - to a smaller pore region. In each of these regions that diffusivities can be calculated using the prevoiusly mentioned equations and corresponding $\phi_p$.

# In[1]:


import math
T=500
MA=16 #g/mol
MB=29 #g/mol
P=120000#𝑝𝑎 
nuA=24.42 
nuB=20.1
phi_p=0.4 #pellet porosity
S=180000  #catalyst surface area
rho=1400  #pellet density
R=8.314   #J/mol/K
tau=4     #tortuosity factor
dp=4*phi_p/(S*rho) #m
DAB=(1.013*pow(10,-2)*pow(T,1.75)*pow((1/MA+1/MB),0.5))/(P*pow(pow(nuA,1/3)+pow(nuB,1/3),2))
#Dk=(dp/3)*pow(((8*R*T)/(3.14*MA*0.01)),0.5)
Dk=48.5*dp*pow((T/(MA)),0.5)
Dpore=pow((1/DAB+1/Dk),-1)
Deff=(phi_p*Dpore)/tau
print("DAB={0:.6f}".format(DAB),'bulk diffusivity [m2/s]')
print("dp={0:.10f}".format(dp),'m')
print("Dk={0:.8f}".format(Dk),'Knudsen diffusivity [m2/s]')
print("Dpore={0:.8f}".format(Dpore),'m2/s')
print("Deff={0:.9f}".format(Deff),'m2/s')
print("DAB/Deff={0:.1f}".format(DAB/Deff))


# ----

# **Q2. Calculation of generalized Thiele modulus and internal effectiveness factor**
# 
# A heterogeneous catalytic reaction follows a Langmuir-Hinshelwood mechanism but under certain conditions the L-H rate law can be simplified to a first-order power law as $r = kC_A$. Specific surface area of catalyst: $100 [\frac{m^2}{g}]$, effective diffusivity: $0.02 [\frac{cm^2}{s}]$ (or $[\frac{cm^3\; fluid}{ cm\; catalyst\;*\;s}]$), intrinsic rate constant: $2.8*10^{-7} [\frac{m^3\;fluid}{m^2\;catalyst\;*\;s}]$  (or $[\frac{m}{s}]$). The reaction is in gas-phase, process $T$ is $500 [K]$, partial pressure of $A$ is $0.24\; [MPa]$, compressibility factor is $0.8$, catalyst bed density is $1.2*10^6 [\frac{g}{m^3\; bed}]$ (independent of the particle size), bed porosity is $0.4$.
# 
# There are no external mass transfer limitations under all conditions of this problem. Assume an isothermal pellet.
# 
# **2a).** Consider $1.5 [mm]$ (radius) spherical catalyst particles. Find generalized Thiele modulus, internal effectiveness factor and observed reaction rate as $[\frac{mol}{m^2\;catalyst\; surface\; area\;*\;s}]$, $[\frac{mol}{m^3\; bed\;*\;s}]$ and as $[\frac{mol}{m^2\;external\; catalyst\; surface\; area\;*\;s}]$.
# 
# 
# **2b).** What should be the particle size to avoid internal mass transfer limitations?

# **Answer to Q2**
# 
# **2a)** There are no external mass transfer limitations under all conditions of this problem. In this case, for a $1^{st}$ order reaction for sphere, $\phi_n$ is (refer to Fogler, Elements of CRE, 5th Ed. page 732):
# 
# \begin{eqnarray}
# \phi_n=\frac{R}{3}(\frac{k}{D_{eff}})^{0.5}
# \end{eqnarray}
# 
# The Thiele modulus,$\phi_n$, will always contain a subscript (e.g., n), which refers to the reaction order. $\phi_n$ is a measure of the ratio of “a” surface reaction rate to “a” rate of diffusion through the catalyst pellet. Here $k$ is a catalyst-volume based constant (intrinsic).
# 
# \begin{eqnarray}
# r [\frac{mol}{m^3_{cat}.s}]=k*C_A[\frac{mol}{m^3\;fluid}]
# \end{eqnarray}
# 
# We are given $k"$ as $[\frac{m^3_{fl}}{m^2_{cat}.s}](=[m/s])$.
# 
# \begin{eqnarray}
# k[\frac{m^3_{fl}}{m^3_{cat}.s}]=k"[\frac{m^3_{fl}}{m^2_{cat}.s}]*SSA[\frac{m^2_{cat}}{g_{cat}}]*\rho_c [\frac{g_{cat}}{m^3_{cat}}]
# \end{eqnarray}
# 
# \begin{eqnarray}
# \rho_c=\frac{\rho_{bed}[\frac{g_{cat}}{m^3_{bed}}]}{1-\Phi}; \;\;\;\Phi= the\; porosity \;of\; the\; bed [\frac{m^3_{cat}}{m^3_{bed}}]
# \end{eqnarray}
# 
# 
# \begin{eqnarray}
# k=2.8*10^{-7}*100*\frac{1.2*10^6}{1-0.4}=56 [\frac{m^3_{fl}}{m^3_{cat}.s}]
# \end{eqnarray}
# 
# \begin{eqnarray}
# \phi_n=Thiele\; modulus=\frac{0.0015}{3}(\frac{56}{0.02*10^{-4}})^{0.5}=2.6
# \end{eqnarray}
# 
# \begin{eqnarray}
# \eta=\frac{observed\;rate}{intrinsic\; rate \; at \; surface\;conditions}=\frac{tanh \phi_n}{\phi_n}=0.37
# \end{eqnarray}
# 
# Since we do not have external limitations, surface conditions $C_{As},T_s=bulk\;conditions\;C_{Ab},T_b$.
# 
# 
# Intrinsic rate: rate from rate law (Assume an isothermal pellet)
# \begin{eqnarray}
# r"=k"*C_A
# \end{eqnarray}
# 
# Gas phase: 
# \begin{eqnarray}
# P*V=Z*N*R*T
# \end{eqnarray}
# \begin{eqnarray}
# C_A=\frac{N_A}{V}=\frac{P_A}{Z*R*T}=\frac{0.24*10^6}{0.8*8.314*500}=72 [\frac{mol}{m^3_{fl}}]
# \end{eqnarray}
# \begin{eqnarray}
# r"=k"*C_A=2.8*10^{-7}[\frac{m^3_{fl}}{m^2_{cat}.s}]*72 \frac{mol}{m^3_{fl}}=2*10^{-5}[\frac{mol}{m^2_{cat}*s}]
# \end{eqnarray}
# 
# Observed rate $=\eta*intrinsic\;rate=0.37*2*10^{-5}=0.74*10^{-5}[\frac{mol}{m^2_{cat}.s}]$
# 
# FYI: same refers to rate constant: observed rate constant$=\eta*intrinsic\;rate\;constant(\;no\; external\; MTL))$ or if there are external MTL/HTL then intrinsic rate constant at $T_s$.
# 
# Observed rate in $[\frac{mol}{m^3_{bed}*s}]$?
# 
# connect to bed volume:
# 
# $r[\frac{mol}{m^3_{bed}*s}]=r"[\frac{mol}{m^2_{cat}.s}]*SSA[\frac{m^2_{cat}}{g_{cat}}]*\phi_{bed}[\frac{g_{cat}}{m^3_{bed}}]=0.74*10^{-5}*100*1.2*10^6=888[\frac{mol}{m^3_{bed}.s}]$
# 
# In units of $[\frac{mol}{m^2_{external\;catalyst\;surface\;area}.s}]$?
# 
# 
# Let's use $m^3_{bed}$ and $a_c=[\frac{m^3_{ext\;surf\;area}}{m^3_{bed}}]$
# 
# \begin{eqnarray}
# a_c(sphere)=\frac{6(1-\phi)}{D_p}
# \end{eqnarray}
# 
# external area $\neq$ total (internal pores $+$ external)
# 
# \begin{eqnarray}
# r"[\frac{mol}{m^2_{ext\;area*s}}]=r[\frac{mol}{m^3_{bed}*s}]/a_c[\frac{m^2_{ext}}{m^3_{bed}}]=\frac{888*0.003}{6*(1-0.4)}=0.74[\frac{mol}{m^2_{external\;catalyst\;surface\;area}.s}]
# \end{eqnarray}

# In[2]:


import math
import numpy as np
T=500
Dp=0.003
R=8.314
z=0.8
kp=2.8*pow(10,-7)
SSA=100
phi=0.4
P_A=0.24*pow(10,6)
rho_bed=1.2*pow(10,6)
rho_c=rho_bed/(1-phi)
r=0.0015
D_eff=0.02*pow(10,-4)
k=kp*SSA*rho_c
phi_n=(r/3)*pow((k/D_eff),0.5)
eta=np.tanh(phi_n)/phi_n
CA=P_A/(z*R*T)
rp=kp*CA
Intr_rate=2*pow(10,-5)
robs=eta*Intr_rate
phi_bed=1.2*pow(10,6)
rate=robs*SSA*phi_bed
a_c=(6.*(1-phi))/Dp
rp2=rate/a_c
print("k={0:.0f}".format(k),'[m^3_{fl}/m^3_{cat}/s]')
print("robs={0:.6f}".format(robs),'[mol/m^2_{cat}/s]')
print("rate={0:.3f}".format(rate),'[mol/m^3_{bed}/s]')
print("phi_n={0:.3f}".format(phi_n))
print("eta={0:.3f}".format(eta))
print("rp2={0:.3f}".format(rp2),'[mol\m^2_{external catalyst surface area}/s]')


# **Q2-2b**.
# 
# To avoid internal MTL, $\eta>=0.95$
# 
# \begin{eqnarray}
# 0.95=\eta=\frac{tanh \Phi_n}{\Phi_n}
# \end{eqnarray}
# 
# $\Phi_n=0.4$
# 
# \begin{eqnarray}
# \Phi_n[sphere\;1^{st} \; order\;rxn]=\frac{R}{3}(\frac{k}{D_{eff}})^{0.5}
# \end{eqnarray}
# 
# $k$ and $D_{eff}$ do not depend on particle size (only on $T^0$)
# 
# $R=0.4*3*(\frac{0.02*10^{-4}}{56})^{0.5}=0.21 mm$
# 
# Particles of this size and smaller are free of internal MTL.

# In[12]:


import numpy as np
from scipy.optimize import fsolve

def f(x):
    return (np.tanh(x)/x)-0.95
    
x = fsolve(f,0.1)
f(x)

print('phi_n=',x)

phi=0.4
Deff=0.02*pow(10,-4)
k=56        #[m^3_{fl}/m^3_{cat}/s]
R=pow(phi*3.*(Deff/k),0.5)

print("R={0:.6f}".format(R), 'm')


# ---
