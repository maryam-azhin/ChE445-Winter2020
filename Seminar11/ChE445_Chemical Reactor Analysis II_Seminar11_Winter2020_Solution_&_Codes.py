#!/usr/bin/env python
# coding: utf-8

# ###### T.A. Maryam Azhin,
# 
# Department of Chemical and Materials Engineering, University of Alberta

# # Seminar 11. Evaluation of internal MTL and conversion in a PBR.

# **Q1. Evaluation of internal MTL and conversion in a PBR.**
# 
# Nitrous oxide ($N_2O$) has a greenhouse gas global warming potential that is almost $300$ times higher than that of carbon dioxide. Its catalytic reduction to harmless $N_2$ can be done using $CH_4$, $NH_3$ or $H_2$.
# 
# Consider gas-phase constant-density $N_2O$ hydrogenation on a $Pt/SiO_2$ catalyst that follows an intrinsic first order to $N_2O$ and apparent $0$ to hydrogen due its large excess. There are no external transfer limitations, assume negligible pressure drop. The catalyst pellet can be considered isothermal. Evaluate internal mass transfer limitations:
# $P=102000\; pa$, Pore diameter in individual catalyst particle, $3*10^{-9}m$, Particle porosity$=0.2$, Tortuosity in a particle$=4$, Intrinsic rate constant (based on $N_2O$), $1.7 \frac{m^3\;fluid}{(kg\; cat*s)}$ at $573\; K$, Catalyst density, $4000kg/m^3$, Reaction activation energy, $120000J/mol$, Volumetric velocity at STP $(0\; ^oC,\; 100000\; Pa),\; 0.8\;m^3/s$, Catalyst mass in PBR, $2kg$
# 
# **a). Assess the particle size effect:** 
# For catalyst particle diameters $= (3,6,12,24,48,96,192,384)\;\mu m$ calculate molecular diffusivity, Knudsen diffusivity of $N_2O$, pore diffusivity of $N_2O$, effective diffusivity of $N_2O$, intrinsic rate constant $[\frac{m^3_{fl}}{m^3_{cat}*s}]$, Thiele modulus, Internal effectiveness factor, Volumetric velosity at reaction conditions, $[m^3_{fl}/s]$, conversion and ideal PBR conversion,\%. 
# 
# write down formulas you used and show the units conversion, when necessary. Sketch a graph “X vs. particle size” for an ideal and real PBR and clearly mark final and initial points. On this graph, circle a part of the real PBR curve where there are no MTL (“kinetic regime”).
# 
# *Remember that internal effectiveness factor is just a coefficient in front of the intrinsic rate law based on external surface concentration.*
# 
# **b). Assess the temperature effect:** 
# For different temperatures $[457,493,533,573,613,653,693,733]$ and catalyst particle diameters of $48\;\mu m$ calculate molecular diffusivity, Knudsen diffusivity of $N_2O$, pore diffusivity of $N_2O$, effective diffusivity of $N_2O$, intrinsic rate constant $[\frac{m^3_{fl}}{m^3_{cat}*s}]$ and Intrinsic rate constant, $[\frac{m^3_{fl}}{m^3_{cat}*s}]$, Thiele modulus, Internal effectiveness factor, Volumetric velosity at reaction conditions, $[\frac{m^3_{fl}}{s}]$, conversion and ideal PBR conversion,\%.
# 
# Write down (an) additional formula(s) you used. Sketch a final “X vs. T” for an ideal and real PBR and clearly mark final and initial points. On this graph, circle a part of the real PBR curve where there are no MTL (“kinetic regime”).
# 
# At what temperatures (lower or higher) the internal MTL become more significant and why? 

# -----

# **Answer to Q1.**
# 
# **a). Molecular diffusivity**
# 
# \begin{eqnarray}
# D_{AB}[m^2/s]=\frac{1.013*10^{-2}*T^{1.75}*(\frac{1}{M_{N2O}}+\frac{1}{M_{H2}})^{0.5}}{P*(\nu_{N2O}^{1/3}+\nu_{H2}^{1/3})^2}
# \end{eqnarray}
# 
# $T=573K$
# 
# $M_A(N_2O)=44\;\frac{g}{mol}$
# 
# $M_B(H_2)=2\;\frac{g}{mol}$
# 
# $P=102000 pa$
# 
# $\nu_A(N_2O)=35.9$, diffusion volume for $N_2O$
# 
# $\nu_B(H_2)=7.07$, diffusion volume for $H_2$
# 
# Binary molecular diffusion in gases $[m2/s]$ can be calculated using the Fuller formula 
# 
# ($P$–pressure in $[Pa]$, $T$ – temperature in $[K]$, $M$ – molar mass in $[g/mol]$, $\nu_i$ is a diffusion volume.):
# 
# \begin{eqnarray}
# D_{AB}=\frac{1.013*10^{-2}*T^{1.75}*(\frac{1}{MA}+\frac{1}{MB})^{0.5}}{P*((\sum \nu_i)_A^{1/3}+(\sum \nu_i)_B^{1/3})^2}=\frac{1.013*10^{-2}*573^{1.75}(\frac{1}{44}+\frac{1}{2})^{0.5}}{102000*(35.9^{1/3}+7.07^{1/3})^2}\\
# D_{AB}=1.76*10^{-4}\frac{m^2}{s}
# \end{eqnarray}
# 
# Knudsen diffusivity $D_K$ for gas species varies with the pore diameter $d_p$ in $[m]$, $T$ the temperature in $[K]$, and the molar mass of the diffusion species, $M$, in $[g/mol]$. 
# 
# $d_p=3*10^{-9}m$
# 
# The equivalent pore diameter is taken as the average diameter of the pores in the pellet, which gives a reasonable approximation provided that there is a relatively narrow pore size distribution. 
# 
# The Knudsen diffusion coefficient, $D_k$ in $[m^2/s]$ in a straight round pore can be calculated from the following formula: ($M=44 \frac{g}{mol}$ molecular mass of the diffusing species.)
# 
# \begin{eqnarray}
# D_k=48.5*d_p*(\frac{T}{M})^{0.5}=48.5*3*10^{-9}*(\frac{573}{44})^{0.5}=5.25*10^{-7}\frac{m^2}{s}
# \end{eqnarray}
# 
# $D_{pore}$ is diffusivity in a pore, which is a combination of molecular and Knudsen diffusivity. In general, the concentration dependence of the pore diffusion coefficient is not large, and in many cases, the following equation (Bosanquet formula) is used to compute the pore diffusion coefficient.
# 
# \begin{eqnarray}
# D_{pore}=(\frac{1}{D_{AB}}+\frac{1}{D_K})^{-1}=5.24*10^{-7}\frac{m^2}{s}
# \end{eqnarray}
# 
# Therefore, it is governed mostly by $D_k$ because of small pore size (width).
# 
# Effective diffusivity in a catalyst for unimodal pore size distribution, $D_{eff}$:
# 
# ($\phi_p$ – pellet porosity, $\tau$ – tortuosity factor, typically $3-4$)
# 
# \begin{eqnarray}
# D_{eff}=\frac{\phi_p*D_{pore}}{\tau}=\frac{0.2*5.24*10^{-7}}{4}=2.618*10^{-8}\frac{m^2}{s}
# \end{eqnarray}
# 
# Intrinsic rate constant $k [\frac{m^3_{fl}}{m^3_{cat}.s}]=k' [\frac{m^3_{fl}}{kgcat.s}].\rho_c [\frac{kg}{m^3 cat}]=17*4000=6800$
# 
# In this case, for a $1^{st}$ order reaction for sphere, $\phi_n$ (Thile modulus) is:
# 
# \begin{eqnarray}
# \phi_n=\frac{d_{part.}}{6}(\frac{k}{D_{eff}})^{0.5}
# \end{eqnarray}
# 
# for example for particle diameter, $d_{part.}=3*10^{-6}$:
# \begin{eqnarray}
# \phi_n=\frac{3*10^{-6}}{6}(\frac{6800}{2.618*10^{-8}})^{0.5}=0.25
# \end{eqnarray}
# 
# \begin{eqnarray}
# \eta=\frac{observed\;rate}{intrinsic\; rate \; at \; surface\;conditions}=\frac{tanh \phi_n}{\phi_n}=0.98
# \end{eqnarray}
# 
# $\phi_n$ and $\eta$ for other particle diameters has been calculated using Python.
# 
# $Q$ at reaction conditions: $Q=Q_0$ (constant density)
# \begin{eqnarray}
# Q_{0,rxn}=\frac{T_{rxn}}{T_{STP}}\frac{P_{STP}}{P_{rxn}}.Q_{0,STP}=1.65 \;\;\frac{m^3}{s}
# \end{eqnarray}
# 
# **Conversion:**
# 
# **PBR MB:**
# 
# \begin{eqnarray}
# F_{A0}\frac{dX}{dW}=-r_A=\eta*k'*\frac{F_{A0}}{Q_{0,rxn}}(1-X)
# \end{eqnarray}
# 
# constant density
# \begin{eqnarray}
# \frac{dX}{1-X}=\frac{\eta*k'*dW}{Q_{0,rxn}}
# \end{eqnarray}
# here
# \begin{eqnarray}
# [k']=[\frac{m^3_{fl}}{kg_{cat}.s}]
# \end{eqnarray}
# \begin{eqnarray}
# ln(\frac{1}{1-X})=\frac{\eta*k'*W}{Q_{0,rxn}}
# \end{eqnarray}
# \begin{eqnarray}
# X=1-exp(\frac{-\eta*k'*W}{Q_{0,rxn}})
# \end{eqnarray}
# 
# Ideal PBR (No MTL): at $\eta=1$, $X=87.3\%$
# 
# based on the following Python calculation for $\eta=0.25$ (real $X$), $X=39.7\%$

# In[1]:


import math
import numpy as np
P=102000          #P𝑎 
T=573             #K
Pstp=pow(10,5)    #Pa
Tstp=273          #K
MA=44             #g/mol
MB=2              #g/mol
nuA=35.9          #diffusion volume N2O
nuB=7.07          #diffusion volume H2
phi_p=0.2         #pellet porosity
tau=4             #tortuosity factor
Q0stp=0.8         #m3/s
W=2               #catalyst mass, kg
rho_c=4000        #catalyst density
kp=1.7            #m3fl/(kgcat.s); Intrinsic kinetic rate
E=120000          #reaction activation energy, J/mol
dp= 3.*pow(10,-9) #m pore size
DAB=(1.013*pow(10,-2)*pow(T,1.75)*pow((1/MA+1/MB),0.5))/(P*pow(pow(nuA,1/3)+pow(nuB,1/3),2))
Dk=48.5*dp*pow((T/(MA)),0.5)
Dpore=pow((1/DAB+1/Dk),-1)
Deff=(phi_p*Dpore)/tau
k=kp*rho_c
Qrxn=(T*Pstp/(Tstp*P))*Q0stp

dpart=[3.*pow(10,-6),6.*pow(10,-6),12.*pow(10,-6),24.*pow(10,-6),48.*pow(10,-6),96.*pow(10,-6),192.*pow(10,-6),384.*pow(10,-6)]
phi_n=np.zeros(len(dpart))
eta=np.zeros(len(dpart))
X=np.zeros(len(dpart))
Xideal=np.zeros(len(dpart))

for i in range(0,len(dpart)):
    
     phi_n[i]=dpart[i]/6*pow((k/Deff),0.5) 
     eta[i]=np.tanh(phi_n[i])/phi_n[i]    
     X[i]=100*(1-math.exp(-1*eta[i]*kp*W/Qrxn))
     Xideal[i]=100*(1-math.exp(-1*kp*W/Qrxn))

print("DAB={0:.6f}".format(DAB),'bulk diffusivity [m2/s]')
print("Dk={0:.8f}".format(Dk),'Knudsen diffusivity [m2/s]')
print("Dpore={0:.8f}".format(Dpore),'m2/s')
print("Deff={0:.9f}".format(Deff),'m2/s')
print("k={0:.3f}".format(k))
print("phi_n",phi_n)  #print("phi_n={0:.2f}".format(phi_n),'m2/s')
print("eta=",eta)     #print("eta={0:.3f}".format(eta))
print("Qrxn={0:.3f}".format(Qrxn),'m3/s')
print("X=",X,'%')
print("Xideal=",Xideal,'%')


# In[2]:


import matplotlib.pyplot as plt

plt.plot(dpart,X,"-ob", label="Real PBR")
plt.plot(dpart,Xideal,"-or", label="Ideal PBR")
plt.legend(loc="lower right")
plt.ylabel('X')
plt.xlabel('Paricle diameter, m')


# **In an ideal PBR:** 
# 
# $X=87.3\%$ (X does not depend on particle size)
# 
# **In a real PBR reactor:**
# 
# at $dp=3 \mu m$, $X=87.3\%$
# 
# at $dp=384 \mu m$, $X=6\%$

# **2b).**
# All formulas are the same as in the section **2a** but now rate constant depends on $T^0$ as 
# \begin{eqnarray}
# k_T=exp(ln(k_{573}+\frac{E}{R}(\frac{1}{573}-\frac{1}{T})))
# \end{eqnarray}
# 
# Values are the same (same $T$ and $D_{particle}$)
# 
# Internal MTL are more significant at higher $T^0$ because $k$ is faster than $D_{eff}$ with $T^0$ increase.

# In[3]:


import math
import numpy as np
R=8.314       #J/mol/K
P=102000          #P𝑎 
#T=573             #K
Pstp=pow(10,5)    #Pa
Tstp=273          #K
MA=44             #g/mol
MB=2              #g/mol
nuA=35.9 
nuB=7.07
phi_p=0.2         #pellet porosity
tau=4             #tortuosity factor
Q0stp=0.8         #m3/s
W=2               #catalyst mass, kg
rho_c=4000        #catalyst density
kp=1.7            #m3fl/(kgcat.s); Intrinsic kinetic rate
E=120000          #reaction activation energy, J/mol
dp=3.*pow(10,-9)  #m pore size

dpart=48.*pow(10,-6)                 #m
T=[453,493,533,573,613,653,693,733]  #K

DAB=np.zeros(len(T))
Dk=np.zeros(len(T))
Dpore=np.zeros(len(T))
Deff=np.zeros(len(T))
kT=np.zeros(len(T))
kT2=np.zeros(len(T))
phi_n=np.zeros(len(T))
eta=np.zeros(len(T))
Qrxn=np.zeros(len(T))
X=np.zeros(len(T))
Xideal=np.zeros(len(T))
for i in range(0,len(T)):
    
     DAB[i]=(1.013*pow(10.,-2)*pow(T[i],1.75)*pow((1/MA+1/MB),0.5))/(P*pow(pow(nuA,1/3)+pow(nuB,1/3),2))
     Dk[i]=48.5*dp*pow((T[i]/(MA)),0.5)
     Dpore[i]=pow((1/DAB[i]+1/Dk[i]),-1)
     Deff[i]=(phi_p*Dpore[i])/tau
     kT[i]=math.exp(math.log(kp)+E/R*(1/573.-1/T[i])) #intrinsic rate constant
     kT2[i]=kT[i]*rho_c
     phi_n[i]=dpart/6*pow((kT2[i]/Deff[i]),0.5)
     eta[i]=np.tanh(phi_n[i])/phi_n[i]   
     Qrxn[i]=(T[i]*Pstp/(Tstp*P))*Q0stp
     X[i]=100*(1-math.exp(-1*eta[i]*kT[i]*W/Qrxn[i]))
     Xideal[i]=100*(1-math.exp(-1*kT[i]*W/Qrxn[i]))

print("DAB=",DAB,'bulk diffusivity [m2/s]')
print("Dk=",Dk,'Knudsen diffusivity [m2/s]')
print("Dpore=",Dpore,'m2/s')
print("Deff=",Deff,'m2/s')
print("kT", kT,'m^3 fl/(kg cat*s)')
print("kT2", kT2,'m^3 fl/(m3 cat*s)')
print("phi_n",phi_n,'m2/s')
print("eta=",eta)     
print("Qrxn=",Qrxn,'m3/s')
print("X=",X,'%')
print("Xideal=",Xideal,'%')


# In[4]:


import matplotlib.pyplot as plt

plt.plot(T,X,"-ob", label="Real PBR")
plt.plot(T,Xideal,"-or", label="Ideal PBR")
plt.legend(loc="upper left")
plt.ylabel('X')
plt.xlabel('Temperature, K')


# Kinetic regime: up to around $453\;K$ when $\eta >=0.95$
# 
# Internal MTL: are more significant at higher $T^0$ because $k$ increases faster than $D_{eff}$ with $T^0$ increase.
