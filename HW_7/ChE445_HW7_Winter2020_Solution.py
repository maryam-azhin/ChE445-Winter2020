#!/usr/bin/env python
# coding: utf-8

# ##### T.A. Maryam Azhin, 
# 
# Department of Chemical and Materials Engineering, University of Alberta

# **HW7. External mass transfer limitations.**
# 
# due Monday, March 30, 2020. Weight: 2%

# Reaction:    $A + O_2 \rightarrow B$ (gas phase)
# 
# Kinetics:    Intrinsically $1^{st}$ order to $A$ and apparent $0^{th}$ order to $O_2$ (excess air and $O_2$)*
# 
# Catalyst:    Spherical, $3\; cm$ diameter, non-porous, $4\; kg$
# 
# Reactor:     tubular PBR, i.d. $0.05\; m$, no pressure drop 
# 
# Bed density and porosity: $500 kg/m^3$, $0.4$
# 
# Other properties at reaction temperature:
#   -Concentration of $A$ in the entering feed $5E-4 mol/m^3$
#   
#   -Diffusivity of $A$ in air $2E-5 m^2/s$
#   
#   -Kinematic viscosity of the fluid $5E-5 m2/s$
#   
#   -Intrinsic rate constant $8 m^3(fluid)/(m3(bed) s)$
#   
# External temperature gradients are negligible at all conditions in this assignment, but there is the possibility of external MTL and we have to find conditions to maximize the production of B. The system can be assumed to be constant density because of the excess oxygen (air).
# 
# **Question 1.** For volumetric flow rates (reaction conditions) of $0.0005$, $0.05$ and $15 m^3/s$, express the surface concentration of $A$ as a function of its bulk concentration. (hint: use $k$ and $kc$). 35 pts
# 
# **Question 2.** Among the three cases, select the flow rate where external MTL can be ignored and calculate $X_A$ and exit molar flow rate of product $B$ (use ideal PBR design equation, and use bed density to have the correct units for the rate constant to be used in the mole balance for the PBR). 10 pts
# 
# **Question 3.** Among these three cases, select the flow rate where external MTL is severe and external diffusion is rate-limiting. Again, calculate $X_A$ and $F_B$ considering that the observed rate is equal to the mass transfer rate. (For bed length, use bed density, catalyst mass and the tube Ac). 15 pts
# 
# **Question 4.** Calculate the conversion and $F_B$ for the third case. (Use PBR design equation with CAs but conversion is related to $C_{Ab}, not $C_{As}. Again, watch out for units of $k$ as in **Q2**). 25 pts
# 
# **Question 5.** Can you explain why the conversion in the cases affected by external MTL is higher than in the case where kinetics are rate determining? If your target were to produce the highest exit $F_B$, which
# case among the three would you prefer? 10 pts
# 
# **Question 6.** Given that diffusivities in liquid phase are lower than those in gas phase, which reactions are more prone to mass-transfer limitations: those in gas phase or in liquid phase? 5 pts
# 
# *Note: an “apparent” order means that the reaction rate at these particular concentrations (pressures) and temperature does not depend on the oxygen pressure. The reaction may still be, for example, an intrinsic
# first order to oxygen:
# 
# \begin{eqnarray}
# -r_A = k·C_A·P_{O2}
# \end{eqnarray}
# 
# But because of the excess oxygen, its partial pressure does not change significantly during the reaction, so it can be lumped into the rate constant (e.g. if there are $100$ times more moles of $O_2$ than of $A$, the pressure of oxygen will drop negligibly even when all $A$ is consumed). So, the rate law (only at these particular conditions) can be simplified to $-r_A = k*·C_A$.

# **Answer to Q1**
# \begin{eqnarray}
# C_{As}=\frac{k_ca_c}{k+k_ca_c}C_{Ab}, \;\;\;1^{st}\;order
# \end{eqnarray}
# 
# \begin{eqnarray}
# k_{c}=\frac{D_{AB}*Sh}{D_p}\\
# Sh=2+0.6*Re_p^{0.5}*Sc^{1/3}\\
# Sc=\frac{\nu}{D_{AB}}=\frac{5*10^{-5}}{2*10^{-5}}=2.5
# \end{eqnarray}
# 
# \begin{eqnarray}
# Re_p=\frac{u*D_p}{\nu}=\frac{Q*D_p}{A_c*\nu}=Q*\frac{0.03}{\pi*0.025^2*5*10^{-5}}=305578*Q
# \end{eqnarray}
# 
# \begin{eqnarray}
# a_c=\frac{6*(1-\phi)}{D_p}=120 \frac{m^2}{m^3_{bed}}
# \end{eqnarray}
# 
# Check units: 
# 
# $[k]=[\frac{m^3_{l}}{m^3_{bed}*s}]$
# 
# $[k_ca_c]=[\frac{m}{s}*\frac{m^2}{m^3_{bed}}]$, match, ok.
# 
# |        |   |      |         |           |                                        |
# |:---------:|-----:|-----:|-----:|----------:|---------------------------------------:|
# |$Q,\;m^3/s$|$Sc$|$Re_p$|$Sh$|$k_c,\;m/s$|$C_{As}=\frac{k_c*a_c}{k+k_c*a_c}*C_{Ab}$|
# |0.0005     |2.5|153   |12.068|$8.04*10^{-3}$   |$C_{As}=0.11*C_{Ab}\;\;\frac{mol}{m^3}$ |
# |0.05       |2.5|15279 |102.68|$6.84*10^{-2}$|$C_{As}=0.51*C_{Ab}\;\;\frac{mol}{m^3}$|
# |15|2.5|4583798|1745.87|1.16|$C_{As}=0.95*C_{Ab}\;\;\frac{mol}{m^3}$|

# In[1]:


import numpy as np

nu=5*pow(10,-5)   #m2/s
DAB=2*pow(10,-5)  #m2/s
Sc=nu/DAB
Dp=0.03           #m
CA0=5*pow(10,-4) #mol/m3
rhob=500           #kg/m3
phi=0.4
k=8              #m3fl/(m3cat.s)
W=4
Id=0.05          #m
Ac=3.14*pow(Id,2)/4
ac=6*(1-phi)/Dp    #1/m

print ("Sc=",Sc)

Q=[0.0005,0.05,15]
Rep=np.zeros(len(Q))
Sh=np.zeros(len(Q))
kc=np.zeros(len(Q))
CAs_coef=np.zeros(len(Q))

for i in range(0,len(Q)):
     Rep[i]=Q[i]*Dp/(nu*Ac)
     Sh[i]=2+0.6*pow(Rep[i],1/2)*pow(Sc,1/3)   
     kc[i]=Sh[i]*DAB/Dp
     CAs_coef[i]=kc[i]*ac/(k+kc[i]*ac)
        

print ("ac=",ac)
print ("Q=",Q)
print ("Rep=",Rep)
print ("Sh=",Sh)
print ("kc=",kc)
print ("CAs=",CAs_coef,'*CAb')


# ----

# **Answer to Q2**
# 
# When $C_{As}=~C_{Ab}$, there are no external MTL, so at $Q=15\;\frac{m^3}{s}$
# 
# **MB PBR:**
# 
# \begin{eqnarray}
# F_{A0}\frac{dX}{dW}=kC_A=k\frac{F_{A0}}{Q}(1-X)
# \end{eqnarray}
# constant $Q$:
# 
# \begin{eqnarray}
# \frac{dX}{1-X}=\frac{k}{Q}dW\\
# ln(\frac{1}{1-X})=\frac{k*W}{Q}
# \end{eqnarray}
# 
# units: from MB: $[\frac{mol}{s*kg_{cat}}]=[k*\frac{mol}{m^3_{fl}}]$
# 
# $[k]=[\frac{m^3_{fl}}{s*kg_{cat}}]$
# 
# $k$ is given as $8\;\;\frac{m^3_{fl}}{m^3_{bed}*s}$ so multiply by $\frac{1}{\rho_b[kg_{cat}/m^3_{bed}]}$
# 
# \begin{eqnarray}
# X=1-exp(\frac{-k*W}{\rho_b*Q})=1-exp(\frac{-8*4}{500*15})=0.4\;%
# \end{eqnarray}
# 
# \begin{eqnarray}
# F_B=F_{A0}*X_A=C_{A0}*Q*X=5*10^{-4}*15*\frac{0.4}{100}=3*10^{-5}\;\;\frac{mol}{s}
# \end{eqnarray}

# In[2]:


import math
X3=1-math.exp((-k*W)/(rhob*Q[2]))
FB3=CA0*Q[2]*X3 

print("X3={0:.4f}".format(X3),'%')
print ("FB3={0:.6f}".format(FB3),'mol/s')


# -----

# **Answer to Q3**
# 
# When $C_{As}=0$ ($<<C_{Ab}$) then there are severe external MTL. This is the case with $0.0005\frac{m^3}{s}$.
# 
# At constant $Q$
# 
# \begin{eqnarray}
# ln(\frac{1}{1-X})=\frac{k_ca_cL}{u}\\
# \frac{L}{u}=\frac{V_{bed}/A_c}{Q/A_c}=\frac{W/\rho_{bed}}{Q}\\
# ln(\frac{1}{1-X})=-8.04*10^{-3}*120*\frac{4/500}{0.0005}\\
# x=100%
# \end{eqnarray}
# 
# \begin{eqnarray}
# F_B=F_{A0}*X_A=C_{A0}*Q*X=5*10^{-4}*0.0005*1=2.5*10^{-7}\;\frac{mol}{s}
# \end{eqnarray}

# In[3]:


X1=1-math.exp(-kc[0]*ac*W/(rhob*Q[0]))
FB1=CA0*Q[0]*X1

print("X1={0:.2f}".format(X1),'%')
print ("FB1={0:.8f}".format(FB1),'mol/s')


# ----

# **Answer to Q4**
# 
# \begin{eqnarray}
# F_{A0}\frac{dX}{dW}=kC_{As}
# \end{eqnarray}
# for $Q=0.05\;m^3/s$:
# 
# \begin{eqnarray}
# C_{As}=C_{Ab}\frac{k_ca_c}{k+k_ca_c}=0.51 C_{Ab}\\
# C_{Ab}=\frac{F_{A0}}{Q}(1-X) 
# \end{eqnarray}
# at $Q$ constant:
# 
# \begin{eqnarray}
# ln(\frac{1}{1-X})=\frac{k*W*0.51}{Q}
# \end{eqnarray}
# 
# Here, $k$ has units $\frac{m^3_{fl}}{s*kg_{cat}}$
# 
# given $k=8 \frac{m^3_{fl}}{m^3_{bed}*s}*\frac{1[m^3_{bed}]}{500[kg_{cat}]}$
# 
# \begin{eqnarray}
# ln(\frac{1}{1-X})=\frac{8*4*0.51}{0.05*500}\\
# X=48
# \end{eqnarray}
# 
# $F_{B}=F_{A0}X=C_{A0}QX=5*10^{-4}*0.05*0.48=1.2*10^{-5} \frac{mol}{s}$
# 

# In[4]:


X2=1-math.exp(-k*W*CAs_coef[1]/(rhob*Q[1]))
FB2=CA0*Q[1]*X2

print("X2={0:.3f}".format(X2),'%')
print ("FB2={0:.7f}".format(FB2),'mol/s')


# ----

# **Answer to Q5**
# 
# **In the cases affected by MTL** 
# 
# External MTL is present when $Q$ and hence $F_{A0}$ are lower values, so the fluid spends more time in the reactor and higher conversion, $X$, is achieved. $k_ca_c<<k$, Slow diffusion and fast reaction. So observed reaction rate = rate of external diffusion.
# 
# **At high $Q$,** when there are no MTL, even at low $X$, $F_B$ is the highest.
# 
# $Q=15 m^3/s$ gives the highest $F_B$, $Q=0.05 m^3/s$ gives a slightly smaller value for $F_B$ while having a larger conversion $X=48%$ So case $Q=0.05 m^3/s$ is perfect to produce the highest exit $F_B$ yet higher conversion.

# ----

# **Answer to Q6**
# 
# Diffusivities in liquid phase are lower (because of higher density), so liquid phase reactions are more prone to MTL. Lower diffusivities give lower $Re$ and hence higher concentration, $C_A$, which means more MTL.
# 
# $kc \propto D_{AB}$ Therefore smaller $D_{AB}$ results in smaller $k_c$ or smaller mass transfer coefficient. The smaller mass transfer coefficient, the more prone the reaction is to mass transfer limitations. Since liquids have low diffusion rate, so they  have higher MTL .
