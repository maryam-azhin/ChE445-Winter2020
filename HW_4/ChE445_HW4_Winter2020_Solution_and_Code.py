#!/usr/bin/env python
# coding: utf-8

# T.A. Maryam Azhin, 
# Department of Chemical and Materials Engineering, University of Alberta

# **Multiple reactions in isothermal PBR with pressure drop, reactor
# descriptors**

# **D** is a target product of conversion of A via the multiple reactions:
# 
# $A\rightarrow 2B$ 
# 
# $3B\rightarrow D$
# 
# The rate laws are:
# 
# $-r_{A1}=k_{A1}C_{A}^{0.5}$
# 
# $-r_{B2}=k_{B2}C_{B}^2$
# 
# The rate constants (at the reaction temperature) are $k_{A1}= 75 kg_{cat}^{-1}mol^{0.5}m_{fluid}^{-0.5}s^{-1}$ , $k_{B2}= 4 m_{fluid}^2 kg_{cat}^{-1}mol^{-1}s^{-1}$.
# 
# A tubular reactor of $0.25\; m$ i.d., packed with $0.2\; kg$ of catalyst with a particle diameter of $125$ $\mu m$, catalyst density $2500$ $kg/m3$ and bed porosity of $0.35$ is used for the reaction. The entering pressure is $800\; kPa$. The feed is pure $A$ with entering molar flow rate 30 $mol/s$. The entering feed dynamic viscosity is $3*10^{-5}\;kg/(ms)$. The molar mass of $A$ is $50\;g/mol$, of $D$ is $40\;g/mol$. The reactor is maintained isothermal at $550\; K$. The gas phase can be assumed to be in ideal gas state. Assume that an ideal plug flow profile is developed in the reactor and no catalyst deactivation occurs. However, pressure drop must be accounted for.

# **Q1.**
# 
# **1a).** Build a model and write code to find molar flow rates of $A$, $B$, and $D$ and pressure drop “$y$” at the reactor exit. Include the calculation of particle Reynolds number in the program. 
# 
# In your submission, include your code and a plot of all molar flow rates vs catalyst weight ($W$). Also include a plot of the pressure drop vs. $W$ and report the overall pressure drop in the reactor. 
# **Each student has to set up their own program.**

# Mole Balance:
# 
# \begin{equation}
# \frac{dF_A}{dW}=r_A\\
# \frac{dF_B}{dW}=r_B\\
# \frac{dF_D}{dW}=r_D
# \end{equation}
# Reactions:
# 
# \begin{equation}
# k_{A1}= 75 \frac{mol^{0.5}}{kg_{cat}m_{fluid}^{0.5}s}\\
# k_{B2}= 4 \frac{m_{fluid}^2}{kg_{cat}mol.s}
# \end{equation}
# 
# Stoich.:
# \begin{equation}
# A\rightarrow 2B;\;\;r_{B1}=-2r_{A1}\\
# 3B\rightarrow D;\;\;r_{D2}=-\frac{1}{3}r_{B2}
# \end{equation}
# 
# \begin{equation}
# r_{A}=r_{A1}\\
# r_{B}=r_{B1}+r_{B2}\\
# r_{D}=r_{D2}
# \end{equation}
# 
# \begin{equation}
# C_{A}=C_{T0}\frac{F_{A}}{F_T}y\\
# C_{B}=C_{T0}\frac{F_{B}}{F_T}y
# \end{equation}
# 
# \begin{equation}
# C_{T0}=C_{A0}=\frac{P_0}{R*T}\\
# F_{T0}=F_{A0}\\
# F_T=F_A+F_B+F_D
# \end{equation}
# 
# \begin{equation}
# Ergun\; equation:\;\;\frac{dy}{dW}=\frac{-\alpha}{2y}\frac{F_T}{F_{T0}}\\
# At\;\;W=1:\;y=1\\
# \alpha=\frac{2\beta_0}{A_c(1-\Phi)\rho_cP_0}\\
# \beta_0=\frac{G(1-\Phi)}{\rho_0g_cD_p\Phi^3}[\frac{150(1-\Phi)\mu}{D_p}+1.75G]\\
# \end{equation}
# 
# \begin{equation}
# A_c=\pi*r^2\;\; [for\; reactor]\\
# Q_0=\frac{F_{T0}}{C_{T0}}=\frac{F_{T0}*R*T}{P_0}=\frac{30[mol/s]*8.314[J/(mol.K)]*550 [K]}{800000 [Pa]}=0.1715 [m^3/s]\\
# \dot{m}=F_{A0}.MW_A\\
# G=\frac{\dot{m}}{A_c}
# \end{equation}
# 
# Particle Re:
# \begin{equation}
# \rho_0[feed]=\frac{\dot{m}}{Q_0}=\frac{F_{A0}MW_A}{Q_0}=\frac{30 [mol/s]*50[g/mol]}{0.1715 [m^3/s]}=8.7475\\
# Re_p=\frac{D_p\rho_0u_s}{\mu}=\frac{125.e-6*8.746*(\frac{0.1715*4}{3.14*(0.25)^2})}{3.e-5}=127.3\\
# u_s=\frac{Q_0}{A_c}=\frac{0.1715 [m^3/s]*4}{3.14*(0.25)^2 [m^2]}=3.49[m/s]\\
# based\; on\; attached\; python\; code\;\;\; y=0.8778\\
# Overall\;Presuure\;Drop=P0-y*P0=800000-0.8778*800000=97775.384 [Pa]
# \end{equation}
# 

# In[1]:


#Q1.a
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

W0=np.linspace(0.,0.2,100)
W=W0
FA=np.zeros(len(W))
FB=np.zeros(len(W))
FD=np.zeros(len(W))
Y=np.zeros(len(W))

def fun(F,W):
    Dr=0.25
    Ac=(np.pi)*(pow(Dr,2))/4
    Wcat=0.2
    Dp=125.e-6
    Rhoc=2500.
    Phi=0.35
    P0=800000.
    FA0=30.
    FT0=FA0
    Mu=0.00003
    MWA=50. #0.050
    MWD=40. #0.040
    T=550.
    R=8.314
    CT0=P0/(R*T)
    Q0=FT0/CT0
    m=FA0*MWA/1000.
    Rho0=m/Q0
    Us=Q0/Ac
    G=m/Ac#Rho0*Us#
    Q0=FT0/CT0
    Rhofeed=m/Q0
    Reb=Dp*Rhofeed*Us/Mu
    FT=F[0]+F[1]+F[2]
    CA0=P0/(R*T)
    CT0=CA0
    CA=CT0*(F[0]*F[3]/FT)
    CB=CT0*(F[1]*F[3]/FT)
    CD=CT0*(F[2]*F[3]/FT)
    Ac=(np.pi/4.)*(pow(Dr,2.))
    beta=(G*(1.-Phi))/(Rho0*Dp*(pow(Phi,3)))*(((150.*(1-Phi)*Mu)/Dp)+(1.75*G))
    alpha=2.*beta/(Ac*(1.-Phi)*Rhoc*P0)
       
    kA1=75.*Ac
    kB2=4.*pow(Ac,2)
    
    dFAdW=-kA1*pow(CA,0.5)
    dFBdW=2.*kA1*pow(CA,0.5)-kB2*pow(CB,2)
    dFDdW=1./3*kB2*pow(CB,2)
    dYdW=(-alpha*FT)/(2.*F[3]*FT0)
    y=np.array([dFAdW,dFBdW,dFDdW,dYdW])
    return y

init_F= [30.,0.,0.,1.]
FS=odeint(fun,init_F,W)

print ("FA={0:.3f}".format(FS[len(FS)-1,0]),'mol/s')
print ("FB={0:.3f}".format(FS[len(FS)-1,1]),'mol/s')
print ("FD={0:.3f}".format(FS[len(FS)-1,2]),'mol/s')

plt.plot(W,FS[:,0],"-b", label="FA")
plt.plot(W,FS[:,1],"-r", label="FB")
plt.plot(W,FS[:,2],"-g", label="FD")
plt.legend(loc="upper left")
plt.ylabel('F, mol/s')
plt.xlabel('W, kg_cat')


# In[2]:


import matplotlib.pyplot as plt
plt.plot(W,FS[:,3])
plt.ylabel('Y')
plt.xlabel('W, kg_cat')


# In[11]:


Dr=0.25
Ac=(np.pi)*(pow(Dr,2))/4
Dp=125.e-6
P0=800000.
FA0=30.
FT0=FA0
Mu=0.00003
MWA=50. #0.050
MWD=40. #0.040
T=550.
R=8.314
CT0=P0/(R*T)
Q0=FT0/CT0
m=FA0*MWA/1000.
Rho0=m/Q0
Us=Q0/Ac
Rhofeed=m/Q0
Rhoc=2500.
Phi=0.35
Wcat=0.2
MWD=40. #0.040
Reb=Dp*Rhofeed*Us/Mu
print ("Re={0:.3f}".format(Reb))
    # Overall_P_Drop
Overall_P_Drop=P0-FS[len(FS)-1,3]*P0
print ("Overall_P_Drop={0:.3f}".format(Overall_P_Drop),'Pa')
Rhofeed


# **1b).** Calculate WHSV $[1/h]$, GHSV $[1/h]$, WTYD in $[kg_D/kg_{cat}h]$, $X_A$ and integral selectivity to $D$.

# \begin{equation}
# \dot{m}_T=\dot{m}_{A0}=F_{A0}*MW_A=30[mol/s]*0.050[kg/mol]=1.5 [kg/s]\\
# WHSV=\frac{\dot{m}_T}{W_{cat}}=\frac{1.5[kg/s]}{0.2 [kg]}*\frac{3600 [s]}{1[hr]}\\
# WHSV=\frac{1.5}{0.2}*\frac{3600 s}{hr}=27000 \frac{1}{hr}\\
# Q_0(STP)=\frac{F_{T0}RT_{STP}}{P_{STP}}=\frac{30*8.314*273.15}{10^5}=0.68 \frac{m^3}{s}\\
# V_{cat.bed}=\frac{m_{bed}}{\rho_{bed}}=\frac{0.2 kg}{2500*(1-0.35)}=0.00012 \;m^3\\
# GHSV=\frac{Q_0(STP)}{V_{cat.bed}}=\frac{0.68 [m^3/s]}{0.00012 [m^3]}=1.9929.e7 [1/hr]\\
# F_{T0}=F_{A0}=30 \frac{mol}{s} \;does\; not\; depend\; on \; P\;\&\;T\\
# GHSV=\frac{Q_0}{V_{bed}}*3600\frac{s}{h}=\frac{0.68}{0.00012}*3600=20400000 \frac{1}{h}\\
# \dot{m}_D=F_{D,exit}*MW_D=1.02\frac{mol}{s}*0.04\frac{kg}{mol}=0.0408\frac{kg}{s}\\
# WTY_D=\frac{\dot{m}_D}{W_{cat}}=\frac{0.0408 \frac{kg}{s}*3600 \frac{s}{hr}}{0.2 kg_{cat}}=734.4 \frac{kg_D}{kg_{cat}.hr}\\
# X_A=\frac{F_{A0}-F_A}{F_{A0}}=1-\frac{21.74}{30}=0.275\\
# F_D \;at\;exit=1.02 [mol/s]\;based\;on\;the\;attahced\;code\\
# S_{D/A}=\frac{mol_D\;formed}{mol_A\;consumed}=\frac{1.02}{30-21.736}=12.3\%\\
# S_D=\frac{F_D}{F_B+F_D}=\frac{1.020 [mol/s]}{13.470+1.020 [mol/s]}=0.07
# \end{equation}
# 

# In[4]:


#Q1.b
WHSV=(m/0.2)*3600
print ("WHSV={0:.1f}".format(WHSV),'1/h')
Q0STP=R*273.15*30/100000
Rhob= Rhoc*(1.-Phi)
Vbed=Wcat/Rhob
GHSV=(Q0STP/Vbed)*3600
print ("GHSV={0:.3f}".format(GHSV))
Wprod=FS[len(FS)-1,2]*MWD/1000.
WTY=(Wprod/Wcat)*3600
print ("WTY={0:.3f}".format(WTY),'kg_D/kg_cat.hr')
X_A=(MWA-FS[len(FS)-1,0])/MWA
print ("X_A={0:.3f}".format(X_A))
S_D=FS[len(FS)-1,2]/(FS[len(FS)-1,1]+FS[len(FS)-1,2])
print ("S_D={0:.3f}".format(S_D))


# **1c).** Based on the particle Reynolds number, is our assumption of ideal plug flow profile reasonable? (assuming no channeling or bypassing).

# Particle $Re=$ corresponds to turbulent flow and ideal PBR conditions might be satisfied.
# 127.32>100
# Assumption of ideal plug flow is not reasonable due to transition flow regime.

# **Q2.**
# 
# **2a).** To minimize pressure drop, the catalyst was pelletized into spheres of $2\; cm$ diameter. Repeat your calculations for this case (all other parameters and conditions remain the same) and report the molar flow rates and pressure drop at the outlet of the reactor.
# 
# **2b).** Calculate $WTY$ of the target product $D$ for the pelletized catalyst. Is there an improvement compared to the powered catalyst of Q1?

# $D_p$=0.02cm
# Based on the attached code
# \begin{equation}
# FA=21.5[mol/s]\\
# FB=13.3[mol/s]\\
# F_D=1.2[mol/s]\\
# \Delta P{overall}=P_0*(1-y)=800000*(1-0.9995)=400[kPa]\\
# WTY_D=\frac{\dot{m}_D}{\dot{catalyst}}=\frac{1.22*40*10^{-3}*3600}{0.2}=878.5[kg_D/(kg_{cat}.hr)]
# \end{equation}
# $WTY_D$ is increased compare to Q1.b. Therefore by increasing the size of the catalyst pellets the result is improved.

# In[10]:


#Q2.a
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

W0=np.linspace(0.,0.2,100)
W=W0
FA=np.zeros(len(W))
FB=np.zeros(len(W))
FD=np.zeros(len(W))
Y=np.zeros(len(W))

def fun(F,W):
    Dr=0.25
    Ac=(np.pi)*(pow(Dr,2))/4
    Wcat=0.2
    Dp=0.02
    Rhoc=2500.
    Phi=0.35
    P0=800000.
    FA0=30.
    FT0=FA0
    Mu=0.00003
    MWA=50. #0.050
    MWD=40. #0.040
    T=550.
    R=8.314
    CT0=P0/(R*T)
    Q0=FT0/CT0
    m=FA0*MWA/1000.
    Rho0=m/Q0
    Us=Q0/Ac
    G=m/Ac#Rho0*Us#
    Q0=FT0/CT0
    Rhofeed=m/Q0
    Reb=Dp*Rhofeed*Us/Mu
    FT=F[0]+F[1]+F[2]
    CA0=P0/(R*T)
    CT0=CA0
    CA=CT0*(F[0]*F[3]/FT)
    CB=CT0*(F[1]*F[3]/FT)
    CD=CT0*(F[2]*F[3]/FT)
    Ac=(np.pi/4.)*(pow(Dr,2.))
    beta=(G*(1.-Phi))/(Rho0*Dp*(pow(Phi,3)))*(((150.*(1-Phi)*Mu)/Dp)+(1.75*G))
    alpha=2.*beta/(Ac*(1.-Phi)*Rhoc*P0)
       
    kA1=75.*Ac
    kB2=4.*pow(Ac,2)
    
    dFAdW=-kA1*pow(CA,0.5)
    dFBdW=2.*kA1*pow(CA,0.5)-kB2*pow(CB,2)
    dFDdW=1./3*kB2*pow(CB,2)
    dYdW=(-alpha*FT)/(2.*F[3]*FT0)
    y=np.array([dFAdW,dFBdW,dFDdW,dYdW])
    return y

init_F= [30.,0.,0.,1.]
FS2=odeint(fun,init_F,W)

print ("FA={0:.3f}".format(FS2[len(FS2)-1,0]),'mol/s')
print ("FB={0:.3f}".format(FS2[len(FS2)-1,1]),'mol/s')
print ("FD={0:.3f}".format(FS2[len(FS2)-1,2]),'mol/s')
print ("y={0:.3f}".format(FS2[len(FS2)-1,3]))


# In[6]:


Dp=0.02
Overall_P_Drop=P0-FS2[len(FS2)-1,3]*P0
print ("Overall_P_Drop={0:.3f}".format(Overall_P_Drop),'Pa')


# In[7]:


#Q2.b
Wprod=FS2[len(FS2)-1,2]*MWD/1000.
WTY=(Wprod/Wcat)*3600
print ("WTY={0:.3f}".format(WTY),'kg_D/kg_cat.hr')

