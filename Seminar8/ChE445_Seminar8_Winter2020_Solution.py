#!/usr/bin/env python
# coding: utf-8

# ###### T.A. Maryam Azhin, 
# 
# Department of Chemical and Materials Engineering, University of Alberta

# # Seminar 9.  Catalyst deactivation in a straight-through transport reactor STTR and in a fluidized CSTR.

# **Q1.**  A gas-phase constant-density isomerization reaction with a solid catalyst is carried out at $400 ^oC$ in a STTR: $A \rightarrow B$
# 
# The rate law is $–r_A = k.P_A^2.\rho$ with $k = 0.15 \frac{mol}{kgcat.s.atm^2}$ at $400 ^oC$. The STTR is $15 m$ high. Pure $A$ enters the reactor at $400 ^oC$ and $12 atm$.  The entering gas velocity is $7 \frac{m}{s}$. $C_{A_0} = 220 mol/m^3$. The fluidized catalyst density $\rho$ is $80 \frac{kgcat}{m^3}$. 
# 
# The catalyst is deactivated by coking. The decay law is $–da/dt = k_d.a.C_{coke}$ with $kd = 0.03 \frac{m^3}{(mol.s)}$ at $400 ^oC$. The coke concentration is described as $C_{coke}=\frac{P_B}{RT}$ with $R = 8.2*10^{-5} \frac{m^3.atm}{K.mol}$.
# 
# Calculate the achieved conversion with the deactivating catalyst and recommend an optimal reactor height.            

# **Answer to Q1:**
# 
# In the STTR the catalyst is fluidized and moves together with the feed at the same velocity. In an ideal case, it can be modeled as an ideal PFR (or, often, the non-idealities can be accounted by tanks-in-series model).
# 
# Although the catalyst deactivates, this is a steady-state reactor because we always maintain the flow of the same regenerated catalyst and at every point of the reactor the deactivation extent and reaction rate will not change with time on stream TOS (i.e., steady-state).

# M.B. 
# \begin{equation}
# \frac{dF_A}{dV}=r_A
# \end{equation}
# \begin{equation}
# dV=A_c.dZ
# \end{equation}
# \begin{equation}
# \frac{dF_A}{dZ}=r_A.A_c
# \end{equation}
# 
# Where $A_c$ is the reactor cross sectional area and $Z$ refers to the height of the reactor.
# 
# \begin{equation}
# X=\frac{F_A-F_{A0}}{F_{A0}}
# \end{equation}
# \begin{equation}
# F_A=F_{A0}-X*F_{A0}
# \end{equation}
# \begin{equation}
# \frac{dF_A}{dX}=-F_{A0}
# \end{equation}
# 
# \begin{equation}
# \frac{dF_A}{dZ}=r_A.A_c
# \end{equation}
# \begin{equation}
# F_{A0}\frac{dX}{dZ}=-r_A.A_c
# \end{equation}
# 
# \begin{equation}
# F_{A0}[mol/s]=u_0 [m/s].A_c [m^2].C_{A0}[mol/m^3]
# \end{equation}
# \begin{equation}
# \frac{dX}{dZ}=\frac{r_A}{U_0.C_{A0}}\;\;\; \textbf{(Eq.1)}
# \end{equation}
# 
# \begin{equation}
# –r_A [mol/(m^3.s)]= k [mol/(kg_{cat}.s.atm^2)].P_A^2 [atm^2].\rho [mol/(m^3)]
# \end{equation}
# 
# This is reaction rate with the ideally stable catalyst. To account for deactivation, a modified rate law is needed. The modified equation is as follows:
# 
# \begin{equation}
# –r_A = k.P_A^2.\rho .a\;\;\; \textbf{(Eq.2)}
# \end{equation}
# \begin{equation}
# –\frac{da}{dt} = k_d.a.C_{coke}\;\;a(t)\;is \;the  \; activity\; drop
# \end{equation}
# \begin{equation}
# C_{coke}=\frac{P_B}{RT}\;\;\; \textbf{(Eq.3)}
# \end{equation}
# 
# We need to get rid of TOS since $\textbf{(Eq.1)}$ has only $X$ and $Z$ as unknowns. $Z$ is connected to $t$, because as the catalyst moves in the reactor together with the feed it deactivates. 
# 
# At a specific time: $u [m/s]=Z [m]/t [s]$.
# \begin{equation}
# u=\frac{dZ}{dt}
# \end{equation}
# \begin{equation}
# –\frac{da}{dZ} = \frac{1}{u}k_d.a.C_{coke}\;\;\; \textbf{(Eq.4)}
# \end{equation}
# 
# To relate $P_A$ to $X$:
# 
# Stoichiometric Table for $A \rightarrow B$
# 
# |  Species   | Initial  | Final         |
# |:---:|---------:|--------------:|
# |A    | $F_{A0}$ | $F_A=F_{A0}(1-X)$ |
# |B    | 0        | $F_{B}=F_{A0}X$     |             
# |Total|$F_{T0}=F_{A0}$  | $F_T=F_{A0}=F_{T0}$|
# 
# \begin{equation}
# P_A=y_A.P_{total}=\frac{F_A}{F_T}P_0=P_0(1-X)\;\;\; \textbf{(Eq.5)}
# \end{equation}
# \begin{equation}
# P_B=y_B.P_{0}=\frac{F_B}{F_T}P_0=P_0X\;\;\; \textbf{(Eq.6)}
# \end{equation}
# \begin{equation}
# P_{total}=P_{0}=2 atm, constant
# \end{equation}

# In[1]:


from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

Z=np.linspace(0.,15,100)
X=np.zeros(len(Z))
a=np.zeros(len(Z))

def fun(F,Z):
    T   = 673.15          #K
    CA0 = 220.            #mol/m3
    k   = 0.15            #mol/(kgcat.s.atm^2)
    kd  = 0.03            #m3/(mol.s)
    P0  = 12.             #atm
    u0  = 7.              #m/s
    rho = 80.             #kgcat/m3
    R   = 8.2*pow(10.,-5) #(m3.atm)/(K.mol)
    
    PA=P0*(1-F[0])                      #Eq.5
    PB=P0*F[0]                          #Eq.6
    Ccoke=PB/(R*T)                      #Eq.3
    dadZ=-1/u0*kd*F[1]*Ccoke            #Eq.4
    dXdZ=k*pow(PA,2)*rho*F[1]/(u0*CA0)  #Eq.1 and Eq.2

    y=np.array([dXdZ,dadZ])
    return y

init_F = [0.,1.]
FS = odeint(fun,init_F,Z)

print ("max X={0:.3f}".format(FS[len(FS)-1,0]))
print ("a={0:.3f}".format(FS[len(FS)-1,1]))

plt.plot(Z,FS[:,0],"-b", label="X")
plt.plot(Z,FS[:,1],"-r", label="a")
plt.legend(loc="upper right")
plt.ylabel('X and a')
plt.xlabel('Z, m')


# As it it illustrated in the above plot, around 8 m catalyst deactivates and no more conversion can be achieved.  

# -----

# **Q2.** **Based on Problem 10-17 and Example 10-6 from the textbook (Fogler’s, 4th Ed.).**
# 
# Consider a first-order gas-phase cracking reaction on a solid catalyst in a fluidized CSTR at $700 K$: $Gas\; oil (g)\rightarrow Products (g)$ or $A\rightarrow B + C$. The feed contains $80\%$ of $A$ and $20\%$ inerts ($I$). The volumetric feed rate is $5000 m3/h$. There are $50$ metric tons of catalyst in the reactor with the bulk density of $500 \frac{kg}{m^3}$. $C_{A0} = 0.8 \frac{mol}{L}$, $C_{T0} = 1.0 \frac{mol}{L}$. At $700 K$, $k =\rho_bk’ = 45 h^{-1}$,  $k_d = 9 \frac{L}{mol.h}$. The gas oil contains sulfur compounds, which poison the catalyst. The rate of catalyst decay is first order in the present activity, and first order in the gas oil concentration. 
# 
# Plot the exiting reactant concentration and activity drop as a function of time on stream.

# **Answer to Q2.**
# 
# In an ideal fluidized CSTR the catalyst is fluidized and perfectly mixed with total fluid volume. It remains in the reactor during operation (as opposed to STTR in **Q1**) and because of its
# deactivation, the CSTR cannot be considered at steady-state:
# 
# **Stoichiometric Table** for $A \rightarrow B+C (+Inerts)$
# 
# |Species| Initial                    | Final                                                        |
# |:-----:|---------------------------:|-------------------------------------------------------------:|
# |A      | $F_{A0}$                   | $F_A=F_{A0}(1-X)$                                            |
# |B      | 0                          | $F_{B}=F_{A0}-F_A=F_{A0}X$                                   |    
# |C      | 0                          | $F_{C}=F_{A0}-F_A=F_{A0}X$                                   | 
# |I      |$F_{A0}/4$                  | $F_{I}=F_{I0}=F_{A0}/4$                                      | 
# |Total  |$F_{T0}=\frac{5}{4}F_{A0}$  | $F_T=\frac{9}{4}F_{A0}-F_A=\frac{5}{4}F_{A0}(1+\frac{4}{5}X)$|
# 
# \begin{equation}
# Q=\frac{Q_0}{F_{T0}}F_{T}=\frac{Q_0(\frac{9}{4}F_{A0}-F_A)}{\frac{5}{4}F_{A0}}
# \end{equation}
# Use $F_A=C_AQ;\;and \;\; F_{A0}=C_{A0}Q_0$ and then divide by $Q_0$.
# \begin{equation}
# \frac{5}{4}QC_{A0}=\frac{9}{4}C_{A0}Q_0-C_AQ
# \end{equation}
# \begin{equation}
# Q=\frac{\frac{9}{4}C_{A0}Q_0}{(\frac{5}{4}C_{A0}+C_A)}\;\;\; \textbf{(Eq.1)}
# \end{equation}
# 
# **Mole balance** on reactant:
# 
# flow in - flow out + rate of generation = rate of accumulation
# \begin{equation}
# Q_0C_{A0}-QC_A+r_A'W=\frac{dN_A}{dt}
# \end{equation}
# For constant volume:
# \begin{equation}
# N_A=C_AV\;\;and\;\;r_A'W=r_AV
# \end{equation}
# \begin{equation}
# Q_0C_{A0}-QC_A+r_AV=V\frac{dC_A}{dt}\;\;\; \textbf{(Eq.2)}
# \end{equation}
# \begin{equation}
# F_A=Q*C_A
# \end{equation}
# 
# **Rate law with decay:**
# 
# \begin{equation}
# -r_A=kC_Aa\;\;\; \textbf{(Eq.3)}
# \end{equation}
# 
# **Decay law:** 
# 
# \begin{equation}
# -\frac{da}{dt}=k_daC_A\;\;\; \textbf{(Eq.4)}
# \end{equation}
# 
# \begin{equation}
# a=\frac{rate\;during\;decay}{rate \;of\;fresh}
# \end{equation}
# 
# By introducing $\textbf{(Eq.1)}$ into $\textbf{(Eq.2)}$ :
# 
# \begin{equation}
# Q_0C_{A0}-\frac{\frac{9}{4}C_{A0}Q_0}{(\frac{5}{4}C_{A0}+C_A)}C_A-kC_AaV=V\frac{dC_A}{dt}\;\;\; \textbf{(Eq.5)}
# \end{equation}
# 
# For a gas phase if $P=P_0$ and $T=T_0$ we also have:
# $\epsilon=y_{A0}*\delta=y_{A0}*(1+1-1)=y_{A0}=\frac{C_{A0}}{C_{T0}}$
# 
# \begin{equation}
# Q/Q_0=\frac{F_T}{F_{T0}}=(1+\epsilon X)
# \end{equation}
# \begin{equation}
# X=1-\frac{F_A}{F_{A0}}=1-\frac{C_AQ}{C_{A0}Q_0}
# \end{equation}
# \begin{equation}
# X=1-\frac{\frac{9}{4}C_{A0}}{(\frac{5}{4}C_{A0}+C_A)}\frac{C_A}{C_{A0}} \;\;\; \textbf{(Eq.6)}
# \end{equation}
# 
# The following Python simulation shows that the catalyst deactivates very fast.

# In[2]:


from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

T=np.linspace(0.,0.5,100)
CA=np.zeros(len(T))
a=np.zeros(len(T))
X=np.zeros(len(T))

def fun(F,Z):
    k   = 45      #h-1
    kd  = 9       #L/mol.h
    CA0 = 0.8     #mol/L
    CT0  = 1.     #mol/L
    W=50000       #kg
    rhob=500      #kg/m3
    V=W/rhob      #m3  
    Q0=5000       #m3/h
        
    dadT=-kd*F[1]*F[0]                                               #Eq.4
    dCAdT=CA0*Q0/V-((9/4*CA0*Q0)/(5/4*CA0+F[0]))*F[0]/V-k*F[0]*F[1]  #Eq.5
         
    y=np.array([dCAdT,dadT])
    return y

init_F = [0.8,1.]
FS = odeint(fun,init_F,T)

CA0= 0.8
for i in range(len(T)):
    X[i]=1-((9./4*CA0*FS[i,0]/(5./4*CA0+FS[i,0])))*(FS[i,0]/CA0)     #Eq.6

print("CA={0:.3f}".format(FS[len(FS)-1,0]))
print("a={0:.3f}".format(FS[len(FS)-1,1]))
print("X={0:.3f}".format(X[len(FS)-1]))
         
plt.plot(T,FS[:,0],"-b", label="CA")
plt.plot(T,FS[:,1],"-r", label="a") 
plt.plot(T,X[:],"-g", label="X")
plt.legend(loc="upper right")
plt.ylabel('CA, X, and a')
plt.xlabel('t(h)')


# As it is illustrated in the above figure (X), there is some conversion while the catalyst is active. After t=0.02 hr, conversion derops as the catalyst deactivates. Moreover, the catalyst deactivates very fast and as it is expected the activity is decreasing over the time and around 0.5 hr it goes almost to zero. 
