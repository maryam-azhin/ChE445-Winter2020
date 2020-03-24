#!/usr/bin/env python
# coding: utf-8

# ###### T.A. Maryam Azhin, 
# 
# Department of Chemical and Materials Engineering, University of Alberta

# # Seminar 9. External mass transfer limitations

# **Q1.** Based on Example 14-4 from the textbook (Fogler's, $5^{th}$ Ed.)
# 
# An externally mass-transfer limited reaction (gas-phase reaction mixture) is carried out in an isothermal packed-bed reactor with no pressure drop at a constant volumetric flow rate. At $400 ^0C$, conversion of $85\%$ is achieved.
# 
# As a known rule of thumb, chemical reaction rate doubles at every $10-$degree temperature increase (the exact increase must be found from Arrhenius law). To ensure full conversion, it was proposed to increase temperature to $430 ^0C$. Analyze the suggestion.
# ![Diffusivity.png](Diffusivity.png)

# **Answer to Q1.**
# 
# **Case1:** $T=400 ^0C$, $X=0.85$
# 
# **Case2:** $T=430 ^0C$, $X=?$ 
# 
# For a constant density system, 
# \begin{equation}
# ln\frac{1}{1-X}=\frac{k_ca_c}{u}L\;\;(external\;MTL)
# \end{equation}
# for $T_{0,1}$ and $T_{0,2}$:
# 
# \begin{equation}
# (*)\frac{ln\frac{1}{1-X_2}}{ln\frac{1}{1-X_1}}=\frac{k_{c2}}{k_{c1}}\frac{u_1}{u_2} \;\;\;(same\;a_c \;and\;L)
# \end{equation}
# 
# $k_c$ and $u$ depends on $T^0$ for the same conditions (except $T^0$):
# $F_{T01}=F_{T02}$ and $P_0=constant$.
# 
# from ideal gas law: 
# 
# \begin{eqnarray}
# F_{T0}=Q_{01}\frac{P_0}{R*T_{0,1}}=Q_{02}\frac{P_0}{R*T_{0,2}}\\
# \end{eqnarray}
# 
# The pressure remains constant so
# \begin{equation}
# \frac{Q_{01}}{T_1}=\frac{Q_{02}}{T_2}
# \end{equation}
# 
# \begin{equation}
# Q_0=A_c[m^2]*u[m/s]
# \end{equation}
# 
# \begin{equation}
# \frac{u_1}{T_1}=\frac{u_2}{T_2}(same\;a_c)
# \end{equation}
# 
# \begin{equation}
# \frac{k_{c2}}{k_{c1}}=?
# \end{equation}
# 
# \begin{equation}
# k_c \propto \frac{u^{1/2}D_{AB}^{2/3}}{D_P^{1/2}\nu^{1/6}}
# \end{equation}
# 
# Taking the ratio of case 2 to case 1 and realizing that the particle diameter is the same for both cases gives us.
# 
# \begin{equation}
# \frac{k_{c2}}{k_{c1}}=(\frac{u_2}{u_1})^{1/2}(\frac{D_{AB2}}{D_{AB1}})^{2/3}(\frac{\nu_1}{\nu_2})^{1/6}
# \end{equation}
# 
# The gas-phase diffusivity is a function of temperature(from Table 11-2)
# \begin{equation}
# D_{AB}\propto T^{1.75}
# \end{equation}
# 
# For most gases, viscosity increases with increasing temperature according to the following relation
# 
# From the ideal gas law : $P*Q=F*R*T$, 
# 
# $\rho[\dot{m}/Q]\propto T^{-1}$
# 
# For gases we assume that the power-law is valid and $\mu \propto T^{2/3}$
# see power-low viscosity law: https://www.cfd-online.com/Wiki/Power-law_viscosity_law
# 
# So $\nu=\frac{\mu}{\rho} \propto T^{5/3}$
# 
# $\nu=\frac{\mu}{\rho} \propto T^{3/2}$
# 
# Combine into $(*)$:
# 
# \begin{equation}
# \frac{ln\frac{1}{1-X_2}}{ln\frac{1}{1-X_1}}=\frac{u_1}{u_2}\frac{k_{c2}}{k_{c1}}=\frac{u_1}{u_2}(\frac{u_2}{u_1})^{1/2}(\frac{D_{AB2}}{D_{AB1}})^{2/3}(\frac{\nu_1}{\nu_2})^{1/6}
# \end{equation}
# 
# \begin{equation}
# =(\frac{T_1}{T_2})^{0.5}*((\frac{T_2}{T_1})^{1.75})^{2/3}*((\frac{T_1}{T_2})^{5/3})^{1/6}=(\frac{T_2}{T_1})^{0.39}=(\frac{703}{673})^{0.39}=1.02
# \end{equation}
# 
# At $T_1=673\;K$, $Ln\frac{1}{1-X_1}=1.897$ where $X_1=85\%$
# 
# so $ln \frac{1}{1-X_2}=1.897*1.02=1.935$
# 
# $X_2=0.855$ practically didn't change with $30^0C$ increase of temperature.
# 
# **Analysis:** Consequantly, conversion in an externally mass-transfer limited (MTL) process is not (almost) affected by $T_0$ as opposed to the intrinsic reaction rate. 
# By increasing the temperature from $400 ^0C$ to $430 ^0C$ the conversion increases by only $2.19\%$ from $0.85$ to $0.87$.

# In[47]:


import math
X1=0.85
T1=273+400
T2=273+430
ratio=pow(T2/T1,0.39)
lninvx1=math.log(1/(1-X1))
X2=1-math.exp(-1*ratio*lninvx1)
print("X2={0:.3f}".format(X2))
print("𝐿𝑛(1/(1−𝑋1))=",lninvx1)
#print ("Increase in conversion by increasing temperature={0:.3f}".format(((0.869-0.85)/0.869)*100))


# **Q2.** Based on Example 14-2 from the textbook (Fogler's, $5^{th}$ Ed.)
# 
# Hydrazine decomposition to nitrogen and hydrogen $N_2H_4$. In a proposed study, a $2\%$ hydrazine in $98\%$ helium mixture is to be passed over a packed bed of cylindrical particles with The volume-average particel diameter of $3.6*10^{-3}$ at a gas-phase velocity of $15 m/s$ and a temperature of $750 K$. The kinematic viscosity of helium at this temperature is $4.5*10^{-4} m^2/s$.
# The hydrazine decomposition reaction is believed to be externally mass transfer-limited under these conditions. If the packed bed is $0.05 m$ in length, what conversion can be expected? Assume isothermal operation.
# 
# $D_{AB}=0.69*10^{-4}\; m^2/s\; at\; 298\;K$
# 
# Bed porosity: $30\%$
# 
# ![S9Q2.png](S9Q2.png)

# **Answer to Q2.** 
# 
# For a constant-density system when the process is limited by external MTL:
# \begin{equation}
# ln\frac{1}{1-X}=\frac{k_ca_c}{u}L
# \end{equation}
# 
# \begin{equation}
# a_c=\frac{6(1-\phi)}{D_p}=\frac{6*0.7}{3.6*10^{-3}}=1167[\frac{m^2}{m^3}]
# \end{equation}
# 
# \begin{equation}
# k_c=\frac{D_{AB}*Sh}{D_p}
# \end{equation}
# 
# $D_{AB}$ at $750$? See the table in Q1.
# 
# $D_{N2H4}$ (at $750)=D_{298}(\frac{750}{298})^{1.75}=3.47*10^{-4}\;\frac{m^2}{s}$
# 
# Frossling correlation:
# \begin{equation}
# Sh=2+0.6*Re^{1/2}Sc^{1/3}
# \end{equation}
# 
# \begin{equation}
# Re(p)=\frac{u_s*D_p*\rho_{fluid}}{\mu_{fluid}}=\frac{u_s*D_p}{\nu_{fluid}}=\frac{15*3.6*10^{-3}}{4.5*10^{-4}}=120
# \end{equation}
# 
# \begin{equation}
# Sc=\frac{\nu}{D_{AB}}=\frac{4.5*10^{-4}}{3.47*10^{-4}}=1.3
# \end{equation}
# 
# \begin{equation}
# Sh=2+0.6*(120)^{0.5}(1.3)^{1/3}=9.2
# \end{equation}
# 
# \begin{equation}
# k_c=\frac{Sh*D_{AB}}{D_P}=\frac{9.2*3.47*10^{-4}}{3.6*10^{-3}}=0.89 [m/s]
# \end{equation}
# 
# Then $X=1-exp[\frac{-k_c*a_c*L}{u}]=0.97$
# where $L=5\; cm$

# In[26]:


import math

T=750
Dp=3.6*pow(10,-3) #m
phi=0.3
ac=6*(1-phi)/Dp
L=0.05  #m
us=15
nu=4.5*pow(10,-4)   #m^2/s
Rep=us*Dp/nu
DAB=0.69*pow(10,-4)*pow(T/298,1.75)
Sc=nu/DAB
Sh=2+0.6*pow(Rep,0.5)*pow(Sc,1/3)
kc=Sh*DAB/Dp
X=1-math.exp(-kc*ac*L/us)

print ("DAB={0:.3f}".format(DAB))
print ("Re={0:.3f}".format(Rep))
print ("Sc={0:.3f}".format(Sc))
print ("Sh={0:.3f}".format(Sh))
print ("kc={0:.3f}".format(kc))
print ("X={0:.3f}".format(X))

