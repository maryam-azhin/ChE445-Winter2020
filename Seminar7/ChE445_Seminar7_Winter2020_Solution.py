#!/usr/bin/env python
# coding: utf-8

# T.A. Maryam Azhin, Department of Chemical and Materials Engineering, University of Alberta

# **Seminar 7. Turnover frequency.**

# **Q1. Supported metal catalyst: turnover frequency.**
# 
# An irreversible reaction ${A + 2B}\rightarrow 2C + D$
# is carried out on a $Pt/Al_2O_3 (2 wt.\% Pt)$ catalyst. $Pt$ dispersion is $25\%$. Initial reaction rate was found from a differential reactor data as $3\frac{molC}{kgcatalyst.min}$. Find TOF with respect to $A$ in $s^{-1}$. 

# Catalytic reaction: ${A + 2B}\rightarrow 2C + D$
# 
# $\frac{r_A}{-1}=\frac{r_C}{2}$ -> $-r_A=\frac{r_C}{2}=\frac{1}{2}\frac{3 molC}{kgcatalyst.min}=1.5 \frac{mol_A}{kg_{cat}.min}$
# 
# PBR: $\frac{dF_A}{dW}=r_A'\frac{mol_A}{kg_{cat}.min}$
# 
# Differential PBR ($X<5\%$): $\frac{\Delta F_A}{W}=r_A'$ or $\frac{\Delta F_A}{mol_{active,sites}}=\frac{mol_A}{mol_{active,sites}.s}=TOF$
# 
# Where an active site=1 surface Pt atom
# 
# The turnover number (or frequency) is defined as the number of molecules that react per active site per unit of time. 
# \begin{eqnarray}
# -r_A=f_{Pt}*D*(\frac{1}{MW_{Pt}})*\frac{\%Pt}{100}\\
# D=dispersion=\frac{mol_{surface,Pt}}{mol_{Pt}}
# \end{eqnarray}
# \begin{eqnarray}
# TOF=f_{Pt}=1.5\frac{mol_A}{kg_{cat}.min}*\frac{1}{0.02}\frac{kg_{cat}}
# {kg_{Pt}}*0.195\frac{kg_{Pt}}{mol_{Pt}}*\frac{1}{0.25}\frac{mol_{Pt}}{mol_{Surface,Pt}}*\frac{1}{60}\frac{min}{s}=0.975s^{-1}
# \end{eqnarray}
# 
# 
# Note: The differential form of the mole balance on PBRs are applied when there is a pressure drop of catalys decay. In the absence of $\Delta P$ or catalyst decay we can integrate the differential form.

# **Q2. (Based on Q1).**
# 
# $80\%$ of the active sites of the catalyst from Q1 became poisoned by sulfur impurities in the feed. What initial reaction rate with respect to $A$ can be expected $[mol_A/(kg_{catalyst}.min)]$ at all other parameters kept constant? 

# TOF of unpoisoned sited remains the same but now there are only $20\%$ of sites left.
# \begin{eqnarray}
# -r_A=f_{Pt}*D*(\frac{1}{MW_{Pt}})*\frac{\%Pt}{100}\\
# f_{Pt} = TOF
# \end{eqnarray}
# 
# \begin{eqnarray}
# -r_A\frac{mol_A}{kg_{catalyst}.min}=0.975\frac{mol_A}{mol_{Surface,Pt}(unpoisoned).s}*0.2\frac{mol_{unpoisoned-sites}}{mol_{all-sites}}*...\\
# ...*0.25\frac{mol_{site}}{mol_{Pt}}*\frac{1}{0.195}*\frac{molPt}{kgPt}*0.02\frac{kgPt}{kgCat}*60\frac{s}{min}\\
# =0.3 \frac{mol_A}{kg_{cat}.min}
# \end{eqnarray}

# Easier solution: 
# 
# $1 kg$ of catalyst has $5$ times less active sites, so rate per kg of catalyst will be $5$ times lower: $\frac{1.5}{5}=0.3$

# **Q3. Zeolites.**
# 
# A zeolite is to be used to catalyze cracking of A with Bronsted acid sites (BAS) being active sites. The zeolite specific surface area (SSA) is $200 m^2/g$, the external surface area is $10 m^2/g$. The BAS density is $2 \frac{H^+}{nm^2_{zeolite,surface}}$. Based on the $TOF$ of $5 s^{-1}$, calculate how much catalyst is required to convert $10 \frac{mol_A}{s}$ in a differential PBR.

# Differential PBR: 
# 
# \begin{equation}
# TOF=\frac{\Delta F_A}{mol_{active,sites}}
# \end{equation}
# 
# Mole of active sites? $(H^+)$
# 
# \begin{equation}
# 2\frac{H^+}{nm^2_{zeolite,surface}}=2*10^{18}\frac{H^+}{m^2_{zeolite,surface}}*\frac{1 mole_{sites}}{N_A(Avogadro's No._{sites})}
# =\frac{2*10^{18}}{6.02*10^{23}}\frac{mol_{sites}}{m^2_{zeolite,surface}}
# \end{equation}
# 
# Mol of active sites required:
# 
# \begin{equation}
# \frac{\Delta F_A}{TOF}=\frac{10}{5}=2
# \end{equation}
# 
# Zeolite surface required:
# \begin{equation}
# \frac{2 mole_{sites}*6.02*10^{23}}{2*10^{18}}=602000 m^2
# \end{equation}
# 
# Mass catalyst required:
# \begin{equation}
# \frac{602000 m^2}{200 m^2/g}=3010 g
# \end{equation}

# **Q4. (Based on Q3).**
# 
# When the catalyst was put into the operation, only $0.5 \frac{mol_A}{s}$ were converted. Suggest a reason assuming no catalyst deactivation and ideal PBR behavior with no pressure drop. 

# Zeolites are crystalline porous materials and their pore size can be smaller than the reactant size, so the reactant cannot access active sites inside the pores. Only active sites on external surface will be accessible indeed.
# 
# \begin{equation}
# \frac{External\; surface}{Total\; surface}=\frac{10}{200}=5\%
# \end{equation}
# 
# The rate is $\frac{0.5 \frac{mol_A}{s}}{10 \frac{mol_A}{s}(from\; Q3)}=5\%$
# 
# So this may explain the observed behavior. To confirm, one needs to know the pore size of the zeolite and size of the reactant molecule.
