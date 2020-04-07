#!/usr/bin/env python
# coding: utf-8

# **Q1. Internal, external and global mass transfer limitations at isothermal conditions. Review of TOF.**
# 
# (Fogler Ed5. Chp 15 & Hayes Ch 5)
# 
# A reaction rate was measured in a PBR without pressure drop at the following isothermal conditions:
#     
# |    |   |   |   |
# |:--:|--:|--:|--:|
# |Run#|Particle radius, mm|Total flow rate,mol/s|Observed rate,$mol/(g_{cat}.h)$|
# |1| 0.1| 600 |0.50|
# |2| 0.2|400 |0.50|
# |3| 0.2| 800| 0.49|
# |4|0.2| 1200 |0.51|
# |5|0.4| 400 | 0.26|
# |6 |0.4|800| 0.24|
# |7 |0.4|1200| 0.25|
# |8|4|400| 0.02|
# |9|4| 800| 0.05|
# |10|4| 1200| 0.05|
# 
# Hint: remember that these are experimental data subject to an experimental error.
# 
# a). What are the internal, external and global effectiveness factors for runs 1-4? Explain your choices and provide relevant calculations, where necessary. 
# 
# b). What are the internal, external and global effectiveness factors for runs 5-7? Explain your choices and provide relevant calculations, where necessary. 
# 
# c). What are the internal, external and global effectiveness factors for runs 9-10? Explain your choices and provide relevant calculations, where necessary. 
# 
# d). What are the internal, external and global effectiveness factors for run 8? Explain your choices and provide relevant calculations, where necessary. 
# 
# e). Calculate TOF if the catalyst contains $5 wt.\%$ Pd of $30\%$ dispersion. Pd molar mass is $106 g/mol$.
# 
# *(Hint: TOF is characteristic of an active site, it reflects intrinsic catalyst activity and must calculated based on the intrinsic rate law).*
# 
# 

# ---

# **Answer to Q1**
# 
# **a.)** For runs 1-4, variables are particle radius, $R_p$ and total flow rate, $F_T$.
# 
# but observed rate is alomost the same $0.50 \pm 0.01 \;\frac{mol}{gr}$ so this is an intrinsic kinetic regime, $\Omega = \eta_{ext}=\eta_{int}=1$

# ---

# **b.)** For runs 5-7 the total flow rate, $F_T$, does not have an effect on the reaction rate, so there are no external MTL $(\eta_{ext}=1)$ 
# 
# but the observed rate$=0.25<0.50$ in kinetic regime (for runs 1-4) , so there are internal MTL.
# 
# $\eta_{int}=\frac{Observed\;rate\;without\;external\;MTL}{Intrinsic\;rate}=\frac{0.25}{0.5}=0.5$ 
# 
# $$\Omega=Overall\;\;effectiveness\;\;factor$$
# 
# \begin{eqnarray}
# \Omega=\frac{Actual/observed\;overal\;rate\;of\;reaction}{Rate\;that\;would\;result\;if\;the\;entire\;surface\;were\;exposed\;to\;the\;bulk\;condition,\;C_{Ab},T_b}
# \end{eqnarray}
# \begin{eqnarray}
# \Omega=\eta_{int}*\eta_{ext}=0.5*1=0.5
# \end{eqnarray}

# ---

# **c.)** For runs 9-10, the total flow rate, $F_T$, does not have an effect on the reaction rate, so there are no external MTL $\eta_{ext}=1$.
# 
# $\eta_{int}=\frac{Observed\;rate\;without\;external\;MTL}{Intrinsic\;rate}=\frac{0.05}{0.50}=0.1$
# 
# $\Omega=\eta_{int}*\eta_{ext}=\eta_{int}*1=0.1$

# ---

# **d.)** For run 8:
# Where both external and internal limitation affect the rate.
# 
# From runs 9 and 10 for $R_{particle}=4mm$, $\eta_{int}=0.1$
# 
# \begin{eqnarray}
# \Omega=\frac{observed\;rate}{Intrinsic\;rate}=\frac{0.02[run 8]}{0.50[run\;1-4]}=0.04\\
# \Omega=\eta_{int}*\eta_{ext}\\
# \eta_{ext}=\frac{0.04}{0.1}=0.4
# \end{eqnarray}

# ---

# **e.)** TOF (turn over frequency) is intrinsic reaction rate as $[\frac{mol_A}{mol_{Active\;sites}*s}]$. 
# 
# Use runs $1-4$, rate $0.5$ to calculate TOF as follows:
# 
# \begin{eqnarray}
# TOF[\frac{mol_A}{mol_{active\;site}*s}]=0.5\frac{mol}{g_{cat}*h}*\frac{1}{0.05}\frac{g_{cat}}{g_{pd}}*106\frac{g_{Pd}}{mol_{Pd}}*\frac{1}{0.3}\frac{mol\;Pd}{mol\;surf(active\;site)\;Pd}*\frac{1}{3600}\frac{h}{s}=0.98s^{-1}
# \end{eqnarray}

# In[ ]:




