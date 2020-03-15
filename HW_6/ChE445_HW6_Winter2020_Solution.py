#!/usr/bin/env python
# coding: utf-8

# ##### T.A. Maryam Azhin, 
# 
# Department of Chemical and Materials Engineering, University of Alberta

# **HW6. Active sites. Turnover frequency.**

# **Q 1. Metal supported catalysts: turnover frequency and reaction rates. 65 points**
# 
# Isomerization catalytic reaction: $A \rightarrow B$
# 
# **Catalyst #1:**    
# 
# $Pd/SiO_2$ (2 wt.% Pd), support specific surface area (SSA) $150 \frac{m^2}{g}$ $Pd$ dispersion is $28\%$
# 
# Rate law: $-r_A = \frac{k K_AC_A }{(1 + K_AC_A)}$
# 
# Activation energy: $60 \frac{kJ}{mol}$
# 
# Heat of adsorption of A: $-40 \frac{kJ}{mol}$
# 
# Rate constant $k$ at $600 K: 2 \frac{mol_A}{(kg_{catalyst}.s)}$
# 
# Adsorption constant of $A$ at $600 K: 0.6 \frac{L}{mol}$
# 
# a). Calculate TOF $[s^{-1}]$ at $500 K$ for $C_A = 0.85 \frac{mol}{L}$ (you can assume that the adsorption follows the Van’t Hoff equation). **25 points**
# 
# b). The catalyst producer offers you other $Pd$ supported catalysts with the following characteristics:
# 
# |Catalyst\# | Catalyst                | Support SSA, $[m2/ g]$|  Pd dispersion, \% |
# |:---------:|------------------------:|----------------------:|-------------------:|
# |2          | $Pd/SiO_2 (7 wt.\% Pd)$ | 150                   |  3\%               |
# |3          | $Pd/SiO_2 (7 wt.\% Pd)$ | 800                   | 40\%               |             
# |4          | $Pd/SiO_2 (2 wt.\% Pd)$ |  10                   |  8\%               |
# |5          |$Pd/SiO_2 (0.1 wt.\% Pd)$| 150                   | 70\%               |
# 
# 
# For each of the catalysts, calculate the reaction rates in $ mol_A/ (kg_{catalyst}.s)$ at $500 K$ for $C_A = 0.85 mol/L$. Show the general formula and report only the final rates for each of the $4$ catalysts listed above. **20 points**
# 
# c). Look at your results and suggest the trends:
# 
# For catalysts of the same surface area (catalysts $1$ and $2$), how does the metal loading affect metal dispersion? Can you suggest a reason for this trend? From the viewpoint of the calculated rates, is it preferable to have a catalyst with a higher or lower metal loading for the same support? **5 points**
# 
# Now include catalyst $5$ to your comparison of catalysts $1$ and $2$. How does your previous conclusion on the higher or lower preferable metal loading change? (The lesson here is that there is an optimal metal loading.) **5 points**
# 
# Now compare catalysts $2$ and $3$, and catalysts $1$ and $4$ (each pair has the same metal loading but different support SSA). From the viewpoint of high metal dispersion and reaction rates, is high SSA preferable or not? **5 points**
# 
# Draw conclusions: which factors affect metal dispersion; what characteristics affect the rate per gram of catalyst? **5 points**

# **Answer to Q1:**
# 
# **1a)** TOS is a way to express reaction rate, which depends on $k$ and $K$. Both constants change with $T$. Rate constant $k$ follows Arrhenius law: 
# 
# \begin{equation}
# k=k_0exp(\frac{-E}{RT})
# \end{equation}
# \begin{equation}
# ln(k_{500})=ln(k_0)-\frac{E}{R*500}
# \end{equation}
# \begin{equation}
# ln(k_{600})=ln(k_0)-\frac{E}{R*600}
# \end{equation}
# By subtracting two above equations we will have the following ralation :
# \begin{equation}
# ln(k_{500}/2)=-\frac{60000}{8.31}(\frac{1}{500}-\frac{1}{600})
# \end{equation}
# \begin{equation}
# k_{500}=0.18\frac{mol_A}{kg_{cat}.s}\;\;\;\;\;\;\;\;\;lower\;than\; k_{600}.
# \end{equation}
# 
# 
# Equilibrium constant $K_A$ follows Van't Hoff's law:
# \begin{equation}
# ln(\frac{K_{2}}{K_{1}})=-\frac{\Delta H_R}{R}(\frac{1}{T_2}-\frac{1}{T_1})
# \end{equation}
# \begin{equation}
# ln(\frac{K_{500}}{0.6})=\frac{40000}{8.31}(\frac{1}{500}-\frac{1}{600})
# \end{equation}
# \begin{equation}
# K_{A,500}=2.99 \frac{L}{mol}
# \end{equation}
# Higher than $k_{600}$ because chemisorption is exothermic ($\Delta H<0$ and $K$ increases by decreasing temperature)
# 
# Rate at $500K$:
# \begin{equation}
# r=\frac{kK_AC_A}{1+K_AC_A}=\frac{0.18*2.99*0.85}{1+2.99*0.85}=0.129 \frac{mol_A}{kg_{cat}.s}
# \end{equation}
# 
# find TOF $\frac{mol_A}{mol_{act.site}.s}$.
# \begin{equation}
# TOF=0.129 [\frac{mol_A}{kg_{cat}.s}]*\frac{1}{0.02}[\frac{kg_{cat}}{kg_{pd}}]*0.106[\frac{kg_{Pd}}{mol_{Pd}}]*\frac{1}{0.28}[\frac{mol_{Pd}}{mol_{act.sites}}]=2.44 [s^{-1}]
# \end{equation}
# 

# **1b)** Rate [$\frac{mol_A}{kg_{cat}.s}$]= TOF [$\frac{mol_A}{mol_{act.site}.s}$]*D[$\frac{mol_{cat.site}}{mol_{Pd}}$]*$\frac{1}{0.106}[\frac{mol_{Pd}}{kg_{Pd}}$]*$\frac{wt\%Pd}{100}\frac{kg_{Pd}}{kg_{cat}}$

# $TOF$ is the same for all catalysts because it reflects activity of each surface atom at the same $C_A$ and $T$

# Rates:
# 
# |Cat.\#|Pd wt\%|SSA|Pd D|Rate |
# |:----:|------:|--:|---:|----:|
# |1     |2      |150|28\%|0.129|
# |2     |7      |150|3\% |0.048|
# |3     |7      |800|40\%|0.644|
# |4     |2      |10 |8\% |0.037|
# |5     |0.1    |150|70\%|0.016|
# 
# **1c)**
# 
# For catalysts $1$ and $2$ (same SSA), higher loading decreases $D$ and rate. 
# At high loadings $Pd$ atoms are located very close to each other on a support and easier form large nanoparticles. Lower loading catalyst (at the same SSA) is preferable for reaction rate. 
# 
# But catalyst $5$, although it shows very high dispersion, the total low loading of Pd reduces the reaction rate per $kg_{cat}$. Thus, there is an optional loading.
# 
# For catalysts $2$ \& $3$ and $1$ \& $4$ high SSA allows for higher metal dispersion and higher rates.
# 
# Conclusions: 
# 
# $D$ is affected by metal loading and SSA of the support. 
# 
# Rates are affected by a combination of metal loading and $D$.

# -----

# **Q2. Zeolites. 10 points**
# 
# A mixture of 1-hexene and 4,4-dimethyl-1-hexene is isomerized catalytically using BAS on metal oxide catalysts. When alumina with an average pore size of $4 nm$ is used as a catalyst containing BAS, the isomerization products of both reactants are obtained. When a zeolite is used as a catalyst, 1-hexene isomerizes selectively, with 4,4-dimethyl-1-hexene remaining unconverted. Explain the reason. (Refer to the presentation “Zeolites.pptx”).

# 1-hexene is a linear molecule 
# ![1-hexene.png](1-hexene.png)
# 
# 4,4 - dimethyl-1-hexene is branched.
# ![4,4-dimethyl-1-hexene.png](4,4-dimethyl-1-hexene.png)
# 
# Both molecules can enter the $4nm$ pores of $Al_2O_3$, but most likely only 1-hexene can enter the zeolite pores, so (B) will not access active Bronsted acid sites (BAS).

# ----

# **Q3. Fischer-Tropsch reaction: a mechanistic insight. 25 pts**
# 
# This Youtube video shows a catalytic cycle of Fischer-Tropsch synthesis on Ru nanoparticles catalyst from $CO$ and $H_2$: “The Fischer-Tropsch reaction” by ICMS https://www.youtube.com/watch?v=44OU4JxEK4k
# It has two parts: in the first part, $Ru$ nanoparticles (composed of green Ru atoms) are small, while they are quite large in the second part. Hydrogen atoms are shown as white, oxygen as red, and carbon as black. Watch the video and follow the fate of the different molecules as they approach, adsorb, dissociate, react on the surface of Ru nanoparticles and desorb as products. (Note that the video encapsulates the results of detailed quantum chemical simulations, and can be taken as an accurate representation of what is truly happening at this scale).
# 
# a). Compare a small $Ru$ nanoparticle (at $22$ seconds) and a larger nanoparticle (at $24$ seconds). Which
# particle has a higher dispersion? **5 pts**
# 
# b). Are the metal nanoparticles porous, i.e. can reactants reach inside of the particles? **5 pts**
# 
# c). Consider the first part of the movie (up to $23 s$) for the small $Ru$ nanoparticle. Is the Fischer-Tropsch
# reaction (i.e., formation of a long-chain hydrocarbon) happening here? What are the main products in this
# part? **5 pts**
# 
# d). At approximately $26$ seconds into the video, there is a worm-like molecule starting to grow on the nanoparticle. What is the approximate chemical composition of this when it desorbs as a product? Is it a desired Fischer-Tropsch product? Note that the Fischer-Tropsch process is used to convert a mixture of carbon monoxide and hydrogen into liquid hydrocarbons to be used as synthetic fuels or lubrication oils or waxes. **5 pts**
# 
# e). Is it a Langmuir-Hinshelwood or Eley-Rideal mechanism? If adsorption occurs, is it dissociative or non-dissociative? **5 pts**

# **a)** Smaller particles (at $22$ s) have higher dispersion ($\frac{surface\; atoms}{total\; atoms}$).
# 
# **b)** No, reactants can access only the surface atoms on the nanoparticles.
# 
# **c)** No. Methane and water form.
# 
# **d)** $C_{10}H_{22}$ (approximately). Yes, it is a desired Fischer-Tropsch product.
# 
# **e)** Langmuir Hinshelwood with dissociative adsorption of both $CO$ and $H_2$
