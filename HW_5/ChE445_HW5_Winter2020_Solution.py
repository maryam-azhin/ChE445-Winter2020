#!/usr/bin/env python
# coding: utf-8

# ###### T.A. Maryam Azhin, 
# 
# Department of Chemical and Materials Engineering, University of Alberta

# ### HW Assignment 5. Mechanisms of heterogeneous catalytic reactions**

# Gas-phase hydrogenation of $A$ to $B$ occurs on a solid catalyst. Two mechanisms were proposed:

# **Langmuir-Hinshelwood (LH):**
# 
# \begin{eqnarray}
# A+S<=>AS\\
# H_2 + 2S <=> 2(HS)\\
# AS + 2(HS) \rightarrow BS + 2S\\
# BS <=> B + S
# \end{eqnarray}
# 
# **Eley-Rideal (ER):**
# 
# \begin{eqnarray}
# H_2 + 2S <=> 2(HS)\\
# A + 2(HS) \rightarrow B + 2S
# \end{eqnarray}
# 
# **1a.** Derive the rate expression for $A$ for the Langmuir-Hinshelwood model. The answer must contain only measurable concentrations of $[A]$, $[H_2]$ and $[B]$, and constants. **25 pts**
# 
# **1b.** Derive the rate expression for $A$ for the Eley-Rideal model (here, with the surface reaction being an RDS). The answer must contain only measurable concentrations of $[A]$, $[H_2]$ and $[B]$, and constants. **25 pts**
# 
# **1c.** Four experiments were conducted in which all three of the components ($A$, $B$ and hydrogen) were present in the feed. The following reaction rates were obtained:
# 
# |Run| $[A]$, $mol/L$ |$[B]$, $mol/L$| $[H2]$, $mol/L$| $-d[A]/dt$, $mol/(L h)$|
# |:--:|---:|--:|--:|--:|
# |1| 0.04| 0.02| 0.06| 0.10|
# |2| 0.04| 0.04| 0.06| 0.11|
# |3| 0.02| 0.02| 0.03| 0.20|
# |4| 0.04| 0.02| 0.03| 0.05|
# 
# Does $[B]$ affect the reaction rate at other parameters kept constant? What is the apparent order to $B$ for these experiments? **5 pts**
# 
# how does $[H_2]$ affect the reaction rate at other parameters kept constant? What is the apparent order to hydrogen for these experiments? **5 pts**
# 
# how does $[A]$ affect the reaction rate at other parameters kept constant? What is the apparent order to $A$ for these experiments? **10 pts**
# 
# **1d.** Go back to the derived LH model in question **1a.** Is it consistent with experimental data and under
# what assumptions (i.e. what terms must be neglected)? What is the most abundant surface intermediate ($A$, $B$ or hydrogen)? **15 pts**
# 
# **1e.** Go back to the derived ER model in question **1b.** Is it consistent with experimental data and under what assumptions (i.e. what terms must be neglected)? **15 pts**
# 

# **Answer to Q1:**
# 
# **1a)** LH: 
# \begin{eqnarray}
# A+S<=>AS;\;\;\;\;\theta_A=K_A[A]\theta_s\\
# H_2 + 2S <=> 2(HS);\;\;\;\;\theta_H=(K_H[H_2])^{0.5}\theta_s\\
# AS + 2(HS) \rightarrow BS + 2S;\;\;\;\;r=k\theta_A\theta_H^2\\
# BS <=> B + S;\;\;\;\;\theta_B=K_B[B]\theta_s
# \end{eqnarray}
# Surface balance:
# \begin{eqnarray}
# \theta_A+\theta_H+\theta_B+\theta_s=1
# \end{eqnarray}
# \begin{eqnarray}
# \theta_s=\frac{1}{K_A[A]+(K_H[H_2])^{0.5}+K_B[B]+1}
# \end{eqnarray}
# \begin{eqnarray}
# r=k\theta_A\theta_H^2=\frac{kK_AK_H[H_2][A]}{(K_A[A]+(K_H[H_2])^{0.5}+K_B[B]+1)^3}
# \end{eqnarray}
# 
# 
# **1b)** ER:
# 
# \begin{eqnarray}
# H_2 + 2S <=> 2(HS)\\
# A + 2(HS) \rightarrow B + 2S
# \end{eqnarray}
# 
# \begin{eqnarray}
# \theta^2_H=K_H[H_2]\theta^2_S\\
# r=k[A]\theta^2_H
# \end{eqnarray}
# 
# Surface balance:
# \begin{eqnarray}
# \theta_H+\theta_s=1\\
# \theta_S=\frac{1}{(K_H[H_2])^{0.5}+1}
# \end{eqnarray}
# 
# \begin{eqnarray}
# r=\frac{k[A]K_H[H_2]}{(K_H[H_2])^{0.5}+1)^2}
# \end{eqnarray}

# **1c)** 
# 
# $[B]$? when $[A]$ and $[H_2]$ are constant: Runs $1$ and $4$.
# 
# $[H_2]$ drops twice, rate drops twice, 
# \begin{equation}
# r \propto [H_2]^1
# \end{equation}
# 
# $[A]$? when $[H_2]$ and $[B]$ are constant: Runs $3$ and $4$.
# 
# \begin{eqnarray}
# r_3 \propto [A]_3^a \\
# r_4 \propto [A]_4^a
# \end{eqnarray}
# 
# \begin{eqnarray}
# ln\frac{r_3}{r_4} \propto  a.ln\frac{[A]_3}{[A]_4}\\
# ln(4)=a.ln(\frac{1}{2})\\
# a=-2
# \end{eqnarray}
# Order to $[A]$ is $-2$.

# **1d)** 
# 
# If $K_A[A]>>((K_H[H_2])^{0.5}+K_B+1)$
# 
# $r=\frac{kK_AK_H[H_2][A]}{K^3_A[A]^3}=\frac{kK_H[H_2]}{K^2_A[A]^2}$
# 
# and the rate law matches the reaction order to $A$, $B$ and $H_2$.
# 
# This assumption means that $A$ is the most abundant surface intermediate ($\theta_A=1$).

# **1e)**
# 
# If $(K_H[H_2])^{0.5}<<1$ then $r=k[A]K_H[H_2]$ matches $1^{st}$ order to $H_2$, $0^{th}$ order to $B$ but doesn't match order to $A$. The rate law and mechanism are not consistent with experimental data and must be rejected.
