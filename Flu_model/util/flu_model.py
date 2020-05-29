# -*- coding: utf-8 -*-
"""
A simple mathematical description of the spread of a flu in a company is the 
so-called the flu model, which divides the (fixed) population of N individuals 
into three "compartments" which may vary as a function of time, t:

V(t) are those vulnerable but not yet infected with the flu;

S(t) is the number of sick individuals;

RI(t) are those individuals who either are immune or have recovered from the 
flu and now have immunity to it.

The model describes the change in the population of each of these compartments 
in terms of two parameters, β and γ. 

β describes the effective contact rate of the flu: an infected individual comes 
into contact with βN other individuals per unit time (of which the fraction 
that are susceptible to contracting the flu is V/N). 

γ is the mean recovery rate: that is, 1/γ is the mean period of time during 
which a sick individual can pass it on.

The differential equations describing this model were first derived by Kermack 
and McKendrick [Proc. R. Soc. A, 115, 772 (1927)]:

dV / dt = -βVS / N
 
dS / dt = βVS / N - γS

dRI / dt = γS

@author: Andres Peñuela
"""
from scipy.integrate import odeint
import numpy as np

class population:
    
    def __init__(self, N = 1000, RI_ratio = 0.1, S_0 = 1, D_0 = 0):
        # Total population, N.
        self.N = N
        # Ratio of population who is immune
        self.RI_ratio = RI_ratio
        # Initial number of sick individuals
        self.S_0 = S_0
        # Initial number of immune individuals
        self.RI_0 = N * RI_ratio
        # Deads, D_0
        self.D_0 = D_0
        # Everyone else, S_0, is susceptible to infection initially.
        self.V_0 = N - S_0 - self.RI_0 - D_0

        
def simulation(t,population, beta, gamma, alpha):
    
    y_0 = population.V_0, population.S_0, population.RI_0, population.D_0
    
    def deriv(y, t, N, beta, gamma,alpha):
        V, S, IR, D = y
        dVdt = -beta * V * S / N
        dSdt = beta * V * S / N - gamma * S - alpha * S * S / N #- sigma * S
        dRIdt = gamma * S
        dDdt = alpha * S * S / N # + sigma * S
        
        return dVdt, dSdt, dRIdt, dDdt
    
    ret = odeint(deriv, y_0, t, args=(population.N, beta, gamma, alpha))
    V, S, RI, D = ret.T
    
    return D, S, RI, V

def function(param,N,t):
    RI_0 = param[0]*N
    S_0 = param[1]
    beta = param[2]
    gamma = 1/param[3]
    alpha = param[4]
    
    D_0 = 0
    V_0 = N - S_0 - RI_0 - D_0
    
    y_0 = V_0, S_0, RI_0, D_0
    
    def deriv(y, t, N, beta, gamma,alpha):
        V, S, IR, D = y
        dVdt = -beta * V * S / N
        dSdt = beta * V * S / N - gamma * S - alpha * S * S / N #- sigma * S
        dRIdt = gamma * S
        dDdt = alpha * S * S / N # + sigma * S  
        
        return dVdt, dSdt, dRIdt, dDdt
    
    ret = odeint(deriv, y_0, t, args=(N, beta, gamma, alpha))
    V, S, RI, D = ret.T
    
#    if output == 0:
#        out = np.array(D[-1])
#    elif output == 1:
#        out = np.array(S[-1])
#    elif output == 2:
#        out = np.array(RI[-1])
#    elif output == 3:
#        out = np.array(V[-1])
    
    out = np.array(np.max(S))
    
    return out