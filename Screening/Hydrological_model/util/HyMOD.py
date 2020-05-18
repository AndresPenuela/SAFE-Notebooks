# -*- coding: utf-8 -*-
"""
This function is a Python implementation of the rainfall-runoff model Hymod

@author: Andres Peñuela andres.penuela-fernandez@bristol.ac.uk, penyuela@gmail.com
"""

# HyMOD model
import numpy as np

class hymod:

    def __init__(self, c_max, beta, alpha, K_s, K_q):
        # Hymod parameters
        self.c_max = c_max                         # Maximum storage capacity
        self.beta  = beta                          # Degree of spatial variability of the c_max ditribution
        self.alpha = alpha                         # Factor distributing slow and quick flows
        self.K_s   = K_s                           # Fractional discharge of the slow release reservoir
        self.K_q   = K_q                           # Fractional discharge of the quick release reservoir
        self.s_max = self.c_max / (1. + self.beta) # Maximum soil moisture

    def simulation(self, P, PET, s_ini = [0,0,0,0,0]):
        
        T = len(P) # number of time steps
        
        self.s_star = s_ini[0]
        self.s_s    = s_ini[1]
        self.s_1    = s_ini[2]
        self.s_2    = s_ini[3]
        self.s_3    = s_ini[4]
        
        ER1  = np.zeros(T) # effective rainfall 1
        ER2  = np.zeros(T) # effective rainfall 2
        Q_s  = np.zeros(T) # slow flow
        Q_q3 = np.zeros(T) # quick flow
        Q    = np.zeros(T) # flow
        e    = np.zeros(T) # actual evaporation
        
        for t in range(T):
            
            c = self.c_max * (1 - (1 - ((self.beta + 1) * self.s_star) / (self.c_max))**(1 / (self.beta+1)))
            ER1[t] = np.max([P[t] - (self.c_max - c), 0.0]) # effective rainfall part 1
            c_star = np.min([P[t] + c, self.c_max])
            s = np.min([(self.c_max / (1 + self.beta)) * (1 - (1 -(c_star / self.c_max))**(self.beta + 1)),self.s_max])
            e[t] = np.min([PET[t] * c_star / self.c_max , s])
            ER2[t] = np.max([(c_star - c) - (s - self.s_star)])
            self.s_star = s - e[t]
            self.s_s = (1 - self.K_s) * self.s_s + (1 - self. K_s) * (1 - self.alpha) * ER2[t]
            Q_s[t] = (self.K_s / (1 - self.K_s)) * self.s_s
            self.s_1 = (1 - self.K_q) * self.s_1 + (1 - self.K_q) * (ER1[t] + self.alpha * ER2[t])
            Q_q1 = (self.K_q / (1 - self.K_q)) * self.s_1
            self.s_2 = (1 - self.K_q) * self.s_2 + (1 - self.K_q) * Q_q1
            Q_q2 = (self.K_q / (1 - self.K_q)) * self.s_2
            self.s_3 = (1 - self.K_q) * self.s_3 + (1 - self.K_q) * Q_q2
            Q_q3[t] = (self.K_q / (1 - self.K_q)) * self.s_3            
            Q[t] = Q_s[t] + Q_q3[t]
            
        return ER2, Q_q3, Q_s, Q

    def get_params(self):
        return self.alpha, self.beta, self.c_max, self.K_q. self.K_s
    
    def get_end_cond(self):
        return self.s_star, self.s_s, self.s_1, self.s_2. self.s_3
        

