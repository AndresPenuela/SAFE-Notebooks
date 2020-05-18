# -*- coding: utf-8 -*-
"""
This function is a Python implementation of the rainfall-runoff model Hymod

@author: Andres Peñuela andres.penuela-fernandez@bristol.ac.uk, penyuela@gmail.com
"""

# HyMOD model
import numpy as np
from util.util import NSE

class hymod_gamma:

    def __init__(self, param):
        
        #######################################################################
        # Recover model parameters
        #######################################################################
        self.Sm    = 0   # Maximum soil moisture 
        self.gamma = 0   # Evaporation rate
        self.alpha = 0 # Factor distributing slow and fast flows
        self.Rf    = 0 # Fast reservoir outflow coefficient (ratio) [1/Dt]
        self.Rs    = 0 # Slow reservoir outflow coefficient (ratio) [1/Dt]
        
    def simulation(self, param, rain, ept):
        
        self.Sm    = param[0]   # Maximum soil moisture 
        self.gamma = param[1]   # Evaporation rate
        self.alpha = 1-param[2] # Factor distributing slow and fast flows
        self.Rf    = 1/param[3] # Fast reservoir outflow coefficient (ratio) [1/Dt]
        self.Rs    = 1/param[4] # Slow reservoir outflow coefficient (ratio) [1/Dt]
        
        T = len(rain)
        
        #######################################################################
        # Initialise variables
        #######################################################################
        Pe = np.zeros((T, )) # Recharge from the soil [mm/Dt]
        Ea = np.zeros((T, )) # Actual Evapotranspiration [mm/Dt]
        sm = np.zeros((T+1, )) # Soil Moisture [mm]
        sL = np.zeros((T+1, )) # Slow reservoir moisture [mm]
        sF1 = np.zeros((T+1, )) # Fast reservoir 1 moisture [mm]
        sF2 = np.zeros((T+1, )) # Fast reservoir 2 moisture [mm]
        sF3 = np.zeros((T+1, )) # Fast reservoir 3 moisture [mm]
        QsL = np.zeros((T, )) # Slow flow [mm/Dt]
        QsF = np.zeros((T, )) # Fast flow [mm/Dt]
        
        beta = 0
        
        for t in range(T):
            ###################################################################
            # Soil moisture dynamics
            ###################################################################           
            F = 1 - (1-sm[t]/self.Sm)**beta
            Pe[t] = F * rain[t] # Compute the value of the outflow
            # (we assumed that this process is faster than evaporation)
            sm_temp = max(min(sm[t] + rain[t] - Pe[t], self.Sm), 0)
            # Compute the water balance with the value of the outflow
            Pe[t] = Pe[t] + max(sm[t] + rain[t] - Pe[t] - self.Sm, 0) + \
                    min(sm[t] + rain[t] - Pe[t], 0)
            # Adjust Pe by an amount equal to the possible negative sm amount or
            # to the possible sm amount above Sm.
    
            W = min(np.abs(sm[t]/self.Sm)*self.gamma, 1) # Correction factor for evaporation (modified equation)
            Ea[t] = W * ept[t] # Compute the evaporation
            sm[t+1] = max(min(sm_temp - Ea[t], self.Sm), 0) # Compute the water balance
            Ea[t] = Ea[t] + max(sm_temp - Ea[t] - self.Sm, 0) + min(sm_temp - Ea[t], 0)
            # Adjust Ea by an amount equal to the possible negative sm amount or to
            # the possible sm amount above Sm.
            
            ###################################################################
            # Groundwater dynamics
            ###################################################################
            # slow flow
            QsL[t] = self.Rs * sL[t]
            sL[t+1] = sL[t] + (1-self.alpha)*Pe[t] - QsL[t]
            # fast flow
            sF1[t+1] = sF1[t] +  self.alpha*Pe[t] - self.Rf*sF1[t]
            sF2[t+1] = sF2[t] +  self.Rf*sF1[t] - self.Rf*sF2[t]
            sF3[t+1] = sF3[t] +  self.Rf*sF2[t] - self.Rf*sF3[t]
            QsF[t] = self.Rf * sF3[t]

        Q_sim = QsL + QsF
        
        self.Pe = Pe # Recharge from the soil [mm/Dt]
        self.Ea = Ea # Actual Evapotranspiration [mm/Dt]
        self.sm = sm # Soil Moisture [mm]
        self.sL = sL # Slow reservoir moisture [mm]
        self.sF1 = sF1 # Fast reservoir 1 moisture [mm]
        self.sF2 = sF2 # Fast reservoir 2 moisture [mm]
        self.sF3 = sF3 # Fast reservoir 3 moisture [mm]
        self.QsL = QsL # Slow flow [mm/Dt]
        self.QsF = QsF # Fast flow [mm/Dt]
            
        return Q_sim

def hymod_gamma_nse(x, rain, ept, flow, warmup):

    """This function runs the rainfall-runoff Hymod model
    and returns the associated Nash-Sutcliffe Efficiency

    y, Q_sim, STATES, FLUXES = HyMod.hymod_nse(x, rain, ept, flow, warmup)

    Input:
         x = 5 elements vector of model parameters         - numpy.ndarray(5, )
             (Smax, gamma, alfa, Rs, Rf)
      rain = time series of rainfall                       - numpy.ndarray(T, )
      ept = time series of potential evaporation          - numpy.ndarray(T, )
      flow = time series of observed flow                  - numpy.ndarray(T, )
    warmup =  number of time steps for model warm-up       - integer

    Output:
         y = Nash-Sutcliffe Efficiency                     - numpy.ndarray(1, )
     Q_sim = time series of simulated flow                 - numpy.ndarray(T, )
    STATES = time series of simulated storages             - numpy.ndarray(T,5)
             (all in mm)
    FLUXES = time series of simulated fluxes               - numpy.ndarray(T,4)
            (all in mm/Dt)

    See also HyMod.hymod_sim about the model parameters, simulated variables,
    and references."""

    ###########################################################################
    # Check inputs
    ###########################################################################
    M = 5 # number of model parameters
    if not isinstance(x, np.ndarray):
        raise ValueError('"x" must be a numpy.array.')
    if x.dtype.kind != 'f' and x.dtype.kind != 'i' and x.dtype.kind != 'u':
        raise ValueError('"x" must contain floats or integers.')
    Nx = x.shape
    if len(Nx) != 1 or len(x) != M:
        raise ValueError('"x" must have shape (5, ).')

    if not isinstance(rain, np.ndarray):
        raise ValueError('"rain" must be a numpy.array.')
    if rain.dtype.kind != 'f' and rain.dtype.kind != 'i' and rain.dtype.kind != 'u':
        raise ValueError('"rain" must contain floats or integers.')
    Nrain = rain.shape
    if len(Nrain) != 1:
        raise ValueError('"rain" must be of shape (T, ).')
    T = Nrain[0]

    if not isinstance(ept, np.ndarray):
        raise ValueError('"ept" must be a numpy.array.')
    if ept.dtype.kind != 'f' and ept.dtype.kind != 'i' and ept.dtype.kind != 'u':
        raise ValueError('"ept" must contain floats or integers.')
    Nept = ept.shape
    if len(Nept) != 1:
        raise ValueError('"ept" must be of shape (T, ).')
    if len(ept) != T:
        raise ValueError('"ept" and "prec" must have the same number of elements.')

    if not isinstance(flow, np.ndarray):
        raise ValueError('"flow" must be a numpy.array.')
    if flow.dtype.kind != 'f' and flow.dtype.kind != 'i' and flow.dtype.kind != 'u':
        raise ValueError('"flow" must contain floats or integers.')
    Nflow = flow.shape
    if len(Nflow) != 1:
        raise ValueError('"flow" must be of shape (T, ).')
    if len(flow) != T:
        raise ValueError('"flow" and "rain" must have the same number of elements.')

    if not isinstance(warmup, (int, np.int8, np.int16, np.int32, np.int64)):
        raise ValueError('"warmup" must be scalar and integer.')
    if warmup < 0 or warmup >= T:
        raise ValueError('"warmup" must be in [0, T).')

    ###########################################################################
    # Simulate HyMod and compute scalar output
    ###########################################################################
    model = hymod_gamma(x)
    Q_sim = model.simulation(x, rain, ept)

    Qs = Q_sim[warmup:len(Q_sim)+1]
    Qo = flow[warmup:len(Q_sim)+1]
    y = np.array(np.max([NSE(Qs, Qo),0])*100)

    return y, Q_sim

def function(param, rain, ept):
    
    #######################################################################
    # Recover model parameters
    #######################################################################
    Sm    = param[0] # Maximum soil moisture 
    gamma = param[1] # Evapotranspiration parameter
    alpha = param[2] # Factor distributing slow and fast flows
    Rs    = param[3] # Slow reservoir outflow coefficient (ratio) [1/Dt]
    Rf    = param[4] # Fast reservoir outflow coefficient (ratio) [1/Dt]

    T = len(rain)
    
    #######################################################################
    # Initialise variables
    #######################################################################
    Pe = np.zeros((T, )) # Recharge from the soil [mm/Dt]
    Ea = np.zeros((T, )) # Actual Evapotranspiration [mm/Dt]
    sm = np.zeros((T+1, )) # Soil Moisture [mm]
    sL = np.zeros((T+1, )) # Slow reservoir moisture [mm]
    sF1 = np.zeros((T+1, )) # Fast reservoir 1 moisture [mm]
    sF2 = np.zeros((T+1, )) # Fast reservoir 2 moisture [mm]
    sF3 = np.zeros((T+1, )) # Fast reservoir 3 moisture [mm]
    QsL = np.zeros((T, )) # Slow flow [mm/Dt]
    QsF = np.zeros((T, )) # Fast flow [mm/Dt]
    
    beta = 0
    
    for t in range(T):
        ###################################################################
        # Soil moisture dynamics
        ###################################################################           
        F = 1 - (1-sm[t]/Sm)**beta
        Pe[t] = F * rain[t] # Compute the value of the outflow
        # (we assumed that this process is faster than evaporation)
        sm_temp = max(min(sm[t] + rain[t] - Pe[t], Sm), 0)
        # Compute the water balance with the value of the outflow
        Pe[t] = Pe[t] + max(sm[t] + rain[t] - Pe[t] - Sm, 0) + \
                min(sm[t] + rain[t] - Pe[t], 0)
        # Adjust Pe by an amount equal to the possible negative sm amount or
        # to the possible sm amount above Sm.
        
        # Ea[t] = ept[t]*(1-np.exp(-6.68*sm[t]-1/Sm))
        W = min(np.abs(sm[t]/Sm)*gamma, 1) # Correction factor for evaporation (modified equation)
        Ea[t] = W * ept[t] # Compute the evaporation
        sm[t+1] = max(min(sm_temp - Ea[t], Sm), 0) # Compute the water balance
        Ea[t] = Ea[t] + max(sm_temp - Ea[t] - Sm, 0) + min(sm_temp - Ea[t], 0)
        # Adjust Ea by an amount equal to the possible negative sm amount or to
        # the possible sm amount above Sm.
        
        ###################################################################
        # Groundwater dynamics
        ###################################################################
        # slow flow
        QsL[t] = Rs * sL[t]
        sL[t+1] = sL[t] + (1-alpha)*Pe[t] - QsL[t]
        # fast flow
        sF1[t+1] = sF1[t] +  alpha*Pe[t] - Rf*sF1[t]
        sF2[t+1] = sF2[t] +  Rf*sF1[t] - Rf*sF2[t]
        sF3[t+1] = sF3[t] +  Rf*sF2[t] - Rf*sF3[t]
        QsF[t] = Rf * sF3[t]

    Q_sim = QsL + QsF
    STATES = np.column_stack((sm, sL, sF1, sF2, sF3))
    FLUXES = np.column_stack((Pe, Ea, QsL, QsF))

    return Q_sim, STATES, FLUXES
