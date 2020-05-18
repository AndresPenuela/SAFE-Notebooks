# -*- coding: utf-8 -*-
"""
This a Python implementation and simplification of the DICE model.

@author: Andres Peñuela andres.penuela-fernandez@bristol.ac.uk, penyuela@gmail.com
"""

# Import necessary libraries
import numpy as np

class simple_DICE:
    
    def __init__(self, pop_param = 0.134, prod_grow_rate = 0.08, ec_grow_rate = 0.08, CO2_temp_factor = 0.13, clim_damage_exp = 2):
        
        ####################
        # Model parameters #
        ####################
        
        self.time_step = 5 # years
        
        ### Population ###
        self.pop_ini   = 7403 # in 2016
        self.pop_param = pop_param   # Shape of the asymptotic curve of population growth 
        self.pop_asym  = 11500 # Maximum population (population stabilizes)
        
        ### Capital ###
        self.saving_rate       = 0.25 # fraction of net output for investment
        self.depreciation_rate = 0.1
        self.capital_share     = 0.3
        
        ### Productivity ###
        self.prod_ini       = 5.12
        self.prod_grow_rate = prod_grow_rate       # Rate of increase of the productivity every 5 years
        
        ### Industrial emissions ###
        self.sigma           = 0.35 # in 2010 and 2015.  (industrial, MTCO2/$1000 2010 US$)
        self.sigma_grow_rate = -0.015
        
        ## Emissions control ###
        self.backstop_price     = 550000 # Price backstop technology (2010 US 000 $ per tCO2)
        self.backstop_grow_rate = -0.025 # per 5 years
        self.ec_exp             = 2.6
        self.ec_rate_ini        = 0.03
        self.ec_grow_rate       = ec_grow_rate
        
        ### Temp-CO2 emissions factor ###
        self.temp_ini        = 0.85 # Atmospheric temperature (degrees Celsius above preindustrial) in 2015
        self.CO2_equil       = 588 # Equil. conc. of CO2 in atmos. (GTC)
        self.CO2_temp_factor = CO2_temp_factor # = 0.13 = 0.101*(3.681/3.1+0.088)
        
        ### Climate change damage ###
        self.clim_damage_coeff = 0.00236
        self.clim_damage_exp   = clim_damage_exp
        
    def simulation(self):
        
        N = 20 # N is the number of time steps (N x time_step = 20 x 5 years = 100 years)
        population       = np.zeros(N)
        capital          = np.zeros(N)
        prod             = np.zeros(N)
        gross_output     = np.zeros(N)
        sigma            = np.zeros(N)
        backstop_price   = np.zeros(N)
        ec_rate          = np.zeros(N)
        ec_cost_coeff    = np.zeros(N)
        ec_frac          = np.zeros(N)
        CO2_emis         = np.zeros(N)
        CO2_conc         = np.zeros(N)
        temp             = np.zeros(N)
        clim_damage_frac = np.zeros(N)
        
        ######################
        # Initial conditions #
        ######################
        
        population[0]       = self.pop_ini
        capital[0]          = 223 # ($trill, 2005$)
        prod[0]             = self.prod_ini
        gross_output[0]     = prod[0] * capital[0]**self.capital_share * (population[0] / 1000) ** (1 - self.capital_share)
        backstop_price[0]   = 550 # Price of backstop technology (1000$ per tC)
        sigma[0]            = 0.35
        ec_rate[0]          = self.ec_rate_ini
        ec_cost_coeff[0]    = backstop_price[0] * sigma[0] / self.ec_exp / 1000
        ec_frac[0]          = ec_cost_coeff[0] * ec_rate[0] ** self.ec_exp
        CO2_emis[0]         = sigma[0] * gross_output[0] * (1 - ec_rate[0])
        CO2_conc[0]         = 851 # Atmospheric concentration of carbon (GTC)
        temp[0]             = self.temp_ini
        clim_damage_frac[0] = self.clim_damage_coeff * temp[0] ** self.clim_damage_exp
        
        ##############
        # Simulation #
        ##############
        
        for n in np.arange(1,N):
            
            ### Population ###
            population[n] = population[n-1] * (self.pop_asym/population[n-1]) ** self.pop_param
            
            ### Capital ###
            net_output = gross_output[n-1] * (1 - clim_damage_frac[n-1] - ec_frac[n-1])
            investment = net_output * self.saving_rate # investment in the previous time-step
            capital[n] = capital[n-1] * (1-self.depreciation_rate)**self.time_step + investment * self.time_step
            
            ### Productivity ###
            prod[n] = prod[n-1] * (1+self.prod_grow_rate)
            
            ### Gross outputs ###
            gross_output[n] = prod[n] * capital[n]**self.capital_share * (population[n] / 1000) ** (1 - self.capital_share) # Before discounting EC_cost and clim_damage
            
            ### Sigma(Industrial emissions) ###
            sigma[n] = sigma[n-1] * np.exp(self.sigma_grow_rate * self.time_step)
            
            ## Emissions control ###
            backstop_price[n] = backstop_price[n-1] * (1+self.backstop_grow_rate)
            ec_rate[n]        = ec_rate[n-1] * (1+self.ec_grow_rate)
            ec_cost_coeff[n]  = backstop_price[n] * sigma[n] / self.ec_exp / 1000
            ec_frac[n]        = ec_cost_coeff[n] * ec_rate[n] ** self.ec_exp 
            
            ### Temp-CO2 emissions factor ###
            CO2_emis[n] = sigma[n] * gross_output[n] * (1 - ec_rate[n])
            CO2_conc[n] = CO2_conc[n-1] + CO2_emis[n] * self.time_step / 3.666 
            
            ### Climate change damage ###
            #print(n)
            #print(ec_rate[n])
            rad_forc_incr       = 3.681 * np.log(CO2_conc[n] / self.CO2_equil)/np.log(2) + 0.5
            temp[n]             = temp[n-1] + 0.101 * rad_forc_incr - (self.CO2_temp_factor * temp[n-1])
            clim_damage_frac[n] = self.clim_damage_coeff * temp[n] ** self.clim_damage_exp
            
        return population, capital, prod, gross_output, sigma, backstop_price, ec_rate, ec_cost_coeff, ec_frac, CO2_emis, CO2_conc, temp, clim_damage_frac
        
            
def function(param):
    
    pop_param       = param[0] 
    prod_grow_rate  = param[1] 
    ec_grow_rate    = param[2] 
    CO2_temp_factor = param[3] 
    clim_damage_exp = param[4] 
    
    ####################
    # Model parameters #
    ####################
    
    time_step = 5 # years
    
    ### Population ###
    pop_ini   = 7403 # in 2016
    pop_param = pop_param   # Shape of the asymptotic curve of population growth 
    pop_asym  = 11500 # Maximum population (population stabilizes)
    
    ### Capital ###
    saving_rate       = 0.25 # fraction of net output for investment
    depreciation_rate = 0.1
    capital_share     = 0.3
    
    ### Productivity ###
    prod_ini       = 5.12
    prod_grow_rate = prod_grow_rate       # Rate of increase of the productivity every 5 years
    
    ### Industrial emissions ###
    sigma           = 0.35 # in 2010 and 2015.  (industrial, MTCO2/$1000 2010 US$)
    sigma_grow_rate = -0.015
    
    ## Emissions control ###
    backstop_price     = 550000 # Price backstop technology (2010 US 000 $ per tCO2)
    backstop_grow_rate = -0.025 # per 5 years
    ec_exp             = 2.6
    ec_rate_ini        = 0.03
    ec_grow_rate       = ec_grow_rate
    
    ### Temp-CO2 emissions factor ###
    temp_ini        = 0.85 # Atmospheric temperature (degrees Celsius above preindustrial) in 2015
    CO2_equil       = 588 # Equil. conc. of CO2 in atmos. (GTC)
    CO2_temp_factor = CO2_temp_factor # = 0.13 = 0.101*(3.681/3.1+0.088)
    
    ### Climate change damage ###
    clim_damage_coeff = 0.00236
    clim_damage_exp   = clim_damage_exp

    N = 20 # N is the number of time steps (N x time_step = 20 x 5 years = 100 years)
    population       = np.zeros(N)
    capital          = np.zeros(N)
    prod             = np.zeros(N)
    gross_output     = np.zeros(N)
    sigma            = np.zeros(N)
    backstop_price   = np.zeros(N)
    ec_rate          = np.zeros(N)
    ec_cost_coeff    = np.zeros(N)
    ec_frac          = np.zeros(N)
    CO2_emis         = np.zeros(N)
    CO2_conc         = np.zeros(N)
    temp             = np.zeros(N)
    clim_damage_frac = np.zeros(N)
    
    ######################
    # Initial conditions #
    ######################
    
    population[0]       = pop_ini
    capital[0]          = 223 # ($trill, 2005$)
    prod[0]             = prod_ini
    gross_output[0]     = prod[0] * capital[0]**capital_share * (population[0] / 1000) ** (1 - capital_share)
    backstop_price[0]   = 550 # Price of backstop technology (1000$ per tC)
    sigma[0]            = 0.35
    ec_rate[0]          = ec_rate_ini
    ec_cost_coeff[0]    = backstop_price[0] * sigma[0] / ec_exp / 1000
    ec_frac[0]          = ec_cost_coeff[0] * ec_rate[0] ** ec_exp
    CO2_emis[0]         = sigma[0] * gross_output[0] * (1 - ec_rate[0])
    CO2_conc[0]         = 851 # Atmospheric concentration of carbon (GTC)
    temp[0]             = temp_ini
    clim_damage_frac[0] = clim_damage_coeff * temp[0] ** clim_damage_exp
    
    ##############
    # Simulation #
    ##############
    
    for n in np.arange(1,N):
        
        ### Population ###
        population[n] = population[n-1] * (pop_asym/population[n-1]) ** pop_param
        
        ### Capital ###
        net_output = gross_output[n-1] * (1 - clim_damage_frac[n-1] - ec_frac[n-1])
        investment = net_output * saving_rate # investment in the previous time-step
        capital[n] = capital[n-1] * (1-depreciation_rate)**time_step + investment * time_step
        
        ### Productivity ###
        prod[n] = prod[n-1] * (1+prod_grow_rate)
        
        ### Gross outputs ###
        gross_output[n] = prod[n] * capital[n]**capital_share * (population[n] / 1000) ** (1 - capital_share) # Before discounting EC_cost and clim_damage
        
        ### Sigma(Industrial emissions) ###
        sigma[n] = sigma[n-1] * np.exp(sigma_grow_rate * time_step)
        
        ## Emissions control ###
        backstop_price[n] = backstop_price[n-1] * (1+backstop_grow_rate)
        ec_rate[n]        = ec_rate[n-1] * (1+ec_grow_rate)
        ec_cost_coeff[n]  = backstop_price[n] * sigma[n] / ec_exp / 1000
        ec_frac[n]        = ec_cost_coeff[n] * ec_rate[n] ** ec_exp 
        
        ### Temp-CO2 emissions factor ###
        CO2_emis[n] = sigma[n] * gross_output[n] * (1 - ec_rate[n])
        CO2_conc[n] = CO2_conc[n-1] + CO2_emis[n] * time_step / 3.666 
        
        ### Climate change damage ###
        #print(n)
        #print(ec_rate[n])
        rad_forc_incr       = 3.681 * np.log(CO2_conc[n] / CO2_equil)/np.log(2) + 0.5
        temp[n]             = temp[n-1] + 0.101 * rad_forc_incr - (CO2_temp_factor * temp[n-1])
        clim_damage_frac[n] = clim_damage_coeff * temp[n] ** clim_damage_exp
        
        ec_cost = ec_frac*gross_output
        clim_damage = clim_damage_frac*gross_output
        rel_clim_damage = (ec_cost+clim_damage)/gross_output *100
        
    return np.array(rel_clim_damage[-1])
    
        