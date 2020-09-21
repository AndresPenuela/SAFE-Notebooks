#!/usr/bin/env python
# coding: utf-8

# # How can we reach the equilibrium in a predator-prey habitat?
# In this Notebook we will try to understand the evolution and interactions in a predator-prey habitat and how to reach the equilibrium of the prey and predator populations.For this purpose we will use a predator-prey model and Global Sensitivity Analysis (GSA).
# 
# # What is a predator-prey habitat?
# In a predator-prey habiat we consider two species. One species, the prey, is the primary source of food source for the other, the predator. The predator-prey food chain persists over time thanks to a **negative feedback loop**, i.e. when the product of a reaction leads to a decrease in that reaction. When predators eat preys, they reduce their population but this drop in the predators food source will soon cause the predator population to decline. This will also reduce the predation rate allowing the prey population to increase again. This cycle creates an up and down wavelike pattern that maintain a long-term equilibrium.
# 
# <left><img src="util/predator_prey_equil.gif" width="700px">
# 
# Now imagine that in a given natural habitat we would like to introduce a predator species to maintain the population of rabbits stable. For this purpose, we can select the initial number and species of predators as well as the initial number of rabbits. We could do this by trial and error by introducing different predator species and then monitoring the evolution of the populations but this will likely require high amounts of time, money and animal lives. But what if we could we simulate the habitat and its evolution on a computer using mathematical equations? 
# 
# A hydrological model is a mathematical model (set of equations) describing the hydrological processes that occur in a catchment as a function of various parameters. These model parameters describe the hydrological characteristics of the catchment, such as the climate and soil characteristics, ultimately enabling the estimation of the river flow at selected river sections.
# 
# # What is a predator-prey model?
# A predator-prey model is a mathematical model (set of equations) describing the interactions that occur in a predator-prey habitat as a function of various parameters. These model parameters describe the characteristics of the habitat, such as the initial number of predators and preys and the characteristics of the predator species, ultimately enabling the estimation of the evolution in time of the predator and prey populations.
# 
# ## Model parameters
# - Initial number of predators
# - Attack rate of predators: number of times that a predator attacks per week
# - Efficiency rate of predators: number of encounters per week that resulted in a kill
# - Death rate of predators: fraction of the predator population that dies per week
# - Initial number of preys
# - Growth rate of preys: number of preys that are born per week
# - Carrying capacity of environment: maximum number of preys that the environment can sustain
# 
# ## Model outputs
# - Daily predator population
# - Daily prey population
# 
# 
# ## Model assumptions
# - The food supply of the predator population depends entirely on the size of the prey population.
# - The rate of change of population is proportional to its size.
# - Predators have limitless appetite.
# - If the prey population grows beyond the carrying capacity of the environment, then their population would be wiped out as all the available food resources would have been consumed.
# - Preys only die as a result of predator attacks.
# 
# Now, using the predator-prey model and changing the model parameters, let's try to **reach the equilibrium** of the predator and prey populations in **less than a year** and that at the equilibrium point there are **at least 5 individuals of each**.

# In[1]:


# # What is Global Sensitivity Analysis? and why shall we use it?
# 
# **Global Sensitivity Analysis** is a set of mathematical techniques which investigate how uncertainty in the output of a numerical model can be attributed to variations of its input factors.
# 
# The benefits of applying GSA are:
# 
# 1. **Better understanding of the model**
#     
# 2. **Test whether the model "behaves well"**
#     
# 3. **Identify the important inputs on which to focus computation time (e.g. for calibration), acquisition of new data, etc.** 
#     
# 4. **Understand the main impacts of input uncertainty on the modelling outcome and thus model-informed decisions**
# 
# ## How Global Sensitivity Analysis works?
# 
# GSA investigates how the uncertainty of the selected model input parameters influences the variability of the model output/performance.
# 
# A '**model parameter**' is any element that can be changed before running the model.
# 
# An '**output**' is any variable that is obtained after the model execution.
# 
# Before executing the model, we will sample the parameters from their ranges of variability and then repeatedly run the model so that for each simulation all the inputs vary simultaneously. After GSA is performed we obtain a set of sensitivity indices for each output. The sensitivity indices measure the relative influence of each input factor on the output (*Refs. 4-5*).
# 
# ## Investigate interactions between model parameters
# In order to investigate the interactions between input factors we plot one input against the other, coloured by the value taken by the output.
# 
# ### a) Model output = predator population at equilibrium

# ### Step 1: Import python modules

# In[2]:


from __future__ import division, absolute_import, print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import plotly.graph_objects as go
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')

from util import PAWN
from util.model_execution import model_execution # module to execute the model
from util.sampling import AAT_sampling# module to perform the input sampling
from util.util import aggregate_boot

import util.predator_prey_model as predator_prey_model

def screening():
    # ### Step 2: Setup the model
    
    # Define: 
    # - the input factors whose influence is to be analysed with GSA, 
    # - their range of variability (choice made by expert judgement, available data or previous studies),
    # - choice of their distributions.
    
    # In[3]:
    
    
    param_names  = ["Predator initial population (x1000)",
                    "Predator attack rate (attacks per week)",
                    "Predator efficiency ratio", 
                    "Predator death rate (deaths per week)", 
                    "Prey initial population (x1000)"]
    M = len(param_names) # number of parameters
    
    # Parameter distributions:
    distr_fun = st.uniform # uniform distribution
    samp_strat = 'lhs' # Latin Hypercube
    # The shape parameters of the uniform distribution are the lower limit and the difference between lower and upper limits:
    distr_par  = [np.nan] * M
    
    # Define output:
    fun_test = predator_prey_model.function
    
    
    # Range of variability
    
    # In[4]:
    
    
    data = [["-",  2,     10 ,param_names[0]],
            ["-",  0.01,   1, param_names[1]],
            ["-",  0.01,  1, param_names[2]],
            ["-",  0.01,  1, param_names[3]],
            ["-",  2,     10, param_names[4]]]
    
    model_param = pd.DataFrame(data, 
                               columns=["Unit", "Min value", "Max value","Description"],
                               index = param_names)
    model_param
    
    
    # In[5]:
    
    
    outputs  = ["equilibrium error"]
    
    model_outputs = pd.DataFrame(outputs, columns = ['model outputs'])
    model_outputs
    
    
    # In[6]:
    
    
    class setup_model:
        def __init__(self, x1, x2, x3, x4, x5):
            # The shape parameters of the uniform distribution are the lower limit and 
            # the difference between lower and upper limits:
            self.xmin = [x1.value[0], x2.value[0], x3.value[0], x4.value[0], x5.value[0]]
            self.xmax = [x1.value[1], x2.value[1], x3.value[1], x4.value[1], x5.value[1]]
            for i in range(M):
                distr_par[i] = [self.xmin[i], self.xmax[i] - self.xmin[i]]
            self.distr_par = distr_par
    
    
    # ### Step 3: Sample inputs space
    
    # The number of model runs ($N$) typically increases proportionally to the number of input factors ($M$) and will depend on the GSA methods chosen as well. 
    
    # In[7]:
    
    
    N = 10000 # number of samples
    class sample_input:
        def __init__(self,distr_par):
            self.X = AAT_sampling(samp_strat, M, distr_fun, distr_par, N)   
    
    
    # ### Step 4: Run the model
    
    # For each sampled input factors combination, we run the model and save the associated model output.
    
    # In[8]:
    
    
    T = 365 # days
    class run_model:
        def __init__(self,X):
            self.Y = model_execution(fun_test, X, T) 
    
    
    # ### Step 5: Apply the PAWN Global Sensitivity Analysis method
    # Let’s now apply Global Sensitivity Analysis: for example the **PAWN** method. 
    # 
    # Its main idea is that the influence of an input factor is proportional to the amount of change in the output distribution produced by fixing that input.
    # <br> The sensitivity of $y$ to $x_{i}$ is measured by the difference between the unconditional CDF of $y$: $F_y( y )$, which is induced by varying all input factors simultaneously, and the conditional CDF that is obtained by varying all inputs but $x_{i}$: $F_{y|x_i}( y | x_i )$.
    
    # ### Step 6: Check model behaviour by visualising input/output samples
    # Scatterplots are plotted to visualise the behaviour of the output over each input factor in turn.
    # 
    # Definition of interactivity
    
    # In[9]:
    
    
    def update_figures(change):
        with fig1.batch_update():
            distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
            X = sample_input(distr_par).X
            Y = run_model(X).Y
            KS_median, _, _, KS_dummy = PAWN.pawn_indices(X, Y, n, dummy = True)
            KS_median, _, _ = PAWN.pawn_indices(X, Y, n, Nboot=1000)
            KS_median_m, KS_median_lb, KS_median_ub = aggregate_boot(KS_median)
    
            k = 0
            for i in range(0,M):
                fig1.data[k].y = np.ones(2)*KS_median_m[i]
                fig1.data[k+1].y = np.ones(2)*KS_median_lb[i]
                fig1.data[k+2].y = np.ones(2)*KS_median_ub[i]
                k = k + 3
                
            fig1.data[k].y = [KS_dummy[0],KS_dummy[0]]
    #         fig1.layout.title = 'model output: <b>'+outputs[output.value]
    
    
    # Definition of the sliders
    
    # In[10]:
    
    
    # output = widgets.IntSlider(value = 0, min = 0, max = 1,
    #                               step = 1, description = 'output',
    #                               style = {'description_width': '300px'} ,layout={'width': '700px'},
    #                               continuous_update=False)
    # output.observe(update_figures,names = 'value')
    
    x1 = widgets.IntRangeSlider(value = [model_param['Min value'][0], model_param['Max value'][0]], 
                                min = model_param['Min value'][0], max = model_param['Max value'][0],
                                  step = 1, description = model_param.index[0], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  continuous_update=False)
    x1.observe(update_figures,names = 'value')
    
    x2 = widgets.FloatRangeSlider(value = [model_param['Min value'][1], model_param['Max value'][1]], 
                                  min = model_param['Min value'][1], max = model_param['Max value'][1], 
                                  step = 0.01, description = model_param.index[1],
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    x2.observe(update_figures,names = 'value')
    
    x3 = widgets.FloatRangeSlider(value = [model_param['Min value'][2], model_param['Max value'][2]], 
                                  min = model_param['Min value'][2], max = model_param['Max value'][2], 
                                  step = 0.01, description = model_param.index[2],
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    x3.observe(update_figures,names = 'value')
    
    x4 = widgets.FloatRangeSlider(value = [model_param['Min value'][3], model_param['Max value'][3]], 
                                  min = model_param['Min value'][3], max = model_param['Max value'][3], 
                                  step = 0.01, description = model_param.index[3], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    x4.observe(update_figures,names = 'value')
    
    x5 = widgets.IntRangeSlider(value = [model_param['Min value'][4], model_param['Max value'][4]], 
                                min = model_param['Min value'][4], max = model_param['Max value'][4], 
                                  step = 1, description = model_param.index[4], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  continuous_update=False)
    x5.observe(update_figures,names = 'value')
    
    # x6 = widgets.FloatRangeSlider(value = [model_param['Min value'][5], model_param['Max value'][5]], 
    #                               min = model_param['Min value'][5], max = model_param['Max value'][5], 
    #                               step = 0.01, description = model_param.index[5], 
    #                               style = {'description_width': '300px'} ,layout={'width': '700px'},
    #                               readout_format = '.1f', continuous_update=False)
    # x6.observe(update_figures,names = 'value')
    
    # x7 = widgets.IntRangeSlider(value = [model_param['Min value'][6], model_param['Max value'][6]], 
    #                             min = model_param['Min value'][6], max = model_param['Max value'][6], 
    #                               step = 1, description = model_param.index[4], 
    #                               style = {'description_width': '300px'} ,layout={'width': '700px'},
    #                               continuous_update=False)
    # x7.observe(update_figures,names = 'value')
    
    #for i in range(M):
    #    
    #    exec('x'+str(i+1)+' = widgets.IntRangeSlider(value = [model_param[\'Min value\']['+str(i)+'], model_param[\'Max value\']['+str(i)+']],min = model_param[\'Min value\']['+str(i)+'], max = model_param[\'Max value\']['+str(i)+'],step = 1, description = model_param.index['+str(i)+'], style = {\'description_width\': \'300px\'} ,layout={\'width\': \'700px\'},continuous_update=False)')
    #    exec('x'+str(i+1)+'.observe(update_figures,names = \'value\')')

    
    
    #### Plot sensitivity indices and identify non-influential parameters
    #
    #The dummy parameter is a numerical artifice, with no influence on the model output, which is used to estimate the threshold for non-influential inputs. 
    #
    #Uninfluential input factors should have zero-valued sensitivity indices, but since sensitivity indices are computed by numerical approximations rather than analytical solutions, an uninfluential factor may still be associated with a non-zero (although small) index value. 
    #
    #Therefore if the index of an input factor is below the value of the sensitivity index of the dummy parameter, then the input factor is deemed uninfluential *(Ref. 6)*.
    # In[11]:
    
    
    n=10
    distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
    X = sample_input(distr_par).X
    Y = run_model(X).Y # NSE*100 = simulation accuracy (%)
    KS_median, _, _, KS_dummy = PAWN.pawn_indices(X, Y, n, dummy = True)
    KS_median, _, _ = PAWN.pawn_indices(X, Y, n, Nboot=1000)
    # KS_median has shape (Nboot, M)
    # Compute mean and confidence intervals of the sensitivity indices across the
    # bootstrap resamples:
    KS_median_m, KS_median_lb, KS_median_ub = aggregate_boot(KS_median) # shape (M,)
    
    
    # In[14]:
    
    
    fig1 = go.FigureWidget(layout = dict(width=500, height=400,showlegend = False,margin=dict(t=10,r=0,l=75)))
    # Soil storage capacity
    fig1.add_trace(go.Scatter(x=[-0.25,0.25], y=np.ones(2)*KS_median_m[0], mode = 'lines', line = dict(color ='rgba(120, 170, 150, 1)')))
    fig1.add_trace(go.Scatter(x=[-0.25,0.25], y=np.ones(2)*KS_median_lb[0], mode = 'lines', line = dict(color ='rgba(120, 170, 150, 0)')))
    fig1.add_trace(go.Scatter(x=[-0.25,0.25], y=np.ones(2)*KS_median_ub[0], mode = 'none', fill='tonexty',fillcolor = 'rgba(120, 170, 150, 0.25)'))
    # Evaporation rate
    fig1.add_trace(go.Scatter(x=[0.75,1.25], y=np.ones(2)*KS_median_m[1], mode = 'lines', line = dict(color ='rgba(250, 0, 0, 1)')))
    fig1.add_trace(go.Scatter(x=[0.75,1.25], y=np.ones(2)*KS_median_lb[1], mode = 'lines', line = dict(color ='rgba(250, 0, 0, 0)')))
    fig1.add_trace(go.Scatter(x=[0.75,1.25], y=np.ones(2)*KS_median_ub[1], mode = 'none', fill='tonexty',fillcolor = 'rgba(250, 0, 0, 0.25)'))
    # Infiltration rate
    fig1.add_trace(go.Scatter(x=[1.75,2.25], y=np.ones(2)*KS_median_m[2], mode = 'lines', line = dict(color ='rgba(105,105,105, 1)')))
    fig1.add_trace(go.Scatter(x=[1.75,2.25], y=np.ones(2)*KS_median_lb[2], mode = 'lines', line = dict(color ='rgba(105,105,105, 0)')))
    fig1.add_trace(go.Scatter(x=[1.75,2.25], y=np.ones(2)*KS_median_ub[2], mode = 'none', fill='tonexty',fillcolor = 'rgba(105,105,105, 0.25)'))
    # Travel time - surface flow (days)
    fig1.add_trace(go.Scatter(x=[2.75,3.25], y=np.ones(2)*KS_median_m[3], mode = 'lines', line = dict(color ='rgba(44,172,206, 1)')))
    fig1.add_trace(go.Scatter(x=[2.75,3.25], y=np.ones(2)*KS_median_lb[3], mode = 'lines', line = dict(color ='rgba(44,172,206, 0)')))
    fig1.add_trace(go.Scatter(x=[2.75,3.25], y=np.ones(2)*KS_median_ub[3], mode = 'none', fill='tonexty',fillcolor = 'rgba(44,172,206, 0.25)'))
    # Travel time - underground flow (days)
    fig1.add_trace(go.Scatter(x=[3.75,4.25], y=np.ones(2)*KS_median_m[4], mode = 'lines', line = dict(color ='rgba(33,76,127, 1)')))
    fig1.add_trace(go.Scatter(x=[3.75,4.25], y=np.ones(2)*KS_median_lb[4], mode = 'lines', line = dict(color ='rgba(33,76,127, 0)')))
    fig1.add_trace(go.Scatter(x=[3.75,4.25], y=np.ones(2)*KS_median_ub[4], mode = 'none', fill='tonexty',fillcolor = 'rgba(33,76,127, 0.25)'))
    # Threshold
    fig1.add_trace(go.Scatter(x=[-1,M],y=[KS_dummy[0],KS_dummy[0]],line=dict(color="black",width=2, dash='dash')))
    
    tick_text  = ["Predator initial population",
                "Predator attack rate",
                "Predator efficiency ratio", 
                "Predator death rate", 
                "Prey initial population"]
    
    fig1.layout.xaxis.range=[-0.5,4.5]
    fig1.update_layout(xaxis = dict(tickmode = 'array',tickvals = [0, 1, 2, 3, 4],ticktext = tick_text))
    fig1.layout.yaxis.range=[0,1]
    fig1.layout.yaxis.title='Sensitivity index'
        
    
    return fig1
    
def mapping():    
    
    
    # ### Step 2: Setup the model
    
    # Define: 
    # - the input factors whose influence is to be analysed with GSA, 
    # - their range of variability (choice made by expert judgement, available data or previous studies),
    # - choice of their distributions.
    
    # In[3]:
    
    
    param_names  = ["Predator initial population (x1000)",
                    "Predator attack rate (attacks per week)",
                    "Predator efficiency ratio", 
                    "Predator death rate (deaths per week)", 
                    "Prey initial population (x1000)"]
    M = len(param_names) # number of parameters
    
    # Parameter distributions:
    distr_fun = st.uniform # uniform distribution
    samp_strat = 'lhs' # Latin Hypercube
    # The shape parameters of the uniform distribution are the lower limit and the difference between lower and upper limits:
    distr_par  = [np.nan] * M
    
    # Define output:
    fun_test = predator_prey_model.function
    
    
    # Range of variability
    
    # In[4]:
    
    
    data = [["-",  2,     10 ,param_names[0]],
            ["-",  0.01,   1, param_names[1]],
            ["-",  0.01,  1, param_names[2]],
            ["-",  0.01,  1, param_names[3]],
            ["-",  2,     10, param_names[4]]]
    
    model_param = pd.DataFrame(data, 
                               columns=["Unit", "Min value", "Max value","Description"],
                               index = param_names)
    model_param
    
    
    # In[5]:
    
    
    outputs  = ["equilibrium error (x1000 individuals)"]
    
    model_outputs = pd.DataFrame(outputs, columns = ['model outputs'])
    model_outputs
    
    
    # In[6]:
    
    
    class setup_model:
        def __init__(self, x1, x2, x3, x4, x5):
            # The shape parameters of the uniform distribution are the lower limit and 
            # the difference between lower and upper limits:
            self.xmin = [x1.value[0], x2.value[0], x3.value[0], x4.value[0], x5.value[0]]
            self.xmax = [x1.value[1], x2.value[1], x3.value[1], x4.value[1], x5.value[1]]
            for i in range(M):
                distr_par[i] = [self.xmin[i], self.xmax[i] - self.xmin[i]]
            self.distr_par = distr_par
    
    
    # ### Step 3: Sample inputs space
    
    # The number of model runs ($N$) typically increases proportionally to the number of input factors ($M$) and will depend on the GSA methods chosen as well. 
    
    # In[7]:
    
    
    N = 10000 # number of samples
    class sample_input:
        def __init__(self,distr_par):
            self.X = AAT_sampling(samp_strat, M, distr_fun, distr_par, N)   
    
    
    # ### Step 4: Run the model
    
    # For each sampled input factors combination, we run the model and save the associated model output.
    
    # In[8]:
    
    
    T = 365 # days
    class run_model:
        def __init__(self,X):
            self.Y = model_execution(fun_test, X, T) 
            
    #### Investigate interactions between input factors
    #    
    #In order to investigate the interactions between input factors we plot one input against the other, coloured by the value taken by the output.
    #
    #Model output = equilibrium error
    
    
    
    # Definition of the sliders
    
    # In[10]:
    
    
    # output = widgets.IntSlider(value = 0, min = 0, max = 1,
    #                               step = 1, description = 'output',
    #                               style = {'description_width': '300px'} ,layout={'width': '700px'},
    #                               continuous_update=False)
    # output.observe(update_figures,names = 'value')
    
    x1 = widgets.IntRangeSlider(value = [model_param['Min value'][0], model_param['Max value'][0]], 
                                min = model_param['Min value'][0], max = model_param['Max value'][0],
                                  step = 1, description = model_param.index[0], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  continuous_update=False)
    
    x2 = widgets.FloatRangeSlider(value = [model_param['Min value'][1], model_param['Max value'][1]], 
                                  min = model_param['Min value'][1], max = model_param['Max value'][1], 
                                  step = 0.01, description = model_param.index[1],
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    
    x3 = widgets.FloatRangeSlider(value = [model_param['Min value'][2], model_param['Max value'][2]], 
                                  min = model_param['Min value'][2], max = model_param['Max value'][2], 
                                  step = 0.01, description = model_param.index[2],
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    
    x4 = widgets.FloatRangeSlider(value = [model_param['Min value'][3], model_param['Max value'][3]], 
                                  min = model_param['Min value'][3], max = model_param['Max value'][3], 
                                  step = 0.01, description = model_param.index[3], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    
    x5 = widgets.IntRangeSlider(value = [model_param['Min value'][4], model_param['Max value'][4]], 
                                min = model_param['Min value'][4], max = model_param['Max value'][4], 
                                  step = 1, description = model_param.index[4], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  continuous_update=False)
    
    # x6 = widgets.FloatRangeSlider(value = [model_param['Min value'][5], model_param['Max value'][5]], 
    #                               min = model_param['Min value'][5], max = model_param['Max value'][5], 
    #                               step = 0.01, description = model_param.index[5], 
    #                               style = {'description_width': '300px'} ,layout={'width': '700px'},
    #                               readout_format = '.1f', continuous_update=False)
    # x6.observe(update_figures,names = 'value')
    
    # x7 = widgets.IntRangeSlider(value = [model_param['Min value'][6], model_param['Max value'][6]], 
    #                             min = model_param['Min value'][6], max = model_param['Max value'][6], 
    #                               step = 1, description = model_param.index[4], 
    #                               style = {'description_width': '300px'} ,layout={'width': '700px'},
    #                               continuous_update=False)
    # x7.observe(update_figures,names = 'value')
    
    #for i in range(M):
    #    
    #    exec('x'+str(i+1)+' = widgets.IntRangeSlider(value = [model_param[\'Min value\']['+str(i)+'], model_param[\'Max value\']['+str(i)+']],min = model_param[\'Min value\']['+str(i)+'], max = model_param[\'Max value\']['+str(i)+'],step = 1, description = model_param.index['+str(i)+'], style = {\'description_width\': \'300px\'} ,layout={\'width\': \'700px\'},continuous_update=False)')
    #    exec('x'+str(i+1)+'.observe(update_figures,names = \'value\')')

    
 
    # In[14]:
    
    x1.value=(2,3)
    x5.value=(2,3)
    distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
    X = sample_input(distr_par).X
    Y = run_model(X).Y
    Y1 = Y.flatten()
    
    Y,X[:, 1] = zip(*sorted(zip(Y1,X[:, 1])))
    Y,X[:, 2] = zip(*sorted(zip(Y1,X[:, 2])))
    Y,X[:, 3] = zip(*sorted(zip(Y1,X[:, 3])))
    Y = np.flipud(Y)
    X = np.flipud(X)
    
    # In[15]:
    
    
    colorscale = 'gist_earth' # black and white plot
    ms=25
    font_size = 12
    
    Labels = param_names
    
    fig = plt.figure(figsize=(10, 10))  
    # Evap rate vs Inf rate
    plt.subplot(2, 2, 1)
    map_plot = plt.scatter(X[:, 1], X[:, 2], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=15)
    plt.xlabel(Labels[1], fontsize = font_size)
    plt.ylabel(Labels[2], fontsize = font_size)
    plt.xlim((np.min(X[:, 1]), np.max(X[:, 1])))
    plt.ylim((np.min(X[:, 2]), np.max(X[:, 2])))
    # Evap rate vs Travel time - surf flow
    plt.subplot(2, 2, 2)
    map_plot = plt.scatter(X[:, 1], X[:, 3], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=15)
    plt.xlabel(Labels[1], fontsize = font_size)
    plt.ylabel(Labels[3], fontsize = font_size)
    plt.xlim((np.min(X[:, 1]), np.max(X[:, 1])))
    plt.ylim((np.min(X[:, 3]), np.max(X[:, 3])))
    # Inf rate vs Travel time - surf flow
    plt.subplot(2, 2, 4)
    map_plot = plt.scatter(X[:, 2], X[:, 3], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=15)
    plt.xlabel(Labels[2], fontsize = font_size)
    plt.ylabel(Labels[3], fontsize = font_size)
    plt.xlim((np.min(X[:, 2]), np.max(X[:, 2])))
    plt.ylim((np.min(X[:, 3]), np.max(X[:, 3])))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.5, wspace=0.5)
    # Create colorbar
    cax = fig.add_axes([0.92, 0.05, 0.02, 0.8]) # Add axes for the colorbar
    cb = plt.colorbar(map_plot, ax=cax, fraction=1, extendfrac=1, extendrect=True)
    cb.set_label(outputs[0], fontsize = font_size)
    cb.ax.tick_params(labelsize=font_size)
    # Make axes of the colorbar invisible
    cax.set_visible(False)  
    
    return fig


# ### References

# 1. [SAFE Website](https://www.safetoolbox.info/)
# 2. [Introductory paper to SAFE - Pianosi et al. (2015)](https://www.sciencedirect.com/science/article/pii/S1364815215001188)
# 3. [PAWN method - Pianosi and Wagener (2018)](https://doi.org/10.1016/j.envsoft.2018.07.019)
# 4. [A review of available methods and workflows for Sensitivity Analysis - Pianosi et al. (2016)](https://www.sciencedirect.com/science/article/pii/S1364815216300287)
# 5. [What has Global Sensitivity Analysis ever done for us? A systematic review to support scientific advancement and to inform policy-making in earth system modelling - Wagener and Pianosi (2019)](https://www.sciencedirect.com/science/article/pii/S0012825218300990)
# 6. [Global Sensitivity Analysis . The primer - Saltelli et al. (2008)](http://www.andreasaltelli.eu/file/repository/A_Saltelli_Marco_Ratto_Terry_Andres_Francesca_Campolongo_Jessica_Cariboni_Debora_Gatelli_Michaela_Saisana_Stefano_Tarantola_Global_Sensitivity_Analysis_The_Primer_Wiley_Interscience_2008_.pdf) 
# 7. [Dummy parameter - Zadeh et al. (2017)](https://www.sciencedirect.com/science/article/pii/S1364815217301159)
