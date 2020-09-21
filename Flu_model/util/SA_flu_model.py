#!/usr/bin/env python
# coding: utf-8

# # Workflow to perform Global Sensitivity Analysis with the PAWN method 

# This document provides:
# * a brief introduction to Global Sensitivity Analysis (GSA);
# * a workflow to apply GSA using the SAFE (Sensitivity Analysis For Everybody) toolbox (*Refs. 1-2*).
# 
# In this example we apply the PAWN GSA method (*Ref. 3*) to the Ishigami-Homma test function.

# ## What is Global Sensitivity Analysis? and why shall we use it?
# 
# **Global Sensitivity Analysis** is a set of mathematical techniques which investigate how uncertainty in the output of a numerical model can be attributed to variations of its input factors.
# 
# The benefits of applying GSA are:
# 
# 1. **Better understanding of the model**: Evaluate the model behaviour beyond default set-up
#     
# 2. **“Sanity check” of the model**: Test whether the model "behaves well" (model validation)
#     
# 3. **Priorities for uncertainty reduction**: Identify the important inputs on which to focus computationally-intensive calibration, acquisition of new data, etc. 
#     
# 4. **More transparent and robust decisions**: Understand the main impacts of input uncertainty on the modelling outcome and thus model-informed decisions
# 
# 

# ## How Global Sensitivity Analysis works
# 
# GSA investigates how the uncertainty of the selected model input factors influences the variability of the model output.
# 
# An '**input factor**' is any element that can be changed before running the model. In general, input factors could be the equations implemented in the model, set-up choices needed for the model execution on a computer, the model's parameters and forcing input data. 
# 
# An '**output**' is any variable that is obtained after the model execution.
# 
# <br>
# 
# The key steps of GSA are summarised in the figure below.
# 
# Before executing the model, we will sample the inputs from their ranges of variability and then repeatedly run the model so that for each simulation all the inputs vary simultaneously (Input Sampling step). For every output of interest a probability distribution is obtained, after which GSA is performed, to obtain a set of sensitivity indices for each output. The sensitivity indices measure the relative influence of each input factor on the output (*Refs. 4-5*).
# <img src="how_GSA_works.png" width="800px">

# ## What can we learn from GSA?
# 
# In general, GSA can have three specific purposes:
# 
# 1. Ranking (or Factor Prioritization) to rank the input factors based on their relative contribution to the output variability. This allows to prioritise resources to reduce uncertainty, so you know on which input factor(s) to focus. 
# 
# 2. Screening (or Factor Fixing) to identify the input factors, if any, which have a negligible influence on the output variability. The input factors which are found to have negligible influence can be fixed to their default values.
# 
# 3. Mapping to determine the regions of the inputs' variability space which produce output values of interest (e.g. extreme values). For example this can be useful when you want to know which values of your inputs produce an output below or above a threshold of interest.

# ### Step 1: Import python modules

# In[1]:


from __future__ import division, absolute_import, print_function

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as st
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')

from util import PAWN
import util.plot_functions as pf # module to visualize the results
from util.model_execution import model_execution # module to execute the model
from util.sampling import AAT_sampling, AAT_sampling_extend # module to perform the input sampling
from util.util import aggregate_boot  # function to aggregate the bootstrap results
from util.util import aggregate_boot

import util.flu_model as flu_model

def screening():
    # ### Step 2: Setup the model
    
    # Define: 
    # - the input factors whose influence is to be analysed with GSA, 
    # - their range of variability (choice made by expert judgement, available data or previous studies),
    # - choice of their distributions.
    
    # In[2]:
    
    
    param_names  = ["Initial number of vaccinated individuals",
                    "Contact rate per day", 
                    "Contagion ratio", 
                    "Recovery time", 
                    "Vaccination rate per day"]
    M = len(param_names) # number of parameters
    
    # Parameter distributions:
    distr_fun = st.uniform # uniform distribution
    samp_strat = 'lhs' # Latin Hypercube
    # The shape parameters of the uniform distribution are the lower limit and the difference between lower and upper limits:
    distr_par  = [np.nan] * M
    
    # Define output:
    fun_test = flu_model.function
    
    
    # Range of variability
    
    # In[3]:
    
    
    data = [["-",     0,   50   ],
            ["-",     0.3,   2 ],
            ["-",     0.3,   1   ],
            ["days",  7,   21  ],
            ["-",     0,   10 ]]
    model_inputs = pd.DataFrame(data, 
                               columns=["Unit", "Min value", "Max value"],
                               index = param_names)
    
    
    # In[4]:
    
    
    outputs  = ["outbreak peak",
                "cost of actions"]
    
#    model_outputs = pd.DataFrame(outputs, columns = ['model outputs'])
    
    
    # In[5]:
    
    
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
    
    # In[6]:
    
    
    N = 5000 # number of samples
    class sample_input:
        def __init__(self,distr_par):
            self.X = AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
    
    
    # ### Step 4: Run the model
    
    # For each sampled input factors combination, we run the model and save the associated model output.
    
    # In[7]:
    
    
    pop = 100 # population
    T = 100 # days
    time_range = np.arange(1,T+1)
    class run_model:
        def __init__(self,X):
            self.Y = model_execution(fun_test, X, time_range, pop, output.value)
    
    
    # ### Step 5: Apply the PAWN Global Sensitivity Analysis method
    # Let’s now apply Global Sensitivity Analysis: for example the **PAWN** method. 
    # 
    # Its main idea is that the influence of an input factor is proportional to the amount of change in the output distribution produced by fixing that input.
    # <br> The sensitivity of $y$ to $x_{i}$ is measured by the difference between the unconditional CDF of $y$: $F_y( y )$, which is induced by varying all input factors simultaneously, and the conditional CDF that is obtained by varying all inputs but $x_{i}$: $F_{y|x_i}( y | x_i )$.
    
    # ### Step 6: Check model behaviour by visualising input/output samples
    # Scatterplots are plotted to visualise the behaviour of the output over each input factor in turn.
    # 
    # Definition of interactivity
    
    # In[8]:
    
    
    def update_figures(change):
        with fig1.batch_update():
            distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
            X = sample_input(distr_par).X
            Y = run_model(X).Y
            KS_median, _, _, KS_dummy = PAWN.pawn_indices(X, Y, n, dummy = True)
            KS_median, _, _ = PAWN.pawn_indices(X, Y, n, Nboot=1000)
            KS_median_m, KS_median_lb, KS_median_ub = aggregate_boot(KS_median)
            
#            k = 0
#            for i in range(0,M):
#                fig1.data[k].y = np.ones(2)*KS_median_m[i]
#                fig1.data[k+1].y = np.ones(2)*KS_median_lb[i]
#                fig1.data[k+2].y = np.ones(2)*KS_median_ub[i]
#                k = k + 3
#                
##            fig1.data[k].y = [KS_dummy[0],KS_dummy[0]]
#            fig1.layout.title = 'model output: <b>'+outputs[output.value]
    
    
    # Definition of the sliders
    
    # In[9]:
    
    
    output = widgets.IntSlider(value = 0, min = 0, max = 1,
                                  step = 1, description = 'output',
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  continuous_update=False)
    output.observe(update_figures,names = 'value')
    
    x1 = widgets.IntRangeSlider(value = [model_inputs['Min value'][0], model_inputs['Max value'][0]], 
                                  min = model_inputs['Min value'][0], max = model_inputs['Max value'][0],
                                  step = 1, description = model_inputs.index[0], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    x1.observe(update_figures,names = 'value')
    
    x2 = widgets.FloatRangeSlider(value = [model_inputs['Min value'][1], model_inputs['Max value'][1]], 
                                min = model_inputs['Min value'][1], max = model_inputs['Max value'][1], 
                                  step = 0.1, description = model_inputs.index[1],
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  continuous_update=False)
    x2.observe(update_figures,names = 'value')
    
    x3 = widgets.FloatRangeSlider(value = [model_inputs['Min value'][2], model_inputs['Max value'][2]],
                                  min = model_inputs['Min value'][2], max = model_inputs['Max value'][2], 
                                  step = 0.1, description = model_inputs.index[2],
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    x3.observe(update_figures,names = 'value')
    
    x4 = widgets.IntRangeSlider(value = [model_inputs['Min value'][3], model_inputs['Max value'][3]], 
                                min = model_inputs['Min value'][3], max = model_inputs['Max value'][3], 
                                  step = 1, description = model_inputs.index[3], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  continuous_update=False)
    x4.observe(update_figures,names = 'value')
    
    x5 = widgets.IntRangeSlider(value = [model_inputs['Min value'][4], model_inputs['Max value'][4]], 
                                  min = model_inputs['Min value'][4], max = model_inputs['Max value'][4], 
                                  step = 1, description = model_inputs.index[4], 
                                  style = {'description_width': '300px'} ,layout={'width': '700px'},
                                  readout_format = '.1f', continuous_update=False)
    x5.observe(update_figures,names = 'value')
    
    
    # ### Step 8: Plot sensitivity indices and identify non-influential parameters
    
    # The dummy parameter is a numerical artifice, with no influence on the model output, which is used to estimate the threshold for non-influential inputs. 
    # 
    # Uninfluential input factors should have zero-valued sensitivity indices, but since sensitivity indices are computed by numerical approximations rather than analytical solutions, an uninfluential factor may still be associated with a non-zero (although small) index value. 
    # 
    # Therefore if the index of an input factor is below the value of the sensitivity index of the dummy parameter, then the input factor is deemed uninfluential *(Ref. 6)*.
    
    # In[10]:
    
    output.value = 0
    n=10
    distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
    X = sample_input(distr_par).X
    Y = run_model(X).Y
#    KS_median, _, _, KS_dummy = PAWN.pawn_indices(X, Y, n, dummy = True)
    KS_median, _, _ = PAWN.pawn_indices(X, Y, n, Nboot=1000)
    KS_median_m, KS_median_lb, KS_median_ub = aggregate_boot(KS_median)
    
    
    # In[21]:
    
    
    fig1 = go.FigureWidget(layout = dict(width=500, height=400,showlegend = False,margin=dict(t=30,r=0,l=75)))
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
    # fig1.add_trace(go.Scatter(x=[-1,M],y=[KS_dummy[0],KS_dummy[0]],line=dict(color="black",width=2, dash='dash')))
    
    tick_text  = param_names
    
    fig1.layout.xaxis.range=[-0.5,4.5]
    fig1.update_layout(xaxis = dict(tickmode = 'array',tickvals = [0, 1, 2, 3, 4],ticktext = tick_text))
    fig1.layout.yaxis.range=[0,0.5]
    fig1.layout.yaxis.title='Sensitivity index'
    fig1.layout.title = 'model output: <b>'+outputs[output.value]
    
    # In[10]:
    
    output.value = 1
    n=10
    distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
    X = sample_input(distr_par).X
    Y = run_model(X).Y
#    KS_median, _, _, KS_dummy = PAWN.pawn_indices(X, Y, n, dummy = True)
    KS_median, _, _ = PAWN.pawn_indices(X, Y, n, Nboot=1000)
    KS_median_m, KS_median_lb, KS_median_ub = aggregate_boot(KS_median)
    
    
    # In[21]:
    
    
    fig2 = go.FigureWidget(layout = dict(width=500, height=400,showlegend = False,margin=dict(t=30,r=0,l=75)))
    # Soil storage capacity
    fig2.add_trace(go.Scatter(x=[-0.25,0.25], y=np.ones(2)*KS_median_m[0], mode = 'lines', line = dict(color ='rgba(120, 170, 150, 1)')))
    fig2.add_trace(go.Scatter(x=[-0.25,0.25], y=np.ones(2)*KS_median_lb[0], mode = 'lines', line = dict(color ='rgba(120, 170, 150, 0)')))
    fig2.add_trace(go.Scatter(x=[-0.25,0.25], y=np.ones(2)*KS_median_ub[0], mode = 'none', fill='tonexty',fillcolor = 'rgba(120, 170, 150, 0.25)'))
    # Evaporation rate
    fig2.add_trace(go.Scatter(x=[0.75,1.25], y=np.ones(2)*KS_median_m[1], mode = 'lines', line = dict(color ='rgba(250, 0, 0, 1)')))
    fig2.add_trace(go.Scatter(x=[0.75,1.25], y=np.ones(2)*KS_median_lb[1], mode = 'lines', line = dict(color ='rgba(250, 0, 0, 0)')))
    fig2.add_trace(go.Scatter(x=[0.75,1.25], y=np.ones(2)*KS_median_ub[1], mode = 'none', fill='tonexty',fillcolor = 'rgba(250, 0, 0, 0.25)'))
    # Infiltration rate
    fig2.add_trace(go.Scatter(x=[1.75,2.25], y=np.ones(2)*KS_median_m[2], mode = 'lines', line = dict(color ='rgba(105,105,105, 1)')))
    fig2.add_trace(go.Scatter(x=[1.75,2.25], y=np.ones(2)*KS_median_lb[2], mode = 'lines', line = dict(color ='rgba(105,105,105, 0)')))
    fig2.add_trace(go.Scatter(x=[1.75,2.25], y=np.ones(2)*KS_median_ub[2], mode = 'none', fill='tonexty',fillcolor = 'rgba(105,105,105, 0.25)'))
    # Travel time - surface flow (days)
    fig2.add_trace(go.Scatter(x=[2.75,3.25], y=np.ones(2)*KS_median_m[3], mode = 'lines', line = dict(color ='rgba(44,172,206, 1)')))
    fig2.add_trace(go.Scatter(x=[2.75,3.25], y=np.ones(2)*KS_median_lb[3], mode = 'lines', line = dict(color ='rgba(44,172,206, 0)')))
    fig2.add_trace(go.Scatter(x=[2.75,3.25], y=np.ones(2)*KS_median_ub[3], mode = 'none', fill='tonexty',fillcolor = 'rgba(44,172,206, 0.25)'))
    # Travel time - underground flow (days)
    fig2.add_trace(go.Scatter(x=[3.75,4.25], y=np.ones(2)*KS_median_m[4], mode = 'lines', line = dict(color ='rgba(33,76,127, 1)')))
    fig2.add_trace(go.Scatter(x=[3.75,4.25], y=np.ones(2)*KS_median_lb[4], mode = 'lines', line = dict(color ='rgba(33,76,127, 0)')))
    fig2.add_trace(go.Scatter(x=[3.75,4.25], y=np.ones(2)*KS_median_ub[4], mode = 'none', fill='tonexty',fillcolor = 'rgba(33,76,127, 0.25)'))
    # Threshold
    # fig2.add_trace(go.Scatter(x=[-1,M],y=[KS_dummy[0],KS_dummy[0]],line=dict(color="black",width=2, dash='dash')))
    
    fig2.layout.xaxis.range=[-0.5,4.5]
    fig2.update_layout(xaxis = dict(tickmode = 'array',tickvals = [0, 1, 2, 3, 4],ticktext = tick_text))
    fig2.layout.yaxis.range=[0,0.5]
    fig2.layout.yaxis.title='Sensitivity index'
    fig2.layout.title = 'model output: <b>'+outputs[output.value]
    
    return fig1, fig2
    # In[22]:
    
    
#    widgets.VBox([widgets.VBox([widgets.VBox([output,x1,x2,x3,x4,x5]),fig1])])
    
    
    # ### Step 11: Investigate interactions between input factors
    
    # In order to investigate the interactions between input factors we plot one input against the other, coloured by the value taken by the output.
    
    # In[13]:
    
    
    # distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
    # X = sample_input(distr_par).X
    # Y = run_model(X).Y
    # Y = Y.flatten()
    
    
    # In[14]:
    
    
    # plt.rcParams['figure.figsize'] = [10, 10]
    # pf.scatter_plots_interaction(X, Y, ms=5)
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.5, wspace=0.5)
    # plt.show()
    
    
    # In[15]:
    
    
    # Y1 = Y
    # Y,X[:, 0] = zip(*sorted(zip(Y1,X[:, 0])))
    # Y,X[:, 1] = zip(*sorted(zip(Y1,X[:, 1])))
    # Y,X[:, 2] = zip(*sorted(zip(Y1,X[:, 2])))
    # Y,X[:, 3] = zip(*sorted(zip(Y1,X[:, 3])))
    # Y,X[:, 4] = zip(*sorted(zip(Y1,X[:, 4])))
    # Y = np.flipud(Y)
    # X = np.flipud(X)
    
    
    # In[16]:
    
    
    # colorscale = 'jet'
    # ms=10
    # font_size = 12
    
    # Labels = param_names
    
    # fig2 = plt.figure(figsize=(2, 9))    
    
    # plt.subplot(4, 1, 1)
    # map_plot = plt.scatter(X[:, 0], X[:, 1], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=30)
    # # plt.xlabel(Labels[0], fontsize = font_size)
    # plt.ylabel(Labels[1], fontsize = font_size)
    # plt.xlim((np.min(X[:, 0]), np.max(X[:, 0])))
    # plt.ylim((np.min(X[:, 1]), np.max(X[:, 1])))
    # # 
    # plt.subplot(4, 1, 2)
    # map_plot = plt.scatter(X[:, 0], X[:, 2], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=30)
    # # plt.xlabel(Labels[0], fontsize = font_size)
    # plt.ylabel(Labels[2], fontsize = font_size)
    # plt.xlim((np.min(X[:, 0]), np.max(X[:, 0])))
    # plt.ylim((np.min(X[:, 2]), np.max(X[:, 2])))
    # # 
    # plt.subplot(4, 1, 3)
    # map_plot = plt.scatter(X[:, 0], X[:, 3], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=30)
    # # plt.xlabel(Labels[0], fontsize = font_size)
    # plt.ylabel(Labels[3], fontsize = font_size)
    # plt.xlim((np.min(X[:, 0]), np.max(X[:, 0])))
    # plt.ylim((np.min(X[:, 3]), np.max(X[:, 3])))
    # #
    # plt.subplot(4, 1, 4)
    # map_plot = plt.scatter(X[:, 0], X[:, 4], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=30)
    # plt.xlabel(Labels[0], fontsize = font_size)
    # plt.ylabel(Labels[4], fontsize = font_size)
    # plt.xlim((np.min(X[:, 0]), np.max(X[:, 0])))
    # plt.ylim((np.min(X[:, 4]), np.max(X[:, 4])))
    
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.2, wspace=5)
    # # Create colorbar
    # cax = fig2.add_axes([1.1, 0.05, 0.1, 0.9]) # Add axes for the colorbar
    # cb = plt.colorbar(map_plot, ax=cax, fraction=1, extendfrac=1, extendrect=True)
    # cb.set_label(outputs[output.value], fontsize = font_size)
    # cb.ax.tick_params(labelsize=font_size)
    # # Make axes of the colorbar invisible
    # cax.set_visible(False)
    
    
    # ### References
    
    # 1. [SAFE Website](https://www.safetoolbox.info/)
    # 2. [Introductory paper to SAFE - Pianosi et al. (2015)](https://www.sciencedirect.com/science/article/pii/S1364815215001188)
    # 3. [PAWN method - Pianosi and Wagener (2018)](https://doi.org/10.1016/j.envsoft.2018.07.019)
    # 4. [A review of available methods and workflows for Sensitivity Analysis - Pianosi et al. (2016)](https://www.sciencedirect.com/science/article/pii/S1364815216300287)
    # 5. [What has Global Sensitivity Analysis ever done for us? A systematic review to support scientific advancement and to inform policy-making in earth system modelling - Wagener and Pianosi (2019)](https://www.sciencedirect.com/science/article/pii/S0012825218300990)
    # 6. [Global Sensitivity Analysis . The primer - Saltelli et al. (2008)](http://www.andreasaltelli.eu/file/repository/A_Saltelli_Marco_Ratto_Terry_Andres_Francesca_Campolongo_Jessica_Cariboni_Debora_Gatelli_Michaela_Saisana_Stefano_Tarantola_Global_Sensitivity_Analysis_The_Primer_Wiley_Interscience_2008_.pdf) 
    # 7. [Dummy parameter - Zadeh et al. (2017)](https://www.sciencedirect.com/science/article/pii/S1364815217301159)
