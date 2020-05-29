#!/usr/bin/env python
# coding: utf-8

# # What parameters are not important in our model?
# In this Notebook we will apply Global Sensitivity Analysis (GSA) to identify the parameters in a model, if any, which have a negligible influence on the model results. In this way, we can reduce the number of parameters and hence the time and money that we need to invest to improve our model. In this example we will apply GSA to a hydrological model.
# 
# In[1]:


# ## What is Global Sensitivity Analysis? and why shall we use it?
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

# ## How Global Sensitivity Analysis works?
# 
# GSA investigates how the uncertainty of the selected model input parameters influences the variability of the model output/performance.
# 
# A '**model parameter**' is any element that can be changed before running the model.
# 
# An '**output**' is any variable that is obtained after the model execution.
# 
# Before executing the model, we will sample the parameters from their ranges of variability and then repeatedly run the model so that for each simulation all the inputs vary simultaneously. After GSA is performed we obtain a set of sensitivity indices for each output. The sensitivity indices measure the relative influence of each input factor on the output (*Refs. 4-5*).

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

from util.HyMOD_gamma import hymod_gamma,hymod_gamma_error

def screening(param_obs):
    # ### Step 2: Setup the model
    
    # Define: 
    # - the forcing input data
    # - the model parameters whose influence is to be analysed with GSA, 
    # - their range of variability (choice made by expert judgement, available data or previous studies),
    # - choice of their distributions.
    
    # #### Model parameters
    
    # In[3]:
    
    
    param_names  = ["Soil_sto" , "Evap_rate" , "Inf_rate", "Time_surf","Time_under"]
    M = len(param_names) # number of parameters
    
    # Parameter distributions:
    distr_fun = st.uniform # uniform distribution
    samp_strat = 'lhs' # Latin Hypercube
    # The shape parameters of the uniform distribution are the lower limit and the difference between lower and upper limits:
    distr_par  = [np.nan] * M
    
    # Define output:
    fun_test = hymod_gamma_error
    
    
    # #### Range of variability
    
    # In[4]:
    
    
    data = [["mm",    10,   90 , "Soil storage capacity"],
            ["-",     0.01, 0.99   , "Evaporation ratio"],
            ["-",     0.01, 0.99, "Infiltration ratio"],
            ["days",  0.8,  2,       "Travel time - surface flow"],
            ["days",  2,    10,      "Travel time - underground flow"]]
    model_param = pd.DataFrame(data, 
                               columns=["Unit", "Min value", "Max value", "Description"],
                               index = param_names)
    model_param
    
    
    # 
    # #### Forcing input data
    
    # In[5]:
    
    
    T = 150 # days
    np.random.seed(1)
    rain = 20 * np.random.random(T)
    np.random.seed(2)
    ept = 5 * np.random.random(T)
    warmup = 31 # days
    
    
    # Observed river flow
    
    # In[6]:
    
    
    model_obs = hymod_gamma(param_obs)
    np.random.seed(3)
    noise = np.random.random(T)*1
    Q_obs = model_obs.simulation(param_obs, rain, ept) + noise
    
    
    # In[7]:
    
    
    outputs  = ["Calibration error"]
    
    model_outputs = pd.DataFrame(outputs, columns = ['model outputs'])
    model_outputs
    
    
    # In[8]:
    
    
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
    
    # In[9]:
    
    
    N = 1000 # number of samples
    class sample_input:
        def __init__(self,distr_par):
            self.X = AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
    
    
    # ### Step 4: Run the model
    
    # For each sampled input factors combination, we run the model and save the associated model output.
    
    # In[10]:
    
    
    class run_model:
        def __init__(self,X):
            self.Y = model_execution(fun_test, X, rain, ept, Q_obs, warmup)

    # ### Step 11: Investigate interactions between input factors
    
    # In order to investigate the interactions between input factors we plot one input against the other, coloured by the value taken by the output.
    
    # In[17]:
    
    
    # ### Step 5: Apply the PAWN Global Sensitivity Analysis method
    # Let’s now apply Global Sensitivity Analysis: for example the **PAWN** method. 
    # 
    # Its main idea is that the influence of an input factor is proportional to the amount of change in the output distribution produced by fixing that input.
    # <br> The sensitivity of $y$ to $x_{i}$ is measured by the difference between the unconditional CDF of $y$: $F_y( y )$, which is induced by varying all input factors simultaneously, and the conditional CDF that is obtained by varying all inputs but $x_{i}$: $F_{y|x_i}( y | x_i )$.
    
    # ### Step 6: Check model behaviour by visualising input/output samples
    # Scatterplots are plotted to visualise the behaviour of the output over each input factor in turn.
    # 
    # Definition of interactivity
    
    # In[11]:
    
    
    def update_figures(change):
        with fig1.batch_update():
            distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
            X = sample_input(distr_par).X
            Y = run_model(X).Y # NSE*100 = simulation accuracy (%)
            KS_median, _, _, KS_dummy = PAWN.pawn_indices(X, Y, n, dummy = True)
    
            for i in range(0,M):
                fig1.data[i].y = [KS_median[i]]
            fig1.data[i+1].y = KS_dummy
    
    
    # Definition of the sliders
    
    # In[12]:
    
    
    x1 = widgets.FloatRangeSlider(value = [model_param['Min value'][0], model_param['Max value'][0]], 
                                  min = model_param['Min value'][0], max = model_param['Max value'][0],
                                  step = 1, description = model_param["Description"][0]+' (mm)', 
                                  style = {'description_width': '250px'} ,layout={'width': '450px'},
                                  readout_format = '.0f', continuous_update=False)
    x1.observe(update_figures,names = 'value')
    
    x2 = widgets.FloatRangeSlider(value = [model_param['Min value'][1], model_param['Max value'][1]], 
                                  min = model_param['Min value'][1], max = model_param['Max value'][1], 
                                  step = 0.01, description = model_param["Description"][1],
                                  style = {'description_width': '250px'} ,layout={'width': '450px'},
                                  readout_format = '.2f', continuous_update=False)
    x2.observe(update_figures,names = 'value')
    
    x3 = widgets.FloatRangeSlider(value = [model_param['Min value'][2], model_param['Max value'][2]], 
                                  min = model_param['Min value'][2], max = model_param['Max value'][2],
                                  step = 0.01, description = model_param["Description"][2],
                                  style = {'description_width': '250px'} ,layout={'width': '450px'},
                                  readout_format = '.2f', continuous_update=False)
    x3.observe(update_figures,names = 'value')
    
    x4 = widgets.FloatRangeSlider(value = [model_param['Min value'][3], model_param['Max value'][3]], 
                                  min = model_param['Min value'][3], max = model_param['Max value'][3],
                                  step = 0.1, description = model_param["Description"][3]+' (days)', 
                                  style = {'description_width': '250px'} ,layout={'width': '450px'},
                                  readout_format = '.1f', continuous_update=False)
    x4.observe(update_figures,names = 'value')
    
    x5 = widgets.FloatRangeSlider(value = [model_param['Min value'][4], model_param['Max value'][4]], 
                                  min = model_param['Min value'][4], max = model_param['Max value'][4], 
                                  step = 0.1, description = model_param["Description"][4]+' (days)', 
                                  style = {'description_width': '250px'} ,layout={'width': '450px'},
                                  readout_format = '.1f', continuous_update=False)
    x5.observe(update_figures,names = 'value')
    
    
    # ### Step 8: Plot sensitivity indices and identify non-influential parameters
    
    # The dummy parameter is a numerical artifice, with no influence on the model output, which is used to estimate the threshold for non-influential inputs. 
    # 
    # Uninfluential input factors should have zero-valued sensitivity indices, but since sensitivity indices are computed by numerical approximations rather than analytical solutions, an uninfluential factor may still be associated with a non-zero (although small) index value. 
    # 
    # Therefore if the index of an input factor is below the value of the sensitivity index of the dummy parameter, then the input factor is deemed uninfluential *(Ref. 6)*.
    
    # In[13]:
    
    
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
    
    
    fig1 = go.FigureWidget(layout = dict(width=500, height=300,showlegend = False,margin=dict(t=10,r=0,l=75)))
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
    fig1.add_shape(dict(type="line",x0=-1,y0=KS_dummy[0],x1=M,y1=KS_dummy[0],line=dict(color="black",width=2, dash='dash')))
    
    fig1.layout.xaxis.range=[-0.5,4.5]
    fig1.update_layout(xaxis = dict(tickmode = 'array',tickvals = [0, 1, 2, 3, 4],ticktext = model_param["Description"]))
    fig1.layout.yaxis.range=[0,1]
    fig1.layout.yaxis.title='Sensitivity index'
    
    
#    trace_dummy = go.Scatter(x=np.arange(M),y=np.array(KS_dummy*M),mode = 'lines')
    
    # In[16]:
    
    return x1,x2,x3,x4,x5,fig1
#    widgets.VBox([widgets.HBox([widgets.VBox([x1,x2,x3,x4,x5]),fig1])])
    
    
def mapping(param_obs,x1,x2,x3,x4,x5):
    # ### Step 2: Setup the model
    
    # Define: 
    # - the forcing input data
    # - the model parameters whose influence is to be analysed with GSA, 
    # - their range of variability (choice made by expert judgement, available data or previous studies),
    # - choice of their distributions.
    
    # #### Model parameters
    
    # In[3]:
    
    
    param_names  = ["Soil_sto" , "Evap_rate" , "Inf_rate", "Time_surf","Time_under"]
    M = len(param_names) # number of parameters
    
    # Parameter distributions:
    distr_fun = st.uniform # uniform distribution
    samp_strat = 'lhs' # Latin Hypercube
    # The shape parameters of the uniform distribution are the lower limit and the difference between lower and upper limits:
    distr_par  = [np.nan] * M
    
    # Define output:
    fun_test = hymod_gamma_error
    
    
    # #### Range of variability
    
    # In[4]:
    
    
    data = [["mm", 10   , 90 , "Soil storage capacity"],
            ["-",  0.01, 0.99   , "Evaporation ratio"],
            ["-",  0.01, 0.99, "Infiltration ratio"],
            ["days",  0.8, 2,       "Travel time - surface flow"],
            ["days",  2, 10,      "Travel time - underground flow"]]
    model_param = pd.DataFrame(data, 
                               columns=["Unit", "Min value", "Max value", "Description"],
                               index = param_names)
    model_param
    
    
    # 
    # #### Forcing input data
    
    # In[5]:
    
    
    T = 150 # days
    np.random.seed(1)
    rain = 20 * np.random.random(T)
    np.random.seed(2)
    ept = 5 * np.random.random(T)
    warmup = 31 # days
    
    
    # Observed river flow
    
    # In[6]:
    
    
    model_obs = hymod_gamma(param_obs)
    np.random.seed(3)
    noise = np.random.random(T)*1
    Q_obs = model_obs.simulation(param_obs, rain, ept) + noise
    
    
    # In[7]:
    
    
    outputs  = ["calibration error"]
    
    model_outputs = pd.DataFrame(outputs, columns = ['model outputs'])
    model_outputs
    
    
    # In[8]:
    
    
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
    
    # In[9]:
    
    
    N = 5000 # number of samples
    class sample_input:
        def __init__(self,distr_par):
            self.X = AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
    
    
    # ### Step 4: Run the model
    
    # For each sampled input factors combination, we run the model and save the associated model output.
    
    # In[10]:
    
    
    class run_model:
        def __init__(self,X):
            self.Y = model_execution(fun_test, X, rain, ept, Q_obs, warmup)

    # ### Step 11: Investigate interactions between input factors
    
    # In order to investigate the interactions between input factors we plot one input against the other, coloured by the value taken by the output.
    
    # In[17]:
    
    
    distr_par = setup_model(x1, x2, x3, x4, x5).distr_par
    X = sample_input(distr_par).X
    Y = run_model(X).Y
    
    
    # In[18]:
#    
#    
#    colorscale = 'gist_earth' # black and white plot
#    ms=15
#    font_size = 12
#    Y = Y.flatten() # shape (N, )
#    Labels = [np.nan]*M
#    for i in range(M):
#        Labels[i] = param_names[i]
#    
#    fig = plt.figure(figsize=(10, 10))    
#    k = 1
#    for i in range(M-1):
#        for j in range(i+1, M, 1):
#            plt.subplot(M-1, M-1, k+i)
#            map_plot = plt.scatter(X[:, i], X[:, j], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=400)
#            plt.xlabel(Labels[i], fontsize = font_size)
#            plt.ylabel(Labels[j], fontsize = font_size)
#            plt.xlim((np.min(X[:, i]), np.max(X[:, i])))
#            plt.ylim((np.min(X[:, j]), np.max(X[:, j])))
#            plt.xticks([])
#            plt.yticks([])
#            k = k + 1
#        k = k + i
#    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.5, wspace=0.5)
#    # Create colorbar
#    cax = fig.add_axes([0.92, 0.05, 0.02, 0.8]) # Add axes for the colorbar
#    cb = plt.colorbar(map_plot, ax=cax, fraction=1, extendfrac=1, extendrect=True)
#    cb.set_label(outputs[0], fontsize = font_size)
#    cb.ax.tick_params(labelsize=font_size)
#    # Make axes of the colorbar invisible
#    cax.set_visible(False)
    # In[19]:
    
    colorscale = 'gist_earth' # black and white plot
    ms=25
    font_size = 12
    Y1 = Y.flatten() # shape (N, )
    
    Y,X[:, 1] = zip(*sorted(zip(Y1,X[:, 1])))
    Y,X[:, 2] = zip(*sorted(zip(Y1,X[:, 2])))
    Y,X[:, 3] = zip(*sorted(zip(Y1,X[:, 3])))
    Y = np.flipud(Y)
    X = np.flipud(X)
    
    Labels = model_param['Description']
    
    fig = plt.figure(figsize=(8, 8))  
    # Evap rate vs Inf rate
    plt.subplot(2, 2, 1)
    map_plot = plt.scatter(X[:, 1], X[:, 2], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=400)
    plt.xlabel(Labels[1], fontsize = font_size)
    plt.ylabel(Labels[2], fontsize = font_size)
    plt.xlim((np.min(X[:, 1]), np.max(X[:, 1])))
    plt.ylim((np.min(X[:, 2]), np.max(X[:, 2])))
    # Evap rate vs Travel time - surf flow
    plt.subplot(2, 2, 2)
    map_plot = plt.scatter(X[:, 1], X[:, 3], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=400)
    plt.xlabel(Labels[1], fontsize = font_size)
    plt.ylabel(Labels[3], fontsize = font_size)
    plt.xlim((np.min(X[:, 1]), np.max(X[:, 1])))
    plt.ylim((np.min(X[:, 3]), np.max(X[:, 3])))
    # Inf rate vs Travel time - surf flow
    plt.subplot(2, 2, 4)
    map_plot = plt.scatter(X[:, 2], X[:, 3], s=ms, c=Y, cmap=colorscale,vmin=0, vmax=400)
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
