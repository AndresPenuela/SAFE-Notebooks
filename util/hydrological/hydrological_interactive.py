#!/usr/bin/env python
# coding: utf-8

# # What is a hydrological model?
# Imagine that we want to simulate the natural water flows draining into a river. Knowing the amount of rainfall that has fallen in the river drainage area (watershed) we can estimate how much of the rainwater flows over (surface flow) and under (subsurface flow) the land surface and finally reach the river. 
# 
# <left><img src="util/Watershed.jpg" width="500px">
# 
# For that purpose, we can use a hydrological model. A hydrological model is a mathematical model (set of equations) describing the hydrological processes that occur in a watershed as a function of various parameters. These model parameters describe the hydrological characteristics of the watershed, such as the climate and soil characteristics, ultimately enabling the estimation of the river flow at selected river sections.
# 
# ## Hydrological processes in a river section
# This diagram represents the main hydrological processes that occur in one of the river sections. 
# 
# <left><img src="util/Rain_runoff_diagram7.gif" width="700px">
#     
# ## Diagram of the hydrological model
# 
# <left><img src="util/HyMOD_diagram_simple4.gif" width="1000px">
# 
# #### Import libraries

# In[1]:


import numpy as np
import pandas as pd
import plotly.graph_objs as go
from ipywidgets import widgets, Layout

from util.hydrological.HyMOD_gamma import hymod_gamma,hymod_gamma_error

def hydrological_model():
    # ## Model parameters
    # To represent the hydrological processes the model use mathematical equations that are as a function of several parameters.
    # 
    # <left><img src="util/Rain_runoff_param.png" width="500px">
    
    # In[2]:
    
    
    data = [["mm", 10   , 90 , "Soil storage capacity"],
            ["-",  0.01, 0.99, "Evaporation ratio"],
            ["-",  0.01, 0.99, "Infiltration ratio"],
            ["-",  0.8, 2,       "Travel time - surface flow"],
            ["-",  2, 10,      "Travel time - underground flow"]]
    model_param = pd.DataFrame(data, 
                               columns=["Unit", "Min value", "Max value", "Description"],
                               index = ["Soil_sto" , "Evap_rate" , "Inf_rate",
                                        "Time_surf","Time_under"])
    #model_param
    
    
    # ## Interactive manual calibration
    # ### Simulation time range
    
    # In[3]:
    
    
    T = 150 # days
    dates = pd.date_range(start = '2000-01-01', periods = T)
    warmup = 31
    
    
    # ### Inputs
    
    # In[4]:
    
    np.random.seed(1)
    rain = 20 * np.random.random(T)
    np.random.seed(2)
    ept = 5 * np.random.random(T)
    
    
    # ### Function to update the Hymod simulation when changing the parameters with the sliders
    
    # In[5]:
    
    
    def update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under):
        param = np.array([Soil_sto.value, Evap_rate.value, Inf_rate.value, Time_surf.value, Time_under.value])
        model = hymod_gamma(param)
        Q_sim = model.simulation(param, rain, ept)[warmup:]
        ER = model.Pe[warmup:]
        ER_q = ER*Inf_rate.value
        ER_s = ER*(1-Inf_rate.value)
        Q_s = model.QsL[warmup:]
        Q_q = model.QsF[warmup:]
        evap = model.Ea[warmup:]
        error = hymod_gamma_error(param, rain, ept, Q_obs, warmup)
        return ER, ER_q, ER_s, Q_q, Q_s, Q_sim, evap, error
    
    
    # ### Function to update the figure when changing the parameters with the sliders
    
    # In[6]:
    
    
    def update_figure(change):
        with fig_hyd.batch_animate(duration=1000):
            fig_hyd.data[0].y = update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under)[5]
            fig_hyd.layout.title = "calibration error = "+\
            str("%.0f" % (update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under)[7]))+' ML'
        with fig_sto.batch_animate(duration=1000):
            fig_sto.data[0].y = np.ones(T+1)*Soil_sto.value
            fig_sto.data[1].line.width = np.mean(update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under)[6])*5
            fig_sto.data[2].marker.size = np.mean(update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under)[6])*10
        with fig_flo.batch_animate(duration=1000):
            fig_flo.data[0].y = Time_surf.value/5 * np.sin(2 * np.pi * (freq_f * x_f + phase_f)) + 8
            fig_flo.data[0].line.width = np.mean(update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under)[3])*2
            fig_flo.data[1].y = Time_under.value/5 * np.sin(2 * np.pi * (freq_s * x_f + phase_s)) + 4
            fig_flo.data[1].line.width = np.mean(update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under)[4])*2
            fig_flo.data[2].marker.size = np.mean(update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under)[3])*6
            fig_flo.data[3].marker.size = np.mean(update_sim(Soil_sto, Evap_rate, Inf_rate, Time_surf, Time_under)[4])*6        
    
    
    # ### Definition of the sliders
    # #### Soil_sto: Maximum soil moisture storage capacity (mm)
    
    # In[7]:
    
    
    Soil_sto = widgets.FloatSlider(min=model_param.loc['Soil_sto','Min value'],
                                max=model_param.loc['Soil_sto','Max value'],
                                value=50, step = 1,
                                description = model_param['Description'][0]+' (mm)',
                                continuous_update=False,
                              style = {'description_width': '250px'} ,layout={'width': '500px'})
    Soil_sto.observe(update_figure,names = 'value')
    
    
    # #### Evap_rate
    
    # In[8]:
    
    
    Evap_rate = widgets.FloatSlider(min=model_param.loc['Evap_rate','Min value'],
                               max=model_param.loc['Evap_rate','Max value'],
                               value=0.5, step = 0.01,
                               description = model_param['Description'][1],
                               continuous_update=False,
                              style = {'description_width': '250px'} ,layout={'width': '500px'})
    Evap_rate.observe(update_figure,names = 'value')
    
    
    # #### Inf_rate
    
    # In[9]:
    
    
    Inf_rate = widgets.FloatSlider(min=model_param.loc['Inf_rate','Min value'],
                                max=model_param.loc['Inf_rate','Max value'],
                                value=0.5, step = 0.01, 
                                description = model_param['Description'][2],
                                continuous_update=False,
                              style = {'description_width': '250px'} ,layout={'width': '500px'})
    Inf_rate.observe(update_figure,names = 'value')
    
    
    # #### Time_surf
    
    # In[10]:
    
    
    Time_surf = widgets.FloatSlider(min=model_param.loc['Time_surf','Min value'],
                              max=model_param.loc['Time_surf','Max value'],
                              value=1.4, step = 0.01,
                              description = model_param['Description'][3]+' (days)',
                              continuous_update=False,
                              style = {'description_width': '250px'} ,layout={'width': '500px'})
    Time_surf.observe(update_figure,names = 'value')
    
    
    # #### Time_under
    
    # In[11]:
    
    
    Time_under = widgets.FloatSlider(min=model_param.loc['Time_under','Min value'],
                              max=model_param.loc['Time_under','Max value'],
                              value= 6, step = 0.1, 
                              description = model_param['Description'][4]+' (days)',
                              continuous_update=False,
                              style = {'description_width': '250px'} ,layout={'width': '500px'})
    Time_under.observe(update_figure,names = 'value')
    
    
    # ### Observed hydrograph

    
    param_obs = [80,0.35,0.30,1.80,5.65]
    model_obs = hymod_gamma(param_obs)
    np.random.seed(3)
    noise = np.random.random(T)*1
    Q_obs = model_obs.simulation(param_obs, rain, ept) + noise
    
    
    # ### Plot the interactive figure
    # #### Initial simulation
    
    # In[15]:
    
    
    param = np.array([Soil_sto.value, Evap_rate.value, Inf_rate.value, Time_surf.value, Time_under.value])
    model = hymod_gamma(param)
    Q_sim = model.simulation(param, rain, ept)
#    ER_q_sim = model.Pe * (1-Inf_rate.value)
#    ER_s_sim = model.Pe * Inf_rate.value
    
    
    # #### Figure: storage capacity
    
    # In[16]:
    
    
    sto_trace = go.Scatter(x=np.linspace(0, 1, num=T+1), y=np.ones(T+1)*Soil_sto.value, name=None, fill="tozeroy", fillcolor = 'rgba(120, 170, 150, 1)',mode='none')
    evap_line = go.Scatter(x=[0.8,0.8], y=[75,95], mode = 'lines', line = dict(color = 'red', width = model.Ea.mean()*5),opacity = 0.9)
    evap_head = go.Scatter(x=[0.8], y=[95], mode = 'markers', opacity = 0.9, marker = dict(symbol = 'triangle-up', size = model.Ea.mean()*10, color = 'red'))
    sto_layout = go.Layout(xaxis = dict(showticklabels=False,range = [0,1], showgrid = False),
                           yaxis = dict(range = [0,100],
                                        showticklabels=False, showgrid = False),
                           width=150, height=175, margin=dict(l=0,r=20,t=50,b=10),plot_bgcolor=None, showlegend=False,
                           annotations = [dict(x = 0.13, y = 6, text = '<b>Soil</b>',showarrow=False,font = dict(size=15)),
                                          dict(x = 0.28, y = 94, text = '<b>Effec rain</b>',showarrow=False,font = dict(size=15)),
                                          dict(x = 1.00, y = 95, text = '<b>Evap</b>',showarrow=False,font = dict(size=15,color = 'red'))])
    fig_sto = go.FigureWidget(data   = [sto_trace,evap_line,evap_head],
                              layout = sto_layout)
    
    
    # #### Figure: quick and slow flow
    
    # In[17]:
    
    
    # fast flow
    amp_f = Time_surf.value/5 # amplitude
    phase_f = 0 # phase
    freq_f = 20 # frequency
    x_f = np.linspace(0,15,150) # x axis from 0 to 15 with a 1/150 step
    y_f = amp_f * np.sin(2 * np.pi * (freq_f * x_f + phase_f)) + 8
    
    # slow flow
    amp_s = Time_under.value/5 # amplitude
    phase_s = 0 # phase
    freq_s = 20 # frequency
    x_s = np.linspace(0,15,150) # x axis from 0 to 15 with a 1/150 step
    y_s = amp_s * np.sin(2 * np.pi * (freq_s * x_s + phase_s)) + 4
    
    
    # In[18]:
    
    
    Q_f_sin_line = go.Scatter(x=x_f, y=y_f, mode = 'lines', line = dict(color = 'rgba(0,176,240, 1)', width = model.QsF.mean()*2))
    Q_s_sin_line = go.Scatter(x=x_s, y=y_s, mode = 'lines', line = dict(color = 'rgba(33,76,127, 1)', width =  model.QsL.mean()*2))
    
    Q_f_head = go.Scatter(x=[15.5],   y=[8], mode = 'markers', marker = dict(symbol = 'triangle-right', size = model.QsF.mean()*6, color = 'rgba(0,176,240, 1)'))
    Q_s_head = go.Scatter(x=[15.5],   y=[4], mode = 'markers', marker = dict(symbol = 'triangle-right', size = model.QsL.mean()*6, color = 'rgba(33,76,127, 1)'))
    
    flo_layout = go.Layout(width=250, height=200, margin=dict(l=0,t=50,b=0,r=0),plot_bgcolor='white',showlegend=False,
                           xaxis = dict(range = [0,18],showticklabels=False),
                           yaxis = dict(range = [0,10],showticklabels=False),
                           annotations = [dict(x = 10, y = 9, text = '<b>surface flow</b>',showarrow=False,font = dict(size=15)),
                                          dict(x = 10, y = 1, text = '<b>underground flow</b>',showarrow=False,font = dict(size=15))])
    
    fig_flo = go.FigureWidget(data   = [Q_f_sin_line,Q_s_sin_line,Q_f_head,Q_s_head],
                              layout = flo_layout)
    
    
    # #### Figure: hydrographs (with two traces: simulated and observed hydrogrphs)
    
    # In[19]:
    
    
    sim_hyd = go.Scatter(x=dates[warmup:], y=Q_sim[warmup:], name='sim hyd',line = dict(color='blue'))
    obs_hyd = go.Scatter(x=dates[warmup:], y=Q_obs[warmup:], name='obs hyd',line = dict(color='darkgrey'))
    error = hymod_gamma_error(param, rain, ept, Q_obs, warmup)
    fig_hyd = go.FigureWidget(data   = [sim_hyd,obs_hyd],
                              layout = go.Layout(xaxis = dict(title = '<b>date</b>'),
                                                 yaxis = dict(range = [0,20],title = '<b>River flow (ML/day)</b>'),
                                                 height = 250,width=500,margin=dict(t=50,r=0,l=0),
                                                 legend = dict(x=0,y=1,bgcolor='rgba(0,0,0,0)'),
                                                 title = "calibration error = "+str("%.0f" % (error))+' ML'))
    
    
    # #### Plot
    
    # In[20]:
    
    
    slider_box_layout = Layout(border='solid',width='650px')
    hbox_layout = Layout(display='inline-flex',height='250px')
    vbox_layout = Layout(align_items = 'flex-start')
#    widgets.VBox([widgets.VBox([Soil_sto,Evap_rate,Inf_rate, Time_surf,Time_under],layout=slider_box_layout),
#                  widgets.HBox([fig_sto,fig_flo,fig_hyd],layout=hbox_layout)],layout = vbox_layout)
    
    
    return Soil_sto,Evap_rate,Inf_rate, Time_surf,Time_under,slider_box_layout,\
           fig_sto,fig_flo,fig_hyd,hbox_layout,vbox_layout,\
           param_obs




