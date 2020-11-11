#!/usr/bin/env python
# coding: utf-8

# # The flu model
# 
# <left><img src="util/flu_image.png" width="400px">
#     
# The flu model is a simple mathematical description of the spread of a flu in a small town. The model divides the population of N individuals into four "compartments" which may vary as a function of time, t:
# 
# - **Vulnerable**: people vulnerable but not yet infected with the flu;
# 
# - **Sick**: individuals who are infected with the flu;
# 
# - **Immune**: individuals have immunity to the flu. This includes individuals who either have recovered from the flu or have been vaccinated.
# 
# ### Model parameters
# 
# The model describes the change in the population of each of these compartments in terms of five parameters: 
# 
# - *Initial number of immune individuals*: people who are immune at the start of the infection because they have been previously vaccinated
# - *Contact rate per day*: number of times that an infected individual comes into contact with vulnerable individuals in a day
# - *Contagion ratio*: proportion of contacts that get infected
# - *Recovery time*: the number of days during which a sick individual can pass it on
# - *Vaccination rate*: number of inviduals who are vaccinated per day# Workflow to perform Global Sensitivity Analysis with the PAWN method 
# 
# #### Import libraries

# In[1]:


import numpy as np
import plotly.graph_objs as go
from ipywidgets import widgets

import util.flu.flu_model as flu_model

def flu_interactive():
    # #### Simulation time range
    
    # In[2]:
    
    
    T = 100 # days
    time_range = np.arange(1,T+1)
    pop = 100
    
    
    # ## Interactive figure
    # ### Function to update the simulation when changing the parameters with the sliders
    
    # In[3]:
    
    
    def update_sim(RI_0, contact, contagion, recovery, vaccination):
        param = np.array([RI_0.value/1000, contact.value, contagion.value, recovery.value, vaccination.value/1000])
        S, RI, V,S_max,cost = flu_model.model(param,time_range, pop)
        return S, RI, V, S_max, cost
    
    
    # ### Function to update the figure when changing the parameters with the sliders
    
    # In[4]:
    
    
    def update_figure(change):
        with fig_flu.batch_animate(duration=0):
            fig_flu.data[0].y = update_sim(RI_0, contact, contagion, recovery, vaccination)[0]*1000
            fig_flu.data[1].y = update_sim(RI_0, contact, contagion, recovery, vaccination)[1]*1000
            fig_flu.data[2].y = update_sim(RI_0, contact, contagion, recovery, vaccination)[2]*1000
            fig_flu.layout.title = 'Peak = '+f"{update_sim(RI_0, contact, contagion, recovery, vaccination)[3]*1000:,.0f}"+        ' sick individuals - Cost = '+f"{update_sim(RI_0, contact, contagion, recovery, vaccination)[4]*1000:,.0f}"+' £'
    
    
    # ### Definition of the sliders
    # #### Ratio of population who is initially immune
    
    # In[5]:
    
    
    # Initial number of immune individuals
    RI_0 = widgets.IntSlider(min = 0,
                            max = 50000,
                            value = 0, step = 1,
                            description = 'Initial number of vaccinated individuals: ',
                            continuous_update=True,
                            style = {'description_width': '300px'} ,
                            layout={'width': '700px'})
    RI_0.observe(update_figure,names = 'value')
    
    
    # #### Contact rate per day
    
    # In[6]:
    
    
    contact = widgets.FloatSlider(min = 0.3,
                               max = 2,
                               value = 2, step = 0.01,
                               description = 'Contact rate per day: ',
                               continuous_update=True,
                               style = {'description_width': '300px'} ,
                               layout={'width': '700px'})
    contact.observe(update_figure,names = 'value')
    
    
    # #### Contagion ratio
    
    # In[7]:
    
    
    contagion = widgets.FloatSlider(min = 0.3,
                               max = 1,
                               value = 1, step = 0.01,
                               description = 'Contagion ratio: ',
                               continuous_update=True,
                               style = {'description_width': '300px'} ,
                               layout={'width': '700px'})
    contagion.observe(update_figure,names = 'value')
    
    
    # #### Recovery time in days
    
    # In[8]:
    
    
    recovery = widgets.IntSlider(min = 7,
                               max = 21,
                               value = 21, step = 1,
                               description = 'Recovery time (days): ',
                               continuous_update=True,
                               style = {'description_width': '300px'} ,
                               layout={'width': '700px'})
    recovery.observe(update_figure,names = 'value')
    
    
    # #### Vaccination rate per day
    
    # In[9]:
    
    
    vaccination = widgets.IntSlider(min = 0,
                               max = 10000,
                               value = 0, step = 1,
                               description = 'Vaccination rate per day: ',
                               continuous_update=True,
                               style = {'description_width': '300px'} ,
                               layout={'width': '700px'})
    vaccination.observe(update_figure,names = 'value')
    
    
    # ### Plot the interactive figure
    # #### Initial simulation
    
    # In[10]:
    
    
    param = np.array([RI_0.value, contact.value, contagion.value, recovery.value, vaccination.value])
    S, RI, V, S_max, cost = flu_model.model(param,time_range, pop)
    
    
    # #### Plot the data on four separate curves for V(t), S(t), RI(t) and D(t)
    
    # In[11]:
    
    
    NHS_capacity = 40*1000
    S_trace  = go.Scatter(x=time_range, y=S*1000,  name='Sick',  mode='lines',line=dict(width=0, color = 'red'),stackgroup='one')
    RI_trace = go.Scatter(x=time_range, y=RI*1000, name='Immune', mode='lines',line=dict(width=0, color = 'green'),stackgroup='one')
    V_trace  = go.Scatter(x=time_range, y=V*1000,  name='Vulnerable', mode='lines',line=dict(width=0, color = 'blue'),stackgroup='one')
    NHS_trace = go.Scatter(x=[0,T+1],y=[NHS_capacity,NHS_capacity], mode = 'lines',name='40% of population',line=dict(color="darkblue",width=1, dash='dash'))
    fig_flu  = go.FigureWidget(data   = [S_trace,RI_trace,V_trace,NHS_trace],
                              layout = go.Layout(xaxis = dict(range = [1,100],title = 'days',rangeslider=dict(visible=False,range=[1,T])),
                                                 yaxis = dict(range = [0,pop*1000],title = 'population'),
                                                 height = 400, title = 'Peak = '+f"{S_max*1000:,.0f}"+' sick individuals - Cost = '+f"{cost*1000:,.0f}"+' £'))
    
    
    # #### Plot
    
    # In[12]:
    
    
#    widgets.VBox([widgets.VBox([RI_0, contact, contagion, recovery, vaccination]),
#                  fig_flu])
    
    
    return RI_0,contact, contagion,recovery,vaccination,fig_flu




