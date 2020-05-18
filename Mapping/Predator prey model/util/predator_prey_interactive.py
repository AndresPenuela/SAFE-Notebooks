#!/usr/bin/env python
# coding: utf-8

# # Predator prey model
# 
# This model simulates the evolution in time of the population of two species. One species, the prey, is the primary source of food source for the other, the predator. The predator-prey food chain persists over time thanks to a **negative feedback loop**, i.e. when the product of a reaction leads to a decrease in that reaction. When predators eat preys, they reduce their population but this drop in the predators food source will soon cause the predator population to decline. This will also reduce the predation rate allowing the prey population to increase again. This cycle creates an up and down wavelike pattern that maintain a long-term equilibrium.
# 
# <left><img src="util/predator_prey_equil.gif" width="700px">
# 
# ## Model assumptions
# - The food supply of the predator population depends entirely on the size of the prey population.
# - The rate of change of population is proportional to its size.
# - Predators have limitless appetite.
# - If the prey population grows beyond the carrying capacity of the environment, then their population would be wiped out as all the available food resources would have been consumed.
# - Preys only die as a result of predator attacks.
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

# #### Import libraries

# In[1]:


import numpy as np # import neccissary libraries for using matplotlib
import plotly.graph_objs as go
from ipywidgets import widgets
import util.predator_prey_model as predator_prey_model

def predator_prey_interactive():
    # #### Simulation time range
    
    # In[2]:
    
    
    T = 3000 # days
    time_range = np.arange(1,T)
    
    
    # ### Function to update the simulation when changing the parameters with the sliders
    
    # In[3]:
    
    
    def update_sim(predator_ini,attack_rate,death_rate,efficiency_rate,prey_ini,growth_rate,carrying_capacity):
        predator = predator_prey_model.predator(predator_ini.value,attack_rate.value,death_rate.value,efficiency_rate.value)
        prey = predator_prey_model.prey(prey_ini.value,growth_rate.value)
        environment = predator_prey_model.environment(carrying_capacity.value)
        predator_pop,prey_pop = predator_prey_model.simulation(T,predator,prey,environment)
        return predator_pop,prey_pop
    
    
    # ### Function to update the figure when changing the parameters with the sliders
    
    # In[4]:
    
    
    def update_figure(change):
        with fig.batch_animate(duration=1000):
            fig.data[0].y = update_sim(predator_ini,attack_rate,death_rate,efficiency_rate,prey_ini,growth_rate,carrying_capacity)[0]
            fig.data[1].y = update_sim(predator_ini,attack_rate,death_rate,efficiency_rate,prey_ini,growth_rate,carrying_capacity)[1]
    
    
    # ### Definition of the sliders
    # #### Initial number of predators
    
    # In[5]:
    
    
    predator_ini = widgets.IntSlider(min=1,
                                max=10,
                                value=2, step = 1,
                                description = 'Predator initial pop: ',
                                continuous_update=False,
                                style = {'description_width': '200px'} ,layout={'width': '500px'})
    predator_ini.observe(update_figure,names = 'value')
    
    
    # #### Attack rate of predators
    
    # In[6]:
    
    
    attack_rate = widgets.FloatSlider(min=0.1,
                               max=1,
                               value=0.5, step = 0.01,
                               description = 'Predator attack rate: ',
                               continuous_update=False,
                                style = {'description_width': '200px'} ,layout={'width': '500px'})
    attack_rate.observe(update_figure,names = 'value')
    
    
    # #### Efficiency rate of predators
    
    # In[7]:
    
    
    efficiency_rate = widgets.FloatSlider(min=0.1,
                              max=1,
                              value=0.5, step = 0.01, 
                              description = 'Predator efficiency: ',
                              continuous_update=False,
                                style = {'description_width': '200px'} ,layout={'width': '500px'})
    efficiency_rate.observe(update_figure,names = 'value')
    
    
    # #### Death rate of predators
    
    # In[8]:
    
    
    death_rate = widgets.FloatSlider(min=0.1,
                                max=1,
                                value=0.5, step = 0.01, 
                                description = 'Predator death rate: ',
                                continuous_update=False,
                                style = {'description_width': '200px'} ,layout={'width': '500px'})
    death_rate.observe(update_figure,names = 'value')
    
    
    # #### Initial number of preys
    
    # In[9]:
    
    
    prey_ini = widgets.IntSlider(min=1,
                                max=10,
                                value=5, step = 1,
                                description = 'Prey initial pop: ',
                                continuous_update=False,
                                style = {'description_width': '200px'} ,layout={'width': '500px'})
    prey_ini.observe(update_figure,names = 'value')
    
    
    # #### Growth rate of preys
    
    # In[10]:
    
    
    growth_rate = widgets.FloatSlider(min=0.001,
                                max=2,
                                value=1, step = 0.1, 
                                description = 'Growth rate of preys: ',
                                continuous_update=False,
                                style = {'description_width': '200px'} ,layout={'width': '500px'})
    growth_rate.observe(update_figure,names = 'value')
    
    
    # #### Carrying capacity of environment
    
    # In[11]:
    
    
    carrying_capacity = widgets.IntSlider(min=1,
                                max=20,
                                value=15, step = 1, 
                                description = 'Carrying capacity of environment: ',
                                continuous_update=False,
                                style = {'description_width': '200px'} ,layout={'width': '500px'})
    carrying_capacity.observe(update_figure,names = 'value')
    
    
    # ### Plot the interactive figure
    # #### Initial simulation
    
    # In[12]:
    
    
    predator_pop,prey_pop = update_sim(predator_ini,attack_rate,death_rate,efficiency_rate,prey_ini,growth_rate,carrying_capacity)
    
    
    # #### Figure with two traces: predator and prey populations ###
    
    # In[15]:
    
    
    predator_trace = go.Scatter(x=time_range, y=predator_pop, name='predator', marker_color = 'red')
    prey_trace = go.Scatter(x=time_range, y=prey_pop, name='prey', marker_color = 'blue')
    fig = go.FigureWidget(data   = [predator_trace,prey_trace],
                          layout = go.Layout(xaxis = dict(title = 'days', range = [1,365],
                                                          rangeslider=dict(visible=True,range=[1,T])),
                                             yaxis = dict(range = [0,carrying_capacity.max+1],
                                                          title = 'population'),
                                             margin=dict(t=30)))
    
    
    # #### Plot
    
    # In[16]:
    
    
#    widgets.VBox([widgets.VBox([predator_ini,attack_rate,efficiency_rate,death_rate,prey_ini,growth_rate,carrying_capacity]),
#                  fig])


    return predator_ini,attack_rate,efficiency_rate,death_rate,prey_ini,growth_rate,carrying_capacity,fig




