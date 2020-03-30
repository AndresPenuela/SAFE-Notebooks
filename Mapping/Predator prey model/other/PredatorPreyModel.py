import numpy as np  # import neccissary libraries for using matplotlib
import matplotlib.pyplot as plt
prey = 1 # initial amount of prey population
preyName = 'sheep' # title of the prey population
predatorName = 'wolves' # tile of the predator population
predator = 1 # initial amount of predator population
period = .001 # time period aka how frequently the addition of population differentials happen
preyGrowthRate = 1.3 # how quickly the prey population grows should be positive
attackRate = 0.5 # a ratio of how much prey the predator population consumes
predatorDeathRate = 0.7 # rate of death for predator populations
predatorEfficiencyRate = 1.6 # how efficiently predators reproduce based on prey populations 
preyArray = [] # array for storing the current prey amounts 
predatorArray = [] # array for storing the current predator amounts 

def preyEquation(P, Q, r, s, h): # equation for determining the differential of prey population
    dP = (r -s*Q)*P*h
    return dP

def predatorEquation(P,Q, u,v,h): # equation for determining the differential of predator populations 
    dQ = (-u + v*P)*Q*h
    return dQ

# # Get user inputs for the different variables of this process
# preyName = input("Enter the type of prey ")
# prey = input("Enter the initial number of prey ")
# prey = float(prey)
# predatorName = input("Enter the type of predator ")
# predator = input("Enter the initial number of predators ")
# predator = float(predator)
# preyGrowthRate = input("Enter the growth rate of the prey species (Must be a positive value)")
# preyGrowthRate = float(preyGrowthRate)
# attackRate = input("Enter the attack rate of the predator population (Must be a positive value)")
# attackRate = float(attackRate)
# predatorDeathRate = input("Enter the mortality rate for the predator population (Must be a positive value)")
# predatorDeathRate = float(predatorDeathRate)
# predatorEfficiencyRate = input("Enter the efficiency rate for the predator population (Must be a positive value)")
# predatorEfficiencyRate = float(predatorEfficiencyRate)
# Initialize the prey and predator amounts into the summing variables 
newPrey = prey
newPredator = predator

# For loop which adds the differentials to the prey and predator populations for a given amount of time
for add in range(10000):
    dP = preyEquation(newPrey,newPredator, preyGrowthRate, attackRate, period) # calculates prey differential
    dQ = predatorEquation(newPrey, newPredator, predatorDeathRate,predatorEfficiencyRate, period) # calculates predator differential
    newPrey = newPrey+dP # adds the prey differential to the prey population
    newPredator = newPredator+dQ # adds the predator differential to the predator population
    preyArray.append(newPrey) # stores the new values in the arrays for plotting
    predatorArray.append(newPredator) 

t = range(10000) # x variable of time
plt.plot(t, preyArray, 'b', t, predatorArray, 'r') # plots preyArray with a blue line and the predatorArray with a red line in respect to time 
plt.xlim(0,10000) # set x limits of graph
plt.ylim(0,10) # set y limits of graph

# provide labels for the graph
plt.xlabel('Period (t)') 
plt.ylabel('Populations')
plt.title('Predator vs Prey Populations of ' +preyName + ' and ' +predatorName)

# Displays graph
plt.show()
