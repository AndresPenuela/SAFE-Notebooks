import numpy as np # import neccissary libraries for using matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# for wolf case
prey = 1 # number of prey
preyName = 'sheep' # name of prey
predatorName = 'wolves' # name of predators
predator = 3 # number of predators
period = 0.1 # time interval for calucating differentials 
carryingCapacity = 5 # carrying capacity of environment
preyGrowthRate = 2.0 # how quickly prey are reproducing
attackRate = 2 # ratio of prey consumed per predator
predatorDeathRate = 1 # the mortality rate of predators 
predatorEfficiencyRate = 0.8 # how efficient predators are
preyArray = [] # array for storing prey population values
predatorArray = [] # array for storing predator values 
t =[] # array used for storing time values 
anim = [] # array for allowing multiple animations to run at once 

def preyEquation(P, Q, K, r, s, h): # equation for calculating prey differential
    dP = (r*(1-P/K) -s*Q)*P*h
    return dP

def predatorEquation(P,Q, u,v,h): # equation for calculating predator differential
    dQ = (-u + v*P)*Q*h
    return dQ

def init(): # function for initializing plotting informaiton
    global newPrey, newPredator # keep these variables global so they can be modified in the function
    # initialize populations
    newPrey = prey 
    newPredator = predator
    # setup the parametrs of the graph
    ax.set_xlim(0, 1000) # set the size of the graph
    ax.set_ylim(0,carryingCapacity+5)
    # provide labels
    ax.set_xlabel('Period (t)')
    ax.set_ylabel('Populations')
    ax.set_title('Predator vs Prey Populations of ' +preyName +' and ' +predatorName)
    ax.legend([ln,ln2],[preyName,predatorName]) # provide a legend
    return ln, ln2 # return the values to line 1 and 2 for FuncAnimation

def update_predator_prey(frame):
    global newPrey
    global newPredator
    dP = preyEquation(newPrey,newPredator, carryingCapacity, preyGrowthRate, attackRate, period) # calculate prey differntial
    dQ = predatorEquation(newPrey, newPredator, predatorDeathRate, predatorEfficiencyRate, period) # calculate predator differential
    newPrey = newPrey+dP # add differential to prey value
    newPredator = newPredator+dQ # add differential to predator value
    preyArray.append(newPrey) # add the values to arrays
    predatorArray.append(newPredator)
    t.append(frame) # add time to t array
    ln.set_data(t,preyArray) # set the prey array to line 1 
    ln2.set_data(t,predatorArray) # set the predator array to line 2
    return ln, ln2 # return line 1 and 2

# get user inputs for variables
preyName = input("Enter the type of prey ")
prey = input("Enter the initial number of prey ")
prey = float(prey)
predatorName = input("Enter the type of predator ")
predator = input("Enter the initial number of predators ")
predator = float(predator)

# plotting animation
fig1, ax = plt.subplots() # setup the plot
ln, = plt.plot([], [], 'b', animated = True) # setup line 1
ln2, = plt.plot([], [], 'r', animated = True) # setup line 2 
anim = animation.FuncAnimation(fig1, update_predator_prey, frames = 100000, init_func = init, interval = 20, blit =True) # animate the plot by calling the update and initialize functions


plt.show()#show animation plot


