from src import agent, environment, default, plotting
import os
from pathlib import Path
import numpy as np

largeChannel = {'width' : 45.0, 'height' : 12.0}
smallChannel = {'width' : 45.0, 'height' : 6.0}
hugeChannel = {'width' : 100.0, 'height' : 100.0}

# PARAMETERS
channel = largeChannel # which channel we're using
initPop = [1] # initial population
dt = 0.1 # time step
expNbr = 1 # experiment number
expDir = default.DATA_DIR + os.sep + f'exp_nbr_{expNbr}' # experiment directory
Path(expDir).mkdir(parents=True, exist_ok=True) # make directory

# initialize the environment
sim = environment.BiDirMM(**channel)

# initialize the population within
#sim.setup_initial_nbr_cells(initPop)
sim.place_init_cells()

# TODO : SAVE PARAMETERS TO A FILE, SAVE RESULTS TO A FILE (or many?)

sim.plot_bacteria() # visualize TODO : got to make sure no overlapping!

# simulate
filenames = []
for i in np.arange(0,1000):

    # could alternate or skip some
    filename = expDir + os.sep + f'{i}.png'
    filenames.append(filename)

    # run sim
    sim.step(dt)
    sim.plot_bacteria(filename)

#while len(sim.bacteriaLst) < np.sum(initPop):
    #sim.step(dt)

# plot end result
plotting.gif_experiment(expDir)
#sim.plot_bacteria()
