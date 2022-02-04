from src import agent, environment, default, plotting
import os
from pathlib import Path
import numpy as np
from tqdm import tqdm

# PARAMETERS
initPop = [1000] # initial population
dt = 1E-4 # time step, minutes
saveTime = 10.0 # every save timestep
expNbr = 31 # experiment number
expDir = default.DATA_DIR + os.sep + f'exp_nbr_{expNbr}' # experiment directory
Path(expDir).mkdir(parents=True, exist_ok=True) # make directory

# initialize the environment
np.random.seed(20) #1
gS = 4.0 # TODO change to max length of cells in simulation
largeChannel = {'width' : 45.0, 'height' : 12.0, 'gridSize' : gS}
smallChannel = {'width' : 45.0, 'height' : 6.0, 'gridSize' : gS}
hugeChannel = {'width' : 100.0, 'height' : 100.0, 'gridSize' : gS}
channel = hugeChannel # which channel we're using
sim = environment.BiDirMM(expDir=expDir, **channel)

# initialize the population within
sim.setup_initial_nbr_cells(initPop)
#sim.place_init_cells()

# TODO : SAVE PARAMETERS TO A FILE, SAVE RESULTS TO A FILE (or many?)

sim.plot_bacteria() # visualize TODO : got to make sure no overlapping!

# simulate
filenames = []; j=0
T = int(24*60/dt)
for i in tqdm(np.arange(0,T)):
    # run sim
    if (i % (saveTime//dt) == 0.0):
        # could alternate or skip some
        filename = expDir + os.sep + f'{j}.png'
        filenames.append(filename)
        sim.plot_bacteria(filename)
        j += 1
    sim.step(dt)

#while len(sim.bacteriaLst) < np.sum(initPop):
    #sim.step(dt)

# plot end result
plotting.gif_experiment(expDir)
#sim.plot_bacteria()
