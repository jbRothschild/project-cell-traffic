from python_src import plotting, analysis
import numpy as np
import os


if __name__ == '__main__':
    sim_nbr = 0
    for exp_nbr in np.arange(14, 15):
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        sim_dir = plotting.plot_simulation(exp_dir, sim_nbr)
        plotting.gif_experiment(sim_dir)
