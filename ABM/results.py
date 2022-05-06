from python_src import plotting
import os


if __name__ == '__main__':
    exp_nbr = 0
    sim_nbr = 1

    exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
    sim_dir = plotting.plot_simulation(exp_dir, sim_nbr)
    plotting.gif_experiment(sim_dir)
