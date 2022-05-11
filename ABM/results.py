from python_src import plotting, analysis
import os


if __name__ == '__main__':
    exp_nbr = 2
    sim_nbr = 0

    exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'

    sim_dir = plotting.plot_simulation(exp_dir, sim_nbr)
    plotting.gif_experiment(sim_dir)
    """
    nbr_traj = 1600
    nbr_spec = 2
    analysis.length_trajectory_plots(exp_dir, nbr_traj, nbr_spec)
    analysis.distribution_extinction(exp_dir, nbr_traj, nbr_spec)
    """
