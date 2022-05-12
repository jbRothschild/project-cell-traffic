from python_src import plotting, analysis
import os


if __name__ == '__main__':
    exp_nbr = 62
    exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'

    nbr_traj = 1600
    nbr_spec = 2
    time = 220
    #analysis.length_trajectory_plots(exp_dir, nbr_traj, nbr_spec, time)
    #analysis.distribution_extinction(exp_dir, nbr_traj, nbr_spec, time)

    sim_nbr = [0, 2, 200, 400, 600, 1800]
    for i in range(len(sim_nbr)):
        sim_dir = plotting.plot_simulation(exp_dir, sim_nbr[i])
        plotting.gif_experiment(sim_dir)
