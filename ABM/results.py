from python_src import plotting, analysis
import os


def full_start(exp_dir, nbr_traj, time, timestep):
    return 0


def initialized_start(exp_dir, nbr_species, nbr_traj, time, timestep, labels):
    analysis.length_trajectory_plots(exp_dir, nbr_traj, nbr_spec, time,
                                     timestep, labels)
    analysis.distribution_extinction(exp_dir, nbr_traj, nbr_spec, time,
                                     timestep, labels)
    return 0


if __name__ == '__main__':
    exp_nbr = 203
    exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
    timestep = 1. / 12.
    labels = [['E. Coli A22 strain A', 'E. Coli A22 strain B.'],
              ['C. Subtillus', 'E. Coli'],
              ['E. Coli strain A', 'E. Coli strain B'],
              ['E. Coli', 'E. Coli A22']]

    nbr_traj = 1000
    nbr_spec = 2
    time = 193
    # initialized_start(exp_dir, nbr_spec, nbr_traj, time, timestep, labels[3])
    # full_start(exp_dir, nbr_traj, time)

    sim_nbr = [0, 2, 200, 400, 600, 1800]
    for i in range(len(sim_nbr)):
        sim_dir = plotting.plot_simulation(exp_dir, sim_nbr[i])
        # sim_dir = plotting.plot_simulation_many_init(exp_dir, sim_nbr[i])
        plotting.gif_experiment(sim_dir)
