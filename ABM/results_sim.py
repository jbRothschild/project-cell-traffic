from python_src import plotting, analysis
import os
import matplotlib.pyplot as plt


def full_start(exp_dir, nbr_traj, time, timestep):
    return 0


def initialized_start(exp_dir, nbr_traj, time, timestep=1. / 12., labels=None):
    analysis.length_trajectory_plots(exp_dir, nbr_traj, time,
                                     timestep, labels=labels)
    analysis.distribution_extinction(exp_dir, nbr_traj, time,
                                     timestep, labels=labels)
    """
    analysis.distribution_extinction(exp_dir, nbr_traj, time,
                                     timestep, labels=labels)
    """

    sim_nbr = [1, 200, 2000]

    for i in range(len(sim_nbr)):
        sim_dir = plotting.plot_simulation(exp_dir, sim_nbr[i])
        # sim_dir = plotting.plot_simulation_many_init(exp_dir, sim_nbr[i])
        plotting.gif_experiment(sim_dir)
    return 0


def plot_gif(exp_dir, sim_nbr):

    sim_dir = plotting.plot_simulation(exp_dir, sim_nbr)
    plotting.gif_experiment(sim_dir)


def richness_N_strains(dir_list, bacteria_type, nbr_traj):
    nbr_init_species = []
    av_fraction_species = []
    std_over_mean_species = []
    for exp_nbr in dir_list:
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        nbr_species, av_rich, std_rich = analysis.richness_traj_plots(exp_dir,
                                                                      nbr_traj)
        nbr_init_species.append(nbr_species)
        av_fraction_species.append(av_rich[-1] / nbr_species)
        std_over_mean_species.append(std_rich[-1] / av_rich[-1])

    fig, ax1 = plt.subplots(1)
    color = 'tab:blue'
    ax1.scatter(nbr_init_species, av_fraction_species, color=color)
    ax1.set_title(r"richness after 12h")
    ax1.set_ylabel(r'$\langle S^* \rangle /S$', color=color)
    ax1.set_ylim(0.0, 1.0)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:red'
    ax2.set_ylabel(
        r'coefficient of variation, SD($S^*\mathrm{)}/\langle S^* \rangle$',
        color=color)  # we already handled the x-label with ax1
    ax2.scatter(nbr_init_species, std_over_mean_species, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    # ax2.set_ylim(0.0, 1.0)

    # plt.xlim([0.0, max_t])
    plt.xlabel(r'number of initial species, $S$')
    fig.tight_layout()
    filename = bacteria_type + "_n_strain_richness"
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".png")
    # plt.show()

    return av_rich, std_rich


if __name__ == '__main__':
    exp_nbr = 42
    exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
    # plot_gif(exp_dir, 5)
    labels = [['E. Coli A22 strain A', 'E. Coli A22 strain B.'],
              ['C. Subtillus', 'E. Coli'],
              ['E. Coli strain A', 'E. Coli strain B'],
              ['E. Coli', 'E. Coli A22']]

    nbr_traj = 3200
    nbr_spec = 2
    max_time = None

    N_strains_a22 = [103, 113, 72, 73, 74, 75, 76]
    N_strains_wt = [102, 112, 62, 63, 64, 65, 66]
    N_strains_a22 = [103, 113]
    N_strains_wt = [102, 112]
    N_strains_wt = [62]
    N_strains_wt = [25]
    N_strains_vs = [82, 83, 84, 85, 86]

    """
    for exp_nbr in N_strains_wt:
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        plot_gif(exp_dir, 1)
        initialized_start(exp_dir, nbr_traj, max_time, labels=labels[3])
    """
    """
    initialized_start(os.getcwd() + os.sep + 'data' + os.sep + 'c_exp_112',
                      3200, max_time, labels=['WT lineage 1', 'WT lineage 2'])
    initialized_start(os.getcwd() + os.sep + 'data' + os.sep + 'c_exp_113',
                      3200, max_time,
                      labels=['A22 lineage 1', 'A22 lineage 2'])
    initialized_start(os.getcwd() + os.sep + 'data' + os.sep + 'c_exp_5',
                      999, max_time, labels=['A22', 'WT'])
    """
    data_folders = [os.getcwd() + os.sep + 'data' + os.sep + 'c_exp_112',
                    os.getcwd() + os.sep + 'data' + os.sep + 'c_exp_113',
                    os.getcwd() + os.sep + 'data' + os.sep + 'c_exp_5'
                    ]
    nbr_simulations = [3200, 3200, 999]
    labels = ['WT', 'A22', 'A22-WT']
    analysis.bar_chart_fixations(os.getcwd() + os.sep + 'data',
                                 data_folders,
                                 nbr_simulations,
                                 max_time=None,
                                 timestep=1. / 12.,
                                 labels=labels)
    """
    for exp_nbr in N_strains_a22:
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        initialized_start(exp_dir, nbr_traj, max_time, labels=labels[3])

    richness_N_strains(N_strains_a22, 'A22', nbr_traj)
    richness_N_strains(N_strains_wt, 'WT', nbr_traj)
    richness_N_strains(N_strains_vs, 'VS', nbr_traj)
    """
