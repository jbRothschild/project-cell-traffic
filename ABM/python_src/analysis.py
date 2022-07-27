import os
# bimport re
import numpy as np
# from pathlib import Path
import matplotlib.pyplot as plt
# import python_src
from python_src.first_passage import MoranFPT, MoranGrowFPT

BACT_COL = {'0': 'r', '1': 'g', '2': 'b', '3': 'c', '4': 'm', '5': 'y', '6': 'k'}


def collect_data_array(data_folder, nbr_simulations, nbr_lines=None,
                       timestep=1. / 12., labels=None):
    with open(data_folder + os.sep + 'sim0_data.txt') as f:
        first_line = f.readline()
    nbr_species = first_line.count(",") + 1
    if nbr_lines is None:
        nbr_lines = sum(1 for line in open(data_folder + os.sep + 'sim0_data.txt'))
    data = np.zeros((nbr_simulations, nbr_species, nbr_lines))
    for i in np.arange(0, nbr_simulations):
        time = 0
        with open(data_folder + os.sep + 'sim' + str(i) + "_data.txt") as f:
            for line in f:
                if time < nbr_lines:
                    strip_line = line.strip("[]\n")
                    for j, spec in enumerate(strip_line.split(",")):
                        data[i, j, time] = float(spec)
                time += 1

    return data


def richness_traj_plots(data_folder, nbr_simulations, max_time=None,
                        timestep=1. / 12.):
    data = collect_data_array(data_folder, nbr_simulations, max_time)
    nbr_species = np.shape(data)[1]
    max_t = np.shape(data)[2] * timestep
    t = np.arange(0., max_t, timestep)
    # dist_extinction =
    # labels = ['C. Bacillus', 'E. Coli Alt']

    richness = np.count_nonzero(data, axis=1)
    av_rich = np.mean(richness, axis=0)
    std_rich = np.std(richness, axis=0)
    fig, ax1 = plt.subplots(1)
    color = 'k'
    ax1.plot(t, av_rich, lw=2, color=color)
    ax1.set_title(r"richness")
    ax1.set_xlabel(r'time, $h$')
    ax1.fill_between(t, av_rich + std_rich,
                     av_rich - std_rich, facecolor=color, alpha=0.5)
    ax1.set_ylabel(r'average number of species, $S^*$', color=color)
    # ax1.set_ylim(0.0, nbr_species)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(r'average fraction of species survival, $S^*/S$', color=color)  # we already handled the x-label with ax1
    ax2.plot(t, av_rich / nbr_species, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylim(0.0, 1.0)

    plt.xlim([0.0, max_t])
    fig.tight_layout()
    plt.savefig(data_folder + os.sep + "richness.pdf")
    plt.savefig(data_folder + os.sep + "richness.png")
    # plt.show()
    return nbr_species, av_rich, std_rich


def length_trajectory_plots(data_folder, nbr_simulations, max_time=None,
                            timestep=1. / 12.,
                            labels=None):
    """
    For simulations of same initial number of species, averages the length of
    the strains.
    """
    data = collect_data_array(data_folder, nbr_simulations, max_time)
    nbr_species = np.shape(data)[1]
    max_t = np.shape(data)[2] * timestep
    t = np.arange(0., max_t, timestep)
    # dist_extinction =
    # labels = ['C. Bacillus', 'E. Coli Alt']

    total_len = np.sum(data, axis=1)
    av_traj = np.mean(data, axis=0)
    std_traj = np.std(data, axis=0)
    tot_av_traj = np.mean(total_len, axis=0)
    tot_std_traj = np.std(total_len, axis=0)
    fig, ax = plt.subplots(1)
    ax.plot(t, tot_av_traj, lw=2, color='k', label='total')
    ax.fill_between(t, tot_av_traj + tot_std_traj,
                    tot_av_traj - tot_std_traj, facecolor='k', alpha=0.5)
    for i in np.arange(0, nbr_species):
        ax.plot(t, av_traj[i], lw=2, color=BACT_COL[str(i)])  # , label=labels[i])
        ax.fill_between(t, av_traj[i] + std_traj[i], av_traj[i] - std_traj[i], facecolor=BACT_COL[str(i)], alpha=0.5)
    ax.set_title(r"Total length of bacteria")
    plt.xlim([0.0, max_t])
    plt.ylabel(r'sum of length of bacteria, $\mu m$')
    plt.xlabel(r'time, $h$')
    # plt.legend()
    plt.savefig(data_folder + os.sep + "length_bact.pdf")
    plt.savefig(data_folder + os.sep + "length_bact.png")
    # plt.show()
    #
    #
    fig, ax = plt.subplots(1)
    ax.hist(data[:, 0, -1] / total_len[:, -1], bins=39, color='green', edgecolor='black', density=True)
    # ax.hist(data[:, 1, -1], bins=30, color='red', edgecolor='black')
    plt.ylabel(r'count')
    plt.xlabel(r'length bacteria')
    # plt.legend()
    # plt.show()


def distribution_extinction(data_folder, nbr_simulations,
                            max_time=None, timestep=1. / 12., labels=None):
    data = collect_data_array(data_folder, nbr_simulations, max_time)
    nbr_species = np.shape(data)[1]
    extinctions = [[] for _ in range(nbr_species)]
    nbr_extinctions = np.zeros((nbr_species))
    for i in np.arange(0, nbr_simulations):
        for j in np.arange(0, nbr_species):
            zeros = np.where(data[i, j, :] == 0.0)[0]
            if zeros != []:
                extinctions[j].append(zeros[0] * timestep)
                nbr_extinctions[j] += 1

    # print(np.sum(nbr_extinctions)/nbr_simulations) fraction sims with fixation
    fig, ax = plt.subplots(1)
    num_bins = 20
    for j in np.arange(0, nbr_species):
        ax.hist(extinctions[j], num_bins, facecolor=BACT_COL[str(j)], alpha=0.5, density=True)  # , label=labels[j])
    ax.set_title(r"distribution extinction times")
    max_t = np.shape(data)[2] * timestep

    # Moran fpt
    times = np.arange(0, 3 * max_t + timestep, timestep)
    moran = MoranFPT(60 * 0.0173, 60, times)
    prob, mfpt = moran.probability_mfpt(30)
    fpt_dist, tot_fpt = moran.fpt_distribution(30)
    plt.plot(times, tot_fpt, 'k', label='moran')
    moran = MoranGrowFPT(60 * 0.0173, 60, times)
    prob, mfpt = moran.probability_mfpt(30, 30)
    fpt_dist, tot_fpt = moran.fpt_distribution(30, 30)
    plt.plot(times, tot_fpt, 'b', label='spatial model')
    print('done')
    # plt.xlim([0.0, max_t])
    # plt.ylim([0.000001, 1])
    plt.yscale('log')
    plt.ylabel(r'count')
    plt.xlabel(r'fixation time, $h$')
    plt.legend()
    plt.savefig(data_folder + os.sep + "extinction.pdf")
    plt.savefig(data_folder + os.sep + "extinctions.png")
    # plt.show()
