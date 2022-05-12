import os
import re
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import python_src

BACT_COL = {'0': 'r', '1': 'g', '2': 'b', '3': 'c', '4': 'm', '5': 'y', '6': 'k'}


def collect_data_array(data_folder, nbr_simulations, nbr_species, nbr_lines=None):
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


labels = ['E. Coli Alt 1', 'E. Coli Alt 2']
labels = ['C. Subtillus', 'E. Coli']
labels = ['E.Coli strain A', 'E. Coli strain B']
timestep = 1. / 12.


def length_trajectory_plots(data_folder, nbr_simulations, nbr_species, max_time=None):
    data = collect_data_array(data_folder, nbr_simulations, nbr_species, max_time)
    max_t = np.shape(data)[2] * 1. / 12.
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
        ax.plot(t, av_traj[i], lw=2, color=BACT_COL[str(i)], label=labels[i])
        ax.fill_between(t, av_traj[i] + std_traj[i], av_traj[i] - std_traj[i], facecolor=BACT_COL[str(i)], alpha=0.5)
    ax.set_title(r"Total length of bacteria")
    plt.xlim([0.0, max_t])
    plt.ylabel(r'sum of length of bacteria, $\mu m$')
    plt.xlabel(r'time, $h$')
    plt.legend()
    plt.show()


def distribution_extinction(data_folder, nbr_simulations, nbr_species, max_time):
    data = collect_data_array(data_folder, nbr_simulations, nbr_species, max_time)
    extinctions = [[] for _ in range(nbr_species)]
    nbr_extinctions = np.zeros((2))
    for i in np.arange(0, nbr_simulations):
        for j in np.arange(0, nbr_species):
            zeros = np.where(data[i, j, :] == 0.0)[0]
            if zeros != []:
                extinctions[j].append(zeros[0] * timestep)
                nbr_extinctions[j] += 1
    fig, ax = plt.subplots(1)
    num_bins = 20
    for j in np.arange(0, nbr_species):
        ax.hist(extinctions[j], num_bins, facecolor=BACT_COL[str(j)], alpha=0.5, label=labels[j])
    ax.set_title(r"distribution extinction times")
    max_t = np.shape(data)[2] * 1. / 12.
    plt.xlim([0.0, max_t])
    plt.ylabel(r'count')
    plt.xlabel(r'fixation time, $h$')
    plt.legend()
    plt.show()
