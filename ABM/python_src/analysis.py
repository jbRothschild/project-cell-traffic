import os
# bimport re
import numpy as np
# from pathlib import Path
import matplotlib.pyplot as plt
# import python_src
# from python_src.first_passage import MoranFPT, MoranGrowFPT

BACT_COL = {'0': 'r', '1': 'g', '2': 'b', '3': 'c', '4': 'm', '5': 'y',
            '6': 'k'}


def collect_data_array(data_folder, nbr_simulations, nbr_lines=None,
                       timestep=1. / 12., labels=None):
    with open(data_folder + os.sep + 'sim0_data.txt') as f:
        first_line = f.readline()
    nbr_species = first_line.count(",") + 1
    if nbr_lines is None:
        nbr_lines = sum(1 for line in open(data_folder
                                           + os.sep + 'sim0_data.txt'))
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
    ax2.set_ylabel(r'average fraction of species survival, $S^*/S$',
                   color=color)  # we already handled the x-label with ax1
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
        ax.plot(t, av_traj[i], lw=2,
                color=BACT_COL[str(i)])  # , label=labels[i])
        ax.fill_between(t, av_traj[i] + std_traj[i], av_traj[i] - std_traj[i],
                        facecolor=BACT_COL[str(i)], alpha=0.5)
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
    ax.hist(data[:, 0, -1] / total_len[:, -1], bins=39, color='green',
            edgecolor='black', density=True)
    # ax.hist(data[:, 1, -1], bins=30, color='red', edgecolor='black')
    plt.ylabel(r'count')
    plt.xlabel(r'length bacteria')
    # plt.legend()
    # plt.show()


def distribution_extinction(data_folder, nbr_simulations,
                            max_time=None, timestep=1. / 12., labels=None):

    data = collect_data_array(data_folder, nbr_simulations, max_time)
    # max_t = np.shape(data)[2] * timestep
    max_t = 12.
    nbr_species = np.shape(data)[1]
    extinctions = [[] for _ in range(nbr_species)]
    nbr_extinctions = np.zeros((nbr_species))
    for i in np.arange(0, nbr_simulations):
        for j in np.arange(0, nbr_species):
            zeros = np.where(data[i, j, :] == 0.0)[0]
            if zeros != []:
                extinctions[j].append(zeros[0] * timestep)
                nbr_extinctions[j] += 1

    # figure for distribution of each species extinction times
    fig, ax = plt.subplots(1)
    num_bins = 20
    for j in np.arange(0, nbr_species):
        ax.hist(extinctions[j], num_bins, facecolor=BACT_COL[str(j)],
                alpha=0.5, density=True, label=labels[j])
        ax.axvline(np.mean(extinctions[j]), color=BACT_COL[str(j)],
                   linestyle='dashed', linewidth=1)
    ax.set_title(r"distribution fixation times")

    # Moran fpt
    """
    times = np.arange(0, 3 * max_t + timestep, timestep)
    moran = MoranFPT(60 * 0.0173, 60, times)
    prob, mfpt = moran.probability_mfpt(30)
    fpt_dist, tot_fpt = moran.fpt_distribution(30)
    plt.plot(times, tot_fpt, 'k', label='moran')
    moran = MoranGrowFPT(60 * 0.0173, 60, times)
    prob, mfpt = moran.probability_mfpt(30, 30)
    fpt_dist, tot_fpt = moran.fpt_distribution(30, 30)
    plt.plot(times, tot_fpt, 'b', label='spatial model')
    # plt.ylim([0.000001, 1])
    """
    plt.xlim([0.0, max_t])
    plt.yscale('log')
    plt.ylabel(r'probability')
    plt.xlabel(r'fixation time, $h$')
    plt.legend()
    plt.savefig(data_folder + os.sep + "fixations.pdf")
    plt.savefig(data_folder + os.sep + "fixations.png")
    # plt.show()

    # figure for fixation vs coexistence
    fig, ax = plt.subplots(1)
    cat_ext = [y for x in extinctions for y in x]
    fig, ax = plt.subplots(1)
    ax.hist(cat_ext, num_bins, facecolor='gray', alpha=0.5, density=True)
    ax.axvline(np.mean(cat_ext), color='gray', linestyle='dashed', linewidth=1)
    ax.set_title(r"distribution fixation times")
    max_t = np.shape(data)[2] * timestep

    plt.xlim([0.0, max_t])
    plt.yscale('log')
    plt.ylabel(r'probability')
    plt.xlabel(r'fixation time, $h$')
    plt.legend()
    plt.savefig(data_folder + os.sep + "fixation_cat.pdf")
    plt.savefig(data_folder + os.sep + "fixation_cat.png")
    # plt.show()

    # figure for fixation vs coexistence

    return


def bar_chart_fixations(sim_dir, data_folders, nbr_simulations,
                        max_time=None, timestep=1. / 12., labels=None):
    # bar chart specifications
    ind = np.flip(np.arange(len(data_folders)))
    width = 0.25

    # find probabilities for all fixations
    fix_prob = []
    fix_coex_prob = []
    exp_fix = [[0.23, 0.08, 0.69], [0.37, 0.23, 0.40], [0.02, 0.7, 0.28]]
    exp_fix_coex = [[0.31, 0.69], [0.6, 0.4], [0.72, 0.28]]
    for i, dfolder in enumerate(data_folders):
        fix = fix_vs_coexi(dfolder, nbr_simulations[i], max_time)
        fix_prob.append(list(fix))
        fix_coex_prob.append([np.sum(fix[:-1]), fix[-1]])

        # placeholder experiments

    # fixation vs coexistence
    fig, ax = plt.subplots(1, 1)

    # sim
    ax.barh(y=ind - width / 2, width=np.array(fix_coex_prob)[:, 0],
            height=width, label='fixation', color='lightgrey', edgecolor=None)
    ax.bar_label(ax.containers[0], label_type='center')
    ax.barh(y=ind - width / 2, width=np.array(fix_coex_prob)[:, 1],
            height=width, left=np.array(fix_coex_prob)[:, 0],
            label='coexistence', color='dimgrey', edgecolor=None)
    ax.bar_label(ax.containers[1], label_type='center')
    # exp
    ax.barh(y=ind + width / 2, width=np.array(exp_fix_coex)[:, 0],
            height=width, color='lightgrey', edgecolor=None)
    ax.bar_label(ax.containers[2], label_type='center')
    ax.barh(y=ind + width / 2, width=np.array(exp_fix_coex)[:, 1],
            height=width, left=np.array(exp_fix_coex)[:, 0], color='dimgrey',
            edgecolor=None)
    ax.bar_label(ax.containers[3], label_type='center')
    ax.barh(y=ind + width / 2, width=[1.0] * len(data_folders), height=width,
            fill=False, edgecolor='k', label='experiment')
    ax.set_xlim([0.0, 1.0])
    plt.setp(ax, yticks=ind, yticklabels=labels)
    plt.xlabel(r'Probability')
    plt.legend()
    plt.savefig(sim_dir + os.sep + "coex_prob.pdf")
    plt.savefig(sim_dir + os.sep + "coex_prob.png")

    fig, ax = plt.subplots(1, 1)
    # sim
    bot = 0
    for j in np.arange(np.shape(np.array(fix_prob))[1] - 1):
        ax.barh(y=ind - width / 2, width=np.array(fix_prob)[:, j],
                height=width, left=bot, color=BACT_COL[str(j)], edgecolor=None)
        bot += np.array(fix_prob)[:, j]
        ax.bar_label(ax.containers[j], label_type='center')
    ax.barh(y=ind - width / 2, width=np.array(fix_prob)[:, -1], height=width,
            color='dimgrey', left=bot, edgecolor=None)
    ax.bar_label(ax.containers[j + 1], label_type='center')
    # exp
    nbr = j + 1
    bot = 0
    for j in np.arange(np.shape(np.array(exp_fix))[1] - 1):
        ax.barh(y=ind + width / 2, width=np.array(exp_fix)[:, j], height=width,
                left=bot, color=BACT_COL[str(j)], edgecolor=None)
        ax.bar_label(ax.containers[nbr + j], label_type='center')
        bot += np.array(exp_fix)[:, j]
    ax.barh(y=ind + width / 2, width=np.array(exp_fix)[:, -1], height=width,
            left=bot, label='coexistence', color='dimgrey', edgecolor=None)
    ax.bar_label(ax.containers[nbr + j + 1], label_type='center')
    ax.barh(y=ind + width / 2, width=[1.0] * len(data_folders), height=width,
            fill=False, edgecolor='k', label='experiment')
    ax.bar_label(ax.containers[nbr + j + 2], label_type='center')

    ax.set_xlim([0.0, 1.0])
    plt.setp(ax, yticks=ind, yticklabels=labels)
    plt.xlabel(r'Probability')
    plt.legend()
    plt.savefig(sim_dir + os.sep + "fix_prob.pdf")
    plt.savefig(sim_dir + os.sep + "fix_prob.png")

    return 0


def fix_vs_coexi(data_folder, nbr_simulations, max_time):
    data = collect_data_array(data_folder, nbr_simulations, max_time)
    nbr_species = np.shape(data)[1]
    nbr_extinctions = np.zeros((nbr_species))
    for i in np.arange(0, nbr_simulations):
        for j in np.arange(0, nbr_species):
            zeros = np.where(data[i, j, :] == 0.0)[0]
            if zeros != []:
                nbr_extinctions[j] += 1.
    nbr_extinctions /= nbr_simulations
    nbr_coexistence = 1. - np.sum(nbr_extinctions)

    return np.append(nbr_extinctions, nbr_coexistence)
