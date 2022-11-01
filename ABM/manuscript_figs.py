from python_src.first_passage import MoranFPT, OneBoundaryFPT, TwoBoundaryFPT,\
    OneBoundaryIntFPT  # , TwoBoundaryIntFPT
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.lines import Line2D

plt.style.use('python_src/custom.mplstyle')


def spatial_vs_moran(r1, r2, N, times, cmap_name, filename, model):
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 1.0, num=len(N))]
    linestyle = [':', '-', '--']

    # colour legend
    col_label = N  # ["N={}".format(nbr) for nbr in N]
    col_lines = [Line2D([0], [0], color=col) for col in colors]
    cl = plt.legend(col_lines, col_label, loc='lower right')

    # Fig1A : P_absorbing
    fig1, ax1 = plt.subplots(figsize=(3.4, 2.5))
    ax1.set(title=r"",
            xlabel=r"Initial species 1 fraction, $x$",
            ylabel=r"Probability fixation, $P_N(x)$")
    ax1.set_xlim([0.0, 1.0])

    legend = ['spatial (ME)', 'spatial (FP)', 'moran (FP)']
    custom_lines = [
        Line2D([0], [0], color='dimgray', linestyle='None', marker='x'),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[1]),
        Line2D([0], [0], color='k', linestyle=linestyle[0])
        ]

    # Fig1B : MFPT function of x
    fig2, ax2 = plt.subplots(figsize=(3.4, 2.5))
    ax2.set_xlim([0.0, 1.0])
    ax2.set_yscale('log')
    ax2.set(title=r"",
            xlabel=r"Initial species 1 fraction, $x$",
            ylabel=r"MFPT, $\tau(x)$")

    # Fig1C : MFPT as function of N
    fig3, ax3 = plt.subplots(figsize=(3.4, 2.5))
    ax3.set_yscale('log')
    ax3.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT, $\tau(x_{max})$")

    moran = []
    space = []
    ME_mfpt_space_N = []
    for i, K in enumerate(N):
        # define models
        moran.append(MoranFPT(r1, r2, K, times))
        space.append(model(r1, r2, K, times))

        # plots
        x = []
        ME_mfpt_space = []
        ME_prob_space = []

        for j in np.linspace(1, K - 1, 9):
            x.append(j / K)
            j = int(j)

            # ME approach
            ME_prob, ME_mfpt = space[i].probability_mfpt(j)
            ME_prob_space.append(ME_prob[1])
            ME_mfpt_space.append(np.dot(ME_prob, ME_mfpt))

        # Fig 1 A
        X, FP_prob_space = space[i].FP_prob()
        ax1.plot(X, FP_prob_space, linestyle=linestyle[1], color=colors[i])
        ax1.scatter(x, ME_prob_space, color=colors[i], marker='x', zorder=2.5)

        # Fig 1 B
        X, FP_mfpt_moran = moran[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_moran, linestyle=linestyle[0], color=colors[i])
        X, FP_mfpt_space = space[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_space, linestyle=linestyle[1], color=colors[i])
        ax2.scatter(x, ME_mfpt_space, color=colors[i], marker='x', zorder=2.5)

        # Fig 1 C
        prob, mfpt = space[i].probability_mfpt(int(K / 2))
        ME_mfpt_space_N.append(np.dot(prob, mfpt))

    # Moran probability
    X, FP_prob_moran = moran[-1].FP_prob()
    ax1.plot(X, FP_prob_moran, linestyle=linestyle[0], color='k', zorder=1.)

    # deterministic approximation times
    x_det = np.linspace(0.52, 1.0, 100)
    t_det = - np.log(2 * x_det - 1.0)
    ax2.plot(x_det, t_det, linestyle=linestyle[2], color='r')
    x_det = np.linspace(0.0, 0.48, 100)
    t_det = - np.log(-2 * x_det + 1.0)
    det = ax2.plot(x_det, t_det, linestyle=linestyle[2], color='r')

    # Fig 1C: MFPT as a function of N at x_max
    N_func = np.linspace(10, N[-1], 100)
    _, FP_mfpt_moran_N = moran[0].FP_mfpt(x=0.5, N=N_func)
    _, FP_mfpt_space_N = space[0].FP_mfpt_x(x_find=0.5, N=N_func)

    ax3.plot(N_func, FP_mfpt_moran_N, color='k', linestyle=linestyle[2])
    ax3.plot(N_func, FP_mfpt_space_N, color='dimgray')

    ax3.scatter(N, ME_mfpt_space_N, label=r'ME Homog.', marker='x',
                color=colors[0: len(N)], zorder=2.5)

    # figure limits and legends
    ax3.set_ylim([1.0, 1000])
    ax2.set_ylim(0.1, 100)
    cl = ax1.legend(col_lines, col_label, loc='lower right', title=r'$N$')
    ax1.legend(custom_lines, legend, title=r'Model')
    ax1.add_artist(cl)
    ax2.legend(det, [r'$\tau_{det}$'])
    # ax3.legend(custom_lines, legend)

    # save figures
    fig1.savefig(filename + '_prob.pdf')
    fig1.savefig(filename + '_prob.png')
    fig2.savefig(filename + '_mfpt.pdf')
    fig2.savefig(filename + '_mfpt.png')
    fig3.savefig(filename + '_N.pdf')
    fig3.savefig(filename + '_N.png')

    return 0


if __name__ == '__main__':
    # mkdir
    dir = 'figures_theory'
    Path(dir).mkdir(parents=True, exist_ok=True)

    # theory parameters
    N = [10, 100, 1000]
    times = np.linspace(0.0, 100.0, 10001)
    cmap_name = 'viridis'
    r1 = 1.0
    r2 = 1.0

    # Figure 1 : spatial vs moran
    fname = dir + os.sep + 'SvM'
    spatial_vs_moran(r1, r2, N, times, cmap_name, fname, OneBoundaryIntFPT)
