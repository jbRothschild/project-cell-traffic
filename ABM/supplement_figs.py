from python_src.first_passage import MoranFPT, Spatial, SpatInv
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.lines import Line2D
import scipy as sp

plt.style.use('python_src/custom.mplstyle')


def moranME_FP(times, filename, model):
    s = [1.0, 1.01, 1.1]
    N = 100

    marker = 'X'
    cmap = matplotlib.cm.get_cmap('plasma')
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 1.0, num=len(s))]
    linestyle = [':', '-', '--']

    # colour legend
    col_label = s
    col_lines = [Line2D([0], [0], color=col) for col in colors]
    cl = plt.legend(col_lines, col_label, loc='lower right')

    # Fig1A : P_absorbing
    fig1, ax1 = plt.subplots(figsize=(3.4, 2.5))
    ax1.set(title=r"",
            xlabel=r"Initial species 1 fraction, $x$",
            ylabel=r"Probability fixation, $P_N(x)$")

    legend = ['moran (ME)', 'moran (FP)']
    custom_lines = [
        Line2D([0], [0], markerfacecolor='dimgray', linestyle='None',
               marker=marker, color='k', linewidth=1),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[0])
        ]

    # Fig1B : MFPT function of x
    fig2, ax2 = plt.subplots(figsize=(3.4, 2.5))
    ax2.set_xlim([0.0, 1.0])
    ax2.set(title=r"",
            xlabel=r"Initial species 1 fraction, $x$",
            ylabel=r"MFPT, $\tau(x)$")

    # Fig1C : MFPT as function of N
    fig3, ax3 = plt.subplots(figsize=(3.4, 2.5))
    ax3.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT, $\tau(x_{max})$")

    moran = []
    for i, fit in enumerate(s):
        # define models
        moran.append(model(fit, 1.0, N, times))

        # plots
        x = []
        ME_mfpt_moran = []
        ME_prob_moran = []

        # Fig 1 A
        for j in np.linspace(1, N - 1, 9):
            x.append(j / N)
            j = int(j)

            # ME approach
            ME_prob, ME_mfpt = moran[i].probability_mfpt(j)
            ME_prob_moran.append(ME_prob[1])
            ME_mfpt_moran.append(np.dot(ME_prob, ME_mfpt))

        X, FP_prob_moran = moran[i].FP_prob()
        ax1.plot(X, FP_prob_moran, linestyle=linestyle[0], color=colors[i],
                 zorder=i)
        ax1.scatter(x, ME_prob_moran, color=colors[i], marker=marker,
                    s=12, edgecolors='k', linewidth=1, zorder=i+.5)

        # Fig 1 B
        X, FP_mfpt_moran = moran[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_moran, linestyle=linestyle[0], color=colors[i],
                 zorder=i)
        ax2.scatter(x, ME_mfpt_moran, color=colors[i], marker=marker,
                    s=12, edgecolors='k', linewidth=1, zorder=i+.5)

        # Fig 1 C
        # prob, mfpt = space[i].probability_mfpt(space[i].nmax)
        # ME_mfpt_space_N.append(np.dot(prob, mfpt))

    # Fig 1C: MFPT as a function of N at x_max
    """
    N_func = np.linspace(10, N[-1], 100)
    _, FP_mfpt_moran_N = moran[0].FP_mfpt(x=moran[0].xmax, N=N_func)
    _, FP_mfpt_space_N = space[0].FP_mfpt_x(x_find=space[0].xmax, N=N_func)

    ax3.plot(N_func, FP_mfpt_moran_N, color='k', linestyle=linestyle[2])
    ax3.plot(N_func, FP_mfpt_space_N, color='dimgray')

    ax3.scatter(N, ME_mfpt_space_N, label=r'ME Homog.', marker='o',
                color=colors[0: len(N)], s=12, edgecolors='k', linewidth=1,
                zorder=i+.5)
    """

    # figure limits and legends
    ax1.set_xlim([0.00, 1.0])
    ax1.set_ylim([0.00, 1.0])
    ax3.set_yscale('log')
    ax2.set_ylim(0.1, 100)
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_ylim([1.0, 1000])
    cl = ax1.legend(col_lines, col_label, loc='lower right', title=r'$s$')
    ax1.legend(custom_lines, legend, title=r'Model', loc='upper right',
               framealpha=0.9)
    ax3.legend(custom_lines, legend, title=r'Model')
    ax1.add_artist(cl)
    # ax2.legend(det, [r'$\tau_{det}$'])
    # ax3.legend(custom_lines, legend)

    # save figures
    fig1.savefig(filename + '_prob.pdf')
    fig1.savefig(filename + '_prob.png')
    fig2.savefig(filename + '_mfpt.pdf')
    fig2.savefig(filename + '_mfpt.png')
    # fig3.savefig(filename + '_N.pdf')
    # fig3.savefig(filename + '_N.png')
    return s


def invasion_check(s, N, cmap_name, filename, model):
    marker = 'o'
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 1.0, num=len(s))]
    linestyle = [':', '-', '--']

    # colour legend
    col_label = s
    col_lines = [Line2D([0], [0], color=col) for col in colors]

    # Fig2A : P_absorbing
    fig1, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(3.4, 5.0))
    plt.subplots_adjust(wspace=0, hspace=0.0)
    # fig1, ax1 = plt.subplots(figsize=(3.4, 2.5))
    ax1.set(title=r"",
            # xlabel=r"Initial species 1 fraction, $x$",
            ylabel=r"Probability fixation, $P_N(x)$")

    # Fig2B : MFPT function of x
    # fig2, ax2 = plt.subplots(figsize=(3.4, 2.5))
    ax2.set_xlim([0.0, 1.0])
    ax2.set(title=r"",
            xlabel=r"Initial species 1 fraction, $x$",
            ylabel=r"MFPT, $\tau(x)$")

    space = []
    # ME_mfpt_space_N = []
    for i, fit in enumerate(s):
        # define models
        space.append(model(fit, 1.0, N, times))

        # plots
        x = []
        ME_mfpt_space = []
        ME_prob_space = []

        # Fig 1 A
        # print(np.linspace(1, N - 1, 9))
        for j in np.linspace(1, N - 1, 9):
            x.append(j / N)
            j = int(j)

            # ME approach
            ME_prob, ME_mfpt = space[i].probability_mfpt(0, space[i].K - j)
            ME_prob_space.append(ME_prob[0])
            ME_mfpt_space.append(np.dot(ME_prob, ME_mfpt))

        ax1.scatter(x, ME_prob_space, color=colors[i], marker=marker,
                    s=12, edgecolors='k', linewidth=1, zorder=i+.5)
        ax2.scatter(x, ME_mfpt_space, color=colors[i], marker=marker,
                    s=12, edgecolors='k', linewidth=1, zorder=i + .5)

    # deterministic approximation times
    x_det = np.linspace(0.2, 1.0, 100)
    S = 100
    t_det = - (np.log(np.abs((S - 1) * x_det**2 + 2 * x_det - 1)) - np.log(S))
    det = ax2.plot(x_det, t_det, linestyle=linestyle[2], color='r', zorder=i+1)
    x_det = np.linspace(0.0, 0.18, 100)
    t_det2 = - (sp.log((S - 1) * x_det**2 + 2 * x_det - 1) - sp.log(-1))
    ax2.plot(x_det, t_det2, linestyle=linestyle[2], color='r', zorder=i+1)

    # figure limits and legends
    # legend
    legend = ['spatial (ME)', 'spatial (FP)', 'moran (FP)']
    custom_lines = [
        Line2D([0], [0], markerfacecolor='dimgray', linestyle='None',
               marker=marker, color='k', linewidth=1),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[1]),
        Line2D([0], [0], color='k', linestyle=linestyle[0])
        ]

    # limits
    ax1.set_xlim([0.00, 1.0])
    ax1.set_ylim([0.00, 1.0])
    cl = ax1.legend(col_lines, col_label, loc='lower right', title=r'$s$')
    ax1.legend(custom_lines, legend, title=r'Model', loc='upper right',
               framealpha=0.8)
    ax1.add_artist(cl)

    ax2.set_yscale('log')
    ax2.set_ylim(0.1, 100)
    ax2.legend(det, [r'$\tau_{\mathrm{det}}$'])

    # save figures
    fig1.savefig(filename + '_prob.pdf')
    fig1.savefig(filename + '_prob.png')
    return 0


if __name__ == '__main__':
    # mkdir
    dir = 'figures_suppl'
    Path(dir).mkdir(parents=True, exist_ok=True)

    # for certain distribution functions
    times = np.linspace(0.0, 100.0, 10001)

    # theory parameters
    model = Spatial
    model2 = SpatInv

    # colormaps
    cmap_name1 = 'plasma'

    # Figure S1 : ME vs FP moran model
    fname1 = dir + os.sep + 'moran'
    moranME_FP(times, fname1, MoranFPT)

    # Figure SX : invasion model check
    s = [1., 10., 100.]
    K = 100
    fnameX = dir + os.sep + 'inv_check'
    invasion_check(s, K, cmap_name1, fnameX, SpatInv)
