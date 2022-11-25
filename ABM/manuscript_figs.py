from python_src.first_passage import MoranFPT, Spatial, SpatInv
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy as sp
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def spatial_vs_moran(N, times, cmap_name, filename, model):
    r1 = 1.0
    r2 = 1.0
    # marker = 'o'
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(N))]
    linestyle = [':', '-', '--']

    # colour legend
    col_label = N  # ["N={}".format(nbr) for nbr in N]
    col_lines = [Line2D([0], [0], color=col) for col in colors]
    cl = plt.legend(col_lines, col_label, loc='lower right')

    # Fig1A : P_absorbing
    fig1, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(3.4, 5.0))
    plt.subplots_adjust(wspace=0, hspace=0)
    # fig1, ax1 = plt.subplots(figsize=(3.4, 2.5))
    ax1.set(title=r"",
            # xlabel=r"Relative abundance, $x$",
            ylabel=r"Probability fixation, $P_N(x)$")
    """
    legend = ['spatial (ME)', 'spatial (FP)', 'moran (FP)']
    custom_lines = [
        Line2D([0], [0], markerfacecolor='dimgray', linestyle='None',
               marker=marker, color='k', linewidth=1),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[1]),
        Line2D([0], [0], color='k', linestyle=linestyle[0])
        ]
    """
    # Fig1B : MFPT function of x
    # fig2, ax2 = plt.subplots(figsize=(3.4, 2.5))
    ax2.set(title=r"",
            xlabel=r"Relative abundance, $x$",
            ylabel=r"MFPT, $\tau(x)$")

    axins2 = inset_axes(ax2, width="35%", height="35%", borderpad=1)

    # Fig1C : MFPT as function of N
    fig3, ax3 = plt.subplots(figsize=(3.4, 2.5))
    ax3.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT, $\tau(x_{\mathrm{max}})$")

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

        # Fig 1 A
        for j in np.linspace(K * (i * 1.75 + 1)
                             / 10, K * (1 - (i * 1.75 + 1) / 10), 9):
            x.append(j / K)
            j = int(j)

            # ME approach
            ME_prob, _ = space[i].probability_mfpt(j)
            ME_prob_space.append(ME_prob[1])

        X, FP_prob_space = space[i].FP_prob()
        ax1.plot(X, FP_prob_space, linestyle=linestyle[1], color=colors[i],
                 zorder=i)
        # ax1.scatter(x, ME_prob_space, color=colors[i], marker=marker,
        #            s=12, edgecolors='k', linewidth=1, zorder=i+.5)

        # Fig 1 B
        x = []
        ME_mfpt_space = []
        for j in np.linspace(1, K - 1, 9):
            x.append(j / K)
            j = int(j)

            # ME approach
            ME_prob, ME_mfpt = space[i].probability_mfpt(j)
            ME_mfpt_space.append(np.dot(ME_prob, ME_mfpt))

        X, FP_mfpt_moran = moran[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_moran, linestyle=linestyle[0], color=colors[i])
        X, FP_mfpt_space = space[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_space, linestyle=linestyle[1], color=colors[i],
                 zorder=i)
        axins2.plot(X, FP_mfpt_space, linestyle=linestyle[1], color=colors[i],
                    zorder=i)
        # ax2.scatter(x, ME_mfpt_space, color=colors[i], marker=marker,
        #            s=12, edgecolors='k', linewidth=1, zorder=i+.5)

        # Fig 1 C
        prob, mfpt = space[i].probability_mfpt(space[i].nmax)
        ME_mfpt_space_N.append(np.dot(prob, mfpt))

    # Moran probability
    X, FP_prob_moran = moran[-1].FP_prob()
    ax1.plot(X, FP_prob_moran, linestyle=linestyle[0], color='k', zorder=1.)

    # deterministic approximation times
    x_det = np.linspace(0.50, 1.0, 1000)
    t_det = - np.log(2 * x_det - 1.0)
    ax2.plot(x_det, t_det, linestyle=linestyle[2], color='k', zorder=i+1)
    x_det = np.linspace(0.0, 0.50, 1000)
    t_det = - np.log(-2 * x_det + 1.0)
    # det=ax2.plot(x_det, t_det, linestyle=linestyle[2], color='k', zorder=i+1)
    ax2.plot(x_det, t_det, linestyle=linestyle[2], color='k', zorder=i+1)

    # Fig 1C: MFPT as a function of N at x_max
    N_func = np.linspace(10, N[-1], 100)
    _, FP_mfpt_moran_N = moran[0].FP_mfpt(x=moran[0].xmax, N=N_func)
    _, FP_mfpt_space_N = space[0].FP_mfpt_x(x_find=space[0].xmax, N=N_func)

    ax3.plot(N_func, FP_mfpt_moran_N, color='k', linestyle=linestyle[0])
    ax3.plot(N_func, FP_mfpt_space_N, color='k')

    # ax3.scatter(N, ME_mfpt_space_N, label=r'ME Homog.', marker=marker,
    #            color=colors[0: len(N)], s=12, edgecolors='k', linewidth=1,
    #            zorder=i+.5)

    # figure limits and legends
    cl = ax1.legend(col_lines, col_label, loc='lower right', title=r'$N$')
    # ax1.legend(custom_lines, legend, title=r'Model')
    ax1.add_artist(cl)
    ax1.set_xlim([0.0, 1.0])

    # ax2.legend(det, [r'$\tau_{\mathrm{det}}$'])
    ax2.set_ylim(0.1, 1000)
    ax2.set_xlim([0.0, 1.0])
    ax2.set_yscale('log')

    axins2.set_ylim(1., 10)
    axins2.set_xlim([0.4, 0.6])
    axins2.set_yscale('log')
    axins2.xaxis.set_minor_formatter(NullFormatter())
    axins2.yaxis.set_minor_formatter(NullFormatter())

    # ax3.legend(custom_lines, legend, title=r'Model')
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_ylim([1.0, 1000])
    ax3.set_xlim([np.min(N_func), np.max(N_func)])
    # ax3.legend(custom_lines, legend)

    # save figures
    fig1.savefig(filename + '_prob.pdf')
    fig1.savefig(filename + '_prob.png')
    # fig2.savefig(filename + '_mfpt.pdf')
    # fig2.savefig(filename + '_mfpt.png')
    fig3.savefig(filename + '_N.pdf')
    fig3.savefig(filename + '_N.png')

    return 0


def fitness_spatial(s, N, times, cmap_name, filename, model):

    marker = 'o'
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(s))]
    linestyle = [':', '-', '--']

    # colour legend
    col_label = s
    col_lines = [Line2D([0], [0], color=col) for col in colors]

    # Fig2A : P_absorbing
    fig1, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(3.4, 5.0))
    plt.subplots_adjust(wspace=0, hspace=0.0)
    # fig1, ax1 = plt.subplots(figsize=(3.4, 2.5))
    ax1.set(title=r"",
            # xlabel=r"Relative abundance, $x$",
            ylabel=r"Probability fixation, $P_N(x)$")

    # Fig2B : MFPT function of x
    # fig2, ax2 = plt.subplots(figsize=(3.4, 2.5))
    ax2.set_xlim([0.0, 1.0])
    ax2.set(title=r"",
            xlabel=r"Relative abundance, $x$",
            ylabel=r"MFPT, $\tau(x)$")

    # Fig2D : MFPT as function of N
    fig4, ax4 = plt.subplots(figsize=(3.4, 2.5))
    ax4.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT, $\tau(x_{\mathrm{max}})$")

    moran = []
    space = []
    # ME_mfpt_space_N = []
    for i, fit in enumerate(s):
        # define models
        moran.append(MoranFPT(fit, 1.0, N, times))
        space.append(model(fit, 1.0, N, times))

        # plots
        x = []
        ME_mfpt_space = []
        ME_prob_space = []

        # Fig 1 A
        for j in np.linspace(1, N - 1, 9):
            x.append(j / N)
            j = int(j)

            # ME approach
            ME_prob, _ = space[i].probability_mfpt(j)
            ME_prob_space.append(ME_prob[1])

        X, FP_prob_space = space[i].FP_prob()
        X, FP_prob_moran = moran[-1].FP_prob()
        ax1.plot(X, FP_prob_moran, linestyle=linestyle[0], color=colors[i],
                 zorder=0)
        ax1.plot(X, FP_prob_space, linestyle=linestyle[1], color=colors[i],
                 zorder=i)
        # ax1.scatter(x, ME_prob_space, color=colors[i], marker=marker,
        #            s=12, edgecolors='k', linewidth=1, zorder=i+.5)

        # Fig 2 D
        x = []
        ME_mfpt_space = []
        for j in np.linspace(1, N - 1, 9):
            x.append(j / N)
            j = int(j)

            # ME approach
            ME_prob, ME_mfpt = space[i].probability_mfpt(j)
            ME_mfpt_space.append(np.dot(ME_prob, ME_mfpt))

        X, FP_mfpt_moran = moran[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_moran, linestyle=linestyle[0], color=colors[i])
        X, FP_mfpt_space = space[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_space, linestyle=linestyle[1], color=colors[i],
                 zorder=i)
        # ax2.scatter(x, ME_mfpt_space, color=colors[i], marker=marker,
        #            s=12, edgecolors='k', linewidth=1, zorder=i + .5)

    # deterministic approximation times
    S = 100
    x_max = 1 / (1 + np.sqrt(S))
    x_det = np.linspace(x_max, 1.0, 1000)
    t_det = - (np.log(np.abs((S - 1) * x_det**2 + 2 * x_det - 1)) - np.log(S))
    ax2.plot(x_det, t_det, linestyle=linestyle[2], color='k', zorder=i+1)
    x_det = np.linspace(0.0, x_max, 1000)
    t_det2 = - (sp.log((S - 1) * x_det**2 + 2 * x_det - 1) - sp.log(-1))
    ax2.plot(x_det, t_det2, linestyle=linestyle[2], color='k', zorder=i+1)

    # Fig 1C: MFPT as a function of N at x_max
    N_func = np.linspace(10, 1000, 100)
    for i, fit in enumerate(s):
        _, FP_mfpt_moran_N = moran[i].FP_mfpt_x(x_find=moran[i].xmax, N=N_func)
        _, FP_mfpt_space_N = space[i].FP_mfpt_x(x_find=space[i].xmax, N=N_func)

        ax4.plot(N_func, FP_mfpt_moran_N, color=colors[i],
                 linestyle=linestyle[0])
        ax4.plot(N_func, FP_mfpt_space_N, color=colors[i],
                 linestyle=linestyle[1])

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
    ax1.set_ylim([0.00, 1.1])
    cl = ax1.legend(col_lines, col_label, loc='lower right', title=r'$s$')
    # ax1.legend(custom_lines, legend, title=r'Model', loc='upper right',
    #           framealpha=0.8)
    ax1.add_artist(cl)

    ax2.set_yscale('log')
    ax2.set_ylim(0.1, 100)
    # ax2.legend(det, [r'$\tau_{\mathrm{det}}$'])

    ax4.set_yscale('log')
    ax4.set_xscale('log')
    ax4.set_ylim([1.0, 100])
    ax4.set_xlim([np.min(N_func), np.max(N_func)])
    c2 = ax4.legend(col_lines, col_label, loc='upper right', title=r'$s$')
    custom_lines.pop(0)
    legend.pop(0)
    # ax4.legend(custom_lines, legend, title=r'Model', loc='upper left',
    #           framealpha=0.8)
    ax4.add_artist(c2)

    # save figures
    fig1.savefig(filename + '_prob.pdf')
    fig1.savefig(filename + '_prob.png')
    # fig2.savefig(filename + '_mfpt.pdf')
    # fig2.savefig(filename + '_mfpt.png')
    fig4.savefig(filename + '_N.pdf')
    fig4.savefig(filename + '_N.png')
    return s


def asymptotics(s, cmap_name, filename, model):

    N = np.linspace(10, 1000, 100)
    K = 100  # a placeholder

    # colour legend
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(s))]
    linestyle = [':', '-', '--']
    col_label = s
    col_lines = [Line2D([0], [0], color=col) for col in colors]

    def x_det_o(s, N):
        a = np.sqrt(np.pi * N * np.sqrt(s)) / ((1 + np.sqrt(s)))
        squ = a * (np.sqrt(s) - 1) + 1 / (2 * np.sqrt(s) * a)
        divP = (s + 1.) + a * (s - 1.)
        divM = (s + 1.) - a * (s - 1.)
        return (1 + squ) / divP, (1 - squ) / divM

    def x_det(s, N):  # from Wolfram
        a = np.sqrt(np.pi * N * np.sqrt(s)) / ((1 + np.sqrt(s)))
        squ = np.sqrt(s * (a ** 2 - 1))
        divP = a * (s - 1.) + (s + 1)
        divM = a * (s - 1.) - (s + 1)
        return (1 + squ - a) / divP, (-1 + squ - a) / divM

    def time_det(s, N):
        x1, x2 = x_det(s, N)
        t1 = -(sp.log((s - 1) * x1 ** 2 + 2 * x1 - 1) - sp.log(s))
        t2 = -(sp.log((s - 1) * x2 ** 2 + 2 * x2 - 1) - sp.log(-1))
        return (t1 + t2) / 2

    # Fig3A : potentials
    fig1, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(3.4, 5.0))
    plt.subplots_adjust(wspace=0, hspace=0.0)
    # fig1, ax1 = plt.subplots(figsize=(3.4, 2.5))
    ax1.set_xlim([0.0, 1.0])
    ax1.set(title=r"",
            # xlabel=r"$x$",
            ylabel=r"$U(x)$")

    # Fig3B : x_t as a function of N
    # fig2, ax2 = plt.subplots(figsize=(3.4, 2.5))
    ax2.set(title=r"",
            xlabel=r"Relative abundance, $x$",
            ylabel=r"Total population size, $N$")

    # Fig3C : MFPT asymptotic function of N
    fig3, ax3 = plt.subplots(figsize=(3.4, 2.5))
    ax3.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT")

    space = []
    for i, fit in enumerate(s):
        space.append(model(fit, 1.0, K, times))
        # Fig 3A - U(x)
        x_pot = np.linspace(0., 1., 100)
        ax1.plot(x_pot, space[i].potential(x_pot),
                 linestyle=linestyle[1], color=colors[i], zorder=i)
        ax1.scatter([space[i].xmax], (space[i].potential(space[i].xmax)
                                      / space[i].K * K),
                    marker='X', color=colors[i], s=20, edgecolors='k',
                    linewidth=1, zorder=i + .5)

        # Fig 3B
        x1, x2 = x_det(fit, N)
        ax2.plot(x1, N, color=colors[i], linestyle=linestyle[2])
        ax2.plot(x2, N, color=colors[i], linestyle=linestyle[2])
        ax2.fill_betweenx(N, x2, x1, color=colors[i], alpha=0.2)
        ax2.axvline(model(fit, 1, K, 0).xmax, linestyle=linestyle[0],
                    color=colors[i])

    # Fig 3C
    for i, fit in enumerate(s):
        space = model(fit, 1, K, 0)
        _, FP_mfpt_space_N = space.FP_mfpt_x(x_find=space.xmax, N=N)

        ax3.plot(N, FP_mfpt_space_N, color=colors[i], linestyle=linestyle[1],
                 zorder=i)

        tdet = time_det(fit, N)
        # tdif = ((-np.log(0.5) + np.log(0.5 + x1))/4 + x1 * np.arctan(x1))
        tdif = np.log(2) / 2
        ax3.plot(N, tdet + tdif, color=colors[i], linestyle=linestyle[2],
                 zorder=i+0.5)

        # approximation
        """
        appr = (np.log(np.pi * N * np.sqrt(fit) * (fit - 1) ** 2)
                / (2 * fit * (1 + np.sqrt(fit)) ** 2 * (np.sqrt(fit) - 1) ** 2)
                / 2.0)
        appr = np.log(2 * np.pi * N) / 2
        ax3.plot(N, appr, color='k', linestyle=linestyle[2],
                 zorder=i+0.5)
        """

    # legend and limits
    cl = ax1.legend(col_lines, col_label, loc='lower right', title=r'$s$')
    ax1.set_yscale('log')
    # ax1.legend([Line2D([0], [0], markerfacecolor='dimgray', linestyle='None',
    #           marker='X', color='k', linewidth=1)],
    #           [r'$x_{\mathrm{max}}$'], loc='upper right', framealpha=0.8)
    ax1.add_artist(cl)

    custom_label1 = [r'$x_{\mathrm{max}}$', r'$x_t$']
    custom_lines1 = [
        Line2D([0], [0], color='dimgray', linestyle=linestyle[0]),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[2])
        ]
    ax2.set_yscale('log')
    ax2.set_xlim(0.0, 1.0)
    ax2.set_ylim(10, 1000)
    ax2.legend(custom_lines1, custom_label1)

    plt.subplots_adjust(wspace=0, hspace=0.0)

    ax3.set_xscale('log')
    ax3.set_yscale('log')
    custom_label2 = [r'$\tau(x_{\mathrm{max}})$',
                     r'$\tau_{a}(x_{\mathrm{max}})$']
    custom_lines2 = [
        Line2D([0], [0], color='dimgray', linestyle=linestyle[1]),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[2])
        ]
    ax3.legend(custom_lines2, custom_label2)
    ax3.set_xlim([np.min(N), np.max(N)])

    fig1.savefig(filename + '_pot.tiff', dpi=300)
    fig1.savefig(filename + '_pot.pdf')
    fig1.savefig(filename + '_pot.png')
    # fig2.savefig(filename + '_xdet.pdf')
    # fig2.savefig(filename + '_xdet.png')
    fig3.savefig(filename + '_tdet.pdf')
    fig3.savefig(filename + '_tdet.png')

    return 0


def invasion(s, N, cmap_name, filename, model):
    marker = ['+', 'x']
    linestyle = [':', '-', '--']
    markeredge = None
    ms = 5
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(s))]

    # colour legend
    col_label = s
    col_lines = [Line2D([0], [0], color=col) for col in colors]

    # Fig4A : P_absorbing
    # fig1, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(3.4, 5.0))
    plt.subplots_adjust(wspace=0, hspace=0.0)
    fig1, ax1 = plt.subplots(figsize=(3.4, 2.5))
    ax1.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"Probability invasion succeeds")

    # Fig4B : MFPT function of x
    fig2, ax2 = plt.subplots(figsize=(3.4, 2.5))
    ax2.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT successful invasion")

    fig3, ax3 = plt.subplots(figsize=(3.4, 2.5))
    ax3.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT failure")

    fig4, ax4 = plt.subplots(figsize=(3.4, 2.5))
    ax4.set(title=r"",
            xlabel=r"Relative abundance, $x$",
            ylabel=r"Probability invasion fixates")

    moran = [[] for fit in s]
    space = [[] for fit in s]
    
    Nplot = 100
    for i, fit in enumerate(s):
        prob_loc = np.zeros(Nplot)
        # Fig 4 calcs
        ME_prob_moran = []
        ME_prob_space = []
        ME_fptS_moran = []
        ME_fptS_space = []
        ME_fptF_moran = []
        ME_fptF_space = []
        for j, K in enumerate(N):
            # define models
            moran[i].append(MoranFPT(fit, 1.0, K, times))
            space[i].append(model(fit, 1.0, K, times))

            # Moran is positionally independent invasion
            prob, mfpt = moran[i][j].probability_mfpt(1)
            ME_prob_moran.append(prob[1])
            ME_fptS_moran.append(mfpt[1])
            ME_fptF_moran.append(mfpt[0])

            # ME approach
            mfpt_tot = np.array([0., 0., 0.])
            prob_tot = np.array([0., 0., 0.])
            for pos in range(K):
                prob, mfpt = space[i][j].probability_mfpt(pos, K - (pos + 1))
                if K == Nplot:
                    prob_loc[pos] = prob[0]
                prob_tot += np.array(prob) / K
                mfpt_tot += np.array(prob) * np.array(mfpt) / K
            ME_prob_space.append(prob_tot[0])
            ME_fptS_space.append(mfpt_tot[0])
            ME_fptF_space.append(mfpt_tot[1] + mfpt_tot[2])

        ax1.plot(N, ME_prob_moran, color=colors[i], marker=marker[0],
                 markersize=ms + 1, markeredgecolor=markeredge, linewidth=1,
                 linestyle=linestyle[0], zorder=i + .5)
        ax1.plot(N, ME_prob_space, color=colors[i], marker=marker[1],
                 markersize=ms, markeredgecolor=markeredge, linewidth=1,
                 linestyle=linestyle[1], zorder=i + .5)
        ax2.plot(N, ME_fptS_moran, color=colors[i], marker=marker[0],
                 markersize=ms + 1, markeredgecolor=markeredge, linewidth=1,
                 linestyle=linestyle[0], zorder=i + .5)
        ax2.plot(N, ME_fptS_space, color=colors[i], marker=marker[1],
                 markersize=ms, markeredgecolor=markeredge, linewidth=1,
                 linestyle=linestyle[1], zorder=i + .5)
        ax3.plot(N, ME_fptF_moran, color=colors[i], marker=marker[0],
                 markersize=ms + 1, markeredgecolor=markeredge, linewidth=1,
                 linestyle=linestyle[0], zorder=i + .5)
        ax3.plot(N, ME_fptF_space, color=colors[i], marker=marker[1],
                 markersize=ms, markeredgecolor=markeredge, linewidth=1,
                 linestyle=linestyle[1], zorder=i + .5)
        ax4.plot((np.arange(Nplot) + 0.5) / Nplot, prob_loc, color=colors[i],
                 marker=marker[1], markersize=ms, markeredgecolor=markeredge,
                 linewidth=1, linestyle=linestyle[1], zorder=i + .5)

    # figure limits and legends
    # legend
    # legend = ['moran (ME)', 'spatial (ME)']
    """
    custom_lines = [
        Line2D([0], [0], markerfacecolor='dimgray', marker=marker[0],
               color='k', linewidth=1),
        Line2D([0], [0], markerfacecolor='dimgray', marker=marker[1],
               color='k', linewidth=1)
        ]
    """

    # limits
    ax1.set_xlim([np.min(N), np.max(N)])
    ax1.set_ylim([0.01, 1.1])
    cl = ax1.legend(col_lines, col_label, loc='lower right', title=r'$s$',
                    framealpha=0.8)
    # ax1.legend(custom_lines, legend, title=r'Model', loc='upper right',
    #           framealpha=0.8)
    ax1.add_artist(cl)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.xaxis.set_minor_formatter(NullFormatter())

    # ax2.set_yscale('log')
    ax2.set_xlim([np.min(N), np.max(N)])
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.xaxis.set_minor_formatter(NullFormatter())

    # ax3.set_yscale('log')
    ax3.set_xlim([np.min(N), np.max(N)])
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.xaxis.set_minor_formatter(NullFormatter())

    # ax4.set_yscale('log')
    ax4.set_xlim([0.0, 1.0])
    ax4.set_ylim([0.0001, 1.0])
    ax4.set_yscale('log')
    # ax4.set_xscale('log')
    # ax4.xaxis.set_minor_formatter(NullFormatter())

    # save figures
    fig1.savefig(filename + '_prob.pdf')
    fig1.savefig(filename + '_prob.png')
    fig2.savefig(filename + '_mfptS.pdf')
    fig2.savefig(filename + '_mfptS.png')
    # fig3.savefig(filename + '_mfptF.pdf')
    # fig3.savefig(filename + '_mfptF.png')
    fig4.savefig(filename + '_probloc.pdf')
    fig4.savefig(filename + '_probloc.png')
    return 0


if __name__ == '__main__':
    plt.style.use('python_src/custom.mplstyle')

    # mkdir
    dir = 'figures_manuscript'
    Path(dir).mkdir(parents=True, exist_ok=True)

    # for certain distribution functions
    times = np.linspace(0.0, 100.0, 10001)

    # Figure 1 : spatial vs moran
    N = [10, 100, 1000]
    fname1 = dir + os.sep + 'SvM'
    cmap_name = 'viridis'
    spatial_vs_moran(N, times, cmap_name, fname1, Spatial)

    # Figure 2 : fitness results
    s = [1., 10., 100.]
    K = 100
    fname2 = dir + os.sep + 'fit'
    cmap_name1 = 'plasma'
    fitness_spatial(s, K, times, cmap_name1, fname2, Spatial)

    # Figure 3 : asymptotics
    s = [1., 10., 100.]
    K = 100
    fname3 = dir + os.sep + 'asy'
    asymptotics(s, cmap_name1, fname3, Spatial)

    # Figure 4 : Invasion
    s = [1., 10, 100.]
    N = np.logspace(np.log10(11), np.log10(101), 10, dtype=int)
    N = np.logspace(1, 2, 10, dtype=int)
    fname4 = dir + os.sep + 'inv'
    invasion(s, N, cmap_name1, fname4, SpatInv)
