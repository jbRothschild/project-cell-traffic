from python_src.first_passage import MoranFPT, OneBoundaryFPT, TwoBoundaryFPT,\
    OneBoundaryIntFPT  # , TwoBoundaryIntFPT
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.lines import Line2D

plt.style.use('python_src/custom.mplstyle')


def one_boundary_comparaison(r1, r2, N, times, cmap_name, filename, model):
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 1.0, num=5)]
    linestyle = [':', '-', '--']

    # Fig1B : P_absorbing
    fig1, ax1 = plt.subplots(figsize=(4, 3))
    ax1.set(title=r"",
            xlabel=r"Fraction initial boundary position, $x$",
            ylabel=r"Probability fixation, $P_N(x)$")
    ax1.set_xlim([0.0, 1.0])

    legend = ['ME Moran', 'ME Homog.', 'FP Moran', 'FP Homog.']
    custom_lines1 = [
        Line2D([0], [0], color='dimgray', linestyle='None', marker='o'),
        Line2D([0], [0], color='dimgray', linestyle='None', marker='x'),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[0]),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[1])
        ]

    # Fig1C : MFPT function of x
    fig2, ax2 = plt.subplots(figsize=(4, 3))
    ax2.set_xlim([0.0, 1.0])
    ax2.set_yscale('log')
    ax2.set(title=r"",
            xlabel=r"Fraction initial boundary position, $x$",
            ylabel=r"MFPT, $\tau(x)$")
    custom_lines2 = custom_lines1.copy()
    legend2 = legend.copy()
    custom_lines2.append(Line2D([0], [0], color='r', linestyle=linestyle[2]))
    legend2.append(r'$\tau_{det}$')

    # Fig1D : MFPT as function of N
    fig3, ax3 = plt.subplots(figsize=(4, 3))
    ax3.set_yscale('log')
    ax3.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT, $\tau(1/2)$")

    # Fig1XX : FPT distribution at x = 1/2
    fig4, ax4 = plt.subplots()
    ax4.set_xlim([1., 1000.])
    ax4.set(title=r"",
            xlabel=r"Fraction initial boundary position, $x$",
            ylabel=r"MFPT, $\tau_N(x)$")
    custom_lines4 = [Line2D([0], [0], color='dimgray', linestyle=linestyle[0],
                            marker='D'),
                     Line2D([0], [0], color='dimgray', linestyle=linestyle[1],
                            marker='o')]

    moran = []
    space = []
    ME_mfpt_moran_N = []
    ME_mfpt_space_N = []
    for i, K in enumerate(N):
        # legend
        # legend.append(r'$N={}$'.format(K))
        # legend2.append(r'$N={}$'.format(K))
        # custom_lines1.append(Line2D([0],
        #     [0], color=colors[i], linestyle=linestyle[1]))
        # custom_lines2.append(Line2D([0],
        #     [0], color=colors[i], linestyle=linestyle[1]))
        # custom_lines4.append(Line2D([0],
        #     [0], color=colors[i], linestyle=linestyle[1]))

        # define models
        moran.append(MoranFPT(r1, r2, K, times))
        space.append(model(r1, r2, K, times))

        # plots
        x = []
        ME_mfpt_moran = []
        ME_mfpt_space = []
        ME_prob_moran = []
        ME_prob_space = []

        for j in np.linspace(1, K - 1, 9):
            x.append(j / K)
            j = int(j)

            # ME approach
            ME_prob, ME_mfpt = moran[i].probability_mfpt(j)
            ME_prob_moran.append(ME_prob[1])
            ME_mfpt_moran.append(np.dot(ME_prob, ME_mfpt))

            ME_prob, ME_mfpt = space[i].probability_mfpt(j)
            ME_prob_space.append(ME_prob[1])
            ME_mfpt_space.append(np.dot(ME_prob, ME_mfpt))

        # Fig 1 B
        X, FP_prob_moran = moran[i].FP_prob()
        ax1.plot(X, FP_prob_moran, linestyle=linestyle[0], color=colors[i])
        X, FP_prob_space = space[i].FP_prob()
        ax1.plot(X, FP_prob_space, linestyle=linestyle[1], color=colors[i])
        ax1.scatter(x, ME_prob_space, color=colors[i], marker='x', zorder=2.5)
        ax1.scatter(x, ME_prob_moran, color=colors[i], marker='o', zorder=2.5)

        # Fig 1 C
        X, FP_mfpt_moran = moran[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_moran, linestyle=linestyle[0], color=colors[i])
        X, FP_mfpt_space = space[i].FP_mfpt()
        ax2.plot(X, FP_mfpt_space, linestyle=linestyle[1], color=colors[i])
        ax2.scatter(x, ME_mfpt_space, color=colors[i], marker='x', zorder=2.5)
        ax2.scatter(x, ME_mfpt_moran, color=colors[i], marker='o', zorder=2.5)

        # Fig 1 D
        prob, mfpt = moran[i].probability_mfpt(int(K / 2))
        ME_mfpt_moran_N.append(np.dot(prob, mfpt))
        prob, mfpt = space[i].probability_mfpt(int(K / 2))
        ME_mfpt_space_N.append(np.dot(prob, mfpt))

        """
        # Fig 1 D
        fpt_moran, _ = moran[i].fpt_distribution(int(K / 2))
        fpt_space, _ = space[i].fpt_distribution(int(K / 2))
        ax4.plot(times, fpt_moran[1], linestyle=linestyle[0], color=colors[i])
        ax4.plot(times, fpt_space[1], linestyle=linestyle[1], color=colors[i])
        """

    # deterministic approximation times
    x_det = np.linspace(0.55, 1.0, 100)
    t_det = - np.log(2 * x_det - 1.0)
    ax2.plot(x_det, t_det, linestyle=linestyle[2], color='r')

    # Fig 1D: MFPT as a function of N at x/2
    N_func = np.linspace(10, 500, 100)
    _, FP_mfpt_moran_N = moran[0].FP_mfpt(x=0.5, N=N_func)
    _, FP_mfpt_space_N = space[0].FP_mfpt_x(x_find=0.5, N=N_func)
    # N_space = (np.pi / 16.) * np.sqrt(np.pi / 4) * np.exp(N_func / 4)
    # / N_func**(10/2)

    ax3.plot(N_func, FP_mfpt_moran_N, label=r'FP Moran', color='dimgray',
             linestyle=linestyle[0])
    ax3.plot(N_func, FP_mfpt_space_N, color='black', label=r'FP Inhomog.')
    ax3.scatter(N, ME_mfpt_moran_N, label=r'ME Moran', marker='o',
                color=colors[0: len(N)], zorder=2.5)
    ax3.scatter(N, ME_mfpt_space_N, label=r'ME Homog.', marker='x',
                color=colors[0: len(N)], zorder=2.5)
    custom_lines3 = [Line2D([0], [0], color='dimgray', linestyle='None',
                            marker='o'),
                     Line2D([0], [0], color='dimgray', linestyle='None',
                            marker='x'),
                     Line2D([0], [0], color='dimgray', linestyle=linestyle[0]),
                     Line2D([0], [0], color='dimgray', linestyle=linestyle[1])
                     ]
    legend3 = legend.copy()
    legend3.pop()
    ax3.set_ylim([1.0, 1000])

    # ax3.plot(N_func, N_space, color='black', label='Asymptotic expansion')

    ax2.set_ylim(0.01, 100)
    ax1.legend(custom_lines1, legend)
    ax2.legend(custom_lines2, legend2)
    ax3.legend(custom_lines3, legend3)
    ax4.legend(custom_lines4, legend)

    # save figures
    fig1.savefig(filename + '_prob.pdf')
    fig1.savefig(filename + '_prob.png')
    fig2.savefig(filename + '_mfpt.pdf')
    fig2.savefig(filename + '_mfpt.png')
    fig3.savefig(filename + '_N.pdf')
    fig3.savefig(filename + '_N.png')

    return 0


def one_boundary_ratio(r1, r2, N, times, cmap_name, filename):
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 1.0, num=5)]
    # linestyle = [':', '-', '--']

    # Fig1B : P_absorbing
    fig1, ax1 = plt.subplots(figsize=(4, 3))
    ax1.set(title=r"",
            xlabel=r"Fraction initial boundary position, $x$",
            ylabel=r"Fixation ratio, $ P_{N}^{direc.}(x)/P_{N}^{homog.}(x)$")
    ax1.set_xlim([0.0, 1.0])
    # ax1.set_yscale('log')

    # Fig1C : MFPT function of x
    fig2, ax2 = plt.subplots(figsize=(4, 3))
    ax2.set_xlim([0.0, 1.0])
    # ax2.set_yscale('log')
    ax2.set(title=r"",
            xlabel=r"Fraction initial boundary position, $x$",
            ylabel=r"MFPT ratio, $\tau_N^{direc.}(x)/\tau_N^{homog.}(x)$")

    # Fig1D : MFPT as function of N
    fig3, ax3 = plt.subplots(figsize=(4, 3))
    # ax3.set_yscale('log')
    ax3.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"Ratio MFPT, $\tau_N^{direc.}(1/2)/\tau_N^{homog.}(1/2)$")

    # Fig1XX : FPT distribution at x = 1/2
    fig4, ax4 = plt.subplots()
    ax4.set_yscale('log')
    ax4.set(title=r"",
            xlabel=r"Fraction initial boundary position, $x$",
            ylabel=r"MFPT ratio, $\tau_N^{direc.}(x)/$")

    homog = []
    direc = []
    ME_mfpt_homog_N = []
    ME_mfpt_direc_N = []
    for i, K in enumerate(N):

        # define models
        homog.append(OneBoundaryFPT(r1, r2, K, times))
        direc.append(OneBoundaryIntFPT(r1, r2, K, times))

        # plots
        x = []
        ME_mfpt_direc = []
        ME_mfpt_homog = []
        ME_prob_direc = []
        ME_prob_homog = []

        for j in np.linspace(1, K - 1, 9):
            x.append(j / K)
            j = int(j)

            # ME approach
            ME_prob, ME_mfpt = homog[i].probability_mfpt(j)
            ME_prob_homog.append(ME_prob[1])
            ME_mfpt_homog.append(np.dot(ME_prob, ME_mfpt))

            ME_prob, ME_mfpt = direc[i].probability_mfpt(j)
            ME_prob_direc.append(ME_prob[1])
            ME_mfpt_direc.append(np.dot(ME_prob, ME_mfpt))

        # Fig 1 B
        ratio_prob = [
            ME_prob_direc[j] / val for j, val in enumerate(ME_prob_homog)]

        ax1.plot(x, ratio_prob, color=colors[i], marker='D')

        # Fig 1 C
        ratio_mfpt = [
            ME_mfpt_direc[j] / val for j, val in enumerate(ME_mfpt_homog)]
        ax2.plot(x, ratio_mfpt, color=colors[i], marker='D')

        # Fig 1 D
        prob, mfpt = homog[i].probability_mfpt(int(K / 2))
        ME_mfpt_homog_N.append(np.dot(prob, mfpt))
        prob, mfpt = direc[i].probability_mfpt(int(K / 2))
        ME_mfpt_direc_N.append(np.dot(prob, mfpt))

        """
        # Fig 1 D
        fpt_moran, _ = moran[i].fpt_distribution(int(K / 2))
        fpt_space, _ = space[i].fpt_distribution(int(K / 2))
        ax4.plot(times, fpt_moran[1], linestyle=linestyle[0], color=colors[i])
        ax4.plot(times, fpt_space[1], linestyle=linestyle[1], color=colors[i])
        """

    # deterministic approximation times
    # x_det = np.linspace(0.55, 1.0, 100)
    # t_det = - np.log(2 * x_det - 1.0)
    # ax2.plot(x_det, t_det, linestyle=linestyle[2], color='k')

    # Fig 1D: MFPT as a function of N at x/2
    # N_func = np.linspace(10, 1000, 100)
    # X, FP_mfpt_moran_N = moran[0].FP_mfpt(x=0.5, N=N_func)
    # X, FP_mfpt_space_N = space[0].FP_mfpt(x=0.5, N=N_func)
    # N_space = (np.pi / 16.) * np.sqrt(np.pi / 4) * ( np.exp(N_func / 4)
    #     / N_func**(10/2) )

    # ax3.scatter(N, ME_mfpt_moran_N, label=r'ME Moran', marker='o',
    # color=colors[0: len(N)])
    # ax3.plot(N_func, FP_mfpt_moran_N, label=r'FP Moran', color='dimgray',)
    ratio_N = [
        ME_mfpt_direc_N[j] / val for j, val in enumerate(ME_mfpt_homog_N)]
    ax3.scatter(N, ratio_N, label=r'ME Homog.', marker='D',
                color=colors[0: len(N)], zorder=2.5)
    # ax3.plot(N_func, FP_mfpt_space_N, color='black', label=r'FP Homog')

    # save figures
    fig1.savefig(filename + '_prob.pdf')
    fig1.savefig(filename + '_prob.png')
    fig2.savefig(filename + '_mfpt.pdf')
    fig2.savefig(filename + '_mfpt.png')
    fig3.savefig(filename + '_N.pdf')
    fig3.savefig(filename + '_N.png')

    return 0


def comparaison_FP_approx(r1, r2, N, times, cmap_name, filename, model):
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 1.0, num=5)]
    linestyle = [':', '-', '--']

    legend = ['FP', 'Asymptotic']
    custom_lines1 = [
        Line2D([0], [0], color='dimgray', linestyle=linestyle[0]),
        Line2D([0], [0], color='dimgray', linestyle=linestyle[1])
        ]

    # Fig2B : MFPT function of x
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.set_xlim([0.0, 1.0])
    ax.set_yscale('log')
    ax.set(title=r"",
           xlabel=r"Fraction initial boundary position, $x$",
           ylabel=r"MFPT, $\tau(x)$")

    # Fig1D : MFPT as function of N
    fig2, ax2 = plt.subplots(figsize=(4, 3))
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT, $\tau_{max}$")

    moran = []
    space = []
    ME_mfpt_moran_N = []
    ME_mfpt_space_N = []
    for i, K in enumerate(N):
        moran.append(MoranFPT(r1, r2, K, times))
        space.append(model(r1, r2, K, times))

        # plots
        x = []
        ME_mfpt_moran = []
        ME_mfpt_space = []
        ME_prob_moran = []
        ME_prob_space = []

        for j in np.linspace(1, K - 1, 9):
            x.append(j / K)
            j = int(j)

            # ME approach
            ME_prob, ME_mfpt = moran[i].probability_mfpt(j)
            ME_prob_moran.append(ME_prob[1])
            ME_mfpt_moran.append(np.dot(ME_prob, ME_mfpt))

            ME_prob, ME_mfpt = space[i].probability_mfpt(j)
            ME_prob_space.append(ME_prob[1])
            ME_mfpt_space.append(np.dot(ME_prob, ME_mfpt))

        # Fig 2 A
        x = np.linspace(0., 1., 101)
        X, FP_mfpt_space = space[i].FP_mfpt()
        FP_asym_space = space[i].mfpt_asymp(x)
        ax.plot(X, FP_mfpt_space, linestyle=linestyle[0], color=colors[i])
        ax.plot(x, FP_asym_space, linestyle=linestyle[1], color=colors[i])

        # Fig 2 B
        prob, mfpt = moran[i].probability_mfpt(int(K / 2))
        ME_mfpt_moran_N.append(np.dot(prob, mfpt))
        prob, mfpt = space[i].probability_mfpt(int(K / 2))
        ME_mfpt_space_N.append(np.dot(prob, mfpt))

    # Fig 2B: MFPT as a function of N at x/2
    N_func = np.linspace(10, 500, 100)
    _, FP_mfpt_moran_N = moran[0].FP_mfpt(x=0.5, N=N_func)
    FP_mfpt_space_N = ((space[0].mfpt_asymp(np.array(space[0].max_pot)))
                       * np.sqrt(N_func) / np.sqrt(space[0].K))

    ax2.plot(N_func, FP_mfpt_moran_N, label=r'FP Moran', color='dimgray',
             linestyle=linestyle[0])
    ax2.plot(N_func, FP_mfpt_space_N, color='black', label=r'FP Inhomog.')
    ax2.scatter(N, ME_mfpt_moran_N, label=r'ME Moran', marker='o',
                color=colors[0: len(N)], zorder=2.5)
    ax2.scatter(N, ME_mfpt_space_N, label=r'ME Homog.', marker='x',
                color=colors[0: len(N)], zorder=2.5)
    custom_lines2 = [Line2D([0], [0], color='dimgray', linestyle='None',
                            marker='o'),
                     Line2D([0], [0], color='dimgray', linestyle='None',
                            marker='x'),
                     Line2D([0], [0], color='dimgray', linestyle=linestyle[0]),
                     Line2D([0], [0], color='dimgray', linestyle=linestyle[1])
                     ]
    legend2 = ['Moran', 'Spatial']
    ax2.set_ylim([1.0, 1000])

    ax.set_ylim(0.01, 100)
    ax.legend(custom_lines1, legend)
    ax2.legend(custom_lines2, legend2)

    # save figures
    fig.savefig(filename + '_mfpt_approx.pdf')
    fig.savefig(filename + '_mfpt_approx.png')
    fig2.savefig(filename + '_mfpt_approx_N.pdf')
    fig2.savefig(filename + '_mfpt_approx_N.png')
    return 0


def two_boundary_comparaison(r1, r2, N, times, cmap_name, filename, model):
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 1.0, num=5)]
    linestyle = [':', '-', '--']

    # Fig1B : P_absorbing
    fig1, ax1 = plt.subplots(figsize=(4, 3))
    ax1.set(title=r"",
            xlabel=r"Fraction initial invasion position, $x$",
            ylabel=r"Probability invasion success, $P_N(x)$")
    ax1.set_xlim([0.0, 1.0])

    legend = ['ME Moran', 'ME Inv. Homog.']
    custom_lines1 = [Line2D([0], [0], color='dimgray', linestyle=linestyle[0]),
                     Line2D([0], [0], color='dimgray', linestyle=linestyle[1],
                            marker='x')]

    # Fig1C : MFPT function of x
    fig2, ax2 = plt.subplots(figsize=(4, 3))
    ax2.set_xlim([0.0, 1.0])
    ax2.set_yscale('log')
    ax2.set(title=r"",
            xlabel=r"Fraction initial invasion position, $x$",
            ylabel=r"MFPT invasion, $\tau_N(x)$")

    # Fig1D : MFPT as function of N
    fig3, ax3 = plt.subplots(figsize=(4, 3))
    # ax3.set_yscale('log')
    ax3.set(title=r"",
            xlabel=r"Total population size, $N$",
            ylabel=r"MFPT, $\tau_N(1/2)$")

    moran = []
    space = []
    ME_mfpt_moran_N = []
    ME_mfpt_space_N = []
    for i, K in enumerate(N):

        # define models
        moran.append(MoranFPT(r1, r2, K, times))
        space.append(model(r1, r2, K, times))

        # plots
        x = []
        ME_mfpt_space = []
        ME_prob_space = []

        for j in np.linspace(0, K - 1, 11):
            x.append(j / (K - 1))
            j = int(np.floor(j))
            k = K - j - 1

            # ME approach
            ME_prob, ME_mfpt = space[i].probability_mfpt(j, k)
            ME_prob_space.append(ME_prob[0])
            ME_mfpt_space.append(ME_mfpt[0])

        ME_prob, ME_mfpt = moran[i].FP_mfpt_N(1 / K)
        ME_mfpt_moran_N.append(ME_mfpt)

        # Fig 3 B
        ax1.plot(x, ME_prob_space, color=colors[i], marker='x',
                 linestyle=linestyle[1])
        ax1.hlines(ME_prob, 0.0, 1.0, color=colors[i],
                   linestyle=linestyle[0])

        # Fig 3 C
        ax2.plot(x, ME_mfpt_space, color=colors[i], marker='x',
                 linestyle=linestyle[1])
        ax2.hlines(ME_mfpt, 0.0, 1.0, color=colors[i], linestyle=linestyle[0])

        # calc Figure 3 D
        prob, mfpt = space[i].probability_mfpt(np.floor(K / 2),
                                               np.floor(K / 2))
        ME_mfpt_space_N.append(mfpt[0])

    ax3.scatter(N, ME_mfpt_moran_N, marker='o', color=colors[0: len(N)])
    ax3.scatter(N, ME_mfpt_space_N, marker='x', color=colors[0: len(N)])
    custom_lines3 = [Line2D([0], [0], color='dimgray', linestyle='None',
                            marker='o'),
                     Line2D([0], [0], color='dimgray', linestyle='None',
                            marker='x')]

    # ax3.plot(N_func, N_space, color='black', label='Asymptotic expansion')

    ax1.legend(custom_lines1, legend)
    ax2.legend(custom_lines1, legend)
    ax3.legend(custom_lines3, legend)

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
    N = [10, 50, 100, 500, 1000]
    # N = [10, 50, 100, 500, 1000]
    times = np.linspace(0.0, 100.0, 10001)
    cmap_name = 'viridis'
    r1 = 1.0
    r2 = 1.0

    fname = dir + os.sep + 'asymp'
    comparaison_FP_approx(r1, r2, N, times, cmap_name, fname, OneBoundaryFPT)
    '''
    fname = dir + os.sep + 'homog'
    one_boundary_comparaison(r1, r2, N, times, cmap_name, fname,
        OneBoundaryFPT)

    N = [10, 50, 100]
    fname = dir + os.sep + 'ratio_lin_line'
    # one_boundary_ratio(r1, r2, N, times, cmap_name, fname)
    '''
    """
    N = [10, 50, 100]
    N = [x + 1 for x in N]
    fname = dir + os.sep + 'inv_homog_100'
    cmap_name = 'plasma'
    two_boundary_comparaison(r1, r2, N, times, cmap_name, fname,
                             TwoBoundaryFPT)
    """
    """
    # single = OneBoundaryIntFPT(r1, r2, N[0], times)
    bound = TwoBoundaryFPT(r1, r2, N[0], times)

    # prob1, mfpt1 = single.probability_mfpt(5)
    prob, mfpt = bound.probability_mfpt(5, 5)
    """

    # print(prob1)
    # print(prob, mfpt)
