import scipy as sp
import numpy as np
from scipy import sparse
import mpmath
import scipy.sparse.linalg as sLA
from scipy.special import erf, erfi  # , dawsn
# import scipy.linalg as LA


class FirstPassage:
    def __init__(self):
        self.index = None

    def probability_mfpt(self, *initialize):
        # initial state of the system
        init_prob = np.zeros(self.nbr_states)
        init_state = self.index(*initialize)
        init_prob[init_state] = 1.0

        # calculation from Iyer-Biswas-Zilman Eq. 47
        prob_absorption = [-np.dot(((self.inverse).T @
                                    self.transition_mat[i, self.remove_abs].T
                                    ).toarray().reshape(-1,),
                           init_prob[self.remove_abs]
                                   ) for i in self.abs_idx
                           ]

        # calculation from Iyer-Biswas-Zilman Eq. 48
        mfpt = [0 if prob_absorption[i] == 0
                else np.dot(((self.inverse @ self.inverse).T @
                             self.transition_mat[state, self.remove_abs].T
                             ).toarray().reshape(-1,),
                            init_prob[self.remove_abs]
                            ) / prob_absorption[i]
                for i, state in enumerate(self.abs_idx)
                ]

        return prob_absorption, mfpt

    def fpt_distribution(self, *initialize):
        init_state = self.index(*initialize)
        init_prob = np.zeros(self.nbr_states)
        init_prob[init_state] = 1.0

        prob_absorption = [-np.dot(((self.inverse).T @
                                    self.transition_mat[i, self.remove_abs].T
                                    ).toarray().reshape(-1,),
                                   init_prob[self.remove_abs]
                                   ) for i in self.abs_idx
                           ]

        # calculation from Iyer-Biswas-Zilman (between Eq. 47 and 48)
        fpt_dist_absorp = [np.matmul(sLA.expm_multiply(self.transition_mat,
                                                       init_prob,
                                                       start=self.times[0],
                                                       stop=self.times[-1],
                                                       num=len(self.times)
                                                       )[:, self.remove_abs],
                                     (
                                     self.transition_mat[state,
                                                         self.remove_abs]
                                     ).toarray().reshape(-1,))
                           for idx, state in enumerate(self.abs_idx)]

        total_fpt = (np.sum(np.vstack(fpt_dist_absorp), axis=0)).tolist()
        fpt_dist_absorp_N = [(x * 0).tolist() if prob_absorption[i] == 0
                             else (x / prob_absorption[i]).tolist()
                             for i, x in enumerate(fpt_dist_absorp)]

        return fpt_dist_absorp_N, total_fpt


class MoranFPT(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times

        def index_moran(arg):
            # simple in Moran model, but is more complicated in general when
            # high dimensional state space
            return arg

        self.index = index_moran

        # total number of states
        self.nbr_states = self.K + 1

        # absorbing points of Moran model
        self.abs_idx = [0, self.K]
        self.abs_idx = [0, self.nbr_states - 1]
        self.remove_abs = np.ones(self.nbr_states, dtype=bool)
        self.remove_abs[self.abs_idx] = False

        # transition matrix
        self.transition_mat = np.zeros((self.nbr_states, self.nbr_states))

        # main matrix elements
        for i in np.arange(1, self.nbr_states - 1):
            self.transition_mat[i - 1, i] = self.r1 * i * (self.K - i)
            self.transition_mat[i, i] = (- (self.r1 + self.r2) * i
                                         * (self.K - i))
            self.transition_mat[i + 1, i] = self.r2 * i * (self.K - i)
        self.transition_mat *= (self.K) / (self.K ** 2)

        # truncate transition matrix to not have singular matrix for inverse
        trunc_trans_mat = np.delete(self.transition_mat, self.abs_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, self.abs_idx, 1)

        # inverse and then make it sparse
        self.inverse = sp.linalg.inv(trunc_trans_mat)
        self.inverse = sparse.csc_matrix(self.inverse)
        self.transition_mat = sparse.csc_matrix(self.transition_mat)

        return

    def FP_prob(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(0, 1, 101)
        P = x
        return x, P

    def FP_mfpt(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(1.0 / N, 1.0 - 1.0 / N, 101)
        mfpt = - N * ((1 - x) * np.log(1 - x) + x * np.log(x)) / self.r1
        return x, mfpt

    def FP_mfpt_N(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(1.0 / N, 1.0 - 1.0 / N, 101)
        mfpt = -N * np.log(1 - x) * (1 - x) / x
        return x, mfpt


class MoranGrowFPT(FirstPassage):

    def __init__(self, growth_rate, carrying_capacity, times):
        self.r = growth_rate
        self.K = carrying_capacity
        self.times = times

        def index_grow_moran(*args):
            # for 2 state system, devised alternate scheme for indexing the
            # states, such that i, j are represented by 1 number
            return int(args[0] * (self.K + 2) - args[0] * (args[0] + 1) / 2
                       + args[1])

        def birth(i, j):
            return self.r * i * (1. - (i + j) / self.K)

        def death(i, j):
            return self.r * 0  # no deaths happen really.

        def switch(i, j):
            return self.r * i * j / self.K

        # for external use outside of function, finding index of a 2 species
        # state.
        self.index = index_grow_moran

        # given that i+j <= self.K, this is the total number of states allowed
        self.nbr_states = int((self.K + 1) * (self.K + 2) / 2)

        # define vectors for constructing sparse transition matrix
        row = np.zeros(5 * self.nbr_states, dtype=int)  # 5 reactions per state
        col = np.zeros(5 * self.nbr_states, dtype=int)
        data = np.zeros(5 * self.nbr_states, dtype=float)

        # list of all absorbing states, gets filled below
        self.abs_idx = []

        # Building transition matrix
        rxn_count = 0
        for i in np.arange(0, self.K + 1):
            for j in np.arange(0, self.K + 1 - i):
                idx = index_grow_moran(i, j)
                if i == 0 or j == 0:
                    (self.abs_idx).append(idx)
                    # self
                row[rxn_count] = idx
                col[rxn_count] = idx
                data[rxn_count] = (- birth(i, j) - birth(j, i)
                                   - 2 * switch(i, j))
                rxn_count += 1

                # birth
                birth_idx = index_grow_moran(i - 1, j)
                if i > 0:
                    row[rxn_count] = idx
                    col[rxn_count] = birth_idx
                    data[rxn_count] = birth(i - 1, j)
                    rxn_count += 1

                # birth
                birth_idx = index_grow_moran(i, j - 1)
                if j > 0:
                    row[rxn_count] = idx
                    col[rxn_count] = birth_idx
                    data[rxn_count] = birth(j - 1, i)
                    rxn_count += 1

                # reaction growth i
                rxn_increase_idx = index_grow_moran(i - 1, j + 1)
                if i > 0 and j < self.K:
                    row[rxn_count] = idx
                    col[rxn_count] = rxn_increase_idx
                    data[rxn_count] = switch(i - 1, j + 1)
                    rxn_count += 1

                # reaction death i
                rxn_decrease_idx = index_grow_moran(i + 1, j - 1)
                if j > 0 and i < self.K:
                    row[rxn_count] = idx
                    col[rxn_count] = rxn_decrease_idx
                    data[rxn_count] = switch(i + 1, j - 1)
                    rxn_count += 1
        self.transition_mat = sparse.csc_matrix((data, (row, col)),
                                                shape=(self.nbr_states,
                                                       self.nbr_states)
                                                )
        # states that are absorbing
        self.remove_abs = np.ones(self.nbr_states, dtype=bool)
        self.remove_abs[self.abs_idx] = False

        # remove rows/cols to make the matrix invertible
        trunc_trans_mat = (self.transition_mat[self.remove_abs]
                           )[:, self.remove_abs]
        self.inverse = sLA.inv(trunc_trans_mat)

        return


class OneBoundaryFPT(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times

        def index(*args):
            # for 2 state system, devised alternate scheme for indexing the
            # states, such that i, j are represented by 1 number
            return args

        # for external use outside of function, finding the index pf the state
        self.index = index

        # total number of states
        self.nbr_states = self.K + 1

        # absorbing points of Moran model
        self.abs_idx = [0, self.nbr_states - 1]
        self.remove_abs = np.ones(self.nbr_states, dtype=bool)
        self.remove_abs[self.abs_idx] = False

        # rate definitions
        def birth(i):
            return self.r1 * i / 2.0

        def death(i):
            return self.r2 * (self.K - i) / 2.0  # no deaths happen really.

        # transition matrix
        self.transition_mat = np.zeros((self.nbr_states, self.nbr_states))

        # main matrix elements
        for i in np.arange(1, self.nbr_states - 1):
            self.transition_mat[i - 1, i] = death(i)
            self.transition_mat[i, i] = - (death(i) + birth(i))
            self.transition_mat[i + 1, i] = birth(i)

        # truncate transition matrix to not have singular matrix for inverse
        trunc_trans_mat = np.delete(self.transition_mat, self.abs_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, self.abs_idx, 1)

        # inverse and then make it sparse
        self.inverse = sp.linalg.inv(trunc_trans_mat)
        self.inverse = sparse.csc_matrix(self.inverse)
        self.transition_mat = sparse.csc_matrix(self.transition_mat)

        return

    def A(self, x):
        return self.r1 * (2. * x - 1.) / 2.

    def timeA(self, x):
        return - np.log(2. * x - 1.) / self.r1

    def B(self, x):
        return self.r1 / 2.

    def curvature(self, x):
        return - 4 * self.K

    def force(self, x):  # - grad potential
        return 2 * self.K * (2 * x - 1)

    def potential(self, x):
        return - 2 * self.K * (x**2 - x)

    def FP_prob(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(0, 1, 101)
        P = (erf(np.sqrt(N / 4.) * (2. * x - 1.)) / (2. * erf(np.sqrt(N / 4.)))
             + 1. / 2.)
        return x, P

    def FP_mfpt(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(1.0 / N, 1.0 - 1.0 / N, 101)

        mfpt = (- 2. * np.pi * erf(np.sqrt(N / 4.) * (1. - 2. * x))
                * erfi(np.sqrt(N / 4.) * (1. - 2. * x))
                + 2. * np.pi * erf(np.sqrt(N / 4.)) * erfi(np.sqrt(N / 4.))
                )
        """
        mfpt = (- 2. * np.pi * erf(np.sqrt(N / 4.) * (1. - 2. * x))
                * (2. / np.sqrt(np.pi)) * np.exp((np.sqrt(N / 4.)
                * (1. - 2. * x)) ** 2)
                * dawsn(np.sqrt(N / 4.) * (1. - 2. * x)) * (1. - 2. * x)
                + 2. * np.pi * erf(np.sqrt(N / 4.)) * (2. / np.sqrt(np.pi))
                * np.exp((np.sqrt(N / 4.)) ** 2) * dawsn(np.sqrt(N / 4.))
                )
        """

        if isinstance(N, np.ndarray):
            mfpt += (N * (1. - 2. * x) ** 2 * np.array(
                    [float(mpmath.hyp2f2(1., 1., 1.5, 2.,
                                         n * (1. - 2. * x) ** 2 / 4.))
                     for n in N]
                                                      )
                     )
            mfpt += (- N
                     * np.array(
                         [float(mpmath.hyp2f2(1., 1., 1.5, 2., n / 4.))
                          for n in N])
                     )
        elif isinstance(x, np.ndarray):
            mfpt += (N * (1. - 2. * x) ** 2
                     * np.array(
                         [float(mpmath.hyp2f2(1., 1., 1.5, 2.,
                                              N * (1. - 2. * X) ** 2 / 4.)
                                ) for X in x]
                                                      )
                     )
            mfpt += (- N * float(mpmath.hyp2f2(1., 1., 1.5, 2., N / 4.))
                     )
        else:
            mfpt += (N * (1. - 2. * x) ** 2
                     * float(mpmath.hyp2f2(1., 1., 1.5, 2.,
                                           N * (1. - 2. * x) ** 2 / 4.))
                     )
            mfpt += (- N * float(mpmath.hyp2f2(1., 1., 1.5, 2., N / 4.))
                     )
        mfpt /= 4 * self.r1
        return x, mfpt


class OneBoundaryIntFPT(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times

        def index(*args):
            # for 2 state system, devised alternate scheme for indexing the
            # states, such that i, j are represented by 1 number
            return args

        # for external use outside of function, finding the index pf the state
        self.index = index

        # total number of states
        self.nbr_states = self.K + 1

        # absorbing points of Moran model
        self.abs_idx = [0, self.nbr_states - 1]
        self.remove_abs = np.ones(self.nbr_states, dtype=bool)
        self.remove_abs[self.abs_idx] = False

        # rate definitions
        def birth(i):
            return self.r1 * i * (i - 1) / (2 * (self.K - 1))

        def death(i):
            return self.r2 * (self.K - i) * ((self.K - i - 1)
                                             / (2 * (self.K - 1)))

        # transition matrix
        self.transition_mat = np.zeros((self.nbr_states, self.nbr_states))

        # main matrix elements
        for i in np.arange(1, self.nbr_states - 1):
            self.transition_mat[i - 1, i] = death(i)
            self.transition_mat[i, i] = - (death(i) + birth(i))
            self.transition_mat[i + 1, i] = birth(i)

        # truncate transition matrix to not have singular matrix for inverse
        trunc_trans_mat = np.delete(self.transition_mat, self.abs_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, self.abs_idx, 1)

        # inverse and then make it sparse
        self.inverse = sp.linalg.inv(trunc_trans_mat)
        self.inverse = sparse.csc_matrix(self.inverse)
        self.transition_mat = sparse.csc_matrix(self.transition_mat)

        return

    def force(self, x):
        return 2 * self.K * (2 * x - 1) / (2 * x**2 - 2 * x + 1)

    def potential(self, x):
        return - 2 * self.K * np.log(2 * x**2 - 2 * x + 1)

    def FP_prob(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(0, 1, 101)
        x = [0]
        P = [0]  # np.zeros(len(x))
        return x, P

    def FP_mfpt(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(1.0 / N, 1.0 - 1.0 / N, 101)
        x = [0]
        mfpt = [0]  # np.zeros(len(x))
        return x, mfpt


class OneBoundaryFitFPT(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times

        def index(*args):
            # for 2 state system, devised alternate scheme for indexing the
            # states, such that i, j are represented by 1 number
            return args

        # for external use outside of function, finding the index pf the state
        self.index = index

        # total number of states
        self.nbr_states = self.K + 1

        # absorbing points of Moran model
        self.abs_idx = [0, self.nbr_states - 1]
        self.remove_abs = np.ones(self.nbr_states, dtype=bool)
        self.remove_abs[self.abs_idx] = False

        # rate definitions
        def birth(i):
            return (self.r1 * i
                    / (2 * (self.r2 * (self.K - i) + self.r1 * i)))

        def death(i):
            return (self.r2 * (self.K - i)
                    / (2 * (self.r2 * (self.K - i) + self.r1 * i)))

        # transition matrix
        self.transition_mat = np.zeros((self.nbr_states, self.nbr_states))

        # main matrix elements
        for i in np.arange(1, self.nbr_states - 1):
            self.transition_mat[i - 1, i] = death(i)
            self.transition_mat[i, i] = - (death(i) + birth(i))
            self.transition_mat[i + 1, i] = birth(i)

        # truncate transition matrix to not have singular matrix for inverse
        trunc_trans_mat = np.delete(self.transition_mat, self.abs_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, self.abs_idx, 1)

        # inverse and then make it sparse
        self.inverse = sp.linalg.inv(trunc_trans_mat)
        self.inverse = sparse.csc_matrix(self.inverse)
        self.transition_mat = sparse.csc_matrix(self.transition_mat)

        return

    def force(self, x):
        # s = self.r1 / self.r2
        return 2 * self.K * (2 * x - 1) / (2 * x**2 - 2 * x + 1)

    def potential(self, x):
        return - 2 * self.K * np.log(2 * x**2 - 2 * x + 1)

    def FP_prob(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(0, 1, 101)
        x = [0]
        P = [0]  # np.zeros(len(x))
        return x, P

    def FP_mfpt(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(1.0 / N, 1.0 - 1.0 / N, 101)
        x = [0]
        mfpt = [0]  # np.zeros(len(x))
        return x, mfpt


class OneBoundaryIntFitFPT(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times

        def index(*args):
            # for 2 state system, devised alternate scheme for indexing the
            # states, such that i, j are represented by 1 number
            return args

        # for external use outside of function, finding the index pf the state
        self.index = index

        # total number of states
        self.nbr_states = self.K + 1

        # absorbing points of Moran model
        self.abs_idx = [0, self.nbr_states - 1]
        self.remove_abs = np.ones(self.nbr_states, dtype=bool)
        self.remove_abs[self.abs_idx] = False

        # rate definitions
        def birth(i):
            return (self.r1 * i * (i - 1)
                    / (2 * (self.K - 1) * (self.r2 * (self.K - i)
                                           + self.r1 * i)))

        def death(i):
            return (self.r2 * (self.K - i) * (self.K - i - 1)
                    / (2 * (self.K - 1) * (self.r2 * (self.K - i)
                                           + self.r1 * i)))

        # transition matrix
        self.transition_mat = np.zeros((self.nbr_states, self.nbr_states))

        # main matrix elements
        for i in np.arange(1, self.nbr_states - 1):
            self.transition_mat[i - 1, i] = death(i)
            self.transition_mat[i, i] = - (death(i) + birth(i))
            self.transition_mat[i + 1, i] = birth(i)

        # truncate transition matrix to not have singular matrix for inverse
        trunc_trans_mat = np.delete(self.transition_mat, self.abs_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, self.abs_idx, 1)

        # inverse and then make it sparse
        self.inverse = sp.linalg.inv(trunc_trans_mat)
        self.inverse = sparse.csc_matrix(self.inverse)
        self.transition_mat = sparse.csc_matrix(self.transition_mat)

        return

    def force(self, x):
        return 2 * self.K * (2 * x - 1) / (2 * x**2 - 2 * x + 1)

    def potential(self, x):
        return - 2 * self.K * np.log(2 * x**2 - 2 * x + 1)

    def FP_prob(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(0, 1, 101)
        x = [0]
        P = [0]  # np.zeros(len(x))
        return x, P

    def FP_mfpt(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(1.0 / N, 1.0 - 1.0 / N, 101)
        x = [0]
        mfpt = [0]  # np.zeros(len(x))
        return x, mfpt


class TwoBoundaryFPT(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times

        def index_grow_moran(*args):
            # for 2 state system, devised alternate scheme for indexing the
            # states, such that i, j are represented by 1 number
            return int(args[0] * (self.K + 2) - args[0] * (args[0] + 1) / 2
                       + args[1])

        # for external use outside of function, finding index of a 2 species
        # state.
        self.index = index_grow_moran

        # rates
        def birth_single(i, r):
            return r * i / 2.0

        def death_single(i, r):
            return r * (self.K - i) / 2.0

        def both_left(i, j):
            return self.r1 * j / 2.0

        def both_right(i, j):
            return self.r1 * i / 2.0

        def one_boundary_move(i, j):
            return self.r2 * (self.K - i - j) / 2.0

        # given that i+j <= self.K, this is the total number of states allowed
        self.nbr_states = int((self.K + 1) * (self.K + 2) / 2)

        # define vectors for constructing sparse transition matrix
        # 5 reactions per state at most
        row = np.zeros(5 * self.nbr_states, dtype=int)
        col = np.zeros(5 * self.nbr_states, dtype=int)
        data = np.zeros(5 * self.nbr_states, dtype=float)

        # list of all absorbing states, gets filled below
        self.abs_idx = []

        # Building transition matrix
        rxn_count = 0
        rate_out = 0.0
        print_idx = self.nbr_states + 1
        for i in np.arange(0, self.K + 1):
            for j in np.arange(0, self.K + 1 - i):
                idx = index_grow_moran(i, j)
                if ((i == 0 and j == 0) or (i == 0 and j == self.K) or
                        (i == self.K and j == 0)):
                    (self.abs_idx).append(idx)

                # 1 boundary only, j moving
                if i == 0:
                    if j > 1:  # increase j, don't change i = 0
                        birth_idx = index_grow_moran(i, j - 1)
                        row[rxn_count] = idx
                        col[rxn_count] = birth_idx
                        data[rxn_count] = birth_single(j - 1, self.r1)
                        rate_out += birth_single(j - 1, self.r1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("1", i, j - 1, birth_single(j - 1, self.r1))

                    if j < self.K - 1:  # decrease j, don't change i = 0
                        death_idx = index_grow_moran(i, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = death_idx
                        data[rxn_count] = death_single(j + 1, self.r2)
                        rate_out += death_single(j + 1, self.r2)
                        rxn_count += 1
                        if idx == print_idx:
                            print("2", i, j + 1, death_single(j + 1, self.r2))

                # 1 boundary only, i moving
                if j == 0:
                    if i > 1:  # increase i, don't change j = 0
                        birth_idx = index_grow_moran(i - 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = birth_idx
                        data[rxn_count] = birth_single(i - 1, self.r2)
                        rate_out += birth_single(i - 1, self.r2)
                        rxn_count += 1
                        if idx == print_idx:
                            print("3", i - 1, j, birth_single(i - 1, self.r2))

                    if i < self.K - 1:  # decrease i, don't change j = 0
                        death_idx = index_grow_moran(i + 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = death_idx
                        data[rxn_count] = death_single(i + 1, self.r1)
                        rate_out += death_single(i + 1, self.r1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("4", i + 1, j, death_single(i + 1, self.r1))

                # 2 boundaries moving

                if i + j < self.K:
                    # both boundaries shift to the right
                    if i > 1:
                        shift_r_idx = index_grow_moran(i - 1, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = shift_r_idx
                        data[rxn_count] = both_right(i - 1, j + 1)
                        rate_out += both_right(i - 1, j + 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("5", i - 1, j + 1, both_right(i - 1, j + 1))

                    # both boundaries shift to the left
                    if j > 1:
                        shift_l_idx = index_grow_moran(i + 1, j - 1)
                        row[rxn_count] = idx
                        col[rxn_count] = shift_l_idx
                        data[rxn_count] = both_left(i + 1, j - 1)
                        rate_out += both_left(i + 1, j - 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("6", i + 1, j - 1, both_left(i + 1, j - 1))

                if i + j < self.K - 1:
                    # left boundary moves to the left
                    if i < self.K - 2 and j != 0:
                        grow_l_idx = index_grow_moran(i + 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = grow_l_idx
                        data[rxn_count] = one_boundary_move(i + 1, j)
                        rate_out += one_boundary_move(i + 1, j)
                        rxn_count += 1
                        if idx == print_idx:
                            print("7", i + 1, j, one_boundary_move(i + 1, j))

                    if j < self.K - 2 and i != 0:
                        # right boundary moves to the right
                        grow_r_idx = index_grow_moran(i, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = grow_r_idx
                        data[rxn_count] = one_boundary_move(i, j + 1)
                        rate_out += one_boundary_move(i, j + 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("8", i, j + 1, one_boundary_move(i, j + 1))

                if rate_out == 0.0:
                    (self.abs_idx).append(idx)

                rate_out = 0.0

        temp_transition_mat = sparse.csc_matrix((data, (row, col)),
                                                shape=(self.nbr_states,
                                                       self.nbr_states)
                                                )
        # diagonal elements, leaving state
        diagonal = (temp_transition_mat).sum(axis=0)
        for i, diag in enumerate(np.asarray(diagonal)[0]):
            row[rxn_count] = i
            col[rxn_count] = i
            data[rxn_count] = - diag
            rxn_count += 1

        self.transition_mat = sparse.csc_matrix((data, (row, col)),
                                                shape=(self.nbr_states,
                                                       self.nbr_states)
                                                )

        # states that are absorbing
        self.remove_abs = np.ones(self.nbr_states, dtype=bool)
        self.remove_abs[self.abs_idx] = False

        # remove rows/cols to make the matrix invertible
        trunc_trans_mat = (self.transition_mat[self.remove_abs]
                           )[:, self.remove_abs]
        self.inverse = sLA.inv(trunc_trans_mat)
        print('--- done inverting ---')

        return


class TwoBoundaryIntFPT(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times

        def index_grow_moran(*args):
            # for 2 state system, devised alternate scheme for indexing the
            # states, such that i, j are represented by 1 number
            return int(args[0] * (self.K + 2) - args[0] * (args[0] + 1) / 2
                       + args[1])

        # for external use outside of function, finding index of a 2 species
        # state.
        self.index = index_grow_moran

        # rates
        def birth_single(i, r):
            return r * i * (i - 1) / (2 * (self.K - 1))

        def death_single(i, r):
            return r * (self.K - i) * (self.K - i - 1) / (2 * (self.K - 1))

        def both_left(i, j):
            return self.r1 * j * (j - 1) / (2 * (self.K - 1))

        def both_right(i, j):
            return self.r1 * i * (i - 1) / (2 * (self.K - 1))

        def one_boundary_move(i, j):
            return self.r2 * (self.K - i - j) / 2.0

        def boundary_move_right(i, j):
            return self.r2 * (self.K - i - j) / 2.0

        def boundary_move_left(i, j):
            return self.r2 * (self.K - i - j) / 2.0

        # given that i+j <= self.K, this is the total number of states allowed
        self.nbr_states = int((self.K + 1) * (self.K + 2) / 2)

        # define vectors for constructing sparse transition matrix
        # 5 reactions per state at most
        row = np.zeros(5 * self.nbr_states, dtype=int)
        col = np.zeros(5 * self.nbr_states, dtype=int)
        data = np.zeros(5 * self.nbr_states, dtype=float)

        # list of all absorbing states, gets filled below
        self.abs_idx = []

        # Building transition matrix
        rxn_count = 0
        rate_out = 0.0
        print_idx = self.nbr_states + 1
        for i in np.arange(0, self.K + 1):
            for j in np.arange(0, self.K + 1 - i):
                idx = index_grow_moran(i, j)
                if ((i == 0 and j == 0) or (i == 0 and j == self.K) or
                        (i == self.K and j == 0)):
                    (self.abs_idx).append(idx)

                # 1 boundary only, j moving
                if i == 0:
                    if j > 1:  # increase j, don't change i = 0
                        birth_idx = index_grow_moran(i, j - 1)
                        row[rxn_count] = idx
                        col[rxn_count] = birth_idx
                        data[rxn_count] = birth_single(j - 1, self.r1)
                        rate_out += birth_single(j - 1, self.r1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("1", i, j - 1, birth_single(j - 1, self.r1))

                    if j < self.K - 1:  # decrease j, don't change i = 0
                        death_idx = index_grow_moran(i, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = death_idx
                        data[rxn_count] = death_single(j + 1, self.r2)
                        rate_out += death_single(j + 1, self.r2)
                        rxn_count += 1
                        if idx == print_idx:
                            print("2", i, j + 1, death_single(j + 1, self.r2))

                # 1 boundary only, i moving
                if j == 0:
                    if i > 1:  # increase i, don't change j = 0
                        birth_idx = index_grow_moran(i - 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = birth_idx
                        data[rxn_count] = birth_single(i - 1, self.r2)
                        rate_out += birth_single(i - 1, self.r2)
                        rxn_count += 1
                        if idx == print_idx:
                            print("3", i - 1, j, birth_single(i - 1, self.r2))

                    if i < self.K - 1:  # decrease i, don't change j = 0
                        death_idx = index_grow_moran(i + 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = death_idx
                        data[rxn_count] = death_single(i + 1, self.r1)
                        rate_out += death_single(i + 1, self.r1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("4", i + 1, j, death_single(i + 1, self.r1))

                # 2 boundaries moving

                if i + j < self.K:
                    # both boundaries shift to the right
                    if i > 1:
                        shift_r_idx = index_grow_moran(i - 1, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = shift_r_idx
                        data[rxn_count] = both_right(i - 1, j + 1)
                        rate_out += both_right(i - 1, j + 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("5", i - 1, j + 1, both_right(i - 1, j + 1))

                    # both boundaries shift to the left
                    if j > 1:
                        shift_l_idx = index_grow_moran(i + 1, j - 1)
                        row[rxn_count] = idx
                        col[rxn_count] = shift_l_idx
                        data[rxn_count] = both_left(i + 1, j - 1)
                        rate_out += both_left(i + 1, j - 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("6", i + 1, j - 1, both_left(i + 1, j - 1))

                if i + j < self.K - 1:
                    # left boundary moves to the left
                    if i < self.K - 2 and j != 0:
                        grow_l_idx = index_grow_moran(i + 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = grow_l_idx
                        data[rxn_count] = one_boundary_move(i + 1, j)
                        rate_out += one_boundary_move(i + 1, j)
                        rxn_count += 1
                        if idx == print_idx:
                            print("7", i + 1, j, one_boundary_move(i + 1, j))

                    if j < self.K - 2 and i != 0:
                        # right boundary moves to the right
                        grow_r_idx = index_grow_moran(i, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = grow_r_idx
                        data[rxn_count] = one_boundary_move(i, j + 1)
                        rate_out += one_boundary_move(i, j + 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("8", i, j + 1, one_boundary_move(i, j + 1))

                if rate_out == 0.0:
                    (self.abs_idx).append(idx)

                rate_out = 0.0

        temp_transition_mat = sparse.csc_matrix((data, (row, col)),
                                                shape=(self.nbr_states,
                                                       self.nbr_states)
                                                )
        # diagonal elements, leaving state
        diagonal = (temp_transition_mat).sum(axis=0)
        for i, diag in enumerate(np.asarray(diagonal)[0]):
            row[rxn_count] = i
            col[rxn_count] = i
            data[rxn_count] = - diag
            rxn_count += 1

        self.transition_mat = sparse.csc_matrix((data, (row, col)),
                                                shape=(self.nbr_states,
                                                       self.nbr_states)
                                                )

        # states that are absorbing
        self.remove_abs = np.ones(self.nbr_states, dtype=bool)
        self.remove_abs[self.abs_idx] = False

        # remove rows/cols to make the matrix invertible
        trunc_trans_mat = (self.transition_mat[self.remove_abs]
                           )[:, self.remove_abs]
        self.inverse = sLA.inv(trunc_trans_mat)
        print('--- done inverting ---')

        return
