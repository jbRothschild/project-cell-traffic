import scipy as sp
import numpy as np
import scipy.linalg as LA
from scipy import sparse
import scipy.sparse.linalg as sLA


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
                                     self.transition_mat[state, self.remove_abs]
                                     ).toarray().reshape(-1,))
                           for idx, state in enumerate(self.abs_idx)]

        total_fpt = (np.sum(np.vstack(fpt_dist_absorp), axis=0)).tolist()
        fpt_dist_absorp_N = [(x * 0).tolist() if prob_absorption[i] == 0
                             else (x / prob_absorption[i]).tolist()
                             for i, x in enumerate(fpt_dist_absorp)]

        return fpt_dist_absorp_N, total_fpt


class MoranFPT(FirstPassage):

    def __init__(self, growth_rate, carrying_capacity, times):
        self.r = growth_rate
        self.K = carrying_capacity
        self.times = times

        def index_moran(arg):
            #print(args)
            # simple in Moran model, but is more complicated in general
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
            self.transition_mat[i - 1, i] = i * (self.K - i)
            self.transition_mat[i, i] = - 2 * i * (self.K - i)
            self.transition_mat[i + 1, i] = i * (self.K - i)
        self.transition_mat *= (self.r * self.K) / (self.K ** 2)

        # truncate transition matrix to not have singular matrix for inverse
        trunc_trans_mat = np.delete(self.transition_mat, self.abs_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, self.abs_idx, 1)

        # inverse and then make it sparse
        self.inverse = sp.linalg.inv(trunc_trans_mat)
        self.inverse = sparse.csc_matrix(self.inverse)
        self.transition_mat = sparse.csc_matrix(self.transition_mat)

        return


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

        # for external use outside of function, finding the index of a 2 species
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
                data[rxn_count] = - birth(i, j) - birth(j, i) - 2 * switch(i, j)
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
