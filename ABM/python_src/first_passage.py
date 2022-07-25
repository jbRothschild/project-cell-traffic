import scipy as sp
import numpy as np
import scipy.linalg as LA
from scipy import sparse
import scipy.sparse.linalg as sLA


class FirstPassage:

    def __init__(self, growth_rate, carrying_capacity, times):
        self.r = growth_rate
        self.K = carrying_capacity
        self.times = times
        self.matrix = []

    def model_moran(self):

        # Unique to Moran model ################
        abs_idx = [0, self.K]
        nbr_states = self.K + 1
        init_state = int(self.K / 2)
        init_prob = np.zeros(nbr_states)
        init_prob[init_state] = 1.0
        transition_mat = np.zeros((nbr_states, nbr_states))
        abs_idx = [0, nbr_states - 1]
        remove_abs = np.ones(nbr_states, dtype=bool)
        remove_abs[abs_idx] = False

        # main matrix elements
        for i in np.arange(1, nbr_states - 1):
            transition_mat[i - 1, i] = i * (self.K - i)
            transition_mat[i, i] = - 2 * i * (self.K - i)
            transition_mat[i + 1, i] = i * (self.K - i)

        transition_mat *= (self.r * self.K) / (self.K ** 2)

        #######################################################
        trunc_trans_mat = np.delete(transition_mat, abs_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, abs_idx, 1)

        inverse = sp.linalg.inv(trunc_trans_mat)

        mfpt = [-(np.matmul(inverse, inverse) / inverse)[i, init_state - 1] for i in [0, self.K - 2]]
        print(mfpt)
        prob_absorption = [-np.dot(np.matmul(np.delete(transition_mat[i, :],
                                                       abs_idx),
                                             inverse),
                                   np.delete(init_prob, abs_idx)
                                   ) for i in abs_idx
                           ]

        fpt_dist_absorp = [[np.dot(np.matmul(np.delete(transition_mat[state, :], abs_idx),
                                             LA.expm(trunc_trans_mat * t)),
                                   np.delete(init_prob, abs_idx)
                                   ) / prob_absorption[idx] for t in self.times
                            ] for idx, state in enumerate(abs_idx)]

        """ remove certain rows in exponentiation... is this right?
        fpt_dist_absorp = [[np.dot(np.matmul(transition_mat[state, :],
                                             LA.expm(transition_mat * t)),
                                   init_prob
                                   ) / prob_absorption[idx] for t in self.times
                            ] for idx, state in enumerate(abs_idx)]
        """
        tot_fpt = [fpt_dist_absorp[0][i] + fpt_dist_absorp[1][i] for i, ele in enumerate(fpt_dist_absorp[1][:])]

        return prob_absorption, fpt_dist_absorp, mfpt, tot_fpt

    def model_grow_moran(self):

        def index(i, j, N):
            return int(i * (N + 2) - i * (i + 1) / 2 + j)

        def birth(i, j):
            return self.r * i * (1. - (i + j) / self.K)

        def death(i, j):
            return self.r * 0  # no deaths happen really.

        def switch(i, j):
            return self.r * i * j / self.K

        init_state_i = 5
        init_state_j = 5
        init_state = index(init_state_i, init_state_j, self.K)

        nbr_states = int((self.K + 1) * (self.K + 2) / 2)
        row = np.zeros(5 * nbr_states, dtype=int)
        col = np.zeros(5 * nbr_states, dtype=int)
        data = np.zeros(5 * nbr_states, dtype=float)

        init_prob = np.zeros(nbr_states)
        init_prob[init_state] = 1.0

        abs_idx = []

        rxn_count = 0
        for i in np.arange(0, self.K + 1):
            for j in np.arange(0, self.K + 1 - i):
                idx = index(i, j, self.K)
                if i == 0 or j == 0:
                    abs_idx.append(idx)
                    # self
                row[rxn_count] = idx
                col[rxn_count] = idx
                data[rxn_count] = - (birth(i, j) + birth(j, i) + 2 * switch(i, j))
                rxn_count += 1

                # birth
                birth_idx = index(i - 1, j, self.K)
                if i > 0:
                    row[rxn_count] = idx
                    col[rxn_count] = birth_idx
                    data[rxn_count] = birth(i - 1, j)
                    rxn_count += 1

                # birth
                birth_idx = index(i, j - 1, self.K)
                if j > 0:
                    row[rxn_count] = idx
                    col[rxn_count] = birth_idx
                    data[rxn_count] = birth(j - 1, i)
                    rxn_count += 1

                # reaction growth i
                rxn_increase_idx = index(i - 1, j + 1, self.K)
                if i > 0 and j < self.K:
                    row[rxn_count] = idx
                    col[rxn_count] = rxn_increase_idx
                    data[rxn_count] = switch(i - 1, j + 1)
                    rxn_count += 1

                # reaction death i
                rxn_decrease_idx = index(i + 1, j - 1, self.K)
                if j > 0 and i < self.K:
                    row[rxn_count] = idx
                    col[rxn_count] = rxn_decrease_idx
                    data[rxn_count] = switch(i + 1, j - 1)
                    rxn_count += 1
        transition_mat = sparse.csc_matrix((data, (row, col)), shape=(nbr_states, nbr_states))
        remove_abs = np.ones(nbr_states, dtype=bool)
        remove_abs[abs_idx] = False
        trunc_trans_mat = (transition_mat[remove_abs])[:, remove_abs]

        # sim_idx = [10, 20, 29, 37, 44, 50, 55, 59, 62, 64, 65]
        # print(transition_mat.A[sim_idx, :][:, sim_idx])

        inverse = sLA.inv(trunc_trans_mat)

        # mfpt = [(np.matmul(inverse, inverse) / inverse)[i, init_state - 1] for i in [0, self.K - 2]]
        mfpt = 0

        prob_absorption = [-np.dot((inverse.T @ transition_mat[i, remove_abs].T).toarray().reshape(-1,),
                           init_prob[remove_abs]
                                   ) for i in abs_idx
                           ]
        """
        fpt_dist_absorp = [np.matmul(sLA.expm_multiply(trunc_trans_mat,
                                                       (transition_mat[state, remove_abs].toarray()).reshape(-1,),
                                                       start=self.times[0],
                                                       stop=self.times[-1],
                                                       num=len(self.times)),
                                     np.delete(init_prob, abs_idx))
                           for idx, state in enumerate(abs_idx)]
                           """
        x = sLA.expm_multiply(transition_mat,
                              (transition_mat[1].toarray()).reshape(-1,),
                              start=self.times[0],
                              stop=self.times[-1],
                              num=120)
        fpt_dist_absorp = [np.matmul(sLA.expm_multiply(transition_mat,
                                                       (transition_mat[state].toarray()).reshape(-1,),
                                                       start=self.times[0],
                                                       stop=self.times[-1],
                                                       num=len(self.times))[:, remove_abs],
                                     np.delete(init_prob, abs_idx))
                           for idx, state in enumerate(abs_idx)]

        # Why not normalized??? seems normalized to 24...
        total_fpt = (np.sum(np.vstack(fpt_dist_absorp), axis=0)).tolist()
        # total_fpt = (np.sum(np.vstack(fpt_dist_absorp), axis=0) / np.sum(np.sum(np.vstack(fpt_dist_absorp), axis=0))).tolist()
        fpt_dist_absorp_N = [(x * 0).tolist() if prob_absorption[i] == 0
                             else (x / prob_absorption[i]).tolist()
                             for i, x in enumerate(fpt_dist_absorp)]
        # transition_mat = np.zeros((nbr_states, nbr_states))

        return prob_absorption, fpt_dist_absorp_N, mfpt, total_fpt
