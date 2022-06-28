import scipy as sp
import numpy as np
import scipy.linalg as LA


class FirstPassage:

    def __init__(self, growth_rate, carrying_capacity, times):
        self.r = growth_rate
        self.K = carrying_capacity
        self.times = times
        self.matrix = []

    def model_moran(self):

        ################# Unique to Moran model ################
        absorption_idx = [0, self.K]
        nbr_states = self.K + 1
        init_state = int(self.K / 2)
        init_prob = np.zeros(nbr_states)
        init_prob[init_state] = 1.0
        transition_mat = np.zeros((nbr_states, nbr_states))

        # main matrix elements
        for i in np.arange(1, nbr_states - 1):
            transition_mat[i - 1, i] = i * (self.K - i)
            transition_mat[i, i] = - 2 * i * (self.K - i)
            transition_mat[i + 1, i] = i * (self.K - i)

        transition_mat *= (self.r * self.K) / (self.K ** 2)

        #######################################################
        trunc_trans_mat = np.delete(transition_mat, absorption_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, absorption_idx, 1)

        inverse = sp.linalg.inv(trunc_trans_mat)

        mfpt = [(np.matmul(inverse, inverse) / inverse)[i, init_state-1] for i in [0, self.K-2]]

        prob_absorption = [-np.dot(np.matmul(np.delete(transition_mat[i, :],
                                                       absorption_idx),
                                             inverse),
                                   np.delete(init_prob, absorption_idx)
                                   ) for i in absorption_idx
                           ]

        fpt_dist_absorp = [[-np.dot(np.matmul(np.delete(transition_mat[i, :], absorption_idx),
                                    LA.expm(trunc_trans_mat * t)),
                                    np.delete(init_prob, absorption_idx)
                                    ) / prob_absorption[idx] for t in self.times
                            ] for idx, state in enumerate(absorption_idx)]
        print(mfpt)
        #print(prob_absorption)
        #print(fpt_dist_absorp)
        return prob_absorption, fpt_dist_absorp[0][:]

    def model_grow_moran(self):

        # nbr_states = (self.K + 1) * (self.K + 1)
        # transition_mat = np.zeros((nbr_states, nbr_states))

        return
