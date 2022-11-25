import scipy as sp
import numpy as np
from scipy import sparse
from scipy.integrate import solve_bvp
import scipy.sparse.linalg as sLA
# import time
# from scipy.special import expi


class FirstPassage:
    def __init__(self):
        self.index = None
        self.K = 100
        self.max_pot = 1. / 2.

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
        prob_abs = [prob * (prob > 0) for prob in prob_absorption]
        # calculation from Iyer-Biswas-Zilman Eq. 48
        """
        mfpt = [0 if prob_absorption[i] == 0
                else np.dot(((self.inv_squ @
                             self.transition_mat[state, self.remove_abs].T
                             ).toarray().reshape(-1,),
                            init_prob[self.remove_abs]
                            ) / prob_absorption[i]
                for i, state in enumerate(self.abs_idx)
                ]
        """
        mfpt = []
        for i, state in enumerate(self.abs_idx):
            if prob_absorption[i] == 0:
                mfpt.append(0.0)
            else:
                A = self.inv_squ  # self.inverse @ self.inverse
                B = self.transition_mat[state, self.remove_abs].T
                C = (A @ B).T
                D = C.dot((init_prob[self.remove_abs]))
                mfpt.append((D / prob_absorption[i])[0])
        """
        mfpt = [0 if prob_absorption[i] == 0
                else (((self.inv_squ @
                        self.transition_mat[state, self.remove_abs].T
                       ).T).dot((init_prob[self.remove_abs]))
                       / prob_absorption[i]
                for i, state in enumerate(self.abs_idx)
                ]
        """
        # print('end')
        return prob_abs, mfpt

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

    def B(self, x):
        return 0

    def A(self, x):
        return 0

    def determ_time(self, a, b):
        return - 0

    def potential(self, x):
        # U = - integral ( 2 A(x) / B(x) )
        return 0

    def force(self, x):
        # F = - U' = 2 A(x) / B(x)
        return 2 * self.A(x) / self.B(x)

    def diffusion(self, x):
        # D = - U'' = ( 2 A(x) / B(x) )'
        return 0

    def FP_prob(self, x_find=None, N=None):
        x = np.linspace(0., 1., 101)
        x_mesh = np.linspace(0, 1, 5)

        if N is None:
            N = self.K

        def fun(x, y):
            return np.vstack((y[1], - N / self.K * self.force(x) * y[1]))

        def bc(ya, yb):
            return np.array([ya[0], yb[0] - 1])

        y_guess = np.zeros((2, x_mesh.size))

        res = solve_bvp(fun, bc, x_mesh, y_guess)

        return x, res.sol(x)[0]

    def FP_mfpt(self, x=None, N=None):
        if x is None:
            x = np.linspace(0., 1., 101)
            x_mesh = np.linspace(0, 1, 5)
        else:
            x_mesh = np.linspace(np.min(x), np.max(x), 5)
            r = 10**(-10)
            x_mesh = np.linspace(0+r, 1-r, 10)
            # x_mesh = np.linspace(0, 1, 5)

        if N is None:
            N = self.K

        def fun(x, y):
            return np.vstack((y[1], (((- 2 * N * self.A(x) * y[1] - 2 * N))
                                     / self.B(x))))

        def bc(ya, yb):
            return np.array([ya[0], yb[0]])

        y_guess = np.zeros((2, x_mesh.size))

        res = solve_bvp(fun, bc, x_mesh, y_guess)

        return x, res.sol(x)[0]

    def FP_mfpt_x(self, x_find, N=None):
        x_mesh = np.linspace(0.001, 0.999, 5)

        if N is None:
            N = self.K

            def fun(x, y):
                return np.vstack((y[1], ((- 2 * N * self.A(x) * y[1]
                                          - 2 * N) / self.B(x))))

            def bc(ya, yb):
                return np.array([ya[0], yb[0]])

            y_guess = np.zeros((2, x_mesh.size))

            res = solve_bvp(fun, bc, x_mesh, y_guess)

            return x_find, res.sol(x_find)[0]  # [np.where(x == x_find)]

        else:
            mfpt_N = np.zeros(len(N))
            for i, n in enumerate(N):
                def fun(x, y):
                    return np.vstack((y[1], (- 2 * n * self.A(x) * y[1]
                                             - 2 * n) / self.B(x)))

                def bc(ya, yb):
                    return np.array([ya[0], yb[0]])

                y_guess = np.zeros((2, x_mesh.size))

                res = solve_bvp(fun, bc, x_mesh, y_guess)
                mfpt_N[i] = res.sol(x_find)[0]

            return x_find, mfpt_N


class MoranFPT(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times
        self.s = self.r1 / self.r2
        if self.s >= 1:
            self.xmax = 1. / (np.sqrt(self.s) + 1)
        else:
            self.xmax = 1. / (1 - np.sqrt(self.s))
        self.nmax = int(self.K * self.xmax)

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
            self.transition_mat[i - 1, i] = (self.r2 * i * (self.K - i)
                                             / (self.r2 * (self.K - i) +
                                                self.r1 * i))
            self.transition_mat[i, i] = (- ((self.r2 + self.r1) * i *
                                            (self.K - i) / (self.r2 *
                                                            (self.K - i) +
                                                            self.r1 * i))
                                         )
            self.transition_mat[i + 1, i] = (self.r1 * i * (self.K - i)
                                             / (self.r2 * (self.K - i) +
                                                self.r1 * i))

        # truncate transition matrix to not have singular matrix for inverse
        trunc_trans_mat = np.delete(self.transition_mat, self.abs_idx, 0)
        trunc_trans_mat = np.delete(trunc_trans_mat, self.abs_idx, 1)

        # inverse and then make it sparse
        self.inverse = sp.linalg.inv(trunc_trans_mat)
        self.inverse = sparse.csc_matrix(self.inverse)
        self.inv_squ = ((self.inverse).dot(self.inverse)).T
        self.transition_mat = sparse.csc_matrix(self.transition_mat)

        return

    def FP_prob(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(0, 1, 101)
        if self.s == 1.:
            P = x
        else:
            R = 2 * N * (self.s - 1) / (self.s + 1)
            P = (1. - np.exp(- R * x)) / (1. - np.exp(- R))
        return x, P

    def FP_mfpt(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(1.0 / N, 1.0 - 1.0 / N, 101)
        if self.s == 1.:
            mfpt = - N * ((1 - x) * np.log(1 - x) + x * np.log(x)) / self.r1
            return x, mfpt
        x, mfpt = super().FP_mfpt(x, N)
        return x, mfpt  # / 2

    def FP_mfpt_N(self, x=None, N=None):
        if N is None:
            N = self.K
        if x is None:
            x = np.linspace(1.0 / N, 1.0 - 1.0 / N, 101)
        if self.s == 1.0:
            mfpt = - N * np.log(1 - x) * (1 - x) / x
            return x, mfpt
        raise SystemExit("Don't have dir. escape time for fit. Moran)")
        return x, mfpt

    def A(self, x):
        return ((self.s - 1) * x * (1. - x)) / (1 + x * (self.s - 1))

    def B(self, x):
        return ((self.s + 1) * x * (1. - x)) / (1 + x * (self.s - 1))

    def potential(self, x):
        # U = - integral ( 2 A(x) / B(x) )
        return - 2 * self.K * (self.s - 1) * x / (self.s + 1)

    def force(self, x):
        # U = - U' = 2 A(x) / B(x)
        # s = self.r1 / self.r2
        return 2 * self.K * (self.s - 1) / (self.s + 1)
        # return 2 * self.K * ((s + 1) * x - 1) / ((s - 1) * x + 1)

    def diffusion(self, x):
        # U = - U'' = ( 2 A(x) / B(x) )'
        raise SystemExit("Haven't calculated diffusion term for fit. Moran")
        return 1. / (4 * self.K)


class Spatial(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times
        self.s = growth_rate1 / growth_rate2
        if self.s >= 1:
            self.xmax = 1. / (np.sqrt(self.s) + 1)
        else:
            self.xmax = 1. / (1 - np.sqrt(self.s))
        self.nmax = int(self.K * self.xmax)

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
            return (self.r1 * i * (i - 1) * self.K
                    / (2 * (self.K - 1) * (self.r2 * (self.K - i)
                                           + self.r1 * i)))

        def death(i):
            return (self.r2 * (self.K - i) * (self.K - i - 1) * self.K
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
        self.inv_squ = ((self.inverse).dot(self.inverse)).T
        self.transition_mat = sparse.csc_matrix(self.transition_mat)

        return

    def A(self, x):
        return (((self.s - 1) * x ** 2 + 2 * x - 1)
                / (2 * (1 + x * (self.s - 1))))

    def B(self, x):
        return (((self.s + 1) * x ** 2 - 2 * x + 1)
                / (2 * (1 + x * (self.s - 1))))

    def potential(self, x):
        # U = - integral ( 2 * self.K * A(x) / B(x) )
        # if self.s == 1.0:
        #    return - np.log(2 * x ** 2 - 2 * x + 1) / 2
        return - (2 * self.s * sp.log((self.s + 1.) * x ** 2 - 2. * x + 1.)
                  - 2 * np.arctan(((self.s + 1) * x - 1.) / np.sqrt(self.s))
                  * np.sqrt(self.s) * (self.s - 1.) + (self.s ** 2 - 1.) * x
                  ) * 2 * self.K / (self.s + 1) ** 2
        # return - 2 * self.K * (((self.s + 1) * (self.s - 1)
        #                        - 2 * self.s * np.log((self.s - 1) * x + 1))
        #                       / (self.s - 1)**2)

    def force(self, x):
        # U = - U' = 2 A(x) / B(x)
        # s = self.r1 / self.r2
        return 2 * self.K * self.A(x) / self.B(x)
        # return 2 * self.K * ((s + 1) * x - 1) / ((s - 1) * x + 1)

    def diffusion(self, x):
        # U = - U'' = ( 2 A(x) / B(x) )'
        raise SystemExit("Haven't calculated diffusion term for spatial")
        return 1. / (4 * self.K)


class SpatInv(FirstPassage):

    def __init__(self, growth_rate1, growth_rate2, carrying_capacity, times):
        # 2 boundaries where cells are aligned as r2 - r1 - r2
        self.r1 = growth_rate1
        self.r2 = growth_rate2
        self.K = carrying_capacity
        self.times = times

        def index_grow_moran(*args):
            # for 2 state system, devised alternate scheme for indexing the
            # states, such that i, j are represented by 1 number
            a = args[0]
            b = args[1]
            return int((a + b) * (a + b + 1) / 2 + a)
            # return int(args[0] * (self.K + 2) - args[0] * (args[0] + 1) / 2
            #           + args[1])

        # for external use outside of function, finding index of a 2 species
        # state.
        self.index = index_grow_moran

        # rates
        def birth(pos, r, i, j):
            # r1 growth of i species
            return (r * pos * (pos - 1) * self.K
                    / (2 * (self.K - 1) * (self.r2 * (i + j)
                                           + self.r1 * (self.K - i - j))))

        # given that i+j <= self.K, this is the total number of states allowed
        self.nbr_states = int((self.K + 1) * (self.K + 2) / 2)

        # define vectors for constructing sparse transition matrix
        # 5 reactions per state at most
        row = np.zeros(5 * self.nbr_states, dtype=int)
        col = np.zeros(5 * self.nbr_states, dtype=int)
        data = np.zeros(5 * self.nbr_states, dtype=float)

        # list of all absorbing states, gets filled below
        self.abs_idx = []
        self.rmv_idx = []

        # Building transition matrix
        rxn_count = 0
        rate_out = 0.0
        print_idx = 1E10  # self.index(2, 18)
        for i in np.arange(0, self.K + 1):
            for j in np.arange(0, self.K + 1 - i):
                idx = index_grow_moran(i, j)
                if ((i == 0 and j == 0) or (i == 0 and j == self.K)
                        or (i == self.K and j == 0)):
                    (self.abs_idx).append(idx)
                    (self.rmv_idx).append(idx)

                # 1 boundary only, j moving
                if i == 0:
                    if j > 1:  # increase j, don't change i = 0
                        birth_idx = index_grow_moran(i, j - 1)
                        row[rxn_count] = idx
                        col[rxn_count] = birth_idx
                        data[rxn_count] = birth(j - 1, self.r2, i, j - 1)
                        rate_out += birth(j - 1, self.r2, i, j - 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("1", i, j - 1, birth(j - 1, self.r2, i,
                                                       j - 1))

                    if j < self.K - 1:  # decrease j, don't change i = 0
                        death_idx = index_grow_moran(i, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = death_idx
                        data[rxn_count] = birth(self.K - (j + 1), self.r1, i,
                                                (j + 1))
                        rate_out += birth(self.K - (j + 1), self.r1, i,
                                          (j + 1))
                        rxn_count += 1
                        if idx == print_idx:
                            print("2", i, j + 1, birth(self.K - (j + 1),
                                                       self.r1, i, (j + 1)))

                # 1 boundary only, i moving
                if j == 0:
                    if i > 1:  # increase i, don't change j = 0
                        birth_idx = index_grow_moran(i - 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = birth_idx
                        data[rxn_count] = birth(i - 1, self.r2, i - 1, j)
                        rate_out += birth(i - 1, self.r1, i - 1, j)
                        rxn_count += 1
                        if idx == print_idx:
                            print("3", i - 1, j, birth(i - 1, self.r2,
                                                       i - 1, j))

                    if i < self.K - 1:  # decrease i, don't change j = 0
                        death_idx = index_grow_moran(i + 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = death_idx
                        data[rxn_count] = birth(self.K - (i + 1), self.r1,
                                                i + 1, j)
                        rate_out += birth(self.K - (i + 1), self.r1, i + 1, j)
                        rxn_count += 1
                        if idx == print_idx:
                            print("4", i + 1, j, birth(self.K - (i + 1),
                                                       self.r2, i + 1, j))

                # 2 boundaries moving
                if i + j < self.K:
                    # both boundaries shift to the right
                    if i > 1:
                        shift_r_idx = index_grow_moran(i - 1, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = shift_r_idx
                        data[rxn_count] = birth(i - 1, self.r2, i - 1, j + 1)
                        rate_out += birth(i - 1, self.r2, i - 1, j + 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("5", i - 1, j + 1, birth(i - 1,
                                                           self.r2, i - 1,
                                                           j + 1))

                    # both boundaries shift to the left
                    if j > 1:
                        shift_l_idx = index_grow_moran(i + 1, j - 1)
                        row[rxn_count] = idx
                        col[rxn_count] = shift_l_idx
                        data[rxn_count] = birth(j - 1, self.r2, i + 1, j - 1)
                        rate_out += birth(j - 1, self.r2, i + 1, j - 1)
                        rxn_count += 1
                        if idx == print_idx:
                            print("6", i + 1, j - 1, birth(j - 1, self.r2,
                                                           i + 1, j - 1))

                if i + j < self.K:
                    # left boundary moves to the left
                    if i < self.K - 1 and j != 0:
                        grow_l_idx = index_grow_moran(i + 1, j)
                        row[rxn_count] = idx
                        col[rxn_count] = grow_l_idx
                        data[rxn_count] = (birth(self.K - (i + 1), self.r1,
                                                 i + 1, j)
                                           - birth(j, self.r1, i + 1, j))
                        rate_out += (birth(self.K - (i + 1), self.r1, i + 1, j)
                                     - birth(j, self.r1, i + 1, j))
                        rxn_count += 1
                        if idx == print_idx:
                            print("7", i + 1, j, (birth(self.K - (i + 1),
                                                        self.r1, i + 1, j)
                                                  - birth(j, self.r1, i + 1, j)
                                                  ))

                    if j < self.K - 1 and i != 0:
                        # right boundary moves to the right
                        grow_r_idx = index_grow_moran(i, j + 1)
                        row[rxn_count] = idx
                        col[rxn_count] = grow_r_idx
                        data[rxn_count] = (birth(self.K - (j + 1), self.r1,
                                                 i, j + 1)
                                           - birth(i, self.r1, i, j + 1))
                        rate_out += (birth(self.K - (j + 1), self.r1, i, j + 1)
                                     - birth(i, self.r1, i, j + 1))
                        rxn_count += 1
                        if idx == print_idx:
                            print("8", i, j + 1, (birth(self.K - (j + 1),
                                                        self.r1, i, j + 1)
                                                  - birth(i, self.r1, i, j + 1)
                                                  ))

                if rate_out == 0.0:
                    (self.rmv_idx).append(idx)
                    # print(i, j)

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
        self.remove_abs[self.rmv_idx] = False

        # remove rows/cols to make the matrix invertible
        trunc_trans_mat = (self.transition_mat[self.remove_abs]
                           )[:, self.remove_abs]
        self.inverse = sLA.inv(trunc_trans_mat)
        self.inv_squ = (((self.inverse).toarray()).dot(
            (self.inverse).toarray())).T
        print('--- done inverting ---')

        return
