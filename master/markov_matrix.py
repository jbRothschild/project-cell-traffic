import numpy as np
from scipy.sparse.linalg import inv as spInv
import scipy as sp

class MasterEquation(N0=1000,N=1000):
    def __init__(self):
        self.min = 0
        self.max = N0
        self.absorbingStateLst = [0, N0]
        self.allAbsorbingStatesLst = self.absorbingStateLst

    def rates(self):
        return 0

    def rate_rcn_left(self, i):

        return 0

    def rate_fcn_right(self, i):

        return 0

    def fpt_state(self, stateNbr):
        if statNbr not in self.allAbsorbingStateLst:
            (self.allAbsorbingStateLst).append(stateNbr)

    def master_equation(self):
        rplus = np.zeros(self.max+2)
        rminus = np.zeros(self.max+2)
        for i, _ in np.arange(0, self.max+2):
            rplus[i]    = self.rate_fcn_right(i)
            rminus[i]   = self.rate_fcn_left(i)

        rplus[-1]  = 0.
        rplus[0]   = 0.
        rminus[0]  = 0.
        rminus[-1] = 0.

        matrix = np.zeros((self.max+2, self.max+2))
        for i in np.arange(1, self.max+1) ):
            right   = np.sum(rplus[:i-1])
            left    = np.sum(rminus[i:])
            matrix[i-1,i] = right
            matrix[i,i]   = left
            matrix[i+1,i] = - ( left + right )

        self.forwardMaster = matrix
        self.backwardMaster = matrix.T

    def remove_absorbing_states(self):
        """
        Removes the states that are the absorbing states of the matrix M
        """
        matrix = np.delete(self.backwardMaster, self.absorbingStateLst, axis=0)
        self.MbRemoved = np.delete(matrix, self.absorbingStateLst, axis=1)

    def inverse(self):
        # Pick one of many ways to do this
        inverseMbRemoved = np.linalg.inv(MbAbsRemoved)
        inverseMbRemoved = spInv(MbAbsRemoved)
        inverseMbRemoved = sp.linalg.inv(MbAbsRemoved)

        self.inverseMbRemoved = inverseMbRemoved

    def vector_rates_into_absorbing_state(self, absorbingState):
        """
        Gives a vector of rates into the absorbing state from other states and removes all absorbing states.
        """
        return np.delete(self.MbRemoved[absorbingState,:],
                                            self.absorbingStateLst)

    def absorption_probability(self, absorbedState):

        rateAbsorbedRemoved = vector_rates_into_absorbing_state(Mb,
                                                    absorbedState, absStateLst)
        probabilityExit = (self.inverseMbRemoved).dot(rateAbssorbedState)

        return probabilityExit

    def mean_first_passage_time(self, absorbedState):
        probabilityExit = self.absorption_probability(absorbedState)

        tau = (self.inverseMbRemoved).dot(probabilityExit)

        mfpt = np.divide(tau, probabilityExit)

        return mfpt, probabilityExit
