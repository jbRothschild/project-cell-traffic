from src import plotting, default
import os

expNbr = 1
expDir = expDir = default.DATA_DIR + os.sep + f'exp_nbr_{expNbr}'

plotting.gif_experiment(expDir, modulus=20)
