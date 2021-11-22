import os, re
import numpy as np
import matplotlib.pyplot as plt
import imageio

def gif_experiment(dir, modulus=1, fileExt=r'.png'):

    #nfind all files in dir with file extension fileExt
    list_of_files = [_ for _ in os.listdir( dir ) if _.endswith(fileExt)]

    def tryint(s):
        try:
            return int(s)
        except:
            return s

    def alphanum_key(s):
        """ Turn a string into a list of string and number chunks.
            "z23a" -> ["z", 23, "a"]
        """
        return [ tryint(c) for c in re.split('([0-9]+)', s) ]

    def sort_nicely(l):
        """ Sort the given list in the way that humans expect.
        """
        l.sort(key=alphanum_key)
        return 0

    sort_nicely(list_of_files)

    # Build GIF
    # could skip every x files?
    with imageio.get_writer(dir + os.sep +'gifExp.gif', mode='I') as writer:
        for i, filename in enumerate( list_of_files ):
            if i % modulus == 0:
                image = imageio.imread(dir + os.sep + filename)
                writer.append_data(image)
    return 0

if __name__ == '__main__':
    simNbr = 1
    simDir = os.getcwd() + os.sep + 'data' + os.sep + f'exp_nbr_{simNbr}'
    gif_experiment(simDir)
