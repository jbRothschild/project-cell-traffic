import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
# from matplotlib.patches import Rectangle
import imageio
# import itertools
from pathlib import Path

BACT_COL = {'0': 'r', '1': 'g', '2': 'b', '3': 'c', '4': 'm', '5': 'y',
            '6': 'k', '7': 'lime'}


def convert_line_2_bacteria(line):
    bacteria = {}
    values = (line.strip()).split(", ")
    bacteria['label'] = values[0]
    bacteria['x'] = float(values[1])
    bacteria['y'] = float(values[2])
    bacteria['angle'] = float(values[3])
    bacteria['length'] = float(values[4])
    bacteria['radius'] = float(values[5])
    bacteria['growth_rate'] = float(values[6])
    bacteria['max_length'] = float(values[7])
    bacteria['split_length'] = float(values[8])
    bacteria['inertia'] = float(values[9])
    bacteria['vel_x'] = float(values[10])
    bacteria['vel_y'] = float(values[11])
    bacteria['vel_angle'] = float(values[12])
    bacteria['acc_x'] = float(values[13])
    bacteria['acc_y'] = float(values[14])
    bacteria['daughter'] = int(values[0][-1])

    edge_length = bacteria['length'] / 2.0 - bacteria['radius']
    p1x = bacteria['x'] + edge_length * np.cos(bacteria['angle'])
    p1y = bacteria['y'] + edge_length * np.sin(bacteria['angle'])
    p2x = bacteria['x'] - edge_length * np.cos(bacteria['angle'])
    p2y = bacteria['y'] - edge_length * np.sin(bacteria['angle'])
    bacteria['p1'] = np.array([p1x, p1y])
    bacteria['p2'] = np.array([p2x, p2y])

    return bacteria


def convert_file_2_environment_dict(params_file):
    environment = {}
    params = open(params_file, 'r')
    names = ((params.readline()).strip()).split(", ")
    values = ((params.readline()).strip()).split(", ")
    for i, parameter in enumerate(names):
        environment[parameter] = float(values[i])
    return environment


def plot_simulation(exp_dir, sim_nbr):
    sim_dir = exp_dir + os.sep + f"sim{sim_nbr}"
    Path(sim_dir).mkdir(parents=True, exist_ok=True)  # make directory

    agents_file = exp_dir + os.sep + f"sim{sim_nbr}.txt"
    params_file = exp_dir + os.sep + "params.txt"
    environment = convert_file_2_environment_dict(params_file)

    agents = open(agents_file, 'r')
    line = agents.readline()

    number_plot = 0
    while line:
        width = environment['CHANNEL_WIDTH']
        height = environment['CHANNEL_HEIGHT']
        resize = 2.0
        plt.figure(figsize=(width / resize, height / resize))
        pixels = 120. / resize
        ax = plt.axes()
        # ax.add_patch(Rectangle((21.0, 5.68), 2, 0.65))
        plt.ylim([0, height])
        plt.xlim([0 - 3.0, width + 3.0])
        while line != "\n":
            bacteria = convert_line_2_bacteria(line)

            plt.plot([bacteria['p1'][0], bacteria['p2'][0]], [bacteria['p1'][1],
                     bacteria['p2'][1]], lw=int(pixels * bacteria['radius']) - 2,
                     solid_capstyle='round',
                     color=BACT_COL[bacteria['label'][0]], zorder=-1,
                     path_effects=[pe.Stroke(linewidth=int(pixels * bacteria['radius']),
                                             foreground='k'), pe.Normal()]
                     )

            line = agents.readline()

        plt.tick_params(axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=True,      # ticks along the bottom edge are off
                        top=True,         # ticks along the top edge are off
                        labelbottom=False)  # labels along the bottom edge are off
        # x.set_frame_on(False)
        right_ax = ax.spines["right"]
        right_ax.set_visible(False)
        left_ax = ax.spines["left"]
        left_ax.set_visible(False)
        plt.vlines(x=0.0, ymin=0, ymax=height, colors='k', ls=':', lw=2)
        plt.vlines(x=width, ymin=0, ymax=height, colors='k', ls=':', lw=2)
        ax.axes.get_yaxis().set_visible(False)  # remove ticks and labels
        save_file = sim_dir + os.sep + f"{number_plot}"
        number_plot += 1
        plt.savefig(save_file + ".png")
        plt.savefig(save_file + ".pdf")
        plt.close()
        line = agents.readline()

    return sim_dir


def plot_simulation_many_init(exp_dir, sim_nbr):
    plt.style.use('seaborn-dark-palette')
    sim_dir = exp_dir + os.sep + f"sim{sim_nbr}"
    Path(sim_dir).mkdir(parents=True, exist_ok=True)  # make directory

    agents_file = exp_dir + os.sep + f"sim{sim_nbr}.txt"
    params_file = exp_dir + os.sep + "params.txt"
    environment = convert_file_2_environment_dict(params_file)

    init_count = open(agents_file, 'r')
    init_line = init_count.readline()
    nbr_species = 0
    while init_line != "\n":
        nbr_species += 1
        # print(init_line)
        init_line = init_count.readline()
    init_count.close()

    agents = open(agents_file, 'r')
    line = agents.readline()

    number_plot = 0
    vals = np.linspace(0, 1, nbr_species)
    np.random.shuffle(vals)
    cmap = plt.cm.colors.ListedColormap(plt.cm.gist_ncar(vals))
    # cmap = plt.cm.gist_ncar
    while line:
        width = environment['CHANNEL_WIDTH']
        height = environment['CHANNEL_HEIGHT']
        resize = 2.0
        pixels = 120. / resize

        plt.figure(figsize=(width / 2., height / 2.))
        ax = plt.axes()
        # ax.set_prop_cycle(plt.rcsetup.cycler('color', cmap())
        # seen = set()
        # colors = list(itertools.takewhile(lambda x: x not in seen and not seen.add(x), (tuple(item['color']) for item in ax._get_lines.prop_cycler)))
        plt.ylim([0, height])
        plt.xlim([0 - 3.0, width + 3.0])
        while line != "\n":
            bacteria = convert_line_2_bacteria(line)
            plt.plot([bacteria['p1'][0], bacteria['p2'][0]], [bacteria['p1'][1],
                     bacteria['p2'][1]], lw=int(pixels * bacteria['radius']) - 2,
                     solid_capstyle='round',
                     # color=colors[int(bacteria['label'][:len(str(nbr_species))])],
                     color=cmap(int(bacteria['label'][:len(str(nbr_species))]) / nbr_species),
                     zorder=-1,
                     path_effects=[pe.Stroke(linewidth=int(pixels * bacteria['radius']),
                                             foreground='k'), pe.Normal()]
                     )

            line = agents.readline()

        plt.tick_params(axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=True,      # ticks along the bottom edge are off
                        top=True,         # ticks along the top edge are off
                        labelbottom=False)  # labels along the bottom edge are off
        # x.set_frame_on(False)
        right_ax = ax.spines["right"]
        right_ax.set_visible(False)
        left_ax = ax.spines["left"]
        left_ax.set_visible(False)
        plt.vlines(x=0.0, ymin=0, ymax=height, colors='k', ls=':', lw=2)
        plt.vlines(x=width, ymin=0, ymax=height, colors='k', ls=':', lw=2)
        ax.axes.get_yaxis().set_visible(False)  # remove ticks and labels
        save_file = sim_dir + os.sep + f"{number_plot}"
        number_plot += 1
        plt.savefig(save_file + ".png")
        plt.savefig(save_file + ".pdf")
        plt.close()
        line = agents.readline()

    return sim_dir


def gif_experiment(dir, modulus=1, fileExt=r'.png'):

    # find all files in dir with file extension fileExt
    list_of_files = [_ for _ in os.listdir(dir) if _.endswith(fileExt)]

    def tryint(s):
        try:
            return int(s)
        except Exception:
            return s

    def alphanum_key(s):
        """ Turn a string into a list of string and number chunks.
            "z23a" -> ["z", 23, "a"]
        """
        return [tryint(c) for c in re.split('([0-9]+)', s)]

    def sort_nicely(list):
        """ Sort the given list in the way that humans expect.
        """
        list.sort(key=alphanum_key)
        return

    sort_nicely(list_of_files)

    # Build GIF
    # could skip every x files?
    with imageio.get_writer(dir + os.sep + 'gifExp.gif', mode='I') as writer:
        for i, filename in enumerate(list_of_files):
            if i % modulus == 0:
                image = imageio.imread(dir + os.sep + filename)
                writer.append_data(image)
    return


if __name__ == '__main__':
    exp_nbr = 0
    sim_nbr = 0
    exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'exp_nbr_{exp_nbr}'
    sim_dir = plot_simulation(exp_dir, sim_nbr)
    gif_experiment(sim_dir)
