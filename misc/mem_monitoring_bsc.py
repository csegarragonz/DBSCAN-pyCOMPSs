#!/usr/bin/env python3

import sys
import time, random
import math
from collections import deque

start = time.time()

class RealtimePlot:
    def __init__(self, axes, max_entries = 1000000):
        self.axis_x = deque(maxlen=max_entries)
        self.axis_y = deque(maxlen=max_entries)
        self.axes = axes

        self.max_entries = max_entries
        self.lineplot, = axes.plot([], [], "ro-")
        self.axes.set_autoscaley_on(False)
        self.axes.set_ylim(0,25)

    def add(self, x, y):
        self.axis_x.append(x)
        self.axis_y.append(y)
        self.lineplot.set_data(self.axis_x, self.axis_y)
        self.axes.set_xlim(self.axis_x[0], self.axis_x[-1] + 1e-15)
        # self.axes.relim(); self.axes.autoscale_view() # rescale the y-axis

    def animate(self, figure, callback, interval = 200):
        import matplotlib.animation as animation
        def wrapper(frame_index):
            self.add(*callback(frame_index))
            self.axes.relim(); self.axes.autoscale_view() # rescale the y-axis
            return self.lineplot
        animation.FuncAnimation(figure, wrapper, interval=interval)

def main(num_nodes):
    from matplotlib import pyplot as plt

    figs = [[] for i in range(num_nodes)]
    axes = [[] for i in range(num_nodes)]
    displays = [[] for i in range(num_nodes)]
    for i in range(num_nodes):
        figs[i], axes[i] = plt.subplots()
        displays[i] = RealtimePlot(axes[i])
        displays[i].animate(figs[i], lambda frame_index: (time.time() - start, random.random() * 25))

    plt.show()

    for i in range(num_nodes):
        figs[i], axes[i] = plt.subplots()
        displays[i] = RealtimePlot(axes[i])

    import fileinput
    filename = 'out_no_color.txt'
    while True:
        # Command to obtain the output without color
        # bsc_load 1376744 | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g" > out_no_color.txt
        # Process
#        with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
#            for line in file:
#                line = line.replace("\t\t", "\t")
#                line = line.replace("\t", " ")
#                line = line.replace(" (%)", "")
#                print(line.replace(" (MB)", ""), end='')
#
#        # Load it
#        df = pd.read_csv("out_no_color.txt", sep=' ', skiprows=[0,1,2,4], index_col=0)
        # Access a value using
        # df.loc['NODE']['param']
        for i in range(num_nodes):
            displays[i].add(time.time() - start, random.random() * 25)
        plt.pause(1)

if __name__ == "__main__": main(int(sys.argv[1]))
