# -*- coding: utf-8 -*-
__author__ = 'qsong'

import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
import matplotlib.ticker as ticker
import numpy as np
from scipy.stats import bernoulli

import scipy.stats as st
import pandas as pd
import glob
from scipy.special import gamma as gamma_f
import scipy.special as ss
FIGSIZE = (8, 6)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
params = {
    'legend.fontsize': 15,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 20,
    'legend.numpoints': 1
}
plt.rcParams.update(params)

if __name__ == "__main__":

    seed  = 2000
    SLOT_NB = 20
    DEVICE_NB = 10
    np.random.seed(seed)

    binomial_p = 0.1

    devices = np.ones((10, 20))

    pkt_index = bernoulli.rvs(binomial_p, size=(DEVICE_NB, SLOT_NB))

    pkt_start_end = []
    for device_index in range(DEVICE_NB):
        pkt_start_end.append([])
        for slot_index in range(SLOT_NB):
            if pkt_index[device_index][slot_index] == 1:
                start_time = np.random.uniform(slot_index, slot_index+1)
                end_time = start_time + 1
                pkt_start_end[device_index].append((start_time, end_time))

    print pkt_start_end


    x_axis = np.linspace(0, SLOT_NB, 1000)

    y_axis = []

    for x_element in x_axis:
        y_element = 0.0
        for device_index, pkt_time_pairs in enumerate(pkt_start_end):
            for pkt_time_pair in pkt_time_pairs:
                if pkt_time_pair[0] <= x_element <= pkt_time_pair[1]:
                    y_element += 1

        y_axis.append(y_element)

    print "y_axis", y_axis

    fadings = np.random.exponential(scale=1.0, size=(DEVICE_NB, SLOT_NB))
    print pkt_index


    slot_cumu_itf = np.sum(pkt_index, axis=0)


    print slot_cumu_itf

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.plot(x_axis, y_axis)

    X_TICKS_NUM = range(SLOT_NB+1)
    X_TICKS_LABEL = [str(element) for element in X_TICKS_NUM]

    Y_TICKS_MAX = int(np.floor(max(y_axis))+1)
    Y_TICKS_NUM = range(Y_TICKS_MAX+1)
    Y_TICKS_LABEL = [str(element) for element in Y_TICKS_NUM]
    plt.xticks(X_TICKS_NUM, X_TICKS_LABEL)
    plt.yticks(Y_TICKS_NUM, Y_TICKS_LABEL)
    axes.set_xlim([0, 10])
    axes.set_ylim([0, np.floor(max(y_axis))+1])
    axes.grid()
    axes.set_xlabel("Time slot")
    axes.set_ylabel("Cumulative interference level")
    plt.savefig("pure_ALOHA_itf_variation.pdf", format='pdf', dpi=300)




    fig2, axes2 = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)
    Y_TICKS_MAX = int(np.floor(max(slot_cumu_itf))+1)
    Y_TICKS_NUM = range(Y_TICKS_MAX+1)
    Y_TICKS_LABEL = [str(element) for element in Y_TICKS_NUM]
    axes2.step(range(0, SLOT_NB), slot_cumu_itf)
    plt.xticks(X_TICKS_NUM, X_TICKS_LABEL)
    plt.yticks(Y_TICKS_NUM, Y_TICKS_LABEL)
    axes2.set_xlim([0, 10])

    axes2.grid()
    axes2.set_xlabel("Time slot")
    axes2.set_ylabel("Cumulative interference level")
    plt.savefig("sloted_ALOHA_itf_variation.pdf", format='pdf', dpi=300)

    plt.show()



    # print devices
