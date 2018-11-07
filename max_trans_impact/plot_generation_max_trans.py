# -*- coding: utf-8 -*-
__author__ = 'qsong'


import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'axes.labelsize': 30,
}
plt.rcParams.update(params)


DATA_FOLDED = '.'
label_fs = 40

plt.figure(1, figsize=(18, 18))
plt.title("The packet loss rate evolution with SINR threshold = 3dB")

plt.subplot(2, 1, 1)
# The case for threshold 3.0 dB
ana_result_1 = os.path.join(DATA_FOLDED, "analytical_result_K=1_threshold=3.0_l=1_m=1.csv")
ana_result_5 = os.path.join(DATA_FOLDED, "analytical_result_K=5_threshold=3.0_l=1_m=1.csv")
ana_result_10 = os.path.join(DATA_FOLDED, "analytical_result_K=10_threshold=3.0_l=1_m=1.csv")

with open(ana_result_1, 'r') as ana_1_f_handler, \
        open(ana_result_5, 'r') as ana_5_f_handler, \
        open(ana_result_10, 'r') as ana_10_f_handler:

    ana_5_csv = list(csv.reader(ana_5_f_handler))
    ana_5_intensity = [float(element[-1]) for element in ana_5_csv]
    ana_5_plr = [float(element[-2]) for element in ana_5_csv]

    ana_1_csv = list(csv.reader(ana_1_f_handler))
    ana_1_intensity = [float(element[-1]) for element in ana_1_csv]
    ana_1_plr = [float(element[-2]) for element in ana_1_csv]

    ana_10_csv = list(csv.reader(ana_10_f_handler))
    ana_10_intensity = [float(element[-1]) for element in ana_10_csv]
    ana_10_plr = [float(element[-2]) for element in ana_10_csv]


    plt.xlabel("Fresh packet arrival intensity")
    plt.ylabel("Packet loss rate")
    plt.plot(ana_10_intensity, ana_10_plr, 'r-', label="K=9")
    plt.plot(ana_1_intensity, ana_1_plr, 'b-', label="K=0")
    plt.plot(ana_5_intensity, ana_5_plr, 'g-', label="K=4")



    plt.xlabel("Fresh packet arrival intensity")
    plt.ylabel("Packet loss rate")
    plt.legend(loc='best', numpoints=2)
    plt.axis([0, 2.05, 0, 1.05])
    plt.grid()

plt.subplot(2, 1, 2)
plt.grid()

plt.semilogy(ana_10_intensity, ana_10_plr, 'r-', label="K=9")
plt.semilogy(ana_1_intensity, ana_1_plr, 'b-', label="K=0")
plt.semilogy(ana_5_intensity, ana_5_plr, 'g-', label="K=4")



plt.legend(loc='lower right', numpoints=2)
plt.ylabel("Packet loss rate")
plt.xlabel("Fresh packet arrival intensity")

plt.axis([0, 2.05, 0, 1.05])

plt.savefig('max_trans_plr_impact.eps', format='eps', dpi=600)
plt.show()
