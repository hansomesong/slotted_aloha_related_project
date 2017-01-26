# -*- coding: utf-8 -*-
__author__ = 'qsong'


import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker
import numpy as np
import scipy.stats as st
import scipy.special as special
import pandas as pd
import glob
from statsmodels.distributions.empirical_distribution import ECDF


from shadowing.ALOHA_analytical_model_shadowing_cf import *

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'axes.labelsize': 30,
}
plt.rcParams.update(params)
LOG_DIR = 'logs'
SUB_DIR = 'shadowing'
DATA_FOLDED = '.'
FIGSIZE = (18, 18)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.2
X_END = 2.01
X_STEP = 0.1
Y_START = 0.0
Y_END = 1.01
Y_STEP = 0.1



# 准备就绪，开始画图



fig = plt.figure(1, figsize = FIGSIZE)
ax = fig.add_subplot(111)    # The big subplot
# ax1 = fig.add_subplot(211)
# ax2 = fig.add_subplot(212)

empirical_cdf = "/Users/qsong/Documents/slotted_aloha_related_project/0.6_1062623330.txt"
traces_paths= \
    [
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_338114394.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_1529364351.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_1572178799.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_1620260940.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_1639001807.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_3483984372.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_3605257699.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_3968231334.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_4044882016.txt",
    "/Users/qsong/Documents/slotted_aloha_related_project/0.6_4221010738.txt"
    ]


x_axis = np.linspace(0, 0.005, 1000)

for empirical_cdf in traces_paths:
    with open(empirical_cdf, 'r') as f_handler:
        ana_csv = list(csv.reader(f_handler))
        for csv_line in ana_csv:
            ecdf = ECDF([float(x) for x in csv_line])
            print csv_line
            print ecdf(x_axis)
            ax.plot(x_axis, ecdf(x_axis))


alpah = 0.6
p = 0.002
y = 1 - special.erf(np.pi**2*alpah*p/(4*np.sqrt(x_axis)))
print y

ax.plot(x_axis, y, 'r', marker="*")
# # Turn off axis lines and ticks of the big subplot
# ax.spines['top'].set_color('none')
# ax.spines['bottom'].set_color('none')
# ax.spines['left'].set_color('none')
# ax.spines['right'].set_color('none')
# ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
#
# # Set common labels
# ax.set_xlabel('Fresh Packet Arrival Intensity')
# ax.set_ylabel('Packet Loss Rate')
#
# # 生成子图1
# ax1.plot(ana_no_intensity, ana_no_plr, color='b',  marker='', linestyle='-', linewidth=2, label=A_P_IDENTIC)
# ax1.plot(ana_more_intensity, ana_more_plr, color='g', marker='', linestyle='-.', linewidth=2, label=A_P_INCREMENT)
# ax1.plot(ana_less_intensity, ana_less_plr, color='r', marker='', linestyle='--', linewidth=2, label=A_P_DECREMENT)
# ax1.plot(sim_no_intensity, sim_no_plr, color='b', marker='*', linestyle="", linewidth=2, label=S_P_IDENTIC)
# ax1.plot(sim_more_intensity, sim_more_plr, color='g', marker='o', linestyle="", linewidth=2, label=S_P_INCREMENT)
# ax1.plot(sim_less_intensity, sim_less_plr, color='r', marker='^', linestyle="", linewidth=2, label=S_P_DECREMENT)
# ax1.legend(loc='best', numpoints=2)
# ax1.set_yticks(np.arange(Y_START, Y_END, Y_STEP))
# ax1.set_xticks(np.arange(X_START, X_END, X_STEP))
# ax1.axis([X_START, X_END, Y_START, Y_END])
# ax1.set_title("SINR Threshold -3dB")
# ax1.grid()
#
# # 生成子图2
# ax2.semilogy(ana_no_intensity, ana_no_plr, color='b',  marker='', linestyle='-', linewidth=2, label=A_P_IDENTIC)
# ax2.semilogy(ana_more_intensity, ana_more_plr, color='g', marker='', linestyle='-.', linewidth=2, label=A_P_INCREMENT)
# ax2.semilogy(ana_less_intensity, ana_less_plr, color='r', marker='', linestyle='--', linewidth=2, label=A_P_DECREMENT)
# ax2.semilogy(sim_no_intensity, sim_no_plr, color='b', marker='*', linestyle="", linewidth=2, label=S_P_IDENTIC)
# ax2.semilogy(sim_more_intensity, sim_more_plr, color='g', marker='o', linestyle="", linewidth=2, label=S_P_INCREMENT)
# ax2.semilogy(sim_less_intensity, sim_less_plr, color='r', marker='^', linestyle="", linewidth=2, label=S_P_DECREMENT)
#
# ax2.legend(loc='best', numpoints=2)
# ax2.set_xticks(np.arange(X_START, X_END, X_STEP))
# ax2.axis([X_START, X_END, Y_START, Y_END])
# ax2.axis([X_START, X_END, Y_START, Y_END])
# ax2.grid()
#
# plt.savefig('ecdf.eps', format='eps', dpi=300)
plt.show()
