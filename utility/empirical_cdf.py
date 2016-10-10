# -*- coding: utf-8 -*-
__author__ = 'qsong'


import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker
import numpy as np
import scipy.stats as st
import pandas as pd
import glob

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

empirical_cdf = "model_validation/empirical_cdf_1.12.csv"

ana_file = "/Users/qsong/Documents/slotted_aloha_related_project/imperfect_power_control/case_-3dB/improved_cf_shadowing_analytical_result_threshold=-3.0dB_l=1_m=1_sigma=1.csv"

x_axis = np.arange(0.0, 6.0, 0.02)

with open(empirical_cdf, 'r') as f_handler:
    ana_csv = list(csv.reader(f_handler))
    for k, csv_line in enumerate(ana_csv):
        ax.plot(x_axis, csv_line, label="Empirical CDF for {0}th retransmission".format(k))
        print csv_line


sim_file = "/Users/qsong/Documents/slotted_aloha_related_project/model_validation/simd=50000_warm=5000_maxtrans=5_N=1120_threshold=-3.0dB_l=1_m=1_backoff=36_alpha=1.12_mufading=0.0_mushadowing=0_sigmashadowing=1.0_tmp=20160912143101.csv"
P = []
# with open(sim_file, 'r') as f_handler:
#     sim_csv = list(csv.reader(f_handler))
#     P = np.array([float(x) for x in csv_line[0:-2]])
# print "P vector:", P

alpha_start = 1.12
SIGMA = 1.0
# with open(ana_file, 'r') as f_handler:
#     ana_csv = list(csv.reader(f_handler))
#     for csv_line in ana_csv:
#         if abs(float(csv_line[-1])- alpha_start) < 0.00001:
#             P = np.array([float(x) for x in csv_line[0:-2]])
#             print "P vector:", P

# P = np.array([1.0, 0.80089337048339304, 0.64729556410596611, 0.53005174472601835, 0.44177612666401311, 0.36937773649993366])
# P = np.array([1.0, 0.63265997290620768, 0.43505458602279068, 0.31922862379472466, 0.24539803968443702, 0.19425850665391664])

P = np.array([1.0, 0.57081199850901454, 0.3367861416828517, 0.20414729366527376, 0.1264198693425931])
x_axis = np.arange(0+0.000000000000000000000000001, 6.0, 0.02)
x_axis_2 = np.array([-10.0*np.log10(x) for x in x_axis])
y = [1-trans_failure_p_2(x, P, alpha_start, 1, 0, SIGMA) for x in x_axis_2]
print 1-trans_failure_p_2(-3.0, P, alpha_start, 1, 0, SIGMA)


ax.legend(loc='best', numpoints=2)
plt.plot(x_axis, y)
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
