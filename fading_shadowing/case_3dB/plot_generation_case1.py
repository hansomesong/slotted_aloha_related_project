# -*- coding: utf-8 -*-
__author__ = 'qsong'

#
# import csv
# import matplotlib.pyplot as plt
# import os
# import numpy as np
# import matplotlib.ticker as ticker
# import numpy as np
# import scipy.stats as st
# import pandas as pd
# import glob

# LOG_DIR = 'logs'
# SUB_DIR = 'fading_shadowing'
# SUBSUB_DIR = "backoff_36"
# CASE_DIR = 'case_3dB'
# NO_DIR = "l_1_m_1_sigma_s_1"
# MORE_DIR = "l_2_m_1_sigma_s_1"
# LESS_DIR = "l_1_m_2_sigma_s_1"
#
# no_all_logs = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, NO_DIR, "*.csv"))
# sim_no_intensity = []
# sim_no_plr = []
# for csv_file in no_all_logs:
#     csv_df = pd.read_csv(csv_file, sep=',', header=None)
#     plr = csv_df.values[:, -2]
#     alpha = csv_df.values[:, -1][0]
#     sim_no_intensity.append(alpha)
#     sim_no_plr.append(np.mean(plr))
#
# more_all_logs = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, MORE_DIR, "*.csv"))
# sim_more_intensity = []
# sim_more_plr = []
# for csv_file in more_all_logs:
#     csv_df = pd.read_csv(csv_file, sep=',', header=None)
#     plr = csv_df.values[:, -2]
#     alpha = csv_df.values[:, -1][0]
#     sim_more_intensity.append(alpha)
#     sim_more_plr.append(np.mean(plr))
#
# less_all_logs = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, LESS_DIR, "*.csv"))
# sim_less_intensity = []
# sim_less_plr = []
# for csv_file in less_all_logs:
#     csv_df = pd.read_csv(csv_file, sep=',', header=None)
#     plr = csv_df.values[:, -2]
#     alpha = csv_df.values[:, -1][0]
#     sim_less_intensity.append(alpha)
#     sim_less_plr.append(np.mean(plr))
#
# params = {
#     'legend.fontsize': 20,
#     'xtick.labelsize': 'large',
#     'ytick.labelsize': 'large',
#     'axes.labelsize': 30,
# }
# plt.rcParams.update(params)
#
# DATA_FOLDED = '.'
# plt.figure(1, figsize=(18, 18))
# plt.title("The packet loss rate evolution with SINR threshold = 3dB")
#
# plt.subplot(2, 1, 1)
# # The case for threshold 3.0 dB
# ana_result_f_no = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=3dB_l=1_m=1_sigma=1.csv")
# ana_result_f_more = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=3dB_l=2_m=1_sigma=1.csv")
# ana_result_f_less = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=3dB_l=1_m=2_sigma=1.csv")
#
# with open(ana_result_f_no, 'r') as ana_no_f_handler, \
#         open(ana_result_f_more, 'r') as ana_more_f_handler, \
#         open(ana_result_f_less, 'r') as ana_less_f_handler:
#
#     ana_more_csv = list(csv.reader(ana_more_f_handler))
#     ana_more_intensity = [float(element[-1]) for element in ana_more_csv]
#     ana_more_plr = [float(element[-2]) for element in ana_more_csv]
#
#     ana_no_csv = list(csv.reader(ana_no_f_handler))
#     ana_no_intensity = [float(element[-1]) for element in ana_no_csv]
#     ana_no_plr = [float(element[-2]) for element in ana_no_csv]
#
#     ana_less_csv = list(csv.reader(ana_less_f_handler))
#     ana_less_intensity = [float(element[-1]) for element in ana_less_csv]
#     ana_less_plr = [float(element[-2]) for element in ana_less_csv]
#
#
#
#
#     plt.xlabel("The fresh packet arrival intensity")
#     plt.ylabel("The packet loss rate")
#     plt.plot(ana_less_intensity, ana_less_plr, 'r-', label="analytical result with power decrement")
#     plt.plot(ana_no_intensity, ana_no_plr, 'b-', label="analytical result case with identical power")
#     plt.plot(ana_more_intensity, ana_more_plr, 'g-', label="analytical result with power increment")
#     plt.plot(sim_no_intensity, sim_no_plr, 'b*', label="simulation result with identical power")
#     plt.plot(sim_more_intensity, sim_more_plr, 'go', label="simulation result with power increment")
#     plt.plot(sim_less_intensity, sim_less_plr, 'r+', label="simulation result with power decrement")
#
#     plt.xlabel("The fresh packe arrival intensity")
#     plt.ylabel("The packet loss rate")
#     plt.legend(loc='best', numpoints=2)
#     plt.axis([0, 2.05, 0, 1.05])
#     plt.grid()
#
# plt.subplot(2, 1, 2)
# plt.grid()
#
# plt.semilogy(ana_less_intensity, ana_less_plr, 'r-', label="analytical result with power decrement")
# plt.semilogy(ana_no_intensity, ana_no_plr, 'b-', label="analytical result case with identical power")
# plt.semilogy(ana_more_intensity, ana_more_plr, 'g-', label="analytical result with power increment")
# plt.semilogy(sim_no_intensity, sim_no_plr, 'b*', label="simulation result with  identical power")
# plt.semilogy(sim_more_intensity, sim_more_plr, 'go', label="simulation result with power increment")
# plt.semilogy(sim_less_intensity, sim_less_plr, 'r+', label="simulation result with power decrement")
#
#
# plt.legend(loc='best', numpoints=2)
# plt.ylabel("The packet loss rate")
# plt.xlabel("The fresh packe arrival intensity")
#
# plt.axis([0, 2.05, 0, 1.05])
#
#
# plt.savefig('fading_case1.eps', format='eps', dpi=300)
# plt.show()


import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.stats as st
import pandas as pd
import glob

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'axes.labelsize': 30,
}
plt.rcParams.update(params)

LOG_DIR = 'logs'
SUB_DIR = 'fading_shadowing'
DATA_FOLDED = '.'
FIGSIZE = (18, 18)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.2
X_END = 1.01
X_STEP = 0.1
Y_START = 0.0
Y_END = 1.01
Y_STEP = 0.1

def sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff):
    SUBSUB_DIR = "backoff_{0}".format(backoff)
    CASE_DIR = 'case_{0}dB'.format(sinr_thrsld)
    POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(l, m, sigma_shadowing)
    logs_dir = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"))
    sim_intensity = []
    sim_plr =  []
    for csv_file in logs_dir:
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        plr = csv_df.values[:, -2]
        alpha = csv_df.values[:, -1][0]
        sim_intensity.append(alpha)
        sim_plr.append((np.sum(plr) - np.max(plr) - np.min(plr))/(len(plr)-2.0))
    return sim_intensity, sim_plr

def analytic_data_process(f_name):
    with open(f_name, 'r') as ana_f_handler:
        ana_csv = list(csv.reader(ana_f_handler))
        ana_intensity = [float(element[-1]) for element in ana_csv]
        ana_plr = [float(element[-2]) for element in ana_csv]
    return ana_intensity, ana_plr


sim_no_intensity, sim_no_plr = sim_data_process(3, 1, 1, 1, 50)
sim_more_intensity, sim_more_plr = sim_data_process(3, 2, 1, 1, 50)
sim_less_intensity, sim_less_plr = sim_data_process(3, 1, 2, 1, 50)

ana_result_f_no = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=3dB_l=1_m=1_sigma=1.csv")
ana_result_f_more = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=3dB_l=2_m=1_sigma=1.csv")
ana_result_f_less = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=3dB_l=1_m=2_sigma=1.csv")

ana_no_intensity, ana_no_plr = analytic_data_process(ana_result_f_no)
ana_more_intensity, ana_more_plr = analytic_data_process(ana_result_f_more)
ana_less_intensity, ana_less_plr = analytic_data_process(ana_result_f_less)

# 准备就绪，开始画图
fig = plt.figure(1, figsize = FIGSIZE)
ax = fig.add_subplot(111)    # The big subplot
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

# Set common labels
ax.set_xlabel('Fresh Packet Arrival Intensity')
ax.set_ylabel('Packet Loss Rate')

# 生成子图1
ax1.plot(ana_no_intensity, ana_no_plr, color='b',  marker='', linestyle='-', linewidth=2, label=A_P_IDENTIC)
ax1.plot(ana_more_intensity, ana_more_plr, color='g', marker='', linestyle='-.', linewidth=2, label=A_P_INCREMENT)
ax1.plot(ana_less_intensity, ana_less_plr, color='r', marker='', linestyle='--', linewidth=2, label=A_P_DECREMENT)
ax1.plot(sim_no_intensity, sim_no_plr, color='b', marker='*', linestyle="", linewidth=2, label=S_P_IDENTIC)
ax1.plot(sim_more_intensity, sim_more_plr, color='g', marker='o', linestyle="", linewidth=2, label=S_P_INCREMENT)
ax1.plot(sim_less_intensity, sim_less_plr, color='r', marker='^', linestyle="", linewidth=2, label=S_P_DECREMENT)
ax1.legend(loc='best', numpoints=2)
ax1.set_yticks(np.arange(Y_START, Y_END, Y_STEP))
ax1.set_xticks(np.arange(X_START, X_END, X_STEP))
ax1.axis([X_START, X_END, Y_START, Y_END])
ax1.set_title("SINR Threshold 3dB")
ax1.grid()

# 生成子图2
ax2.semilogy(ana_no_intensity, ana_no_plr, color='b',  marker='', linestyle='-', linewidth=2, label=A_P_IDENTIC)
ax2.semilogy(ana_more_intensity, ana_more_plr, color='g', marker='', linestyle='-.', linewidth=2, label=A_P_INCREMENT)
ax2.semilogy(ana_less_intensity, ana_less_plr, color='r', marker='', linestyle='--', linewidth=2, label=A_P_DECREMENT)
ax2.semilogy(sim_no_intensity, sim_no_plr, color='b', marker='*', linestyle="", linewidth=2, label=S_P_IDENTIC)
ax2.semilogy(sim_more_intensity, sim_more_plr, color='g', marker='o', linestyle="", linewidth=2, label=S_P_INCREMENT)
ax2.semilogy(sim_less_intensity, sim_less_plr, color='r', marker='^', linestyle="", linewidth=2, label=S_P_DECREMENT)

ax2.legend(loc='best', numpoints=2)
ax2.set_xticks(np.arange(X_START, X_END, X_STEP))
ax2.axis([X_START, X_END, Y_START, Y_END])
ax2.axis([X_START, X_END, Y_START, Y_END])
ax2.grid()

plt.savefig('fading_case1.eps', format='eps', dpi=300)
plt.show()
