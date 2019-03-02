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

params = {
    'legend.fontsize': 15,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 30,
}
plt.rcParams.update(params)

LOG_DIR = 'logs'
SUB_DIR = 'shadowing'
DATA_FOLDED = '.'
FIGSIZE = (15, 6)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.3
X_END = 1.01
X_STEP = 0.1
Y_START = 0.0
Y_END = 1.01
Y_STEP = 0.1

MAX_TRANS = 5

def sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff):
    SUBSUB_DIR = "backoff_{0}".format(backoff)
    CASE_DIR = 'case_{0}dB'.format(sinr_thrsld)
    POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(l, m, sigma_shadowing)
    logs_dir = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"))
    sim_intensity = []
    sim_plr =  []
    sim_thrpt = []
    sim_ee = []
    sim_avg_nb = []
    for csv_file in logs_dir:
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        plr = csv_df.values[:, -2]
        # print plr
        alpha = csv_df.values[:, -1][0]
        ee_l  = csv_df.values[:, 6]
        prob_vector_df = csv_df.iloc[:, MAX_TRANS+2:-1]
        # print prob_vector_df
        avg_nb = prob_vector_df.apply(lambda x: sum(x.values[0:-1]), axis=1).mean()
        sim_avg_nb.append(avg_nb)

        sim_intensity.append(alpha)
        # plt = 0.0
        if len(plr) > 10:
            plt = (np.sum(plr) - np.max(plr) - np.min(plr))/(len(plr)-2.0)
            ee = (np.sum(ee_l) - np.max(ee_l) - np.min(ee_l))/(len(ee_l)-2.0)
        else:
            plt = np.sum(plr)/len(plr)
            ee =  np.sum(ee_l)/len(ee_l)

        sim_plr.append(plt)
        sim_thrpt.append(alpha*(1-plt))
        sim_ee.append(ee)

        print POWER_DIR, alpha, plt, alpha*(1-plt), ee, avg_nb
    return sim_intensity, sim_plr, sim_thrpt, sim_ee, sim_avg_nb

def analytic_data_process(f_name, l, m):
    with open(f_name, 'r') as ana_f_handler:
        csv_df = pd.read_csv(f_name, sep=',', header=None)
        # The input csv file has a row format as follows:
        # 1.0,0.9607,0.9229,0.8866,0.8518,0.81830014947865548,0.7000000000000005
        # The last element is the fresh packet arrival intensity
        # The last but one is the packet loss rate
        ana_intensity = csv_df.values[:, -1]
        ana_plr = csv_df.values[:, -2]
        ana_thrpt = csv_df.apply(lambda x: x.values[-1]*(1-x.values[-2]), axis=1)
        power_levels = pd.Series([np.power(l, k)*np.power(m, MAX_TRANS-1-k) for k in range(MAX_TRANS)])
        ana_ee = csv_df.apply(lambda x: (x.iloc[0:MAX_TRANS]*power_levels).sum()/(1-x.values[-2]), axis=1)
        ana_avg_nb = csv_df.apply(lambda x: x.iloc[0:MAX_TRANS].sum(), axis=1)

    return ana_intensity, ana_plr, ana_thrpt, ana_ee, ana_avg_nb


sim_no_intensity, sim_no_plr, sim_no_thrpt, sim_no_ee, sim_avg_no_nb = sim_data_process(3, 1, 1, 1, 36)
sim_more_intensity, sim_more_plr, sim_more_thrpt, sim_more_ee, sim_avg_more_nb = sim_data_process(3, 2, 1, 1, 36)
sim_less_intensity, sim_less_plr, sim_less_thrpt, sim_less_ee, sim_avg_less_nb = sim_data_process(3, 1, 2, 1, 36)

ana_result_f_no = os.path.join(DATA_FOLDED, "shadowing_analytical_result_threshold=3dB_l=1_m=1_sigma=1.0.csv")
ana_result_f_more = os.path.join(DATA_FOLDED, "improved_cf_shadowing_analytical_result_threshold=3.0dB_l=2_m=1_sigma=1.csv")
ana_result_f_less = os.path.join(DATA_FOLDED, "improved_cf_shadowing_analytical_result_threshold=3.0dB_l=1_m=2_sigma=1.csv")

ana_no_intensity, ana_no_plr, ana_no_thrpt, ana_no_ee, ana_avg_no_nb = analytic_data_process(ana_result_f_no, 1, 1)
ana_more_intensity, ana_more_plr, ana_more_thrpt, ana_more_ee, ana_avg_more_nb = analytic_data_process(ana_result_f_more, 2, 1)
ana_less_intensity, ana_less_plr, ana_less_thrpt, ana_less_ee, ana_avg_less_nb = analytic_data_process(ana_result_f_less, 1, 2)

# 准备就绪，开始画图
fig = plt.figure(1, figsize = FIGSIZE)
ax1 = fig.add_subplot(141)
ax2 = fig.add_subplot(142)
ax3 = fig.add_subplot(143)
ax4 = fig.add_subplot(144)

# 生成子图1
ax1.plot(ana_no_intensity, ana_no_thrpt, color='b',  marker='', linestyle='-', linewidth=2, label=A_P_IDENTIC)
ax1.plot(ana_more_intensity, ana_more_thrpt, color='g', marker='', linestyle='-.', linewidth=2, label=A_P_INCREMENT)
ax1.plot(ana_less_intensity, ana_less_thrpt, color='r', marker='', linestyle='--', linewidth=2, label=A_P_DECREMENT)
ax1.plot(sim_no_intensity, sim_no_thrpt, color='b', marker='*', linestyle="", linewidth=2, label=S_P_IDENTIC)
ax1.plot(sim_more_intensity, sim_more_thrpt, color='g', marker='o', linestyle="", linewidth=2, label=S_P_INCREMENT)
ax1.plot(sim_less_intensity, sim_less_thrpt, color='r', marker='^', linestyle="", linewidth=2, label=S_P_DECREMENT)
# ax1.legend(loc='best', numpoints=2)
ax1.set_yticks(np.arange(Y_START, Y_END, Y_STEP))
ax1.set_xticks(np.arange(X_START, X_END, X_STEP))
ax1.axis([X_START, X_END, Y_START, Y_END])
ax1.set_title("Throughput")
ax1.grid()

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
ax2.set_title("Packet loss rate")

ax2.grid()

ax3.plot(ana_no_intensity, ana_no_ee, color='b',  marker='', linestyle='-', linewidth=2, label=A_P_IDENTIC)
ax3.plot(ana_more_intensity, ana_more_ee, color='g', marker='', linestyle='-.', linewidth=2, label=A_P_INCREMENT)
ax3.plot(ana_less_intensity, ana_less_ee, color='r', marker='', linestyle='--', linewidth=2, label=A_P_DECREMENT)
ax3.plot(sim_no_intensity, sim_no_ee, color='b', marker='*', linestyle="", linewidth=2, label=S_P_IDENTIC)
ax3.plot(sim_more_intensity, sim_more_ee, color='g', marker='o', linestyle="", linewidth=2, label=S_P_INCREMENT)
ax3.plot(sim_less_intensity, sim_less_ee, color='r', marker='^', linestyle="", linewidth=2, label=S_P_DECREMENT)
# ax1.legend(loc='best', numpoints=2)
# ax3.set_yticks(np.arange(0, 10, 1))
# ax3.set_xticks(np.arange(X_START, X_END, X_STEP))
ax3.axis([X_START, X_END, 0, 60])
ax3.set_title("Energy efficiency")
ax3.grid()

ax4.plot(ana_no_intensity, ana_avg_no_nb, color='b',  marker='', linestyle='-', linewidth=2, label=A_P_IDENTIC)
ax4.plot(ana_more_intensity, ana_avg_more_nb, color='g', marker='', linestyle='-.', linewidth=2, label=A_P_INCREMENT)
ax4.plot(ana_less_intensity, ana_avg_less_nb, color='r', marker='', linestyle='--', linewidth=2, label=A_P_DECREMENT)
ax4.plot(sim_no_intensity, sim_avg_no_nb, color='b', marker='*', linestyle="", linewidth=2, label=S_P_IDENTIC)
ax4.plot(sim_more_intensity, sim_avg_more_nb, color='g', marker='o', linestyle="", linewidth=2, label=S_P_INCREMENT)
ax4.plot(sim_less_intensity, sim_avg_less_nb, color='r', marker='^', linestyle="", linewidth=2, label=S_P_DECREMENT)
# ax1.legend(loc='best', numpoints=2)
# ax3.set_yticks(np.arange(0, 10, 1))
# ax3.set_xticks(np.arange(X_START, X_END, X_STEP))
ax4.axis([X_START, X_END, 0, 5])
ax4.set_title("Expected number of transmissions")
ax4.grid()

plt.savefig('throughput_shadowing_case1.eps', format='eps', dpi=300)
plt.show()
