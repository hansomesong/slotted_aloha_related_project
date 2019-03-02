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
SUB_DIR = 'perfect'
DATA_FOLDED = '.'
FIGSIZE = (15, 12)

A_P_DECREMENT ="analytical, $v=0.5$"
A_P_IDENTIC = "analytical, $v=1.0$"
A_P_INCREMENT = "analytical, $v=2.0$"
S_P_IDENTIC = "simulation, $v=1.0$"
S_P_INCREMENT = "simulation, $v=2.0$"
S_P_DECREMENT = "simulation, $v=0.5$"

X_START = 0.4
X_END = 1.5
X_STEP = 0.1
Y_START = 0.0
Y_END = 1.01
Y_STEP = 0.1

MAX_TRANS = 5

def cal_avg_trans_nb(row):
    max_trans = len(row) - 1
    tmp = [i+1 for i in range(max_trans)]
    tmp.append(tmp[-1])


def sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff):
    SUBSUB_DIR = "backoff_{0}".format(backoff)
    CASE_DIR = 'case_{0}dB'.format(sinr_thrsld)
    POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(l, m, sigma_shadowing)
    logs_dir = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"))
    sim_intensity = []
    sim_plr =  []
    sim_thrpt = []
    sim_ee = []
    # Average number of transmission trials per packet delivered
    sim_nb = []
    for csv_file in logs_dir:
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        plr = csv_df.values[:, -2]
        alpha = csv_df.values[:, -1][0]
        ee_l  = csv_df.values[:, 6]
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
        print POWER_DIR, alpha, plt, alpha*(1-plt), ee
    return sim_intensity, sim_plr, sim_thrpt, sim_ee

def analytic_data_process(f_name, l, m):
    with open(f_name, 'r') as ana_f_handler:
        ana_csv = list(csv.reader(ana_f_handler))
        ana_intensity = [float(element[-1]) for element in ana_csv]
        ana_plr = [float(element[-2]) for element in ana_csv]
        ana_thrpt = [float(element[-1])*(1-float(element[-2])) for element in ana_csv]
        ana_aa = []

        for line in ana_csv:
            # Do not take the last element in line
            line_length = len(line[0:-2])  #5
            alpha = float(line[-1])
            tmp_plr = float(line[-2]) # packet loss rate
            prob = [float(line[i]) for i in range(line_length)]
            # print alpha, prob, sum([e*l**(k)*m**(4-k) for k, e in enumerate(prob)])
            energies = sum([e*l**(k)*m**(4-k) for k, e in enumerate(prob)])/(1.0-tmp_plr)

            # cumu_powers = [sum([l**(k)*m**(4-k) for k in range(i)]) for i in range(1, 6)]
            # prob = [float(line[i]) - float(line[i+1]) if i < line_length-1 else float(line[i]) for i in range(line_length)]
            # # print alpha, prob, sum([e*l**(k)*m**(4-k) for k, e in enumerate(prob)])
            # energies = sum([e*cumu_powers[k] for k, e in enumerate(prob)])/(1.0-tmp_plr)
            ana_aa.append(energies)

    return ana_intensity, ana_plr, ana_thrpt, ana_aa


sim_no_intensity, sim_no_plr, sim_no_thrpt, sim_no_ee = sim_data_process(0, 1, 1, 0, 150)
sim_more_intensity, sim_more_plr, sim_more_thrpt, sim_more_ee = sim_data_process(0, 2, 1, 0, 150)
sim_less_intensity, sim_less_plr, sim_less_thrpt, sim_less_ee = sim_data_process(0, 1, 2, 0, 150)

ana_result_f_no = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.0_l=1_m=1.csv")
ana_result_f_more = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.0_l=2_m=1.csv")
ana_result_f_less = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.0_l=1_m=2.csv")

ana_no_intensity, ana_no_plr, ana_no_thrpt, ana_no_ee = analytic_data_process(ana_result_f_no, 1, 1)
ana_more_intensity, ana_more_plr, ana_more_thrpt, ana_more_ee = analytic_data_process(ana_result_f_more, 2, 1)
ana_less_intensity, ana_less_plr, ana_less_thrpt, ana_less_ee = analytic_data_process(ana_result_f_less, 1, 2)

# 准备就绪，开始画图
fig = plt.figure(1, figsize = FIGSIZE)
# ax = fig.add_subplot(111)    # The big subplot
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
# ax4 = fig.add_subplot(224)


# Turn off axis lines and ticks of the big subplot
# ax.spines['top'].set_color('none')
# ax.spines['bottom'].set_color('none')
# ax.spines['left'].set_color('none')
# ax.spines['right'].set_color('none')
# ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
#
# # Set common labels
# ax.set_xlabel('Fresh Packet Arrival Intensity')
# ax.set_ylabel('Packet Loss Rate')

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

# 生成子图2
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

print "sim_no_ee", sim_no_ee

plt.savefig('throughput_case2.eps', format='eps', dpi=300)
plt.show()
