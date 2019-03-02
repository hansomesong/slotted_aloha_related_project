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
    'legend.fontsize': 20,
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'axes.labelsize': 30,
}
plt.rcParams.update(params)

LOG_DIR = 'log2'
DATA_FOLDED = '.'
FIGSIZE = (18, 18)

A_K_1 = "analytical, Max.Retrans=0"
A_K_2 = "analytical, Max.Retrans=1"
A_K_3 ="analytical, Max.Retrans=2"
A_K_4 ="analytical, Max.Retrans=3"
A_K_5 ="analytical, Max.Retrans=4"
A_K_6 ="analytical, Max.Retrans=5"

S_K_1 = "simulation, Max.Retrans=0"
S_K_2 = "simulation, Max.Retrans=1"
S_K_3 = "simulation, Max.Retrans=2"
S_K_4 = "simulation, Max.Retrans=3"
S_K_5 = "simulation, Max.Retrans=4"
S_K_6 = "simulation, Max.Retrans=5"


X_START = 0.2
X_END = 2.01
X_STEP = 0.1
Y_START = 0.0
Y_END = 1.01
Y_STEP = 0.1

def sim_data_process(sinr_thrsld, k, l, m, sigma_shadowing, backoff):
    SUBSUB_DIR = "backoff_{0}".format(backoff)
    CASE_DIR = 'case_{0}dB'.format(sinr_thrsld)
    POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(l, m, sigma_shadowing)
    K_DIR = "K={0}".format(k)
    logs_dir = glob.glob(os.path.join(CASE_DIR, SUBSUB_DIR, POWER_DIR, K_DIR, "*.csv"))
    sim_intensity = []
    sim_plr =  []
    for csv_file in logs_dir:
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        plr = csv_df.values[:, -2]
        alpha = csv_df.values[:, -1][0]
        sim_intensity.append(alpha)
        # sim_plr.append((np.sum(plr) - np.max(plr) - np.min(plr))/(len(plr)-2.0))
        sim_plr.append(np.mean(plr))
    return sim_intensity, sim_plr

def analytic_data_process(f_name):
    with open(f_name, 'r') as ana_f_handler:
        ana_csv = list(csv.reader(ana_f_handler))
        ana_intensity = [float(element[-1]) for element in ana_csv]
        ana_plr = [float(element[-2]) for element in ana_csv]
    return ana_intensity, ana_plr


sim_intensity_k_1, sim_plr_k_1 = sim_data_process(-3, 1, 1, 1, 1, 36)
sim_intensity_k_2, sim_plr_k_2 = sim_data_process(-3, 2, 1, 1, 1, 36)
sim_intensity_k_3, sim_plr_k_3 = sim_data_process(-3, 3, 1, 1, 1, 36)
sim_intensity_k_4, sim_plr_k_4 = sim_data_process(-3, 4, 1, 1, 1, 36)
sim_intensity_k_5, sim_plr_k_5 = sim_data_process(-3, 5, 1, 1, 1, 36)
sim_intensity_k_6, sim_plr_k_6 = sim_data_process(-3, 6, 1, 1, 1, 36)



ana_result_f_k_1 = os.path.join(DATA_FOLDED, "improved_cf_shadowing_analytical_result_threshold=-3.0dB_l=1_m=1_sigma=1.0_max_trans=1.csv")
ana_result_f_k_2 = os.path.join(DATA_FOLDED, "improved_cf_shadowing_analytical_result_threshold=-3.0dB_l=1_m=1_sigma=1.0_max_trans=2.csv")
ana_result_f_k_3 = os.path.join(DATA_FOLDED, "improved_cf_shadowing_analytical_result_threshold=-3.0dB_l=1_m=1_sigma=1.0_max_trans=3.csv")
ana_result_f_k_4 = os.path.join(DATA_FOLDED, "improved_cf_shadowing_analytical_result_threshold=-3.0dB_l=1_m=1_sigma=1.0_max_trans=4.csv")
ana_result_f_k_5 = os.path.join(DATA_FOLDED, "improved_cf_shadowing_analytical_result_threshold=-3.0dB_l=1_m=1_sigma=1.0_max_trans=5.csv")
ana_result_f_k_6 = os.path.join(DATA_FOLDED, "improved_cf_shadowing_analytical_result_threshold=-3.0dB_l=1_m=1_sigma=1.0_max_trans=6.csv")


ana_intensity_k_1, ana_plr_k_1 = analytic_data_process(ana_result_f_k_1)
ana_intensity_k_2, ana_plr_k_2 = analytic_data_process(ana_result_f_k_2)
ana_intensity_k_3, ana_plr_k_3 = analytic_data_process(ana_result_f_k_3)
ana_intensity_k_4, ana_plr_k_4 = analytic_data_process(ana_result_f_k_4)
ana_intensity_k_5, ana_plr_k_5 = analytic_data_process(ana_result_f_k_5)
ana_intensity_k_6, ana_plr_k_6 = analytic_data_process(ana_result_f_k_6)


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
ax1.plot(ana_intensity_k_1, ana_plr_k_1, color='b',  marker='', linestyle='-', linewidth=2, label=A_K_1)
ax1.plot(ana_intensity_k_2, ana_plr_k_2, color='g', marker='', linestyle='-.', linewidth=2, label=A_K_2)
ax1.plot(ana_intensity_k_3, ana_plr_k_3, color='r', marker='', linestyle='--', linewidth=2, label=A_K_3)
ax1.plot(ana_intensity_k_4, ana_plr_k_4, color='m', marker='', linestyle='-', linewidth=2, label=A_K_4)
ax1.plot(ana_intensity_k_5, ana_plr_k_5, color='k', marker='', linestyle='-', linewidth=2, label=A_K_5)
ax1.plot(ana_intensity_k_6, ana_plr_k_6, color='c', marker='', linestyle='-', linewidth=2, label=A_K_6)


ax1.plot(sim_intensity_k_1, sim_plr_k_1, color='b', marker='*', linestyle="", linewidth=2, label=S_K_1)
ax1.plot(sim_intensity_k_2, sim_plr_k_2, color='g', marker='o', linestyle="", linewidth=2, label=S_K_2)
ax1.plot(sim_intensity_k_3, sim_plr_k_3, color='r', marker='^', linestyle="", linewidth=2, label=S_K_3)
ax1.plot(sim_intensity_k_4, sim_plr_k_4, color='m', marker='D', linestyle="", linewidth=2, label=S_K_4)
ax1.plot(sim_intensity_k_5, sim_plr_k_5, color='k', marker='s', linestyle="", linewidth=2, label=S_K_5)
ax1.plot(sim_intensity_k_6, sim_plr_k_6, color='c', marker='p', linestyle="", linewidth=2, label=S_K_6)

ax1.legend(loc='best', numpoints=2)
ax1.set_yticks(np.arange(Y_START, Y_END, Y_STEP))
ax1.set_xticks(np.arange(X_START, X_END, X_STEP))
ax1.axis([X_START, X_END, Y_START, Y_END])
ax1.set_title("SINR Threshold -3dB, shadowing variance 1dB")
ax1.grid()

# 生成子图2
ax2.semilogy(ana_intensity_k_1, ana_plr_k_1, color='b',  marker='', linestyle='-', linewidth=2, label=A_K_1)
ax2.semilogy(ana_intensity_k_2, ana_plr_k_2, color='g', marker='', linestyle='-.', linewidth=2, label=A_K_2)
ax2.semilogy(ana_intensity_k_3, ana_plr_k_3, color='r', marker='', linestyle='--', linewidth=2, label=A_K_3)
ax2.semilogy(ana_intensity_k_4, ana_plr_k_4, color='m', marker='', linestyle='-', linewidth=2, label=A_K_4)
ax2.semilogy(ana_intensity_k_5, ana_plr_k_5, color='k', marker='', linestyle='-', linewidth=2, label=A_K_5)
ax2.semilogy(ana_intensity_k_6, ana_plr_k_6, color='c', marker='', linestyle='-', linewidth=2, label=A_K_6)

ax2.semilogy(sim_intensity_k_1, sim_plr_k_1, color='b', marker='*', linestyle="", linewidth=2, label=S_K_1)
ax2.semilogy(sim_intensity_k_2, sim_plr_k_2, color='g', marker='o', linestyle="", linewidth=2, label=S_K_2)
ax2.semilogy(sim_intensity_k_3, sim_plr_k_3, color='r', marker='^', linestyle="", linewidth=2, label=S_K_3)
ax2.semilogy(sim_intensity_k_4, sim_plr_k_4, color='m', marker='D', linestyle="", linewidth=2, label=S_K_4)
ax2.semilogy(sim_intensity_k_5, sim_plr_k_5, color='k', marker='s', linestyle="", linewidth=2, label=S_K_5)
ax2.semilogy(sim_intensity_k_6, sim_plr_k_6, color='c', marker='p', linestyle="", linewidth=2, label=S_K_6)


ax2.legend(loc='best', numpoints=2)
ax2.set_xticks(np.arange(X_START, X_END, X_STEP))
ax2.axis([X_START, X_END, Y_START, Y_END])
ax2.axis([X_START, X_END, Y_START, Y_END])
ax2.grid()

plt.savefig('shadowing_case1.eps', format='eps', dpi=300)
plt.show()
