# -*- coding: utf-8 -*-
__author__ = 'qsong'

import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import numpy as np
import matplotlib.ticker as ticker
import numpy as np
import scipy.stats as st
import pandas as pd
import glob

params = {
    'legend.fontsize': 17,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'axes.labelsize': 20,
}
plt.rcParams.update(params)

# mpl.rcParams['text.usetex'] = True
# mpl.rcParams.update({'figure.autolayout': True})
FIGSIZE = (21, 18)

A_P_IDENTIC = "analytical, $v=1.0$"
A_P_INCREMENT = "analytical, $v=2.0$"
A_P_DECREMENT ="analytical, $v=0.5$"

S_P_IDENTIC = "simulation, $v=1.0$"
S_P_INCREMENT = "simulation, $v=2.0$"
S_P_DECREMENT = "simulation, $v=0.5$"

MAX_TRANS = 5
LOG_DIR = 'logs'
SUB_DIR = 'shadowing'
ANA_DIR = 'analytical'
SINR_THLD = [3.0, 0.0, -3.0]
BACKOFF = [36, 36, 36]
MU_FADING = 0.0
SIGMA_S = 3.0
LINEWIDTH = 3
METRICS = ['Packet Loss Rate', 'Throughput', 'Energy Efficiency', 'Expected Nb. of Transmissions']
# X and Y Range for perfect case
# X_RANGE = [0.2, 1.2]
# Y_RANGE = [[1e-4, 1.001], [0.1, 1.01], [0, 40], [1, 5]]

# X and Y range for shadowing case
X_RANGE = [0.2, 1.4]
Y_RANGE = [[1e-4, 1.001], [0.1, 1.1], [0, 40], [1, 5]]
FIG_NAME = "{0}_performance.eps".format(SUB_DIR)

def sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff):
    SUBSUB_DIR = "backoff_{0}".format(backoff)
    CASE_DIR = 'case_{0}dB'.format(int(sinr_thrsld))
    POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(l, m, sigma_shadowing)
    logs_dir = glob.glob(os.path.join('.', LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"))
    sim_intensity = []
    sim_plr =  []
    sim_thrpt = []
    sim_ee = []
    sim_avg_nb = []
    for csv_file in logs_dir:
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        plr = csv_df.values[:, -2]
        alpha = csv_df.values[:, -1][0]
        ee_l  = csv_df.values[:, 6]
        prob_vector_df = csv_df.iloc[:, MAX_TRANS+2:-1]
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

        # print POWER_DIR, alpha, plt, alpha*(1-plt), ee, avg_nb
    return sim_intensity, sim_plr, sim_thrpt, sim_ee, sim_avg_nb

def analytic_data_process(f_name, l, m):
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

fig, axes = plt.subplots(3, 4, figsize = FIGSIZE, sharex=False)
# fig, axes = plt.subplots(3, 4, sharex=False)

for case_nb in range(len(SINR_THLD)):
    sinr = int(SINR_THLD[case_nb])
    CASE_DIR = "case_{0}dB".format(sinr)
    ANA_DATA_FOLDED = os.path.join(LOG_DIR, ANA_DIR, SUB_DIR, CASE_DIR)
    sim_no_intensity, sim_no_plr, sim_no_thrpt, sim_no_ee, sim_avg_no_nb = sim_data_process(sinr, 1, 1, int(SIGMA_S), BACKOFF[case_nb])
    sim_no_metrics = [sim_no_plr, sim_no_thrpt, sim_no_ee, sim_avg_no_nb]
    sim_more_intensity, sim_more_plr, sim_more_thrpt, sim_more_ee, sim_avg_more_nb = sim_data_process(sinr, 2, 1, int(SIGMA_S), BACKOFF[case_nb])
    sim_more_metrics = [sim_more_plr, sim_more_thrpt, sim_more_ee, sim_avg_more_nb]
    sim_less_intensity, sim_less_plr, sim_less_thrpt, sim_less_ee, sim_avg_less_nb = sim_data_process(sinr, 1, 2, int(SIGMA_S), BACKOFF[case_nb])
    sim_less_metrics = [sim_less_plr, sim_less_thrpt, sim_less_ee, sim_avg_less_nb]

    ana_result_f_no = os.path.join(
        ANA_DATA_FOLDED,
        "analytical_result_K={0}_threshold={1}dB_l=1_m=1_mufading={2}_sigma={3}.csv".format(
            MAX_TRANS,
            SINR_THLD[case_nb],
            int(MU_FADING),
            int(SIGMA_S)
            )
    )

    ana_result_f_more = os.path.join(
        ANA_DATA_FOLDED,
        "analytical_result_K={0}_threshold={1}dB_l=2_m=1_mufading={2}_sigma={3}.csv".format(
            MAX_TRANS,
            SINR_THLD[case_nb],
            int(MU_FADING),
            int(SIGMA_S)
            )
    )

    ana_result_f_less = os.path.join(
        ANA_DATA_FOLDED,
        "analytical_result_K={0}_threshold={1}dB_l=1_m=2_mufading={2}_sigma={3}.csv".format(
            MAX_TRANS,
            SINR_THLD[case_nb],
            int(MU_FADING),
            int(SIGMA_S)
            )
    )

    ana_no_intensity, ana_no_plr, ana_no_thrpt, ana_no_ee, ana_avg_no_nb = analytic_data_process(ana_result_f_no, 1, 1)
    ana_no_metrics = [ana_no_plr, ana_no_thrpt, ana_no_ee, ana_avg_no_nb]
    ana_more_intensity, ana_more_plr, ana_more_thrpt, ana_more_ee, ana_avg_more_nb = analytic_data_process(ana_result_f_more, 2, 1)
    ana_more_metrics = [ana_more_plr, ana_more_thrpt, ana_more_ee, ana_avg_more_nb]
    ana_less_intensity, ana_less_plr, ana_less_thrpt, ana_less_ee, ana_avg_less_nb = analytic_data_process(ana_result_f_less, 1, 2)
    ana_less_metrics = [ana_less_plr, ana_less_thrpt, ana_less_ee, ana_avg_less_nb]

    for metric_nb in range(len(METRICS)):
        curr_ax = axes[case_nb, metric_nb]

        if metric_nb == 0:
            curr_ax.semilogy(ana_no_intensity, ana_no_metrics[metric_nb], color='b',  marker='', linestyle='-', linewidth=LINEWIDTH, label=A_P_IDENTIC)
            curr_ax.semilogy(ana_more_intensity, ana_more_metrics[metric_nb], color='g', marker='', linestyle='-.', linewidth=LINEWIDTH, label=A_P_INCREMENT)
            curr_ax.semilogy(ana_less_intensity, ana_less_metrics[metric_nb], color='r', marker='', linestyle='--', linewidth=LINEWIDTH, label=A_P_DECREMENT)

            curr_ax.semilogy(sim_no_intensity, sim_no_metrics[metric_nb], color='b', marker='*', linestyle="", linewidth=LINEWIDTH, label=S_P_IDENTIC)
            curr_ax.semilogy(sim_more_intensity, sim_more_metrics[metric_nb], color='g', marker='o', linestyle="", linewidth=LINEWIDTH, label=S_P_INCREMENT)
            curr_ax.semilogy(sim_less_intensity, sim_less_metrics[metric_nb], color='r', marker='^', linestyle="", linewidth=LINEWIDTH, label=S_P_DECREMENT)
        else:

            curr_ax.plot(ana_no_intensity, ana_no_metrics[metric_nb], color='b',  marker='', linestyle='-', linewidth=LINEWIDTH, label=A_P_IDENTIC)
            curr_ax.plot(ana_more_intensity, ana_more_metrics[metric_nb], color='g', marker='', linestyle='-.', linewidth=LINEWIDTH, label=A_P_INCREMENT)
            curr_ax.plot(ana_less_intensity, ana_less_metrics[metric_nb], color='r', marker='', linestyle='--', linewidth=LINEWIDTH, label=A_P_DECREMENT)

            curr_ax.plot(sim_no_intensity, sim_no_metrics[metric_nb], color='b', marker='*', linestyle="", linewidth=LINEWIDTH, label=S_P_IDENTIC)
            curr_ax.plot(sim_more_intensity, sim_more_metrics[metric_nb], color='g', marker='o', linestyle="", linewidth=LINEWIDTH, label=S_P_INCREMENT)
            curr_ax.plot(sim_less_intensity, sim_less_metrics[metric_nb], color='r', marker='^', linestyle="", linewidth=LINEWIDTH, label=S_P_DECREMENT)
        curr_ax.set_xlim(X_RANGE)
        curr_ax.set_ylim(Y_RANGE[metric_nb])
        curr_ax.grid()

rows = ['{}'.format(row) for row in ['SINR Threshold, 3dB', 'SINR Threshold, 0dB', 'SINR Threshold, -3dB']]
for ax, row in zip(axes[:, 0], rows):
    ax.set_ylabel(row, rotation=90, size=15)

cols = ['{}'.format(col) for col in METRICS]
for ax, col in zip(axes[0], cols):
    ax.set_title(col)

plt.legend(bbox_to_anchor=(0.12, 0.01, 0.79, 1), loc=8, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
plt.savefig(os.path.join("figures", FIG_NAME), format='eps', dpi=300, bbox_inches='tight', transparent=True)
plt.show()
