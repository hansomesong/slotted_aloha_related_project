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
    'legend.fontsize': 15,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'axes.labelsize': 20,
}
plt.rcParams.update(params)

# mpl.rcParams['text.usetex'] = True
# mpl.rcParams.update({'figure.autolayout': True})
FIGSIZE = (21, 6)

A_P_IDENTIC = "analytical, $v=1.0$"
A_P_INCREMENT = "analytical, $v=2.0$"
A_P_DECREMENT ="analytical, $v=0.5$"

S_P_IDENTIC = "simulation, $v=1.0$,\n$95\%$ confidence \ninterval"
S_P_INCREMENT = "simulation, $v=2.0$,\n$95\%$ confidence \ninterval"
S_P_DECREMENT = "simulation, $v=0.5$,\n$95\%$ confidence \ninterval"

MAX_TRANS = 5
LOG_DIR = 'logs'
SUB_DIR = ''
# SUB_DIR_LIST = ['perfect', 'shadowing']
SUB_DIR_LIST = ['shadowing']

ANA_DIR = 'analytical'
SINR_THLD = [3.0, 0.0, -3.0]
BACKOFF = [36, 36, 36]
MU_FADING = 0.0
# SIGMA_S = [0, 1.0]
SIGMA_S = [1.0]
LINEWIDTH = 2
METRICS = ['Packet Loss Rate']
# X and Y Range for perfect case
# X_RANGE = [0.2, 1.2]
# Y_RANGE = [[1e-4, 1.001], [0.1, 1.01], [0, 40], [1, 5]]

# X and Y range for shadowing case
X_RANGE = [[0.2, 0.8], [0.4, 1.0], [0.6, 1.3]]
X_RANGE_STEP = 0.1
Y_RANGE = [1e-4, 1.0]
FIG_NAME = "packet_loss_rate_ci.eps"


def sim_data_process(max_trans, traces_dir):
    '''
        @:param
        traces_dir: of type string, indicate the path to the folder saving the simulation traces, for example:
                    ./logs/shadowing/case_-3dB/backoff_36/l_1_m_1_sigma_s_1
        @:return
        result: of type list of list, the first element is the list of fresh packet arrival intensity, the second is the
        list packet loss rate, the third is about througphut, the forth is energy efficiency, the last one is the average
        number of transmisions
    '''
    sim_intensity = []
    sim_plr =  []
    sim_thrpt = []
    sim_ee = []
    sim_avg_nb = []
    print traces_dir
    for csv_file in glob.glob(os.path.join(traces_dir, "*.csv")):
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        # One csv file contains a column with identical alpha(i.e., fresh packet arrival rate)
        alpha = csv_df.values[:, -1][0]
        plr_df = csv_df.iloc[:, -2]
        ee_df  = csv_df.iloc[:, 6]
        prob_vector_df = csv_df.iloc[:, max_trans+2:-1]
        nb_df = prob_vector_df.apply(lambda x: sum(x.values[0:-1]), axis=1)

        avg_plr = plr_df.mean()
        # Calculates the standard error of the mean (or standard error of measurement)
        # Here we have to be careful to keep all y values positive:
        min_plr = max(avg_plr-1.96*st.sem(plr_df), 1e-7)
        max_plr = avg_plr + 1.96*st.sem(plr_df)
        print alpha, avg_plr, min_plr, max_plr
        avg_thrpt = alpha*(1-avg_plr)
        min_thrpt = alpha*(1-min_plr)
        max_thrpt = alpha*(1-max_plr)

        avg_ee = ee_df.mean()
        min_ee = max(avg_ee-1.96*st.sem(ee_df), 1e-7)
        max_ee = avg_ee + 1.96*st.sem(ee_df)

        avg_nb = nb_df.mean()
        min_nb = max(avg_nb-1.96*st.sem(nb_df), 1e-7)
        max_nb = avg_nb + 1.96*st.sem(nb_df)

        # Append corresponding list
        sim_intensity.append(alpha)
        sim_plr.append([min_plr, avg_plr, max_plr])
        sim_thrpt.append([min_thrpt, avg_thrpt, max_thrpt])
        sim_ee.append([min_ee, avg_ee, max_ee])
        sim_avg_nb.append([min_nb, avg_nb, max_nb])

    return [sim_intensity, sim_plr, sim_thrpt, sim_ee, sim_avg_nb]

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

fig, axes = plt.subplots(len(SUB_DIR_LIST), len(SINR_THLD), figsize = FIGSIZE, sharey=True)
# fig, axes = plt.subplots(3, 4, sharex=False)
for case_nb in range(len(SUB_DIR_LIST)):
    for sinr_index in range(len(SINR_THLD)):
        sinr = int(SINR_THLD[sinr_index]) # such as 3, 0, -3
        CASE_DIR = "case_{0}dB".format(sinr)
        SUB_DIR = SUB_DIR_LIST[case_nb]
        ANA_DATA_FOLDED = os.path.join(LOG_DIR, ANA_DIR, SUB_DIR, CASE_DIR)
        sigma_s = int(SIGMA_S[case_nb])
        backoff = BACKOFF[sinr_index]
        SUBSUB_DIR = "backoff_{0}".format(backoff)
        CASE_DIR = 'case_{0}dB'.format(sinr)
        POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(1, 1, sigma_s)
        traces_dir = os.path.join('.', LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR)
        sim_no_result = sim_data_process(MAX_TRANS, traces_dir)
        POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(2, 1, sigma_s)
        traces_dir = os.path.join('.', LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR)
        sim_more_result = sim_data_process(MAX_TRANS, traces_dir)
        POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(1, 2, sigma_s)
        traces_dir = os.path.join('.', LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR)
        sim_less_result = sim_data_process(MAX_TRANS, traces_dir)

        ana_result_f_no = os.path.join(
            ANA_DATA_FOLDED,
            "analytical_result_K={0}_threshold={1}dB_l=1_m=1_mufading={2}_sigma={3}.csv".format(
                MAX_TRANS,
                SINR_THLD[sinr_index],
                int(MU_FADING),
                sigma_s
                )
        )

        ana_result_f_more = os.path.join(
            ANA_DATA_FOLDED,
            "analytical_result_K={0}_threshold={1}dB_l=2_m=1_mufading={2}_sigma={3}.csv".format(
                MAX_TRANS,
                SINR_THLD[sinr_index],
                int(MU_FADING),
                sigma_s
                )
        )

        ana_result_f_less = os.path.join(
            ANA_DATA_FOLDED,
            "analytical_result_K={0}_threshold={1}dB_l=1_m=2_mufading={2}_sigma={3}.csv".format(
                MAX_TRANS,
                SINR_THLD[sinr_index],
                int(MU_FADING),
                sigma_s
                )
        )

        ana_no_intensity, ana_no_plr, ana_no_thrpt, ana_no_ee, ana_avg_no_nb = analytic_data_process(ana_result_f_no, 1, 1)
        ana_no_metrics = [ana_no_plr, ana_no_thrpt, ana_no_ee, ana_avg_no_nb]
        ana_more_intensity, ana_more_plr, ana_more_thrpt, ana_more_ee, ana_avg_more_nb = analytic_data_process(ana_result_f_more, 2, 1)
        ana_more_metrics = [ana_more_plr, ana_more_thrpt, ana_more_ee, ana_avg_more_nb]
        ana_less_intensity, ana_less_plr, ana_less_thrpt, ana_less_ee, ana_avg_less_nb = analytic_data_process(ana_result_f_less, 1, 2)
        ana_less_metrics = [ana_less_plr, ana_less_thrpt, ana_less_ee, ana_avg_less_nb]

        sim_no_plr = [element[1] for element in sim_no_result[1]]
        sim_no_plr_lower = [element[1]-element[0] for element in sim_no_result[1]]
        sim_no_plr_upper = [element[2]-element[1] for element in sim_no_result[1]]
        sim_more_plr = [element[1] for element in sim_more_result[1]]
        sim_more_plr_lower = [element[1]-element[0] for element in sim_more_result[1]]
        sim_more_plr_upper = [element[2]-element[1] for element in sim_more_result[1]]
        sim_less_plr = [element[1] for element in sim_less_result[1]]
        sim_less_plr_lower = [element[1]-element[0] for element in sim_less_result[1]]
        sim_less_plr_upper = [element[2]-element[1] for element in sim_less_result[1]]


        curr_ax = axes[sinr_index]
        curr_ax.set_yscale('log') # Y axis is at log scale
        curr_ax.errorbar(sim_no_result[0], sim_no_plr, yerr=[sim_no_plr_lower, sim_no_plr_upper], fmt='*', ecolor='b', capthick=2, label=S_P_IDENTIC)
        curr_ax.errorbar(sim_more_result[0], sim_more_plr, yerr=[sim_more_plr_lower, sim_more_plr_upper], fmt='o', ecolor='g', capthick=2, label=S_P_INCREMENT)
        curr_ax.errorbar(sim_less_result[0], sim_less_plr, yerr=[sim_less_plr_lower, sim_less_plr_upper], fmt='^', ecolor='r', capthick=2, label=S_P_DECREMENT)

        curr_ax.plot(ana_no_intensity, ana_no_metrics[0], color='b',  marker='', linestyle='-', linewidth=LINEWIDTH, label=A_P_IDENTIC)
        curr_ax.plot(ana_more_intensity, ana_more_metrics[0], color='g', marker='', linestyle='-.', linewidth=LINEWIDTH, label=A_P_INCREMENT)
        curr_ax.plot(ana_less_intensity, ana_less_metrics[0], color='r', marker='', linestyle='--', linewidth=LINEWIDTH, label=A_P_DECREMENT)

        print X_RANGE[sinr_index]
        curr_ax.set_xlim(X_RANGE[sinr_index])
        curr_ax.set_xticks(np.arange(X_RANGE[sinr_index][0], X_RANGE[sinr_index][1]+0.0001, X_RANGE_STEP))

        curr_ax.set_ylim(Y_RANGE)
        curr_ax.grid()

# rows = ['{}'.format(row) for row in ["Power Control Error $\sigma=1.0dB$"]]
# for ax, row in zip(axes[0], rows):
axes[0].set_ylabel(r"Packet Loss Rate, $\sigma=1.0$ dB", rotation=90)
#
cols = ['{}'.format(col) for col in ['(a) SINR Threshold, 3dB', '(b) SINR Threshold, 0dB', '(c) SINR Threshold, -3dB']]
for ax, col in zip(axes, cols):
    ax.set_title(col, size=20)
    ax.set_xlabel("Fresh Packet Arrival Intensity")


# Separately place legend into two sub-figures.
handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles[0:3:1], labels[0:3:1], loc='best', numpoints=1, fancybox=True, framealpha=0.5)

handles, labels = axes[1].get_legend_handles_labels()
axes[1].legend(handles[3:4:1], labels[3:4:1], loc='lower right', numpoints=1, fancybox=True, framealpha=0.5)

handles, labels = axes[2].get_legend_handles_labels()
axes[2].legend(handles[4::1], labels[4::1], loc='best', numpoints=1, fancybox=True, framealpha=0.5)

plt.subplots_adjust(wspace=0.06, hspace=0.02)
# plt.legend(bbox_to_anchor=(0.12, -0.04, 0.79, 1), loc=8, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
plt.savefig(os.path.join("figures", FIG_NAME), format='eps', dpi=300, bbox_inches='tight', transparent=True)
# plt.show()
