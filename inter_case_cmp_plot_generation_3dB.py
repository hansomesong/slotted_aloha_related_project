# -*- coding: utf-8 -*-
__author__ = 'qsong'
# In this script, we just plot the analytical curve,
# since the effectiveness of our analytical model has been proved
# by a figure where both analytical and simulation (95% confidence interval) results are drawn in the same figure.

# When this script

import matplotlib.pyplot as plt
import os
import numpy as np
import numpy as np
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
FIGSIZE = (21, 6.5)

A_P_IDENTIC = "$v=1.0$, $\sigma=1.0$"
A_P_INCREMENT = "$v=2.0$, $\sigma=1.0$"
A_P_DECREMENT ="$v=0.5$, $\sigma=1.0$"

P_A_P_IDENTIC = "$v=1.0$, \n$\sigma=0.0$"
P_A_P_INCREMENT = "$v=2.0$, \n$\sigma=0.0$"
P_A_P_DECREMENT = "$v=0.5$, \n$\sigma=0.0$"

MAX_TRANS = 5
LOG_DIR = 'logs'
SUB_DIR = 'shadowing'
ANA_DIR = 'analytical'
SINR_THLD = [3.0]

SETTING_0 ={

    'MAX_TRANS': 5,
    'LOG_DIR':  'logs',
    'SUB_DIR':  'shadowing',
    'ANA_DIR':  'analytical',
    'SINR_THLD': 3.0,
    'BACKOFF': 36,
    'MU_FADING': 0.0,
    'SIGMA_S':  0.0,
    'LINEWIDTH': 2,
    'LM_PAIR':  [(1, 1), (2, 1), (1, 2)],
    'LINE_COLOR': ['b', 'g', 'r'],
    'LINE_TYPE': ['-', '-.', '--'],
    'MARKER_T': ['', '', ''],
    'MARKER_EVERY': [20, 20, 20],
    'LABEL':  ["$v=1.0$, $\sigma=0.0$", "$v=2.0$, $\sigma=0.0$", "$v=0.5$, $\sigma=0.0$"],
}

SETTING_1 ={

    'MAX_TRANS': 5,
    'LOG_DIR':  'logs',
    'SUB_DIR':  'shadowing',
    'ANA_DIR':  'analytical',
    'SINR_THLD': 3.0,
    'BACKOFF': 36,
    'MU_FADING': 0.0,
    'SIGMA_S':  1.0,
    'LINEWIDTH': 2,
    'LM_PAIR':  [(1, 1), (2, 1), (1, 2)],
    'LINE_COLOR': ['b', 'g', 'r'],
    'LINE_TYPE': ['-', '-.', '--'],
    'MARKER_T': ['*', 'o', '^'],
    'MARKER_EVERY': [20, 10, 23],
    'LABEL':  ["$v=1.0$, $\sigma=1.0$", "$v=2.0$, $\sigma=1.0$", "$v=0.5$, $\sigma=1.0$"],
}

SETTING_3 ={
    'MAX_TRANS': 5,
    'LOG_DIR':  'logs',
    'SUB_DIR':  'shadowing',
    'ANA_DIR':  'analytical',
    'SINR_THLD': 1.0,
    'BACKOFF': 36,
    'MU_FADING': 0.0,
    'SIGMA_S':  3.0,
    'LINEWIDTH': 2,
    'LM_PAIR':  [(1, 1), (2, 1), (1, 2)],
    'LINE_COLOR': ['k', 'c', 'm'],
    'LINE_TYPE': ['-', '-.', '--'],
    'MARKER_T': ['s', 'D', 'v'],
    'MARKER_EVERY': [20, 10, 23],
    'LABEL':  ["$v=1.0$, $\sigma=3.0$", "$v=2.0$, $\sigma=3.0$", "$v=0.5$, $\sigma=3.0$"],
}
#
SETTING_0['SINR_THLD'] = SINR_THLD[0]
SETTING_1['SINR_THLD'] = SINR_THLD[0]
SETTING_3['SINR_THLD'] = SINR_THLD[0]

SETTING = [SETTING_0, SETTING_1, SETTING_3]

MU_FADING = 0.0
SIGMA_S = 1.0
LINEWIDTH = 2
LM_Pair = [(1, 1), (2, 1), (1, 2)]


METRICS = ['Packet Loss Rate', 'Throughput', 'Energy Efficiency']
Y_SCALE = ['log', 'linear', 'linear']
SUBFIG_ENUM = ['a', 'b', 'c']

# METRICS = ['Packet Loss Rate', 'Throughput', 'Energy Efficiency', 'Expected Nb. of Transmissions']
# Y_SCALE = ['log', 'linear', 'linear', 'linear']
# SUBFIG_ENUM = ['a', 'b', 'c', 'd']


# X and Y Range for capture ratio 0dB
X_RANGE = [[0.3, 0.65], [0.3, 0.8], [0.3, 0.8]]
Y_RANGE = [[1e-3, 1.001], [0.25, 0.6], [0, 0.8], [1, 5]]
FIG_NAME = "{0}_performance_case_{1}.eps".format(SUB_DIR, SINR_THLD[0])

def sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff):
    SUBSUB_DIR = "backoff_{0}".format(backoff)
    CASE_DIR = 'case_{0}dB'.format(int(sinr_thrsld))
    POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(l, m, sigma_shadowing)
    logs_dir = glob.glob(os.path.join('.', LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"))
    #TODO: to be deleted
    print os.path.join('.', LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"), logs_dir

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
    print  l, m, sigma_shadowing, backoff
    print sim_intensity
    print sim_plr
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
    # ana_ee = csv_df.apply(lambda x: (x.iloc[0:MAX_TRANS]*power_levels).sum()/(1-x.values[-2]), axis=1)
    ana_ee = csv_df.apply(lambda x: (1-x.values[-2])/(x.iloc[0:MAX_TRANS]*power_levels).sum(), axis=1)
    ana_avg_nb = csv_df.apply(lambda x: x.iloc[0:MAX_TRANS].sum(), axis=1)
    return ana_intensity, ana_plr, ana_thrpt, ana_ee, ana_avg_nb

fig, axes = plt.subplots(1, 3, figsize=FIGSIZE, sharex=False)

for setting in SETTING:
    # First is perfect power control case, Second is sigma 1, Third is 3
    sinr = int(setting['SINR_THLD'])
    CASE_DIR = "case_{0}dB".format(sinr)
    ANA_DATA_FOLDED = os.path.join(LOG_DIR, ANA_DIR, SUB_DIR, CASE_DIR)
    for LM_comb, color, marker, marker_every, l_type, label in zip(
            setting['LM_PAIR'],
            setting['LINE_COLOR'],
            setting['MARKER_T'],
            setting['MARKER_EVERY'],
            setting['LINE_TYPE'],
            setting['LABEL']
    ):
        ana_result_f = os.path.join(
            ANA_DATA_FOLDED,
            "analytical_result_K={0}_threshold={1}dB_l={2}_m={3}_mufading={4}_sigma={5}.csv".format(
            MAX_TRANS,
            setting['SINR_THLD'],
            LM_comb[0], LM_comb[1],
            int(setting['MU_FADING']),
            int(setting['SIGMA_S'])
            )
        )
        ana_intensity, ana_plr, ana_thrpt, ana_ee, ana_avg_nb = analytic_data_process(
            ana_result_f,
            LM_comb[0],
            LM_comb[1]
        )
        ana_metrics = [ana_plr, ana_thrpt, ana_ee, ana_avg_nb]
        for i, (ax, title) in enumerate(zip(axes, METRICS)):
            # Iterate to plot for each metric
            ax.plot(
                ana_intensity,
                ana_metrics[i],
                color=color,
                marker=marker,
                markevery=marker_every,
                linestyle=l_type,
                linewidth=setting['LINEWIDTH'],
                label=label
            )
            ax.set_xlim(X_RANGE[i])
            ax.set_ylim(Y_RANGE[i])
            ax.set_title("({0}).".format(SUBFIG_ENUM[i])+METRICS[i], size=20)
            ax.set_yscale(Y_SCALE[i])
            # ax.set_xlabel("({0})\nFresh Packet Arrival Rate".format(SUBFIG_ENUM[i]))
            ax.set_xlabel("Fresh Packet Arrival Rate")

            ax.grid()

# Separately place legend into two sub-figures.
handles, labels = axes[1].get_legend_handles_labels()
# axes[0].legend(handles[0:3:1], labels[0:3:1], loc='best', numpoints=1, fancybox=True, framealpha=0.5)
# # axes[1].legend(handles[3:6:1], labels[3:6:1], loc='best', numpoints=1, fancybox=True, framealpha=0.5)
# axes[2].legend(handles[3::1], labels[3::1], loc='best', numpoints=1, fancybox=True, framealpha=0.5)


# handles, labels = axes[0].get_legend_handles_labels()

# axes[2, 0].annotate('cross point', xy=(1.08, 0.5), xycoords='data',
#             xytext=(1.05, 1e-3), textcoords='cross point',
#             arrowprops=dict(facecolor='black', shrink=0.05),
#             horizontalalignment='right', verticalalignment='top',
# )
# plt.subplots_adjust(hspace=0.05)
# ordered_handles = [handles[1:3:1], handles[3::1], handles[3]]
# labels = []
# plt.legend(bbox_to_anchor=(0.12, -0.12, 0.79, 1), loc=8, ncol=5, mode="expand", bbox_transform=plt.gcf().transFigure)
# # Figurs of type eps does not support transparencies natively.
# plt.subplots_adjust(wspace=0.08)


sinr_thrsld, l, m, sigma_shadowing, backoff= 3.0, 1, 1, 1, 36
sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff)

# plt.savefig(os.path.join("figures", FIG_NAME), format='eps', dpi=300, bbox_inches='tight', transparent=True)
# plt.show()
