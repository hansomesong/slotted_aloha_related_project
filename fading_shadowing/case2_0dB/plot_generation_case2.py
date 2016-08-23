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

LOG_DIR = 'logs'
SUB_DIR = 'fading_shadowing_diff_ipc'
SUBSUB_DIR = "backoff_36"
CASE_DIR = 'case_0dB'
NO_DIR = "l_1_m_1_sigma_s_1"
MORE_DIR = "l_2_m_1_sigma_s_1"
LESS_DIR = "l_1_m_2_sigma_s_1"

no_all_logs = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, NO_DIR, "*.csv"))
sim_no_intensity = []
sim_no_plr = []
for csv_file in no_all_logs:
    csv_df = pd.read_csv(csv_file, sep=',', header=None)
    plr = csv_df.values[:, -2]
    alpha = csv_df.values[:, -1][0]
    sim_no_intensity.append(alpha)
    sim_no_plr.append((np.sum(plr) -np.max(plr) - np.min(plr))/(len(plr)-2.0) )

# more_all_logs = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, MORE_DIR, "*.csv"))
# sim_more_intensity = []
# sim_more_plr = []
# for csv_file in more_all_logs:
#     csv_df = pd.read_csv(csv_file, sep=',', header=None)
#     plr = csv_df.values[:, -2]
#     alpha = csv_df.values[:, -1][0]
#     sim_more_intensity.append(alpha)
#     sim_more_plr.append(np.mean(plr))

# less_all_logs = glob.glob(os.path.join("..", "..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, LESS_DIR, "*.csv"))
# sim_less_intensity = []
# sim_less_plr = []
# for csv_file in less_all_logs:
#     csv_df = pd.read_csv(csv_file, sep=',', header=None)
#     plr = csv_df.values[:, -2]
#     alpha = csv_df.values[:, -1][0]
#     sim_less_intensity.append(alpha)
#     sim_less_plr.append(np.mean(plr))

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'axes.labelsize': 30,
}
plt.rcParams.update(params)

DATA_FOLDED = '.'
plt.figure(1, figsize=(18, 18))
plt.title("The packet loss rate evolution with SINR threshold = 0dB")

plt.subplot(2, 1, 1)
# The case for threshold 3.0 dB
ana_result_f_no = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=0dB_l=1_m=1_sigma=1.csv")
# ana_result_f_more = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=0dB_l=2_m=1_sigma=1.csv")
ana_result_f_less = os.path.join(DATA_FOLDED, "fading_shadowing_analytical_result_threshold=0dB_l=1_m=2_sigma=1.csv")

# with open(ana_result_f_no, 'r') as ana_no_f_handler, \
#         open(ana_result_f_more, 'r') as ana_more_f_handler, \
#         open(ana_result_f_less, 'r') as ana_less_f_handler:

with open(ana_result_f_no, 'r') as ana_no_f_handler:

    # ana_more_csv = list(csv.reader(ana_more_f_handler))
    # ana_more_intensity = [float(element[-1]) for element in ana_more_csv]
    # ana_more_plr = [float(element[-2]) for element in ana_more_csv]

    ana_no_csv = list(csv.reader(ana_no_f_handler))
    ana_no_intensity = [float(element[-1]) for element in ana_no_csv]
    ana_no_plr = [float(element[-2]) for element in ana_no_csv]

    # ana_less_csv = list(csv.reader(ana_less_f_handler))
    # ana_less_intensity = [float(element[-1]) for element in ana_less_csv]
    # ana_less_plr = [float(element[-2]) for element in ana_less_csv]

    plt.xlabel("The fresh packet arrival intensity")
    plt.ylabel("The packet loss rate")
    # plt.plot(ana_less_intensity, ana_less_plr, 'r-', label="analytical result with power decrement")
    plt.plot(ana_no_intensity, ana_no_plr, 'b-', label="analytical result case with identical power")
    # plt.plot(ana_more_intensity, ana_more_plr, 'g-', label="analytical result with power increment")
    plt.plot(sim_no_intensity, sim_no_plr, 'b*', label="simulation result with identical power")
    # plt.plot(sim_more_intensity, sim_more_plr, 'go', label="simulation result with power increment")
    # plt.plot(sim_less_intensity, sim_less_plr, 'r+', label="simulation result with power decrement")

    plt.xlabel("The fresh packet arrival intensity")
    plt.ylabel("The packet loss rate")
    plt.legend(loc='best', numpoints=2)
    plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.xticks(
        [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5],
        [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    )

    plt.axis([0, 1.55, 0, 1.05])
    plt.grid()

plt.subplot(2, 1, 2)
plt.grid()

# plt.semilogy(ana_less_intensity, ana_less_plr, 'r-', label="analytical result with power decrement")
plt.semilogy(ana_no_intensity, ana_no_plr, 'b-', label="analytical result case with identical power")
# plt.semilogy(ana_more_intensity, ana_more_plr, 'g-', label="analytical result with power increment")
plt.semilogy(sim_no_intensity, sim_no_plr, 'b*', label="simulation result with  identical power")
# plt.semilogy(sim_more_intensity, sim_more_plr, 'go', label="simulation result with power increment")
# plt.semilogy(sim_less_intensity, sim_less_plr, 'r+', label="simulation result with power decrement")


plt.legend(loc='best', numpoints=2)
plt.ylabel("The packet loss rate")
plt.xlabel("The fresh packet arrival intensity")

plt.axis([0, 1.55, 0, 1.05])


plt.savefig('fading_case2.eps', format='eps', dpi=300)
plt.show()
