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

LOG_DIR = 'logs'
SUB_DIR = 'fading_shadowing_diff_ipc'
SUBSUB_DIR = "backoff_50"
CASE_DIR = 'case_3dB'
NO_DIR = "l_1_m_1_sigma_s_1"
MORE_DIR = "l_2_m_1_sigma_s_1"
LESS_DIR = "l_1_m_2_sigma_s_1"

diff_all_logs = glob.glob(os.path.join(LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, NO_DIR, "*.csv"))
sim_diff_intensity = []
sim_diff_plr = []
for csv_file in diff_all_logs:
    csv_df = pd.read_csv(csv_file, sep=',', header=None)
    plr = csv_df.values[:, -2]
    alpha = csv_df.values[:, -1][0]
    sim_diff_intensity.append(alpha)
    sim_diff_plr.append(np.mean(plr))

SUB_DIR = 'fading_shadowing_same_ipc'
same_all_logs = glob.glob(os.path.join(LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, NO_DIR, "*.csv"))
sim_same_intensity = []
sim_same_plr = []
for csv_file in same_all_logs:
    csv_df = pd.read_csv(csv_file, sep=',', header=None)
    plr = csv_df.values[:, -2]
    alpha = csv_df.values[:, -1][0]
    sim_same_intensity.append(alpha)
    sim_same_plr.append(np.mean(plr))
    
    
plt.figure(1, figsize=(18, 18))
plt.title("The simulation result comparison in terms of packet loss rate, SINR threshold = 3dB")
plt.subplot(2, 1, 1)

plt.plot(sim_diff_intensity, sim_diff_plr, 'b+', label="Same imperfect power control error for each transmission")
# plt.plot(ana_more_intensity, ana_more_plr, 'g-', label="analytical result with power incremnt")
plt.plot(sim_same_intensity, sim_same_plr, 'r*', label="Different imperfect power control error for each transmission")

plt.legend(loc='best', numpoints=2)

plt.subplot(2, 1, 2)
plt.semilogy(sim_diff_intensity, sim_diff_plr, 'b+', label="Same imperfect power control error for each transmission")
# plt.plot(ana_more_intensity, ana_more_plr, 'g-', label="analytical result with power incremnt")
plt.semilogy(sim_same_intensity, sim_same_plr, 'r*', label="Different imperfect power control error for each transmission")
plt.legend(loc='best', numpoints=2)

plt.savefig('diff_same_comp.eps', format='eps', dpi=300)

plt.show()