# -*- coding: utf-8 -*-
__author__ = 'qsong'

import csv
import matplotlib.pyplot as plt
from matplotlib  import cm

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
SUB_DIR = 'multiple_reception'
CASE_DIR = 'case_3dB'
SUB_CASE_DIR = "bs_0.01_p_0.002"
DATA_FOLDED = '.'
FIGSIZE = (15, 6)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.1
X_END = 1.0
X_STEP = 0.1
Y_START = 10e-5
Y_END = 1.01
Y_STEP = 0.1

MAX_TRANS = 1
LINEWIDTH = 2


def sim_data_process(folder_dir):
    sim_intensity = []
    sim_plr =  []
    for csv_file in logs_dir:
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        plr_df = csv_df.values[:, -2]
        # sort the obtained to remove the max and min
        plr_df.sort()
        # remove the max and min
        plr_df = plr_df[1:-1]
        avg_plr = plr_df.mean()
        # Calculates the standard error of the mean (or standard error of measurement)
        # Here we have to be careful to keep all y values positive:
        ci_min_plr = max(avg_plr-1.96*st.sem(plr_df), 1e-7)
        ci_max_plr = avg_plr + 1.96*st.sem(plr_df)
        alpha = csv_df.values[:, -1][0]
        sim_intensity.append(alpha)
        sim_plr.append([ci_min_plr, avg_plr, ci_max_plr])
    return sim_intensity, sim_plr

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

if __name__ == '__main__':
    # logs_dir = glob.glob(os.path.join("..",  LOG_DIR, SUB_DIR, CASE_DIR, SUB_CASE_DIR , "*.csv"))
    # sim_intensity, sim_plr_list = sim_data_process(logs_dir)
    # sim_plr = [element[1] for element in sim_plr_list]
    # sim_plr_lower = [element[1]-element[0] for element in sim_plr_list]
    # sim_plr_upper = [element[2]-element[1] for element in sim_plr_list]
    # print sim_intensity
    # print sim_plr

    width  = 40
    intensity = 0.4
    density_bs = 0.04
    bs_nb = int(np.random.poisson(density_bs*np.power(width, 2)))
    device_nb = int(np.random.poisson(intensity*np.power(width, 2)))
    print "BS_STATION NUMBER:", bs_nb, "DEVICE_NUMBER:", device_nb


    # device_x_array = np.random.uniform(-width/2.0, width/2.0, device_nb)
    # device_y_array = np.random.uniform(-width/2.0, width/2.0, device_nb)
    # bs_x_array = np.random.uniform(-width/2.0, width/2.0, bs_nb)
    # bs_y_array = np.random.uniform(-width/2.0, width/2.0, bs_nb)
    # coordinates_devices_array = zip(device_x_array, device_y_array)
    # coordinates_bs_array = zip(bs_x_array, bs_y_array)

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    # ax.set_title("One realization of device and Bs spatial distribution", fontsize=18)
    ax.set_xlabel("X-axis", fontsize=12)
    ax.set_ylabel("Y-axis", fontsize=12)
    ax.grid(True, linestyle='-', color='0.75')

    # scatter with colormap mapping to z value
    # ax.scatter(device_x_array, device_y_array, marker ='o', cmap = cm.jet)
    # ax.scatter(bs_x_array, bs_y_array, marker ='*', color="r")
    # plt.show()

    radius = 20
    # rho = radius*np.sqrt(np.random.uniform(0, 1, device_nb))
    rho = np.concatenate(([0.0], radius*np.sqrt(np.random.uniform(0, 1, device_nb-1))))

    arguments = np.random.uniform(-np.pi-10e-4, np.pi+10e-4, device_nb)
    device_x_array_2 = rho*np.cos(arguments)
    device_y_array_2 = rho*np.sin(arguments)

    bs_rho = radius*np.sqrt(np.random.uniform(0, 1, bs_nb))
    bs_arguments = np.random.uniform(-np.pi-10e-4, np.pi+10e-4, bs_nb)
    bs_x_array_2 = bs_rho*np.cos(bs_arguments)
    bs_y_array_2 = bs_rho*np.sin(bs_arguments)
    ax.scatter(device_x_array_2, device_y_array_2, marker='o', cmap=cm.jet, label="Device")
    ax.scatter(bs_x_array_2, bs_y_array_2, marker='*', color='r', s=80, cmap=cm.jet, label="BS")
    ax.legend(loc='best', numpoints=1)

    # plt.savefig('throughput_shadowing_case1.eps', format='eps', dpi=300)
    plt.savefig("nodes_spatial_dist.eps", format='eps', dpi=300, bbox_inches='tight', transparent=True)
    plt.show()
