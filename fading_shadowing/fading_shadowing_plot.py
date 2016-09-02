# -*- coding: utf-8 -*-
__author__ = 'qsong'

import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.stats as st
import pandas as pd
import glob
import json

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'axes.labelsize': 30,
}
plt.rcParams.update(params)

LOG_DIR = 'logs'
SUB_DIR = 'fading_shadowing'
# DATA_FOLDED = 'case_0dB'
FIGSIZE = (18, 18)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.0
X_END = 2.01
X_STEP = 0.1
Y_START = 0.0
Y_END = 1.01
Y_STEP = 0.1


CSV_PREFIX = "fading_shadowing_analytical_result"


def sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff):
    SUBSUB_DIR = "backoff_{0}".format(backoff)
    CASE_DIR = 'case_{0}dB'.format(sinr_thrsld)
    POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(l, m, sigma_shadowing)
    logs_dir = glob.glob(os.path.join("..", LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"))
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


def sim_plot_line(target_axes, x, y, line_properties, mode):
    if mode != "SEMILOG":
        target_axes.plot(
            x, y,
            color=line_properties["color"],
            linestyle=line_properties["linestyle"],
            linewidth=line_properties["linewidth"],
            label=line_properties["sim_label"]
        )
    else:
         target_axes.semilogy(
            x, y,
            color=line_properties["color"],
            marker=line_properties["marker"],
            label=line_properties["sim_label"]
        )

def ana_plot_line(target_axes, x, y, line_properties, mode):

    if mode != "SEMILOGY":
        target_axes.plot(
            x, y,
            color=line_properties["color"],
            linestyle=line_properties["linestyle"],
            linewidth=line_properties["linewidth"],
            label=line_properties["sim_label"]
        )
    else:
        target_axes.semilogy(
            x, y,
            color=line_properties["color"],
            marker=line_properties["marker"],
            label=line_properties["ana_label"]
        )




def generate_plot(plot_config_f):

    figs= []

    with open(plot_config_f) as json_file:
        json_config = json.load(json_file)

    for i in range(1, 4, 1):

        figs.append(plt.figure(i, figsize=FIGSIZE))

    for x, i in enumerate(["1"]):
        fig = figs[x]
        ax = fig.add_subplot(111)    # The big subplot
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        # Turn off axis lines and ticks of the big subplot
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
        case_config = json_config[i]
        SIGMA_SHADOWING = case_config["SIGMA_SHADOWING"]
        SINR_THRESLD = case_config["SINR_THRSLD"]
        BACKOFF = case_config["BACKOFF"]
        for j in ["IDENTIC", "INCREMENT", "DECREMENT"]:
            l = json_config["POWER_CHANGE"][j]["l"]
            m = json_config["POWER_CHANGE"][j]["m"]
            sim_intensity, sim_plr = sim_data_process(SINR_THRESLD, l, m, SIGMA_SHADOWING, BACKOFF)
            ana_result_f = os.path.join(
                "{0}_threshold={1}dB_l={2}_m={3}_sigma={4}.csv".format(
                    CSV_PREFIX, SINR_THRESLD, l, m, SIGMA_SHADOWING)
            )
            ana_intensity, ana_plr = analytic_data_process(ana_result_f)
            sim_plot_line(ax1, sim_intensity, sim_plr, json_config["POWER_CHANGE"][j], "NORMAL")
            ana_plot_line(ax1, ana_intensity, ana_plr, json_config["POWER_CHANGE"][j], "NORMAL")
            ax1.legend(loc='best', numpoints=2)
            ax1.set_yticks(np.arange(case_config["Y_AXIS"][0], case_config["Y_AXIS"][1], case_config["Y_AXIS"][2]))
            ax1.set_xticks(np.arange(case_config["X_AXIS"][0], case_config["X_AXIS"][1], case_config["X_AXIS"][2]))
            ax1.axis(
                [case_config["X_AXIS"][0],
                 case_config["X_AXIS"][1],
                 case_config["X_AXIS"][0],
                 case_config["Y_AXIS"][1]
                ]
            )
            ax1.grid()

            sim_plot_line(ax2, sim_intensity, sim_plr, json_config["POWER_CHANGE"][j], "SEMILOGY")
            ana_plot_line(ax2, ana_intensity, ana_plr, json_config["POWER_CHANGE"][j], "SEMILOGY")
            ax2.legend(loc='best', numpoints=2)
            ax2.set_xticks(np.arange(case_config["X_AXIS"][0], case_config["X_AXIS"][1], case_config["X_AXIS"][2]))
            ax2.axis(
                [case_config["X_AXIS"][0],
                 case_config["X_AXIS"][1],
                 case_config["X_AXIS"][0],
                 case_config["Y_AXIS"][1]
                ]
            )
            ax2.grid()

        fig.savefig('fading_case{0}.eps'.format(x), format='eps', dpi=300)

if __name__ == "__main__":

        plot_config_f = "plot_cfg.json"
        generate_plot(plot_config_f)
        plt.show()



