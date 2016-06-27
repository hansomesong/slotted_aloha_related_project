# -*- coding: utf-8 -*-
__author__ = 'qsong'


import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

plt.figure(1)
# The case for threshold 1.0
sim_result_f_with = "sim_result_sim=1_N=500_threshold=1.0_l=2_m=1.csv"
ana_result_f_with = "analytical_result_threshold=1.0_l=2_m=1.csv"

sim_result_f_no = "sim_result_sim=1_N=500_threshold=1.0_l=1_m=1.csv"
ana_result_f_no = "analytical_result_threshold=1.0_l=1_m=1.csv"

with open(ana_result_f_with, 'r') as ana_with_f_handler, \
        open(sim_result_f_with, 'r') as sim_with_f_handler, \
        open(ana_result_f_no, 'r') as ana_no_f_handler, \
        open(sim_result_f_no, 'r') as sim_no_f_handler:
    ana_with_csv = list(csv.reader(ana_with_f_handler))
    sim_with_csv = list(csv.reader(sim_with_f_handler))
    ana_no_csv = list(csv.reader(ana_no_f_handler))
    sim_no_csv = list(csv.reader(sim_no_f_handler))

    ana_with_intensity = [float(element[-1]) for element in ana_with_csv]
    ana_with_plr = [float(element[-2]) for element in ana_with_csv]
    sim_with_intensity = [float(e[-1]) for e in sim_with_csv]
    sim_with_plr = [float(e[-2]) for e in sim_with_csv]

    ana_no_intensity = [float(element[-1]) for element in ana_no_csv]
    ana_no_plr = [float(element[-2]) for element in ana_no_csv]
    sim_no_intensity = [float(e[-1]) for e in sim_no_csv]
    sim_no_plr = [float(e[-2]) for e in sim_no_csv]

    # plt.xticks(np.arange(0, max(ana_with_intensity), 0.1))

    plt.plot(ana_with_intensity, ana_with_plr,
             sim_with_intensity, sim_with_plr, 'g^',
             ana_no_intensity, ana_no_plr, 'r-',
             sim_no_intensity, sim_no_plr, 'm*',
    )
    plt.grid()


############################################################################################
# The case threshould: -3dB
plt.figure(2)
sim_result_f_with = "sim_result_sim=1_N=500_threshold=0.5_l=2_m=1.csv"
ana_result_f_with = "analytical_result_threshold=0.5_l=2_m=1.csv"

sim_result_f_no = "sim_result_sim=1_N=500_threshold=0.5_l=1_m=1.csv"
ana_result_f_no = "analytical_result_threshold=0.5_l=1_m=1.csv"

with open(ana_result_f_with, 'r') as ana_with_f_handler, \
        open(sim_result_f_with, 'r') as sim_with_f_handler, \
        open(ana_result_f_no, 'r') as ana_no_f_handler, \
        open(sim_result_f_no, 'r') as sim_no_f_handler:
    ana_with_csv = list(csv.reader(ana_with_f_handler))
    sim_with_csv = list(csv.reader(sim_with_f_handler))
    ana_no_csv = list(csv.reader(ana_no_f_handler))
    sim_no_csv = list(csv.reader(sim_no_f_handler))

    ana_with_intensity = [float(element[-1]) for element in ana_with_csv]
    ana_with_plr = [float(element[-2]) for element in ana_with_csv]
    sim_with_intensity = [float(e[-1]) for e in sim_with_csv]
    sim_with_plr = [float(e[-2]) for e in sim_with_csv]

    ana_no_intensity = [float(element[-1]) for element in ana_no_csv]
    ana_no_plr = [float(element[-2]) for element in ana_no_csv]
    sim_no_intensity = [float(e[-1]) for e in sim_no_csv]
    sim_no_plr = [float(e[-2]) for e in sim_no_csv]

    # plt.xticks(np.arange(0, max(ana_with_intensity), 0.1))

    plt.plot(ana_with_intensity, ana_with_plr,
             sim_with_intensity, sim_with_plr, 'g^',
             ana_no_intensity, ana_no_plr, 'r-',
             sim_no_intensity, sim_no_plr, 'm*',
    )
    plt.grid()
plt.show()