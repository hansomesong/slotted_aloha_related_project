# -*- coding: utf-8 -*-
__author__ = 'qsong'


import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker



sim_result_f = "sim_result_sim_1_N_500_threshold_0.1_l=1_m=1.csv"
ana_result_f = "analytical_result_threshold_0.1_l=1_m=1.csv"


plt.figure(1)
with open(ana_result_f, 'r') as ana_f_handler, open(sim_result_f, 'r') as sim_f_handler:
    ana_csv = list(csv.reader(ana_f_handler))
    sim_csv = list(csv.reader(sim_f_handler))

    ana_intensity = [float(element[-1]) for element in ana_csv]
    ana_plr = [float(element[-2]) for element in ana_csv]


    sim_intensity = [float(e[-1]) for e in sim_csv]
    sim_plr = [float(e[-2]) for e in sim_csv]
    print len(ana_intensity)
    print len(ana_plr)

    print [e[-1] for e in sim_csv]
    print [e[0] for e in sim_csv]

    plt.xticks(np.arange(0, max(ana_intensity), 0.1))

    plt.plot(ana_intensity, ana_plr, sim_intensity, sim_plr, 'g^')
    plt.grid()
    plt.show()