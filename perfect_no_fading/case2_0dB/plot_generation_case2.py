# -*- coding: utf-8 -*-
__author__ = 'qsong'


import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker

params = {
    'legend.fontsize': 25,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 30,
}
plt.rcParams.update(params)


label_fs = 40
DATA_FOLDED = '.'
plt.figure(1, figsize=(15, 12))
plt.title("The packet loss rate evolution with SINR threshold = 3dB")

plt.subplot(1, 1, 1)
# The case for threshold 3.0 dB
ana_result_f_no = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.0_l=1_m=1.csv")
ana_result_f_more = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.0_l=2_m=1.csv")
ana_result_f_less = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.0_l=1_m=2.csv")

with open(ana_result_f_no, 'r') as ana_no_f_handler, \
        open(ana_result_f_more, 'r') as ana_more_f_handler, \
        open(ana_result_f_less, 'r') as ana_less_f_handler:

    ana_more_csv = list(csv.reader(ana_more_f_handler))
    ana_more_intensity = [float(element[-1]) for element in ana_more_csv]
    ana_more_plr = [float(element[-2]) for element in ana_more_csv]

    ana_no_csv = list(csv.reader(ana_no_f_handler))
    ana_no_intensity = [float(element[-1]) for element in ana_no_csv]
    ana_no_plr = [float(element[-2]) for element in ana_no_csv]

    ana_less_csv = list(csv.reader(ana_less_f_handler))
    ana_less_intensity = [float(element[-1]) for element in ana_less_csv]
    ana_less_plr = [float(element[-2]) for element in ana_less_csv]

sim_result_f_less = os.path.join(
    DATA_FOLDED,
    "simd=50000_N=500_threshold=0.0_l=1_m=2_backoff=150_start=0.26_end=1.2_simstep=0.02_112343.csv"
)

with open(sim_result_f_less, 'r') as sim_less_f_handler:
    sim_less_csv = list(csv.reader(sim_less_f_handler))
    sim_less_intensity = [float(element[-1]) for element in sim_less_csv]
    sim_less_plr = [float(element[-2]) for element in sim_less_csv]

sim_result_f_more = os.path.join(
    DATA_FOLDED, "simd=50000_N=500_threshold=0.0_l=2_m=1_backoff=150_start=0.8_end=1.2_simstep=0.02_073345.csv"
)
with open(sim_result_f_more, 'r') as sim_more_f_handler:
    sim_more_csv = list(csv.reader(sim_more_f_handler))
    sim_more_intensity = [float(element[-1]) for element in sim_more_csv]
    sim_more_plr = [float(element[-2]) for element in sim_more_csv]
#     # plt.xticks(np.arange(0, max(ana_with_intensity), 0.1))

sim_result_f_no = os.path.join(
    DATA_FOLDED, "simd=50000_N=500_threshold=0.0_l=1_m=1_backoff=150_start=0.3_end=1.5_simstep=0.02_164805.csv"
)
with open(sim_result_f_no, 'r') as sim_no_f_handler:
    sim_no_csv = list(csv.reader(sim_no_f_handler))
    sim_no_intensity = [float(element[-1]) for element in sim_no_csv]
    sim_no_plr = [float(element[-2]) for element in sim_no_csv]
#     # plt.xticks(np.arange(0, max(ana_with_intensity), 0.1))



#     plt.plot(ana_less_intensity, ana_less_plr, 'r-', label="analytical result with power decrement")
#     plt.plot(ana_no_intensity, ana_no_plr, 'b-', label="analytical result with identical power")
#     plt.plot(ana_more_intensity, ana_more_plr, 'g-', label="analytical result with power increment")
#     plt.plot(sim_less_intensity, sim_less_plr, 'r+', label="simulation result with power decrement")
#     plt.plot(sim_more_intensity, sim_more_plr, 'go', label="simulation result with power increment")
#     plt.plot(sim_no_intensity, sim_no_plr, 'b*', label="simulation result with identical power")
#
#     plt.xlabel("Fresh packet arrival intensity")
#     plt.ylabel("Packet loss rate")
#     plt.legend(loc='best', numpoints=2)
#     plt.axis([0, 2.05, 0, 1.05])
#     plt.grid()
#
# plt.subplot(2, 1, 2)
plt.grid()

plt.semilogy(ana_less_intensity, ana_less_plr, 'r-', label="analytical result with power decrement")
plt.semilogy(ana_no_intensity, ana_no_plr, 'b-', label="analytical result case with identical power")
plt.semilogy(ana_more_intensity, ana_more_plr, 'g-', label="analytical result with power increment")
plt.semilogy(sim_less_intensity, sim_less_plr, 'r+', label="simulation result with power decrement")
plt.semilogy(sim_more_intensity, sim_more_plr, 'go', label="simulation result with power increment")
plt.semilogy(sim_no_intensity, sim_no_plr, 'b*', label="simulation result with identical power")

plt.legend(loc='lower right', numpoints=2)
plt.ylabel("Packet loss rate")
plt.xlabel("Fresh packet arrival intensity")

plt.axis([0, 2.05, 0, 1.05])
plt.savefig('logarithme_case2.eps', format='eps', dpi=300)
plt.show()
