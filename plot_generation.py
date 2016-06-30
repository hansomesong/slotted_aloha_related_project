# -*- coding: utf-8 -*-
__author__ = 'qsong'


import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker


DATA_FOLDED = os.path.join('.')
# DATA_FOLDED = os.path.join('data_files', '20160627')
plt.figure(1)
plt.subplot(1, 2, 1)
# The case for threshold 1.0
sim_result_f_no = os.path.join(DATA_FOLDED, "sim_result_simd=100000_N=500_threshold=1.0_l=1_m=1_PID=6075.csv")
DATA_FOLDED = os.path.join('.')
ana_result_f_no = os.path.join(DATA_FOLDED, "analytical_result_threshold=1.0_l=1_m=1.csv")

DATA_FOLDED = os.path.join('data_files', '20160627')
sim_result_f_with = os.path.join(DATA_FOLDED, 'sim_result_sim=1_N=500_threshold=1.0_l=2_m=1.csv')
ana_result_f_with = os.path.join(DATA_FOLDED, "analytical_result_threshold=1.0_l=2_m=1.csv")


with open(ana_result_f_with, 'r') as ana_with_f_handler, \
        open(sim_result_f_with, 'r') as sim_with_f_handler, \
        open(ana_result_f_no, 'r') as ana_no_f_handler, \
        open(sim_result_f_no, 'r') as sim_no_f_handler:
    ana_with_csv = list(csv.reader(ana_with_f_handler))
    sim_with_csv = list(csv.reader(sim_with_f_handler))

    ana_no_csv = list(csv.reader(ana_no_f_handler))
    sim_no_csv = list(csv.reader(sim_no_f_handler))
    print len(sim_no_csv)

    ana_with_intensity = [float(element[-1]) for element in ana_with_csv]
    ana_with_plr = [float(element[-2]) for element in ana_with_csv]
    sim_with_intensity = [float(e[-1]) for e in sim_with_csv]
    sim_with_plr = [float(e[-2]) for e in sim_with_csv]

    ana_no_intensity = [float(element[-1]) for element in ana_no_csv]
    ana_no_plr = [float(element[-2]) for element in ana_no_csv]
    sim_no_intensity = [float(e[-1]) for e in sim_no_csv]
    sim_no_plr = [float(e[-2]) for e in sim_no_csv]

    # plt.xticks(np.arange(0, max(ana_with_intensity), 0.1))
    plt.title("The packe losss rate evolutaion with SINR = 0 dB, 500 devices")
    plt.xlabel("The fresh packe arrival intensity")
    plt.ylabel("The packet loss rate")
    plt.plot(ana_with_intensity, ana_with_plr, 'b-', label="analytical result with power increment")
    plt.plot(sim_with_intensity, sim_with_plr, 'g^', label="simulation result with power increment")
    plt.plot(ana_no_intensity, ana_no_plr, 'r--', label="analytical result case without power increment")
    plt.plot(sim_no_intensity, sim_no_plr, 'm*', label="simulation result case without power increment")
    plt.xlabel("The fresh packe arrival intensity")
    plt.ylabel("The packet loss rate")
    plt.legend(loc='best', numpoints=2)
    plt.axis([0, 1.05, 0, 1.05])
    plt.grid()

plt.subplot(1, 2, 2)
plt.grid()
plt.xlabel("The fresh packe arrival intensity")
plt.ylabel("The packet loss rate")
plt.semilogy(sim_no_intensity, sim_no_plr, 'g^')
plt.semilogy(ana_no_intensity, ana_no_plr, 'b-')


############################################################################################
# The case threshould: -3dB
DATA_FOLDED = os.path.join('data_files', '20160627')
plt.figure(2)
plt.subplot(1, 2, 1)
sim_result_f_with = os.path.join(DATA_FOLDED, "sim_result_sim=1_N=500_threshold=0.5_l=2_m=1.csv")
ana_result_f_with = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.5_l=2_m=1.csv")

sim_result_f_decre = os.path.join(DATA_FOLDED, "sim_result_sim=1_N=500_threshold=0.5_l=1_m=2.csv")
ana_result_f_decre = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.5_l=1_m=2.csv")

sim_result_f_no = os.path.join(DATA_FOLDED, "sim_result_sim=1_N=500_threshold=0.5_l=1_m=1.csv")
ana_result_f_no = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.5_l=1_m=1.csv")

with open(ana_result_f_with, 'r') as ana_with_f_handler, \
        open(sim_result_f_with, 'r') as sim_with_f_handler, \
        open(ana_result_f_no, 'r') as ana_no_f_handler, \
        open(sim_result_f_no, 'r') as sim_no_f_handler,\
        open(ana_result_f_decre, 'r') as ana_decre_f_handler, \
        open(sim_result_f_decre, 'r') as sim_decre_f_handler:
    ana_with_csv = list(csv.reader(ana_with_f_handler))
    sim_with_csv = list(csv.reader(sim_with_f_handler))
    ana_no_csv = list(csv.reader(ana_no_f_handler))
    sim_no_csv = list(csv.reader(sim_no_f_handler))
    ana_decre_csv = list(csv.reader(ana_decre_f_handler))
    sim_decre_csv = list(csv.reader(sim_decre_f_handler))

    ana_with_intensity = [float(element[-1]) for element in ana_with_csv]
    ana_with_plr = [float(element[-2]) for element in ana_with_csv]
    sim_with_intensity = [float(e[-1]) for e in sim_with_csv]
    sim_with_plr = [float(e[-2]) for e in sim_with_csv]
    ana_decre_intensity = [float(element[-1]) for element in ana_decre_csv]
    ana_decre_plr = [float(element[-2]) for element in ana_decre_csv]
    sim_decre_intensity = [float(e[-1]) for e in sim_decre_csv]
    sim_decre_plr = [float(e[-2]) for e in sim_decre_csv]

    ana_no_intensity = [float(element[-1]) for element in ana_no_csv]
    ana_no_plr = [float(element[-2]) for element in ana_no_csv]
    sim_no_intensity = [float(e[-1]) for e in sim_no_csv]
    sim_no_plr = [float(e[-2]) for e in sim_no_csv]

    plt.title("The packe losss rate evolutaion with SINR = -3 dB, 500 devices")
    plt.xlabel("The fresh packe arrival intensity")
    plt.ylabel("The packet loss rate")
    plt.plot(ana_with_intensity, ana_with_plr, 'b-', label="analytical result with power increment")
    plt.plot(sim_with_intensity, sim_with_plr, 'g^', label="simulation result with power increment")
    plt.plot(ana_no_intensity, ana_no_plr, 'r--', label="analytical result case without power increment")
    plt.plot(sim_no_intensity, sim_no_plr, 'm*', label="simulation result case without power increment")
    plt.plot(ana_decre_intensity, ana_decre_plr, 'm:', label="analytical result with power decrement")
    plt.plot(sim_decre_intensity, sim_decre_plr, 'ko', label="simulation result with power decrement")
    plt.xlabel("The fresh packe arrival intensity")
    plt.ylabel("The packet loss rate")
    plt.legend(loc='best', numpoints=2)
    plt.axis([0, 2.05, 0, 1.05])
    plt.grid()
plt.subplot(1, 2, 2)
plt.semilogy(sim_with_intensity, sim_with_plr, 'g^')
plt.semilogy(ana_with_intensity, ana_with_plr, 'b-')

# The case 3: threshold: -10dB
plt.figure(3)
DATA_FOLDED = 'data_files/20160627'
sim_result_f_with = os.path.join(DATA_FOLDED, 'sim_result_sim=1_N=2000_threshold=0.1_l=2_m=1.csv')
ana_result_f_with = os.path.join(DATA_FOLDED, 'analytical_result_threshold=0.1_l=2_m=1.csv')
DATA_FOLDED = 'data_files/20160627'
sim_result_f_no = os.path.join(DATA_FOLDED, "sim_result_sim=1_N=2000_threshold=0.1_l=1_m=1.csv")
ana_result_f_no = os.path.join(DATA_FOLDED, "analytical_result_threshold=0.1_l=1_m=1.csv")

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

    plt.title("The packe losss rate evolutaion with SINR = -10 dB, 2000 devices")
    plt.plot(ana_with_intensity, ana_with_plr, 'b-', label="analytical result with power increment")
    plt.plot(sim_with_intensity, sim_with_plr, 'g^', label="simulation result with power increment")
    plt.plot(ana_no_intensity, ana_no_plr, 'r--', label="analytical result case without power increment")
    plt.plot(sim_no_intensity, sim_no_plr, 'm*', label="simulation result case without power increment")
    plt.xlabel("The fresh packe arrival intensity")
    plt.ylabel("The packet loss rate")
    plt.legend(loc='best', numpoints=2)
    plt.grid()
    plt.axis([4.5, 8.0, 0, 1.05])


############################################################################################
# The case threshould: -3dB, but with power decrement
# plt.figure(4)
# sim_result_f_decre = "data_files/sim_result_sim=1_N=500_threshold=0.5_l=1_m=2.csv"
# ana_result_f_decre = "data_files/analytical_result_threshold=0.5_l=1_m=2.csv"
#
# sim_result_f_no = "data_files/sim_result_sim=1_N=500_threshold=0.5_l=1_m=1.csv"
# ana_result_f_no = "data_files/analytical_result_threshold=0.5_l=1_m=1.csv"
#
# with open(ana_result_f_with, 'r') as ana_with_f_handler, \
#         open(sim_result_f_with, 'r') as sim_with_f_handler, \
#         open(ana_result_f_no, 'r') as ana_no_f_handler, \
#         open(sim_result_f_no, 'r') as sim_no_f_handler:
#     ana_with_csv = list(csv.reader(ana_with_f_handler))
#     sim_with_csv = list(csv.reader(sim_with_f_handler))
#     ana_no_csv = list(csv.reader(ana_no_f_handler))
#     sim_no_csv = list(csv.reader(sim_no_f_handler))
#
#     ana_with_intensity = [float(element[-1]) for element in ana_with_csv]
#     ana_with_plr = [float(element[-2]) for element in ana_with_csv]
#     sim_with_intensity = [float(e[-1]) for e in sim_with_csv]
#     sim_with_plr = [float(e[-2]) for e in sim_with_csv]
#
#     ana_no_intensity = [float(element[-1]) for element in ana_no_csv]
#     ana_no_plr = [float(element[-2]) for element in ana_no_csv]
#     sim_no_intensity = [float(e[-1]) for e in sim_no_csv]
#     sim_no_plr = [float(e[-2]) for e in sim_no_csv]
#
#     plt.title("The packe losss rate evolutaion with SINR = -3 dB, 500 devices")
#     plt.xlabel("The fresh packe arrival intensity")
#     plt.ylabel("The packet loss rate")
#     plt.plot(ana_with_intensity, ana_with_plr, 'b-', label="analytical result with power decrement")
#     plt.plot(sim_with_intensity, sim_with_plr, 'g^', label="simulation result with power decrement")
#     plt.plot(ana_no_intensity, ana_no_plr, 'r--', label="analytical result case without power decrement")
#     plt.plot(sim_no_intensity, sim_no_plr, 'm*', label="simulation result case without power decrement")
#     plt.xlabel("The fresh packe arrival intensity")
#     plt.ylabel("The packet loss rate")
#     plt.legend(loc='best', numpoints=2)
#     plt.axis([0, 2.05, 0, 1.05])
#     plt.grid()

plt.show()