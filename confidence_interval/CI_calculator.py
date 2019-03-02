__author__ = 'qsong'

# -*- coding: utf-8 -*-
__author__ = 'qsong'

import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import os
import pandas as pd

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'axes.labelsize': 30,
}
plt.rcParams.update(params)


DATASET = 'logs'

csv_file = os.path.join(DATASET, "simd=10000_N=400_threshold=3.0_l=1_m=1_backoff=150_start=0.46_end=0.46_simstep=0.02_045446.csv")

csv_df = pd.read_csv(csv_file, sep=',', header=None)


intensity = csv_df.values[:, 6]
print intensity

plr = csv_df.values[:, 5]
ci =st.t.interval(0.95, len(plr)-1, loc=np.mean(plr), scale=st.sem(plr))
print plr
print ci

# target_df.set_index(intensity)
# target_df.to_csv('result.csv', spe=',')
#
# print target_df
#
# ci_results = []
# for n in range(target_df.shape[0]):
#     test = target_df.iloc[n, :].values
#     ci =[]
#     ci.append(np.mean(test)-2*st.sem(test))
#     # ci =st.t.interval(0.95, len(test)-1, loc=np.mean(test), scale=st.sem(test))
#     ci.append(np.mean(test)+2*st.sem(test))
#     ci_results.append(ci)
#     # print ci
#
# ci_low = [abs(e[0]) for e in ci_results]
# ci_high = [abs(e[1]) for e in ci_results]




intensity = [0.2, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44]
ci_list = [
    (0.00057309854841668009, 0.00082032273503805956),
    (0.00090782709992968057, 0.001263821596166926),
    (0.0014510867222395905, 0.001895169682176955),
    (0.0029988938645252763, 0.0035825841126334757),
    (0.0050126127128484015, 0.0056074709087366847),
    (0.0074241434238066367, 0.0085921424525662107),
    (0.012596555732940157, 0.013942639203443337),
    (0.019423532015275961, 0.021393914044568016),
    (0.029650392640485725, 0.032737002589868509),
    (0.047785142536643767, 0.053066901056349015),
    (0.077394621815116152, 0.083538855592149117),
    (0.11329860696265433, 0.12364207273915266),
    (0.15497381193853241, 0.17434131819783036),
    # (0.21095106246169074, 0.23408939973738499)
]

plt.figure(1)
DATASET = 'confidence_interval'
ana_result_f_with = os.path.join(DATASET, "analytical_result_K=5_threshold=3.0_l=1_m=1.csv")

with open(ana_result_f_with, 'r') as ana_with_f_handler:
        ana_with_csv = list(csv.reader(ana_with_f_handler))

ana_with_intensity = [float(element[-1]) for element in ana_with_csv]
ana_with_plr = [float(element[-2]) for element in ana_with_csv]

# plt.title("The packe losss rate evolutaion with SINR = -3 dB, 500 devices")
# plt.xlabel("The fresh packe arrival intensity")
# plt.ylabel("The packet loss rate")
# plt.plot(ana_with_intensity, ana_with_plr, 'b-', label="analytical result with power increment")
# plt.xlabel("The fresh packe arrival intensity")
# plt.ylabel("The packet loss rate")
# plt.legend(loc='best', numpoints=2)
# plt.axis([0.2, 0.6])
plt.grid()

plt.semilogy(ana_with_intensity, ana_with_plr, 'b-', label ="analytical result")
plt.semilogy(intensity, [e[0] for e in ci_list], 'r--', label = "lower bound of confidence interval 95%")
plt.semilogy(intensity, [e[1] for e in ci_list], 'k--', label = "upper bound of confidence interval 95%")
plt.legend(loc='lower right', numpoints=2)
plt.ylabel("Packet loss rate")
plt.xlabel("Fresh packet arrival intensity")
plt.savefig('ci.eps', format='eps', dpi=300)
plt.show()
############################################################################################






