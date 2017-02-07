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
from scipy.special import gamma as gamma_f
import scipy.special as ss

from analytical_model import sgam


params = {
    'legend.fontsize': 15,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 30,
}
plt.rcParams.update(params)

LOG_DIR = 'logs'
SUB_DIR = 'multiple_reception'
CASE_DIR = 'fading'
SUB_CASE_DIR = "bs_0.01_p_0.002_R_40|100"
DATA_FOLDED = '.'
FIGSIZE = (15, 6)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.0
X_END = 1.5
X_STEP = 0.1
Y_START = 1e-3
Y_END = 0.6
Y_STEP = 0.1

MAX_TRANS = 1
LINEWIDTH = 2

SCALE = ["log", "linear"]



def sim_data_process(folder_dir):
    sim_intensity = []
    sim_plr =  []
    for csv_file in folder_dir:
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        plr_df = csv_df.values[:, -2]
        # sort the obtained to remove the max and min
        plr_df.sort()
        # remove the max and min
        # plr_df = plr_df[1:-1]
        avg_plr = plr_df.mean()
        # Calculates the standard error of the mean (or standard error of measurement)
        # Here we have to be careful to keep all y values positive:
        ci_min_plr = max(avg_plr-1.96*st.sem(plr_df), 1e-7)
        ci_max_plr = avg_plr + 1.96*st.sem(plr_df)
        alpha = csv_df.values[:, -1][0]
        sim_intensity.append(alpha)
        sim_plr.append([ci_min_plr, avg_plr, ci_max_plr])
    return sim_intensity, sim_plr

def sim_data_process_2(folder_dir):
    sim_intensity = []
    sim_plr =  []
    for csv_file in folder_dir:
        csv_df = pd.read_csv(csv_file, sep=',').drop(["1nb", "2nb", "1pr"], axis=1)
        for name, group in csv_df.groupby(["STAT_PERCENT"]):
            print name, group["2pr"].mean(), group["2pr"].std(), st.sem(group["2pr"])
            print type(group["2pr"])
        # csv_df
        # .groupby(["STAT_PERCENT"]).mean()
        # print csv_df
        # plr_df = csv_df.values[:, -2]
        # # sort the obtained to remove the max and min
        # plr_df.sort()
        # # remove the max and min
        # # plr_df = plr_df[1:-1]
        # avg_plr = plr_df.mean()
        # # Calculates the standard error of the mean (or standard error of measurement)
        # # Here we have to be careful to keep all y values positive:
        # ci_min_plr = max(avg_plr-1.96*st.sem(plr_df), 1e-7)
        # ci_max_plr = avg_plr + 1.96*st.sem(plr_df)
        # alpha = csv_df.values[:, -1][0]
        # sim_intensity.append(alpha)
        # sim_plr.append([ci_min_plr, avg_plr, ci_max_plr])
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


def bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB):
    """
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    variant_1 = p*lambda_m/lambda_b
    constant_B = variant_1*gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)*np.exp(0.5*(2.0*sigma/gamma)**2)*np.power(THETA, 2.0/gamma)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)
    return 1.0 - 1.0/(1+np.power(constant_B, B_power))


def bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB):
    """
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    variant_1 = p*lambda_m/lambda_b
    constant_B = variant_1*gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)*np.exp(0.5*(2.0*sigma/gamma)**2)*np.power(THETA, 2.0/gamma)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)
    return np.exp(1.0/(1+np.power(constant_B, B_power))-np.exp(0.5*(2.0*sigma/gamma)**2)/constant_B)



if __name__ == '__main__':
    #logs_dir_1 = glob.glob(os.path.join("..",  LOG_DIR, SUB_DIR, CASE_DIR, SUB_CASE_DIR, "*.csv"))
    logs_dir_1 = glob.glob(os.path.join("..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_RX_DIVERS", "sigma_0dB_50per_p_0.002", "*.csv"))


    # logs_dir_2 = glob.glob(os.path.join("..",  LOG_DIR, SUB_DIR, CASE_DIR, "bs_0.01_p_0.002_R_25|100_nearest", "*.csv"))
    logs_dir_2 = glob.glob(os.path.join("..",  LOG_DIR, SUB_DIR, "fading_shadowing", "bs_0.01_p_0.002_R_25|100", "*.csv"))


    logs_dir_3 = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_RX_DIVERS", "sigma_8dB_50per_p_0.002_bs_0.01", "*.csv")
    )

    logs_dir_4 = glob.glob(os.path.join("..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_NST_ATT", "sigma_0dB_50per", "*.csv"))

    bs_nst_att_8dB = glob.glob(os.path.join("..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_NST_ATT", "sigma_8dB_50per", "*.csv"))


    sim_intensity, sim_plr_list = sim_data_process(logs_dir_1)
    sim_intensity_2, sim_plr_list_2 = sim_data_process(logs_dir_2)
    sim_intensity_3, sim_plr_list_3 = sim_data_process(logs_dir_3)
    sim_intensity_4, sim_plr_list_4 = sim_data_process(logs_dir_4)
    sim_intensity_5, sim_plr_list_5 = sim_data_process(bs_nst_att_8dB)

    sim_plr = [element[1] for element in sim_plr_list]
    sim_plr_lower = [element[1]-element[0] for element in sim_plr_list]
    sim_plr_upper = [element[2]-element[1] for element in sim_plr_list]


    sim_plr_2 = [element[1] for element in sim_plr_list_2]
    sim_plr_lower_2 = [element[1]-element[0] for element in sim_plr_list_2]
    sim_plr_upper_2 = [element[2]-element[1] for element in sim_plr_list_2]


    sim_plr_3 = [element[1] for element in sim_plr_list_3]
    sim_plr_lower_3 = [element[1]-element[0] for element in sim_plr_list_3]
    sim_plr_upper_3 = [element[2]-element[1] for element in sim_plr_list_3]

    sim_plr_4 = [element[1] for element in sim_plr_list_4]
    sim_plr_lower_4 = [element[1]-element[0] for element in sim_plr_list_4]
    sim_plr_upper_4 = [element[2]-element[1] for element in sim_plr_list_4]

    sim_plr_5 = [element[1] for element in sim_plr_list_5]
    sim_plr_lower_5 = [element[1]-element[0] for element in sim_plr_list_5]
    sim_plr_upper_5 = [element[2]-element[1] for element in sim_plr_list_5]

    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 190)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.01
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.002
    thetha_dB = 3.0 # unit dB
    theta = np.power(10, 3.0/10)
    mu_shadow = 0.0
    # shadowing effect 的标准差，单位分贝
    sigma_dB = 12
    BETA = np.log(10.0)/10.0
    # 需要对标准差做个转换，因为对数基底的不同
    sigma_G = BETA*sigma_dB
    sigma_X = 2.0*sigma_G/gamma

    constant_A = gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)*theta**(2.0/gamma)

    # Define p_f_2 as the outage probability over infinite plane
    p_f_rx_div_0 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0)


    p_f_bs_nst_att_12 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB)
    p_f_bs_nst_att_0 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0)
    p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8)

    #p_f_22 refers to case no shadowing, only fading.


    logs_dir_x = glob.glob(os.path.join("..",  LOG_DIR, "test", "*.csv"))

    sim_data_process_2(logs_dir_x)

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, sharey=False)

    for i in range(len(axes)):
        axes[i].set_yscale(SCALE[i])
        axes[i].plot(lambda_m, p_f_rx_div_0, color='k',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_RX_DIVERS ANA")
        axes[i].plot(lambda_m, p_f_bs_nst_att_8, color='g',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_NST_ATT, 8dB")
        axes[i].plot(lambda_m, p_f_bs_nst_att_0, color='g',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_NST_ATT, 0dB")


        axes[i].errorbar(
            sim_intensity,
            sim_plr,
            yerr=[sim_plr_lower, sim_plr_upper],
            fmt='*',
            ecolor='r',
            capthick=2,
            label="BS_RX_DIVERS SIM,0dB")
        axes[i].errorbar(
            sim_intensity_3,
            sim_plr_3,
            yerr=[sim_plr_lower_3, sim_plr_upper_3],
            fmt='d',
            ecolor='m',
            capthick=2,
            label="BS_RX_DIVERS SIM,8dB"
        )
        axes[i].errorbar(
            sim_intensity_4,
            sim_plr_4,
            yerr=[sim_plr_lower_4, sim_plr_upper_4],
            fmt='d',
            ecolor='m',
            capthick=2,
            label="BS_NST_ATT SIM,0dB"
        )
        axes[i].errorbar(
            sim_intensity_5,
            sim_plr_5,
            yerr=[sim_plr_lower_5, sim_plr_upper_5],
            fmt='d',
            ecolor='m',
            capthick=2,
            label="BS_NST_ATT SIM,8dB"
        )

        axes[i].grid()
        axes[i].axis([X_START, X_END, Y_START, Y_END])
        axes[i].set_xticks(np.arange(X_START, X_END, X_STEP))
        axes[i].set_title("Packet loss rate")
        axes[i].legend(loc='best', numpoints=2)
        axes[i].set_xlabel(r"Spatial density of devices")



    plt.savefig('packet_loss_rate_mpr.eps', format='eps', dpi=300)

    plt.show()
