# -*- coding: utf-8 -*-
# This script collects all involved method for maximum ratio combining and curve fitting.

__author__ = 'qsong'

from analytical_model import sgam
import numpy as np
from scipy.special import gamma as gamma_f
from scipy.special import erf as erf # Import error function
from scipy.special import erfc as erfc # Import error function
from scipy.special import erfcinv as erfcinv # Import error inverse function
from scipy.special import binom as binom
from scipy.optimize import leastsq  # 引入最小二乘法算法

import pandas as pd
import json
import glob
import os
import re
import matplotlib.pyplot as plt
import pprint



SIM_FILE_PATH = "/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuAvecFittageDifTetas"


def cumu_sir_lt(s, load, gamma, pure=False, itf_mean=True):
    """
    This method is used to calculate the value of Laplace transform, given a complex value s, for the output SIR when
    MRC technique is applied.
    The closed-form expression: L[s] = exp(-[\Gamma(1 - 2.0/\gamma)]/[A L] s^{2/\gamma})
    :param s:
    :param load:
    :param gamma:
    :param pure:
    :param itf_mean:

    :return:
    """
    A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    if pure:
        if itf_mean:
            A = 2*gamma/(2+gamma) * A
        else:
            A = 2.0 * A

    return np.exp(-gamma_f(1-2.0/gamma)*np.power(s, 2.0/gamma)/A/load)
    # return np.exp(-gamma_f(1-2.0/gamma)*(2.0/gamma)**(1.0/5)*np.power(s, 2.0/gamma)/A/load)



def lt2plr(thetha_dB, load, gamma, pure=False, itf_mean=True):
    '''
    This method is used to numerically calculate packet loss rate from laplace transform.
    The underlying formula originiates from Eq.(35) of reference: "Unified Stochastic Geometry Modeling and Analysis of
    Cellular Networks in LOS/NLOS and Shadowed Fading".

    :param thetha_dB:       scalar, in unit of dB, capture ratio threshold
    :param load:
    :param gamma:
    :param pure:
    :param itf_mean:
    :return:
    '''
    T = np.power(10.0, thetha_dB/10.0)
    # To avoid any name collision with our namespace. For parameters used in above mentioned is added with append _num
    A_NUM = 18.4
    M_NUM = 11
    N_NUM = 15

    def real_part_lt_of_cdf_theta(x, load, gamma, pure, itf_mean):
        return np.real(
            cumu_sir_lt(x, load, gamma, pure, itf_mean)/x
        )

    binom_coeff = np.array([binom(M_NUM, k) for k in range(M_NUM + 1)])
    s_n_k_t_list = []
    for k in range(M_NUM+1):
        a_coeff = np.exp(0.5*A_NUM) * np.array(
            [
                np.power(2*T, -1.0)
                    if l == 0
                    else
                np.power(T, -1.0)
                    for l in range(N_NUM + k + 1)
            ]
        )

        b_coeff = np.array(
            [
                np.power(-1.0, l)*
                real_part_lt_of_cdf_theta(
                    (A_NUM+2j*l*np.pi)/2/T,
                    load, gamma, pure, itf_mean
                )
                for l in range(N_NUM + k + 1)
            ]
        )

        s_n_k_t_list.append(np.sum(a_coeff * b_coeff))

    s_n_k_t = np.array(s_n_k_t_list)
    result = np.sum(binom_coeff * s_n_k_t) * np.power(2.0, -1.0*M_NUM)

    return result



def numercal_invert_plr2load(plr, thetha_dB, gamma, pure=False, itf_mean=True):
    """
    This method is used to numerically invert a given packet loss rate to a normalized load.
    :return:
    """

    load_min = 0.0
    load_max = 10.0
    tmp_plr = 0.0
    load = 0.0

    N_MAX = 500
    N = 0
    precision = 0.0000001

    while np.abs(tmp_plr-plr)/plr > precision and N < N_MAX:
        load = (load_min + load_max) / 2.0
        tmp_plr = lt2plr(thetha_dB, load, gamma, pure, itf_mean)
        if tmp_plr > plr:
            # Packet loss rate should be monotone increasing with respect to normalized load
            # if tmp_plr is greater than the target packet loss rate
            # uppuer bound of seaching interval should be updated
            load_max = load
        else:
            load_min = load

        N += 1


    return load, N



# Methods for curve fitting for Maximum Ratio Combining techniques
## The target function to fit.
## Given that
## 1) when path-loss exponent \gamma is 4.0, the packet loss rate when MRC is applied is a function of erfc function
## 2) Laplace transform of output SIR in case of MRC is in terms of exponential function and depends
##      on A*L*theta^{2/\gamma}
#
def func(para, x):
    k, b = para
    return k*x + b

## 偏差函数：x,y都是列表:这里的x,y更上面的Xi,Yi中是一一对应的
def error(para, x, y):
    return func(para, x) - y

def curve_fitting(thetha_dB, gamma, pure, itf_mean, plr_min=0.001, plr_max=0.1):

    plrs = np.linspace(plr_min, plr_max, 50)
    Y_ref = np.power(erfcinv(plrs), -1.0)
    loads = []

    for plr in plrs:
        load = numercal_invert_plr2load(plr, thetha_dB, gamma, pure, itf_mean)[0]
        loads.append(load)

    loads = np.array(loads)
    X_ref = loads*np.power(np.power(10.0, thetha_dB/10.0), 2.0/gamma)


    fit_para = leastsq(error, (1.0, 1.0), args=(X_ref, Y_ref))

    # print fit_para
    K, B = fit_para[0][0], fit_para[0][1]

    def fit_func(loads, thetha_dB, gamma):
        return erfc(1.0/(K*loads*np.power(np.power(10.0, thetha_dB/10.0), 2.0/gamma)+B))
    plr_fit = fit_func(loads, thetha_dB, gamma)

    fit_std = np.sqrt(np.sum(np.power((plrs - plr_fit)/plrs, 2.0))/loads.size)

    return loads, plrs, X_ref, Y_ref, K, B, fit_func, fit_std


def sim_curve_fitting(file_name, thetha_dB, sigma_dB, gamma, pure, itf_mean, plr_min=0.001, plr_max=0.1):

    print "Target simulation log file: ", file_name

    if pure == False:
        MRC = 12
    else:
        # Case of pure ALOHA
        if itf_mean == True:
            # case of pure ALOHA with mean interference
            MRC = 13
        else:
            MRC = 14

    plr_df = pd.read_csv(file_name, sep=";", decimal=',', skiprows=9, nrows=15).dropna(axis=1, how='all')
    #remove the first column and we get a Sereis
    plrs = plr_df.loc[MRC, plr_df.columns != 'Charge']
    plrs = plrs[plrs > plr_min]
    # Choose those columns that value is greater than a plr_min, e.g, 0.001
    loads =  pd.to_numeric(plrs.index.str.replace(',', '.'), errors='coerce')

    Y_ref = np.power(erfcinv(plrs), -1.0)

    X_ref = loads*np.power(np.power(10.0, thetha_dB/10.0), 2.0/gamma)

    fit_para = leastsq(error, (1.0, 1.0), args=(X_ref, Y_ref))

    # print fit_para
    K, B = fit_para[0][0], fit_para[0][1]

    def fit_func(loads, thetha_dB, gamma):
        return erfc(1.0/(K*loads*np.power(np.power(10.0, thetha_dB/10.0), 2.0/gamma)+B))
    plr_fit = fit_func(loads, thetha_dB, gamma)

    fit_std = np.sqrt(np.sum(np.power((plrs - plr_fit)/plrs, 2.0))/loads.size)

    return loads, plrs, X_ref, Y_ref, K, B, fit_func, fit_std


def fitted_function(thetha_dB, load, gamma, pure=False, itf_mean=True):
    thetha = np.power(10.0, thetha_dB/10.0)
    ## The following part is estimated from gamma [3.5, 4.5]
    # B = 0.1473*gamma -0.5968
    # if pure == False:
    #     # slotted case
    #     K = -0.1745*gamma + 2.4787,
    # else:
    #     if itf_mean == True:
    #         K = -0.0324*gamma + 2.4989
    #     else:
    #         K = -0.3489*gamma + 4.9574
    # # print "curve fitting fucntion parameter:", K, B
    # result =list(1 - erf(np.power(K*load*thetha**(2.0/gamma)+B, -1.0)))
    # result = [0 if element > 1 else element for element in result]
    # return np.array(result)


    # B = 0.1183*gamma - 0.5156
    #
    # if pure == False:
    #     # slotted case
    #     K = -0.1412*gamma + 2.3843
    #
    # else:
    #     if itf_mean == True:
    #         K = -0.01693*gamma + 2.4650
    #     else:
    #         K = -0.2824*gamma + 4.7685
    # # print "curve fitting fucntion parameter:", K, B
    # result =list(1 - erf(np.power(K*load*thetha**(2.0/gamma)+B, -1.0)))
    # result = [0 if element > 1 else element for element in result]
    # return np.array(result)

    # Polynomial fit
    B = 0.02488149*np.power(gamma, 3) - 0.39097563 * np.power(gamma, 2) + 2.0897039*np.power(gamma, 1) - 3.69313914

    if pure == False:
        # slotted case
        K = -0.02689164*np.power(gamma, 3) + 0.42429377*np.power(gamma, 2) - 2.2875156*np.power(gamma, 1) + 5.85222392


    else:
        if itf_mean == True:
            K = -0.0235603*np.power(gamma, 3) + 0.3591914*np.power(gamma, 2) - 1.78448733*np.power(gamma, 1) + 5.2594168

        else:
            K = -0.05378328*np.power(gamma, 3) + 0.84858754*np.power(gamma, 2) -4.57503119 *np.power(gamma, 1) + 11.70444784
    # print "curve fitting fucntion parameter:", K, B
    result =list(1 - erf(np.power(K*load*thetha**(2.0/gamma)+B, -1.0)))
    result = [0 if element > 1 else element for element in result]
    return np.array(result)


def sim_fitted_function(sim_fit_log, thetha_dB, load, gamma, pure, itf_mean):
    """ fitted function from simulations...
    :param thetha_dB:
    :param load:
    :param gamma:
    :param pure:
    :param itf_mean:
    :return:
    """
    sim_table_content = pd.read_csv(sim_fit_log)
    thetha = np.power(10.0, thetha_dB/10.0)
    ## linear regression
    # if pure == False:
    #     # slotted case
    #     K_label  = 'K_slot'
    #     B_label  = 'B_slot'
    # else:
    #     if itf_mean == True:
    #         K_label = 'K_pure_mean'
    #         B_label  = 'B_pure_mean'
    #     else:
    #         K_label = 'K_pure_max'
    #         B_label  = 'B_pure_max'
    #
    # K = sim_table_content.loc[sim_table_content.gamma == gamma, K_label].values[0]
    # B = sim_table_content.loc[sim_table_content.gamma == gamma, B_label].values[0]
    # print "Gamma:", gamma
    # print "K:", K
    # print "B:", B
    # # print "curve fitting fucntion parameter:", K, B
    # result =list(1 - erf(np.power(K*load*thetha**(2.0/gamma)+B, -1.0)))
    # result = [0 if element > 1 else element for element in result]
    # return np.array(result)

    # Polynomial fit

    if pure == False:
        # slotted case
        K = 0.0188573*np.power(gamma, 3) - 0.2517638*np.power(gamma, 2) + 1.0200586*np.power(gamma, 1) + 0.47600459
        B = -3.70488722e-04*np.power(gamma, 3) - 2.75401227e-02 * np.power(gamma, 2) + 3.94416600e-01*np.power(gamma, 1) - 1.06071310e+00

    else:
        if itf_mean == True:
            K = 0.01758788*np.power(gamma, 3) -0.23578095*np.power(gamma, 2) + 0.99829322*np.power(gamma, 1) + 0.9940888
            B = 0.00648887*np.power(gamma, 3) -0.13081056 * np.power(gamma, 2) + 0.91538891*np.power(gamma, 1) -1.93890062

        else:
            K = -0.00870517*np.power(gamma, 3) + 0.16127219*np.power(gamma, 2) -1.10255796  *np.power(gamma, 1) + 5.54168573
            B = 0.01945308*np.power(gamma, 3) -0.31958942 * np.power(gamma, 2) + 1.82832705*np.power(gamma, 1) -3.3912209

    # print "curve fitting fucntion parameter:", K, B
    result =list(1 - erf(np.power(K*load*thetha**(2.0/gamma)+B, -1.0)))
    result = [0 if element > 1 else element for element in result]
    return np.array(result)

def empirical_plr_mrc(thetha_dB, lambda_m, lambda_b, gamma, p, pure, itf_mean):
    """
        This formula has nothing to do with shadowing value
    """
    T = np.power(10.0, thetha_dB/10.0)

    L = p*lambda_m / lambda_b

    A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    if pure:
        if itf_mean:
            A = 2*gamma/(2+gamma) * A
        else:
            A = 2.0 * A

    return 1 - erf(gamma_f(1.0-2.0/gamma)/(2*A*L*T**(2.0/gamma)))

def sim_parser(sim_log_dir, pure, itf_mean, itf_independence):
    # UPDATE: 2018-01-22
    # The output data contains two sets of data for SC and MRC macro diversity:
    # with and without interference independence assumption
    SIM_CONFIG_PATH = os.path.join(sim_log_dir, 'sim_config.json')
    with open(SIM_CONFIG_PATH) as json_file:
        sim_config_dict = json.load(json_file)

    SKIP_NB_AVG = sim_config_dict['avg']
    SKIP_NB_CI = sim_config_dict['ci']
    # The number of entries in a table.
    CASE_NB = sim_config_dict['case_nb']

    # The dict for csv file row index for basic-based simulation, i.e., dependence interference, more close to reality
    B_SIM_MAP = sim_config_dict['dept_itf_case']
    # The dict for csv file row index for model-based simulation, i.e., independence interference
    M_SIM_MAP = sim_config_dict['idept_itf_case']


    SC = {'B_SIM':0, 'M_SIM':0}
    MRC = {'B_SIM': 0, 'M_SIM': 0}

    if pure == False:
        # Index number in csv file minus 1 is the index in pandas table
        SC['B_SIM'] = B_SIM_MAP['slotted_sc_with_shadowing'] - 1
        SC['M_SIM'] = M_SIM_MAP['slotted_sc_with_shadowing'] - 1 # actually the row with number 1 in csv file
        MRC['B_SIM'] = B_SIM_MAP['slotted_mrc_with_shadowing'] - 1
        MRC['M_SIM'] = M_SIM_MAP['slotted_mrc_with_shadowing'] - 1 # actually the row with number 13 in csv file

    else:
        # Case of pure ALOHA
        if itf_mean == True:
            # case of pure ALOHA with mean interference
            SC['M_SIM'] = M_SIM_MAP['pure_sc_with_shadowing_avg_itf'] - 1
            MRC['M_SIM'] = M_SIM_MAP['pure_mrc_with_shadowing_avg_itf'] - 1
            SC['B_SIM'] = B_SIM_MAP['pure_sc_with_shadowing_avg_itf'] - 1
            MRC['B_SIM'] = B_SIM_MAP['pure_mrc_with_shadowing_avg_itf'] - 1
        else:
            SC['M_SIM'] = M_SIM_MAP['pure_sc_with_shadowing_max_itf'] - 1
            MRC['M_SIM'] = M_SIM_MAP['pure_mrc_with_shadowing_max_itf'] - 1
            SC['B_SIM'] = B_SIM_MAP['pure_sc_with_shadowing_max_itf'] - 1
            MRC['B_SIM'] = B_SIM_MAP['pure_mrc_with_shadowing_max_itf'] - 1

    sim_plr_sc_divers ={'B_SIM':{}, 'M_SIM':{}}
    sim_plr_sc_divers_semi_ci ={'B_SIM':{}, 'M_SIM':{}}
    sim_plr_mrc_divers ={'B_SIM':{}, 'M_SIM':{}}
    sim_plr_mrc_divers_semi_ci ={'B_SIM':{}, 'M_SIM':{}}
    # Simulation results from Xavier
    sim_intensity = {}

    for filename in glob.iglob(os.path.join(sim_log_dir, 'Aloha*.csv')):
        # print filename
        pattern = re.compile(r"AlohaMultg(\d+)s(\d+)t(\d+)")
        simu_setting = re.findall(pattern, filename)[0]
        label = simu_setting[0]
        print "The extracted gamma is:", label
        # gamma = float(simu_setting[0])/10.0
        # sigma_dB = float(simu_setting[1])
        # theta_dB = float(simu_setting[2])
        # print gamma, sigma_dB, theta_dB
        plr_df = \
            pd.read_csv(filename, sep=";", decimal=',', skiprows=SKIP_NB_AVG, nrows=CASE_NB).dropna(axis=1, how='all')
        semi_ci = \
            pd.read_csv(filename, sep=";", decimal=',', skiprows=SKIP_NB_CI, nrows=CASE_NB).dropna(axis=1, how='all')

        sim_intensity[label] = [float(x.replace(',', '.')) for x in plr_df.columns[1:].values.tolist()]

        for key in SC.keys():
            sim_plr_sc_divers[key][label] = plr_df.loc[SC[key], plr_df.columns != 'Charge'].values
            sim_plr_sc_divers_semi_ci[key][label] = semi_ci.loc[SC[key], semi_ci.columns != 'Charge'].values
        for key in MRC.keys():
            sim_plr_mrc_divers[key][label] = plr_df.loc[MRC[key], plr_df.columns != 'Charge'].values
            sim_plr_mrc_divers_semi_ci[key][label] = semi_ci.loc[MRC[key], semi_ci.columns != 'Charge'].values

    return sim_intensity, sim_plr_sc_divers, sim_plr_sc_divers_semi_ci, sim_plr_mrc_divers, sim_plr_mrc_divers_semi_ci


def sc_mrc_anayltical_parser(lambda_m, lambda_b, p, thetha_dB, gammas, pure, itf_mean):
    # Analytical results for SC macro diversity
    p_f_rx_div ={}
    # Analytical results for MRC macro diversity
    p_f_mrc_div = {}
    fit_p_f_mrc_div = {}
    empirical_p_f_mrc_div = {}

    L = p*lambda_m/lambda_b

    for gamma in gammas:
        label = str(int(gamma*10))
        p_f_rx_div[label] = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8.0, pure, itf_mean)
        fit_p_f_mrc_div[label] = fitted_function(thetha_dB, L, gamma, pure, itf_mean)
        empirical_p_f_mrc_div[label] = empirical_plr_mrc(thetha_dB, lambda_m, lambda_b, gamma, p, pure, itf_mean)
        p_f_mrc_div[label] = []

        for l in L:
            p_f_mrc_div[label].append(lt2plr(thetha_dB, l, gamma, pure, itf_mean))

        p_f_mrc_div[label] = np.array(p_f_mrc_div[label])

        if pure == True and itf_mean == True:
            print 'thetha_dB', thetha_dB, 'gamma:', gamma, "pure", pure, "avg. interference", itf_mean
            print "load and p_f_mrc_div", pprint.pprint(zip(L, p_f_mrc_div[label]))

    return p_f_rx_div, p_f_mrc_div, fit_p_f_mrc_div, empirical_p_f_mrc_div



if __name__ == "__main__":

    gamma = 3.7
    thetha_dB = 6.0

    pure = [False, True, True]
    itf_mean = [True, True, False]
    title = ['Slotted ALOHA', 'Pure ALOHA, mean', 'Pure ALOHA, max']
    x_interval = [[0.15, 0.3], [0.12, 0.22], [0.08, 0.15]]

    print 2*gamma_f(1+2.0/gamma), 2*2*gamma_f(1+2.0/gamma)

    # fig, axes = plt.subplots(1, 3, figsize=(8, 6))
    #
    # for i, ax in enumerate(fig.axes):
    #     loads, plrs, X_ref, Y_ref, K, B, fit_func, fit_std = curve_fitting(thetha_dB, gamma, pure[i], itf_mean[i], plr_min=0.001, plr_max=0.1)
    #     print i, K, B, fit_std
    #
    #     # ax.scatter(loads, plrs, c='r')
    #     # ax.plot(loads, fit_func(loads, thetha_dB, gamma), 'b')
    #     ax.scatter(X_ref, Y_ref, c='r')
    #     ax.plot(X_ref, K*X_ref+B, 'b')
    #     ax.set_yscale('log')
    #     ax.set_title(title[i])
    #     # ax.set_xlim(x_interval[i])
    #     # ax.set_ylim([0.001, 0.1])
    #     ax.grid()


    fig, axes = plt.subplots(1, 2, figsize=(8, 6))

    TITLES = ['Curve fitting for parameter K','Curve fitting for parameter B']


    gammas = np.array([3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5])
    g_f = gamma_f(1+2.0/gammas)
    K_list_slot = np.array([1.88265070202, 1.85529426, 1.83105926055, 1.80943155043, 1.79000474, 1.77245383, 1.75651528, 1.74197337, 1.72864972, 1.7163953, 1.70508502])
    K_list_pure_mean = np.array([2.39610084685, 2.38537835, 2.37716473221, 2.370979243, 2.36644699, 2.36327175, 2.36121726, 2.36009297, 2.35974399, 2.36004355, 2.36088691])
    K_list_pure_max = np.array([3.76530140568, 3.71058852536, 3.66211852114, 3.61886310213, 3.58000949, 3.54490770, 3.51303055,3.48394674, 3.45729944, 3.43279061, 3.41017003])
    B_list = np.array([-0.0944713161168, -0.07071986, -0.0498277441762, -0.0313163457645, -0.01480779, 0, 0.0133516, 0.02544702, 0.03645134,  0.04650218, 0.05571503])
    # B_list = np.exp(B_list)
    def func(para, x):
        k, b = para
        return k*x + b
        ## 偏差函数：x,y都是列表:这里的x,y更上面的Xi,Yi中是一一对应的
    def error(para, x, y):
        return func(para, x) - y

    K_fit_para_slot = leastsq(error, (1.0, 1.0), args=(gammas, K_list_slot))
    K_fit_para_pure_mean = leastsq(error, (1.0, 1.0), args=(gammas, K_list_pure_mean))
    K_fit_para_pure_max = leastsq(error, (1.0, 1.0), args=(gammas, K_list_pure_max))

    B_fit_para = leastsq(error, (1.0, 1.0), args=(gammas, B_list))



    K_fit_slot = K_fit_para_slot[0][0]*gammas + K_fit_para_slot[0][1]
    K_fit_pure_mean = K_fit_para_pure_mean[0][0]*gammas + K_fit_para_pure_mean[0][1]
    K_fit_pure_max = K_fit_para_pure_max[0][0]*gammas + K_fit_para_pure_max[0][1]


    B_fit = B_fit_para[0][0]*gammas + B_fit_para[0][1]

    fit_std1 = np.sqrt(np.sum(np.power((K_list_slot- K_fit_slot)/K_list_slot, 2.0))/K_list_slot.size)
    fit_std2 = np.sqrt(np.sum(np.power((B_list- B_fit)/B_list, 2.0))/B_list.size)

    print "slot ", K_fit_para_slot, fit_std1
    print "pure mean ", K_fit_para_pure_mean, fit_std1
    print "pure max ", K_fit_para_pure_max, fit_std1
    print B_fit_para, fit_std1
    # for i, ax in enumerate(fig.axes):
    # ax.scatter(gammas, K_list_slot)
    # ax.plot(gammas, K_fit_slot)
    Y_LABELS =['K', 'B']
    for i, ax in enumerate(fig.axes):
        if i == 0:
            ax.scatter(gammas, K_list_slot)
            ax.plot(gammas, K_fit_slot)
        if i == 1:
            ax.scatter(gammas, B_list)
            ax.plot(gammas, B_fit)

        ax.set_xlabel(r"gamma")
        ax.set_ylabel(Y_LABELS[i])
        ax.set_title(TITLES[i])



    plt.savefig()
    plt.show()






