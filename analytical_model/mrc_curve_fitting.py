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
from scipy.stats import linregress

import pandas as pd
import glob
import os
import re
import matplotlib.pyplot as plt


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
    K,  B = fit_para[0][0], fit_para[0][1]

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
    This function use simulation results to fit the relationship between load and packet loss rate.
    The hard coded parameters are obtained from script mrc_curve_fitting_plot.py. Should consider to modify
    """
    sim_table_content = pd.read_csv(sim_fit_log)
    thetha = np.power(10.0, thetha_dB/10.0)

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

def sim_fitted_from_plr_load_function(tgt_plr, thetha_dB, gamma, pure, itf_mean):
    """ fitted function from simulations...
    :param tgt_plr: target packet loss rate, we need to find the corresponding load
    :param thetha_dB:
    :param load:
    :param gamma:
    :param pure:
    :param itf_mean:
    :return:
    """
    thetha = np.power(10.0, thetha_dB/10.0)

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

    # update: 25/09/2018: add the function from target loss rate to load
    # Target loss rate is temporaily set as 0.01,
    tgt_load = (np.power(erfcinv(tgt_plr), -1) - B)/K/thetha**(2.0/gamma)
    return tgt_load

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

def sim_parser(sim_log_dir, pure, itf_mean):
    if pure == False:
        SC = 0 # actually the row with number 1 in csv file
        MRC = 12
    else:
        # Case of pure ALOHA
        if itf_mean == True:
            # case of pure ALOHA with mean interference
            SC = 4 # actually the row with number 5 in csv file
            MRC = 13
        else:
            SC = 9
            MRC = 14

    sim_plr_sc_divers ={}
    sim_plr_sc_divers_semi_ci ={}
    sim_plr_mrc_divers ={}
    sim_plr_mrc_divers_semi_ci ={}
    # Simulation results from Xavier
    sim_intensity = {}

    for filename in glob.iglob(os.path.join(sim_log_dir, 'Aloha*.csv')):
        # print filename
        pattern = re.compile(r"AlohaMultg(\d+)s(\d+)t(\d+)")
        simu_setting = re.findall(pattern, filename)[0]
        label = simu_setting[0]
        # gamma = float(simu_setting[0])/10.0
        # sigma_dB = float(simu_setting[1])
        # theta_dB = float(simu_setting[2])
        # print gamma, sigma_dB, theta_dB
        plr_df = pd.read_csv(filename, sep=";", decimal=',', skiprows=9, nrows=15).dropna(axis=1, how='all')
        semi_ci = pd.read_csv(filename, sep=";", decimal=',', skiprows=26, nrows=15).dropna(axis=1, how='all')

        sim_intensity[label] = [float(x.replace(',', '.')) for x in plr_df.columns[1:].values.tolist()]

        sim_plr_sc_divers[label] = plr_df.loc[SC, plr_df.columns != 'Charge'].values
        sim_plr_sc_divers_semi_ci[label] = semi_ci.loc[SC, semi_ci.columns != 'Charge'].values

        sim_plr_mrc_divers[label] = plr_df.loc[MRC, plr_df.columns != 'Charge'].values
        sim_plr_mrc_divers_semi_ci[label] = semi_ci.loc[MRC, semi_ci.columns != 'Charge'].values

    return sim_intensity, sim_plr_sc_divers, sim_plr_sc_divers_semi_ci, sim_plr_mrc_divers, sim_plr_mrc_divers_semi_ci


def sim_get_load_interval(sim_log_dir, tgt_plr, pure, itf_mean):
    """
        Given a target packet loss rate and simulation result log file, this method returns the minimum load interval
        that covers the load corresponding to the given target packet loss rate.
        For example, the given target packet loss rate is 10%, the simulation results around the 10% may be 9.5% (with
        loda 0.051) and 10.1% (0.062), the method is able to retrieve the load interval [0.051, 0.062].
    :return:
    """
    if pure == False:
        SC = 0 # actually the row with number 1 in csv file
        MRC = 12
    else:
        # Case of pure ALOHA
        if itf_mean == True:
            # case of pure ALOHA with mean interference
            SC = 4 # actually the row with number 5 in csv file
            MRC = 13
        else:
            SC = 9
            MRC = 14

    sim_plr_sc_divers ={}
    sim_plr_mrc_divers ={}
    # Simulation results from Xavier
    sim_intensity = {}
    mrc_load_interval = {}
    mrc_plr_interval = {}
    tgt_fitted_loads = {}
    for filename in glob.iglob(os.path.join(sim_log_dir, 'Aloha*.csv')):
        # print filename
        pattern = re.compile(r"AlohaMultg(\d+)s(\d+)t(\d+)")
        simu_setting = re.findall(pattern, filename)[0]
        label = simu_setting[0]
        # gamma = float(simu_setting[0])/10.0
        # sigma_dB = float(simu_setting[1])
        # theta_dB = float(simu_setting[2])
        # print gamma, sigma_dB, theta_dB
        plr_df = pd.read_csv(filename, sep=";", decimal=',', skiprows=9, nrows=15).dropna(axis=1, how='all')

        sim_intensity[label] = [float(x.replace(',', '.')) for x in plr_df.columns[1:].values.tolist()]

        sim_plr_sc_divers[label] = plr_df.loc[SC, plr_df.columns != 'Charge'].values

        sim_plr_mrc_divers[label] = plr_df.loc[MRC, plr_df.columns != 'Charge'].values

        # The obtained sim_plr_sc_divers[label], of type numpy ndarray, will be surely one dimension.

        # We try to find the minimum interval that covers the target packet loss rate. We need to identify the left and
        # right index
        idx_left = np.max(np.nonzero(sim_plr_mrc_divers[label] < tgt_plr))
        idx_right = np.min(np.nonzero(sim_plr_mrc_divers[label] >= tgt_plr))

        mrc_load_interval[label] = np.array([sim_intensity[label][idx_left], sim_intensity[label][idx_right]])
        mrc_plr_interval[label] = np.take(sim_plr_mrc_divers[label], [idx_left, idx_right])


        slop = linregress( mrc_plr_interval[label], mrc_load_interval[label])[0]
        intercept = linregress( mrc_plr_interval[label], mrc_load_interval[label])[1]

        tgt_fitted_loads[label] = slop*tgt_plr + intercept

        # print label, mrc_load_interval[label], mrc_plr_interval[label], tgt_fitted_loads[label]
        print label, slop, intercept

    return tgt_fitted_loads


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
        print 'gamma:', gamma, empirical_p_f_mrc_div[label]
        p_f_mrc_div[label] = []

        for l in L:
            p_f_mrc_div[label].append(lt2plr(thetha_dB, l, gamma, pure, itf_mean))

        p_f_mrc_div[label] = np.array(p_f_mrc_div[label])

    return p_f_rx_div, p_f_mrc_div, fit_p_f_mrc_div, empirical_p_f_mrc_div


if __name__ == '__main__':
    # To test method: sim_get_load_interval

    thetha_dB = 6.0

    gammas = [3.0, 3.3, 3.5, 3.7, 4.0, 4.2, 4.5, 4.7, 5.0, 6.0]
    SIM_LOG_DIR = '/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuNov23PlusdePoints'
    tgt_plr = 0.1
    pure, itf_mean = True, True
    sim_dict = sim_get_load_interval(SIM_LOG_DIR, tgt_plr, pure, itf_mean)

    fitted_dict = {}
    for gamma in gammas:
        label = str(int(gamma*10))
        fitted_dict[label] = sim_fitted_from_plr_load_function(tgt_plr, thetha_dB, gamma, pure, itf_mean)
        print gamma, fitted_dict[label]


    print fitted_dict


    final =  pd.DataFrame([fitted_dict, sim_dict], index=["fitted", "simulation"], columns = ["33", "35", "37", "40", "45"]).transpose()

    err_per = np.abs(100*(final["fitted"] -final["simulation"])/final["simulation"])

    final["error(%)"] = err_per


    decimal = pd.Series([4,4,2], index=["fitted", "simulation", "error(%)"])
    print final.round(decimal)


    # Update 03-02-2019
    print "Update 03-02-2019..."
    SIM_LOG_DIR_2 = '/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuApril23'
    SIM_LOG_DIR_3 = '/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuNov23PlusdePoints'

    # sim_get_load_interval(sim_log_dir=SIM_LOG_DIR_3, tgt_plr=0.1, pure=True, itf_mean=True)

    # 03-02-2019: why the parameters for linear regression obtained are different from the contents in file: sim_fit_result_theta_3.csv?
    loads, plrs, X_ref, Y_ref, K, B, fit_func, fit_std = sim_curve_fitting("/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuNov23PlusdePoints/AlohaMultg37s8t3.csv",
                      thetha_dB=3.0, sigma_dB=8, gamma=3.7, pure=True, itf_mean=True, plr_min=0.001, plr_max=0.1)

    print K, B
    SIM_LOG = os.path.join(SIM_LOG_DIR_2, "AlohaMultg40s8t3b20.csv")

    plr_df = pd.read_csv(SIM_LOG, sep=";", decimal=',', skiprows=9, nrows=15).dropna(axis=1, how='all')

    # print plr_df












