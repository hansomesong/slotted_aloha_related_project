# -*- coding: utf-8 -*-
# The name of this script, sga is short for "stochastic geometry analytical model".
# A short name is convinient for import usage elsewhere
__author__ = 'qsong'

import numpy as np
import pandas as pd
from scipy.special import gamma as gamma_f
# from scipy.sepcial import erf as erf
import scipy.stats as st


# Some common constant declaration
EPS = np.finfo(np.float64).eps
BETA = np.log(10.0)/10.0

def bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where each device should attach to the nearest
    BS(i.e., Base Station) before transmission. The path-loss function is in form of simple law and has a singularity
    at the origin. The channel randomness effect we take into account fading and large-scale shadowing.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent, if gamma=0, no shadowing effect
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor

    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    constant_A = p*lambda_m*k*fm_fading*fm_shadowing*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)

    return 1.0 - 1.0/(1+np.power(constant_B, B_power))

def bs_best_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where each device should attach to the best
    BS(i.e., Base Station) before transmission. The path-loss function is in form of power-law and has a singularity
    at the origin. The channel randomness effect we take into account Rayleigh fading and large-scale shadowing.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar, path-loss exponent, if gamma=0, no shadowing effect
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor

    fm_fading = gamma_f(1+2.0/gamma)
    constant_A = p*lambda_m*k*fm_fading*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)

    return 1.0 - 1.0/(1+np.power(constant_B, B_power))

def bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where each device employs a 'fire and forget'
    emission strategy. All the BSes are the potential receiver for the considered device.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    variant_1 = p*lambda_m/lambda_b
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor

    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    # In this case, shadowing effect related term is not present in the constant_A...
    constant_A = p*lambda_m*k*fm_fading*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    # B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)
    # return np.exp(1.0/(1+np.power(constant_B, B_power))-np.exp(0.5*(2.0*sigma/gamma)**2)/constant_B)
    return np.exp(-1.0/constant_B)


def bs_rx_div_mrc_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where
    1) each device employs a 'fire and forget' emission strategy;
    2) Maximum Ratio Combining is used.
    All the BSes are the potential receiver for the considered device.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma



    variant_1 = p*lambda_m/lambda_b
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor

    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    # In this case, shadowing effect related term is not present in the constant_A...
    constant_A = p*lambda_m*k*fm_fading*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    # B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)
    # return np.exp(1.0/(1+np.power(constant_B, B_power))-np.exp(0.5*(2.0*sigma/gamma)**2)/constant_B)
    return np.exp(-1.0/constant_B)

def bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where each device should attach to the nearest
    BS(i.e., Base Station) before transmission. The path-loss function is in form of simple law and has a singularity
    at the origin. The channel randomness effect we take into account fading and large-scale shadowing.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent, if gamma=0, no shadowing effect
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor

    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    constant_A = p*lambda_m*k*fm_fading*fm_shadowing*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)

    return p*lambda_m*(1.0/(1+np.power(constant_B, B_power)))/lambda_b


def bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):

    return p*lambda_m*(1-bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure, itf_mean))/lambda_b

def sim_data_process(folder_dir):
    sim_intensity = []
    sim_plr = []
    for csv_file in folder_dir:
        csv_df = pd.read_csv(csv_file, sep=',', skiprows=[0], header=None, dtype=np.float64)
        # Hard code is not good...
        plr_df = csv_df.values[1:, -2]
        # sort the obtained to remove the max and min
        plr_df.sort()
        # remove the max and min
        # plr_df = plr_df[1:-1]
        avg_plr = plr_df.mean()
        # Calculates the standard error of the mean (or standard error of measurement)
        # Here we have to be careful to keep all y values positive:
        ci_min_plr = max(avg_plr-1.96*st.sem(plr_df), 1e-7)
        ci_max_plr = avg_plr + 1.96*st.sem(plr_df)
        alpha = csv_df.values[:, 1][0]
        sim_intensity.append(alpha)
        sim_plr.append([ci_min_plr, avg_plr, ci_max_plr])
    return np.array(sim_intensity), np.array(sim_plr)

def div_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=True):
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor
    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)
    return np.power(k*fm_fading*np.power(THETA, 2.0/gamma)*np.log(1.0/p_max)/np.pi, -1.0)


def nst_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=True):
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor
    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)
    # print "fm_shadowing", fm_shadowing
    B_power = np.power(1+np.pi*sigma_X**2/8.0, 0.5)
    return np.power(k*fm_fading*fm_shadowing*np.power(THETA, 2.0/gamma)/np.pi, -1.0)*np.power(p_max/(1.0-p_max), B_power)

def best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=True):
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor
    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, 0.5)
    return np.power(k*fm_fading*np.power(THETA, 2.0/gamma)/np.pi, -1.0)*np.power(p_max/(1.0-p_max), B_power)


def macro_div_gain(p_f_max, target, sigma_dB, gamma, trans_rep_nb, max_trans_nb):
    # Note that p_f_max cannot be zero, but to facilitate the declaration of ndarray. we add a very small positive
    # constant.
    # p_f_max here refers to the maximum allowed packet loss rate
    # p_f_once referst to the packet transmission failure probability 
    p_f_once = np.power(p_f_max, 1.0/max_trans_nb)
    p_f_max += EPS
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    C = np.power(1+0.5*np.pi*sigma**2/gamma**2, -0.5)

    avg_trans_nb = (1-p_f_max)/(1-p_f_once)

    trans_nb_ratio = trans_rep_nb/avg_trans_nb

    gain = np.power((1-p_f_once)/p_f_once, C**(-1))/np.log(1.0/p_f_once)/trans_nb_ratio
    if target == "BEST":
        return gain
    elif target == "NEAREST":
        return np.exp(2*sigma**2/gamma**2)*gain

def num_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    '''
        Numerically calculate the packet loss rate in case of Maximum Ratio Combining.
    '''
    thrld = np.power(10, thetha_dB/10.0)
    omega_end = 50
    N = 2000
    omega = np.linspace(0, omega_end, N)
    eta = 1.0
    # Why I should add an extremly small positive value to omega?
    EPS = np.finfo(float).eps
    shift_omega = 1j*eta + omega + EPS

    result = []

    for element in lambda_m:
        Y = ((1+EPS*1j)/(eta-1j*omega))*np.exp(-1j*omega*thrld)*total_sinr_cf(shift_omega, element, lambda_b, gamma, p)
        result.append(np.exp(eta*1.0*thrld)*np.real(np.trapz(y=Y, x=omega))/np.pi)
    # print "proba:", 1-result

    return np.array(result)

def total_sinr_cf(omega, lambda_m, lambda_b, gamma, p, pure=False, itf_mean=True):
    # A matrix of dimension len(lambda_m) * len(omega) will be returned.
    return np.exp(-1*lambda_b*np.pi*frac_moment_calculator(lambda_m, p, gamma)*gamma_f(1-2/gamma)*np.exp(-1j*np.pi/gamma)*np.power(omega, 2/gamma))


def frac_moment_calculator(lambda_m, p, gamma):
    # TODO: The following is just for the special case where gamma = 4. Not a general formula
    return 2*np.power(np.pi, -2)*np.power(p*lambda_m, -1)





if __name__ == "__main__":
    X_START = 0.0
    X_END = 0.3
    X_STEP = 0.002
    Y_START = 1e-3
    Y_END = 0.6
    Y_STEP = 0.1

    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 100)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.004
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.04
    thetha_dB = 3.0 # unit dB

    p_f_bs_nst_att_0 = bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0)
    p_f_bs_bst_att_8 = bs_best_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 10)

    # print lambda_m
    # print p_f_bs_nst_att_0
    # print p_f_bs_bst_att_8


    print "Slotted Aloha"
    sigma_dB = 0.0
    pure=False
    print div_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=False)
    print nst_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=True)
    print best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=True)
    print "Non slotted Aloha, avg interference"
    print div_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=True)
    print nst_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=True)
    print best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=True)
    print "Non slotted Aloha, max interference"
    print div_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=False)
    print nst_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=False)
    print best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=False)



