# -*- coding: utf-8 -*-
# The name of this script, sga is short for "stochastic geometry analytical model".
# A short name is convinient for import usage elsewhere
__author__ = 'qsong'

import numpy as np
import pandas as pd
from scipy.special import gamma as gamma_f
from scipy.special import erf as erf
from scipy.special import binom as binom
import scipy.stats as st
import scipy.integrate as integrate



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

def bs_nearest_atch_op2(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
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

    A = fm_fading * fm_shadowing * np.power(THETA, 2.0/gamma)
    B = A*k*p*lambda_m/(lambda_b*np.pi)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)

    return 1.0 - 1.0/(1+np.power(B, B_power))

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


def bs_best_atch_op_with_noise(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, noise_power, pure=False, itf_mean=True):
    """
    This method is used to calculate the packet loss rate for the case where each device should attach to the best
    BS(i.e., Base Station) before transmission. Background noise is taken into account.
    The path-loss function is in form of power-law and has a singularity at the origin.
    The channel randomness effect we take into account Rayleigh fading and large-scale shadowing.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar, path-loss exponent, if gamma=0, no shadowing effect
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param noise_power: scalar, the normalized background noise power level,unit dBm.
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    noise_power = np.power(10, (noise_power - 30.0)/10)

    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)*factor

    U = THETA * noise_power
    V = fm_shadowing*np.pi*(p*lambda_m*A*np.power(THETA, 2.0/gamma) + lambda_b)

    return 1.0 - np.pi * lambda_b * fm_shadowing * 0.5 * np.sqrt(np.pi/U) * np.exp(V**2/4/U) * (1- erf(V/2/np.sqrt(U)))

def bs_best_atch_op_with_noise_trap_rule(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, noise_power, pure=False, itf_mean=True):
    """
    This method is used to calculate the packet loss rate for the case where each device should attach to the best
    BS(i.e., Base Station) before transmission. Background noise is taken into account.
    The path-loss function is in form of power-law and has a singularity at the origin.
    The channel randomness effect we take into account Rayleigh fading and large-scale shadowing.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar, path-loss exponent, if gamma=0, no shadowing effect
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param noise_power: scalar, the normalized background noise power level,unit dBm.
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    noise_power = np.power(10, (noise_power - 30.0)/10)

    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)*factor

    U = THETA * noise_power
    print "V", fm_shadowing*np.pi*(p*0*A*np.power(THETA, 2.0/gamma) + lambda_b)
    result = []
    for m in lambda_m:

        V = fm_shadowing*np.pi*(p*m*A*np.power(THETA, 2.0/gamma) + lambda_b)

        integral = 1 - np.pi*lambda_b*fm_shadowing*integrate.quad(lambda r: np.exp(-U*r**2-V*r), 0, 10000)[0]

        result.append(integral)

    return np.array(result)


def new_bs_rx_div_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where each device should attach to the best
    BS(i.e., Base Station) before transmission. The path-loss function is in form of power-law and has a singularity
    at the origin. The channel randomness effect we take into account Rayleigh fading and large-scale shadowing.
    This formula is obtained from reference "optimal SINR-based coverage in Poisson Cellular networks with Power Density
    Constraints"
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


def another_bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
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
    fm_shadowing = np.exp(0.5*np.power(sigma_X, 2))

    # In this case, shadowing effect related term is not present in the constant_A...
    constant_A = p*lambda_m*k*fm_fading*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)
    # return np.exp(1.0/(1+np.power(constant_B, B_power))-np.exp(0.5*(2.0*sigma/gamma)**2)/constant_B)
    return np.exp(1.0/(1+np.power(constant_B*fm_shadowing, B_power))-1.0/constant_B)

def bs_rx_div_mrc_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where
    1) each device employs a 'fire and forget' emission strategy;
    2) Maximum Ratio Combining is used.
    3) This method is only valid for gamma = 4.0!!!
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
    if gamma != 4.0:
        exit("This method is only valide for gamma = 4.0")
    THETA = np.power(10, thetha_dB/10)

    if not pure:
        A = 1.0
    else:
        if itf_mean:
            A = 4.0/3
        else:
            A = 2.0

    return 1-erf(np.power(A*p*lambda_m/lambda_b*np.sqrt(np.pi*THETA), -1))

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

def bs_rx_div_mrc_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):

    return p*lambda_m*(1-bs_rx_div_mrc_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure, itf_mean))/lambda_b

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

def div_max_load(gamma, thetha_dB, p_max=0.1, pure=False, itf_mean=True):
    THETA = np.power(10, thetha_dB/10.0)
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
    # sigma_X = 2.0*sigma/gamma
    sigma_X = 0

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



def mrc_bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    '''
        Numerically calculate the packet loss rate in case of Maximum Ratio Combining. A general numerical framework
    '''
    theta = np.power(10.0, thetha_dB/10.0)
    omega_end = 100.0
    N = 20000
    omega = np.linspace(0, omega_end, N)
    eta = 1.0
    shift_omega = 1j*eta + omega

    result = []

    for element in lambda_m:
        Y = np.exp(-1j*omega*theta)*total_sinr_cf(
            shift_omega, element, lambda_b, gamma, p, sigma_dB, pure, itf_mean
        )/(eta-1j*omega)

        plr = np.exp(eta*1.0*theta)*np.real(np.trapz(y=Y, x=omega))/np.pi
        print "pure ALOHA:", pure, ", lambda_m:", element, ", packet loss rate:", plr, \
            ", ref:", 1-erf(0.75*np.power(p*element/lambda_b*np.sqrt(np.pi*theta), -1))
        result.append(plr)

    return np.array(result)


def total_sinr_cf(omega, lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean):
    """
    :param omega:
    :param lambda_m:
    :param lambda_b:                    scalar, spatial BS density
    :param gamma:                       scalar, path-loss exponent
    :param p:                           scalar, transmission probability
    :param pure:
    :param itf_mean:
    :return:
    """
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    C = np.pi*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)
    return np.exp(
        -lambda_b*fm_shadowing*C*itf_frac_moment_calculator(lambda_m, p, sigma_dB, gamma, pure, itf_mean)*np.exp(-1j*np.pi/gamma)*np.power(omega, 2/gamma)
    )

    # return np.exp(
    #     -lambda_b*fm_shadowing*C*itf_frac_moment_calculator(lambda_m, p, sigma_dB, gamma, pure, itf_mean)*np.power(-1j*omega, 2/gamma)
    # )



def frac_moment_calculator(lambda_m, p, gamma):
    # TODO: The following is just for the special case where gamma = 4. Not a general formula
    return 2*np.power(np.pi, -2)*np.power(p*lambda_m, -1)


def itf_frac_moment_calculator(lambda_m, p, sigma_dB, gamma, pure, itf_mean):
    # TODO: The following is just for the special case where gamma = 4. Not a general formula
    # 2017-10-12. I have found the general expression for fractional moment calculation problem.
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB

    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma

    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)
    A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    if pure:
        if itf_mean:
            A = 2*gamma/(2+gamma) * A
        else:
            A = 2.0 * A

    return gamma*np.power(2*gamma_f(2.0/gamma)*p*lambda_m*np.pi*fm_shadowing*A, -1.0)

def cumu_itf_cf(omega, p, lambda_m, gamma, sigma_dB, pure, itf_mean):

    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma

    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)
    A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    if pure:
        if itf_mean:
            A = 2*gamma/(2+gamma) * A
        else:
            A = 2.0 * A

    return np.exp(
        -p*lambda_m*np.pi*A*fm_shadowing*np.exp(-1j*np.pi/gamma) * np.power(omega, 2.0/gamma)
    )


def cumu_sir_lt(s, lambda_m, lambda_b, gamma, p, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the value of Laplace transform, given a complex value s
    :param s:
    :param lambda_m:
    :param lambda_b:
    :param gamma:
    :param p:
    :param sigma_dB:
    :param pure:
    :param itf_mean:
    :return:
    """
    # Calculate the fractional moment of cumulative interference
    fm_cumu_itf = itf_frac_moment_calculator(lambda_m, p, sigma_dB, gamma, pure, itf_mean)

    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma

    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    C = np.pi*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    return np.exp(
        -lambda_b * fm_shadowing * C * fm_cumu_itf * np.power(s, 2.0/gamma)
    )


def lt2plr(thetha_dB, lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean):
    '''
    This method is used to numerically calculate packet loss rate from laplace transform.
    The underlying formula originiates from Eq.(35) of reference: "Unified Stochastic Geometry Modeling and Analysis of
    Cellular Networks in LOS/NLOS and Shadowed Fading".

    :param thetha_dB:       scalar, in unit of dB, capture ratio threshold
    :param lambda_m:        scalar,
    :param lambda_b:
    :param gamma:
    :param p:
    :param sigma_dB:
    :param pure:
    :param itf_mean:
    :return:
    '''
    T = np.power(10.0, thetha_dB/10.0)
    #
    # To avoid any name collision with our namespace. For parameters used in above mentioned is added with append _num
    A_num = 18.4
    m_num = 11
    n_num = 15

    def real_part_lt_of_cdf_theta(x, lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean):
        return np.real(
            cumu_sir_lt(x, lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean)/x
        )

    binom_coeff = np.array([binom(m_num, k) for k in range(m_num + 1)])
    s_n_k_t_list = []
    for k in range(m_num+1):
        a_coeff = np.exp(0.5*A_num) * np.array(
            [
                np.power(2*T, -1)
                    if l == 0
                    else
                np.power(T, -1)
                    for l in range(n_num + k + 1)
            ]
        )

        b_coeff = np.array(
            [
                np.power(-1.0, l)*
                real_part_lt_of_cdf_theta(
                    (A_num+2j*l*np.pi)/2/T,
                    lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean
                )
                for l in range(n_num + k + 1)
            ]
        )

        s_n_k_t_list.append(np.sum(a_coeff * b_coeff))

    s_n_k_t = np.array(s_n_k_t_list)
    result = np.sum(binom_coeff * s_n_k_t) * np.power(2.0, -1.0*m_num)

    return result


# def itf_cf2cdf(x, p, lambda_m, gamma, sigma_dB, pure, itf_mean):
#     """
#     :param x:                       ndarray,    list of interference values
#     :param p:                       scalar,     transmission probability
#     :param lambda_m:                scalar,     devices spatial density
#     :param gamma:                   scalar,     path-loss exponent
#     :param thetha_dB:
#     :param sigma_dB:
#     :param pure:
#     :param itf_mean:
#     :return:
#     """
#     EPS = np.finfo(float).eps
#     omega_end = 500.0
#     N = 5000
#     omega = np.linspace(0.0, omega_end, N)
#     eta = 1.0
#
#     shift_omega = 1j*eta + omega
#     result = []
#     for element in x:
#         Y = (1.0/(eta-1j*omega))*np.exp(-1j*omega*element)*cumu_itf_cf(shift_omega, p, lambda_m, gamma, sigma_dB, pure, itf_mean)
#         result.append(np.exp(eta*element)*np.real(np.trapz(y=Y, x=omega))/np.pi)
#
#     result = [0.0 if element<0.00001 else element for element in result]
#     result[0] = 0.0
#     print "CDF of CUMU I:", np.array(result)
#
#     BETA = np.log(10.0)/10.0
#     sigma = BETA*sigma_dB
#     A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
#     if pure:
#         if itf_mean:
#             A = 2*gamma/(2+gamma) * A
#         else:
#             A = 2.0 * A
#     com_f = 1-erf(p*lambda_m*np.pi*A*np.exp(sigma**2/8)/2*np.sqrt(1.0/x))
#     print "CDF of CUMU I according to formula:", com_f
#     return result


# The following function is depreciated. The quantization error is so huge!
# def neg_frac_moment_calculator(lambda_m, p, gamma, sigma_dB, pure, itf_mean):
#     # TODO: The following is just for the special case where gamma = 4. Not a general formula
#     EPS = np.finfo(float).eps
#     # The upper bound of intergral is 10.0
#     X_END = 15.0
#     # The number of SAMPLE_N also has impact to the finally exactness of framework.
#     SAMPLE_N = 3000
#     X = np.linspace(EPS, X_END, SAMPLE_N)
#     Y =  np.power(X, -2.0/gamma - 1.0) * itf_cf2cdf(X, p, lambda_m, gamma, sigma_dB, pure, itf_mean)
#     neg_frac_I = 2.0/gamma * np.trapz(y=Y, x=X) - 1.5
#     print "Numerically obtained neg.frac. I:", neg_frac_I, \
#         "Analytically obtained neg.frac. I:", itf_frac_moment_calculator(lambda_m, p, sigma_dB, 4.0, pure, itf_mean)
#     return neg_frac_I


def min_tx_power_best_case(p_outage, p_f_target, lambda_b, gamma, thetha_dB, sigma_dB):
    """
        This method is used to calculate the minimum transmit power in noise-limited system with outage
        probability target. Interference from other nodes is neglected.
        Thus it has nothing to do with "p", lambda_m.
        The background noise is normalized to 1.
        p_f_target: the target packet loss rate, for example 10%.
        p_outage: the portion of devices whose actual packet loss rate bypass p_f_target.
    """
    THETA = np.power(10, thetha_dB/10.0)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    term1 = np.pi * lambda_b * fm_shadowing
    print "term1", term1

    term2 = THETA * np.power(np.log(1.0/p_outage)/term1, gamma/2) * np.power(np.log(1.0/(1-p_f_target)), -1)

    return term2

def min_bs_intensity_best_case(p_outage, p_f_target, normalized_n, gamma, thetha_dB, sigma_dB):
    """
    This method is used to calculate the minimum required BS intensity, if the network wants to keep a QoS level
    (measured by outage probability). Note that we assume that the transmit power is always 1.0, and the background
    is normalized.

    Note that for a noise limited system, the slotted or non-slotted ALOHA has no impact to the performance.
    :param p_outage:
    :param p_f_target:
    :param lambda_b:
    :param gamma:
    :param thetha_dB:
    :param sigma_dB:
    :return:
    """
    THETA = np.power(10, thetha_dB/10.0)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    return np.power(np.pi*fm_shadowing, -1.0) * np.log(1.0/p_outage) * \
           np.power(np.log(1.0/(1-p_f_target))/THETA/normalized_n, -2.0/gamma)








if __name__ == "__main__":
    X_START = 0.0
    X_END = 3.0
    X_STEP = 0.002
    Y_START = 1e-3
    Y_END = 0.6
    Y_STEP = 0.1

    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 100)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.08
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.008
    thetha_dB = 3.0 # unit dB

    p_f_bs_nst_att_0 = bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0)
    p_f_bs_bst_att_8 = bs_best_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 10)

    # print lambda_m
    # print p_f_bs_nst_att_0
    # print p_f_bs_bst_att_8



    # print "Slotted Aloha"
    # sigma_dB = 0.0
    # pure=False
    # print div_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=False)
    # print nst_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=True)
    # print best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=False, itf_mean=True)
    # print "Non slotted Aloha, avg interference"
    # print div_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=True)
    # print nst_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=True)
    # print best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=True)
    # print "Non slotted Aloha, max interference"
    # print div_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=False)
    # print nst_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=False)
    # print best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max=0.1, pure=True, itf_mean=False)


    #============
    print "Test A general negative moment calculation method..."
    gamma = 6.0
    lambda_m = 2.4
    sigma_dB = 8.0
    pure = True
    max_itf = False
    avg_itf = True
    # print "Numerically obtained neg.frac. I:", neg_frac_moment_calculator(lambda_m, p, gamma, sigma_dB, pure=pure, itf_mean=avg_itf)
    # print "Analytically obtained neg.frac. I:", itf_frac_moment_calculator(lambda_m, p, sigma_dB, 4.0, pure=pure, itf_mean=avg_itf)
    print mrc_bs_rx_div_op([lambda_m], lambda_b, gamma, p, thetha_dB, sigma_dB, pure=pure, itf_mean=avg_itf)
    print "xxxx:", lt2plr(thetha_dB, lambda_m, lambda_b, gamma, p, sigma_dB, pure=pure, itf_mean=avg_itf)
    theta = np.power(10.0, 3.0/10)
    # print 1-erf(0.5*np.power(p*lambda_m/lambda_b*np.sqrt(np.pi*theta), -1))


    print gamma_f(1 - 2.0/3.0), gamma_f(1 - 2.0/3.3), gamma_f(1 - 2.0/3.5), gamma_f(1 - 2.0/3.8), gamma_f(1 - 2.0/4.0)


    # print num_op(EPS, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True)
