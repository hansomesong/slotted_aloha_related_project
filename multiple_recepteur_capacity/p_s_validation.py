# -*- coding: utf-8 -*-
# 这个Python脚本是为了验证，关于successful transmission probability的解析解是否正确

__author__ = 'qsong'

import numpy as np
from scipy.special import lambertw
from scipy.special import gamma as gamma_f
import matplotlib.pyplot as plt


if __name__ == "__main__":
    r = 15
    p = 0.002
    lambda_m = 0.6
    theta = np.power(10, 3.0/10) # unit dB
    gamma = 4.0
    sigma_dB = 12.0
    BETA = np.log(10)/10.0
    sigma_G = BETA*sigma_dB
    sigma_X = 2.0*sigma_G/gamma

    A_1 = p*lambda_m*np.pi*gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)*np.exp(np.power((1.0/np.sqrt(2))*(2.0/gamma-1)*BETA*sigma_dB/gamma, 2))*np.power(theta, 2.0/gamma)
    A_2 = p*lambda_m*np.pi*gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)*np.power(theta, 2.0/gamma)

    p_s_1 = np.exp(-A_2*r**2)

    l_solution = lambertw(A_1*r**2*sigma_X**2)
    p_s_2 = 1.0/np.sqrt(1+l_solution)*np.exp(-1*(l_solution**2 + 2*l_solution)/(2*sigma_X**2))

    # p_s_3 = np.exp(-A_1*r**2)

    print l_solution
    print p*lambda_m*np.pi*gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)
    print np.power(theta, 2.0/gamma)
    print np.exp(np.power(np.sqrt(2)*BETA*sigma_dB/gamma, 2))
    print "p_s_1", p_s_1
    print "p_s_2", p_s_2
    # print "p_s_3", p_s_3

    # 我想验证下 E[(XY)^{a}] 和 E【X^a】*E[Y^a]
    # 我们采用大量的数据求和来估计 随机变量的期望
    N = 100000
    X = np.random.exponential(1, N)
    Y = np.random.lognormal(0, sigma_G, N)
    left_value = np.mean(np.power(X*Y, 2.0/gamma))
    right_value = np.mean(np.power(X, 2.0/gamma))*np.mean(np.power(Y, 2.0/gamma))
    formula = gamma_f(1+2.0/gamma)*np.exp(0.5*(2.0/gamma-1)**2*sigma_G**2)
    print left_value, right_value, formula

    print np.mean(np.power(X, 2.0/gamma)), gamma_f(1+2.0/gamma)

    print np.mean(np.power(Y, 2.0/gamma)), np.exp(0.5*(2.0/gamma)**2*sigma_G**2)

    constant_B = gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)*np.exp(0.5*(2.0*sigma_G/gamma)**2)*theta**(2.0/gamma)


    constant_A = gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)*theta**(2.0/gamma)


    # 生成 lambda_m 的 ndarray
    lambda_m = 0.1
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.01
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.002
    thetha_dB = 3.0 # unit dB
    theta = np.power(10, 3.0/10)
    mu_shadow = 0.0
    # shadowing effect 的标准差，单位分贝
    sigma_dB = np.linspace(0, 12, 100)
    BETA = np.log(10.0)/10.0
    # 需要对标准差做个转换，因为对数基底的不同
    sigma_G = BETA*sigma_dB
    sigma_X = 2.0*sigma_G/gamma

    B_power = np.power(1+np.pi*sigma_X**2/1.0, -0.5)

    print "B_power", B_power







