# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats

N = 5


mu = 0.0
sigma = 1.0
beta = 0.1*np.log(10.0)
rho = 1.0
v = 1.0
K = 4

# First order moment of Y
fm_Y = N*rho*np.exp(0.5*np.power(beta, 2)*np.power(sigma, 2))
# Second order moment of Y



sm_Y = intensity*sum(vector_p)*np.exp(2*np.power(beta*sigma, 2)) \
       + np.power(intensity*sum(vector_p), 2)*np.exp(np.power(beta*sigma, 2))

sm_Y_plus = sm_Y/(1-np.exp(-sum(intensities)))

mu_eq = 2*np.log(fm_Y_plus) - 0.5*np.log(sm_Y_plus)
sigma_eq = np.sqrt(np.log(sm_Y_plus) - 2*np.log(fm_Y_plus))

print fm_Y
print sm_Y

print "estimated equivalent mean", mu_eq
print "estimated equivalent stand deviation:", sigma_eq


plt.figure(1)
mu, sigma = 0.0, 1.0 # mean and standard deviation
intensity = sum(vector_p)


nb = 4000000
z1 = [sum(np.random.lognormal(mu, sigma, max([np.random.poisson(sum(intensities)), 1]))) for i in range(nb)]
log_z = np.log(z1)
log_z_mean = np.mean(log_z)
log_z_std = np.std(log_z)
print "MIN:", min(z1)
print "MAX:", max(z1)
print "mean:", log_z_mean
print "stand error:", log_z_std

print "0:", sum([1 for x in z1 if x == 0])
# count, bins, ignored = plt.hist(z, 100, normed=True, align='mid')
# sns.distplot(z1, bins=[-0.1+0.05*x for x in range(600)]).set(xlim=(0, 30))
# sns.distplot(z1, 600, hist=True)


weights = np.ones_like(np.array(log_z))/float(len(z1))
count, bins, ignored = plt.hist(log_z, 100, normed=True, color='w')
plt.title(r"Histogram of $ln^{Z}$, with packet arrival rate 0.8", fontsize=30)

xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 1000)
mu, std = mu_eq, sigma_eq
pdf = scipy.stats.norm.pdf(x, mu, std)
plt.plot(x, pdf, linewidth=2, color='r', label=r"pdf of log(Z) with Mean:{:0.3f}, STD:{:1.3f}".format(mu_eq, sigma_eq))
plt.legend(loc='upper left')

# weights=weights,

plt.xlabel("Value of Z")
plt.ylabel("Frequency")
plt.show()


