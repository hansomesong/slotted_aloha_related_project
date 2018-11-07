# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats

intensity = 0.4
vector_p = [1.0, 0.60, 0.36, 0.22, 0.13]
# vector_p = [1.0, 0.00, 0.00, 0.00, 0.00]

intensities = [intensity*p for p in vector_p]

print sum(intensities)

mu = 0.0
sigma = 1.0
beta = 0.1*np.log(10.0)
# beta = 1.0
rho = 1.0
v = 1.0
K = 4

# First order moment of Y
fm_Y = intensity*sum(vector_p)*rho*np.exp(0.5*np.power(beta*sigma, 2))
fm_Y_plus = fm_Y/(1-np.exp(-sum(intensities)))
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

# z1 = [sum([sum(np.exp(beta*np.random.normal(0, sigma, np.random.poisson(x)))) for x in intensities]) for i in range(nb)]
# z1 = [y for y in z1 if y != 0.0]


z1 = [sum(np.random.lognormal(mu, sigma, max([np.random.poisson(sum(intensities)), 1]))) for i in range(nb)]
log_z = np.log(z1)
log_z_mean = np.mean(log_z)
log_z_std = np.std(log_z)
print "mean:", log_z_mean
print "stand error:", log_z_std
# count, bins, ignored = plt.hist(z, 100, normed=True, align='mid')
# sns.distplot(z1, bins=[-0.1+0.05*x for x in range(600)]).set(xlim=(0, 30))
# sns.distplot(z1, 600, hist=True)


weights = 100*np.ones_like(np.array(log_z))/float(len(z1))
count, bins, ignored = plt.hist(log_z,  bins=400, normed=True, color='w')
print "Intergral:", (count * np.diff(bins)).sum()
print "Count:", count
print "Bins:", bins
plt.title(r"Histogram of $ln^{0}$, with packet arrival rate {1}".format('Z', sum(intensities)), fontsize=30)

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


