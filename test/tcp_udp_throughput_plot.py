# -*- coding: utf-8 -*-
# This script is used to generate the figure, which illustrates the constatnt level
# of interference during a given time slot, in my PhD defense.
__author__ = 'qsong'

import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
import matplotlib.ticker as ticker
import numpy as np
from scipy.stats import bernoulli


from scipy import arange
x_axis = arange(0, 40.0, 0.2)
y_udp = 27*[0.0]
y_udp.append(0.34)
y_udp.extend(98*[0.4])
y_udp.extend(5*[0.0])
y_udp.append(0.38)
y_udp.extend(68*[.04])

y_tcp = 56*[0.0]
y_tcp.append(2.74432)
y_tcp.append(3.23568)
y_tcp.extend(68*[0.4])
y_tcp.extend(15*[0.0])
y_tcp.extend([0.00144, 1.1792, 3.13024, 2.08912])
y_tcp.extend(55*[0.4])

print len(x_axis), len(y_udp), len(y_tcp)


