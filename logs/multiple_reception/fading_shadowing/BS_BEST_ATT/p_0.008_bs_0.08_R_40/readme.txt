Finally, I think I have identified the cause why there exists always some gap between simulation and analytical
results for the nearest BS attach method, especially when wireless channel is influenced by log-normal shadowing.
The root cause is that the selected spatial device intensity is still too small.
Since the simulation is run on a limited region, which is actually characterized by a Bononimal Point Process
(BPP), however our analytical model is derived on PPP. When intensity spatial device is enough large, BPP can be
well approximated by a PPP.

In the pass, the basic setting for simulation is:
p = 0.008, => OK, situable for modeling rare events
lambda_m = 0.01 : 0.01 : 0.05
lambda_b = 0.004

Although normalized load, which is defined as p*lambda_m/lamdba_b is in interval [0.02, 0.10]
However, a spatial device intensity such as 0.02 is too small to be modeled as a PPP!

In the past, the reasons