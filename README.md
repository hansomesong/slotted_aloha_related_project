# slotted_aloha_related_project

This python-written small project is used to simulate slotted ALOHA protocol, which allows to do some simulations to evaluate packet loss rate.

The basic idea is to modelize two classes: wireless channel and M2M device. 
As the fresh packet arrival process is assumed to be type of Poisson, we use a finite population of devices, each of whom generate a new
packet with a probability in each slot (namely Bernoulli distribution), to approximate the required Poission process in slotted ALOHA.

In each slot, a device first try to get the packet to be transmitted then generate a new packet with probability. For the selection of 
packet to be transmitted, the device first check the buffer for the backloggged packets. If one packet in the buffer is scheduled to be
transmitted at this slot, the device send this packet.Then the device decrement one for backoff timers for the remaining backlogged packets 

Note that packet loss rate is measured by the ratio between the number of dropped packets (still failed after the maximum allowed transmissions)
and the number of sent packets (delivered or dropped).

Some thinking:
    With large number of devices, the simulation result is more approach to the analytical one. For example, for fading and shadowing analysis, with 0dB as SINR threshold, around 1.0 alpha,
    1000 devices with 20,000 observation slot, the packet loss rate is greater than a simulation with 500 devices and 5000 observed slots.

    Due to list comprehension optimization, this simulation script is more resistant to the increase of device number.
