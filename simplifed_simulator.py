# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
from scipy.stats import bernoulli
from time import time

NB_DEVICE, SIM_DURATION = 1000, 500000
history = np.zeros((SIM_DURATION, NB_DEVICE), dtype=np.int)
PKT_BUFFERE = [[] for i in range(NB_DEVICE)]
BACKOFF = 100
THERSHOLD = 1.0
L = 2
M = 1
MAX_TRANS = 5
INTENSITY = 0.78

WARM_UP = 0
def run_simulation(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t,  *args, **kwargs):
    start_t = int(time())
    sim_history = np.zeros((sim_duration, device_nb), dtype=np.int)
    for slot in range(sim_history.shape[0]-1):
        total_p = sum([l**(k-1)*m**(max_trans-k) for k in sim_history[slot] if k != 0])
        if total_p != 0:
            # that means the current slot is not idle
            for devic_id, x in enumerate(sim_history[slot]):
                if total_p > L**(x-1)*M**(max_trans-x) / threshold:
                    if x != 0 and x != max_trans:
                        # max_trans has not been reached. Execute backoff procedure
                        # The new slot index is the sum of current one and backoff length
                        new_slot = int(np.random.exponential(scale=backoff)) + 1 + slot
                        # print new_slot
                        if new_slot <= sim_duration-1:
                        # Take care that the selected new slot should not be out of range.
                        # Also we should note that selected new slot has not yet scheduled for another retransmission
                        # transmission trial should be incremented by 1
                            sim_history[(new_slot, devic_id)] = x+1
                        # Do not forget to 清零 for this slot.
                        sim_history[(slot, devic_id)] = 0
                    elif x == max_trans:
                        # The case where max_trans has been reached. Failure of this packet transmission
                        sim_history[(slot, devic_id)] = max_trans + 1

        # Second, generate packet and schedule the transmission in the next slot for each device
        new_pkts = [bernoulli.rvs(alpha/device_nb) for i in range(device_nb)]
        sim_history[(slot+1, devic_id)] = [new_pkts[devic_id] if k == 0 else k for devic_id, k in enumerate(sim_history[slot+1])]


    statistics = np.bincount(
        sim_history[WARM_UP:SIM_DURATION, ::].reshape(1, device_nb*(sim_duration-warm_t))[0]
    )[1:]
    vector_p = [sum(statistics[i:])*1.0/sum(statistics) for i in range(len(statistics))]
    # The above code is equivalent to the following code:
    # truncated_history = sim_history[WARM_UP:SIM_DURATION, ::]
    # statistics = np.bincount(truncated_history.reshape(1, NB_DEVICE*SIM_DURATION)[0])[1:]
    return alpha, statistics, vector_p
    # print history
    end = int(time())
    time_elapsed = float(end-start)/60.0
    print "Execution time: ", time_elapsed

start = int(time())
print history.shape
for slot in range(history.shape[0]-1):
    # First judge which transmission is successful, which one should be backlogged.
    total_p = sum([L**(k-1)*M**(MAX_TRANS-k) for k in history[slot] if k != 0])
    if total_p != 0:
    # that means the current slot is not idle
        for devic_id, x in enumerate(history[slot]):
            if total_p > L**(x-1)*M**(MAX_TRANS-x) / THERSHOLD:
                if x != 0 and x != MAX_TRANS:
                    # MAX_TRANS has not been reached. Execute backoff procedure
                    # The new slot index is the sum of current one and backoff length
                    new_slot = int(np.random.exponential(scale=BACKOFF)) + 1 + slot
                    # print new_slot
                    if new_slot <= SIM_DURATION-1:
                    # Take care that the selected new slot should not be out of range.
                    # Also we should note that selected new slot has not yet scheduled for another retransmission
                    # transmission trial should be incremented by 1
                        history[(new_slot, devic_id)] = x+1
                    # Do not forget to 清零 for this slot.
                    history[(slot, devic_id)] = 0
                elif x == MAX_TRANS:
                    # The case where MAX_TRANS has been reached. Failure of this packet transmission
                    history[(slot, devic_id)] = MAX_TRANS + 1

    # Second, generate packet and schedule the transmission in the next slot for each device
    new_pkts = [bernoulli.rvs(INTENSITY/NB_DEVICE) for i in range(NB_DEVICE)]
    # print new_pkts
    history[(slot+1, devic_id)] = [new_pkts[devic_id] if k == 0 else k for devic_id, k in enumerate(history[slot+1])]
    for devic_id, k in enumerate(history[slot+1]):
        if k == 0:
            history[(slot+1, devic_id)] = new_pkts[devic_id]


partial_history = history[WARM_UP:SIM_DURATION, ::]

statistics = np.bincount(partial_history.reshape(1, NB_DEVICE*SIM_DURATION)[0])[1:]
vector_p = [sum(statistics[i:])*1.0/sum(statistics) for i in range(len(statistics))]
print statistics, vector_p
# print history
end = int(time())

time_elapsed = float(end-start)/60.0

print "Execution time: ", time_elapsed
