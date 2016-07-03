# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
from scipy.stats import bernoulli
from time import time


def run_simulation(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t,  *args, **kwargs):
    start_t = int(time())
    sim_history = np.zeros((sim_duration, device_nb), dtype=np.int)
    for slot in range(sim_history.shape[0]-1):
        total_p = sum([l**(k-1)*m**(max_trans-k) for k in sim_history[slot] if k != 0])
        if total_p != 0:
            # that means the current slot is not idle
            for device_id, x in enumerate(sim_history[slot]):
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
                            sim_history[(new_slot, device_id)] = x+1
                        # Do not forget to 清零 for this slot.
                        sim_history[(slot, device_id)] = 0
                    elif x == max_trans:
                        # The case where max_trans has been reached. Failure of this packet transmission
                        sim_history[(slot, device_id)] = max_trans + 1

        # Second, generate packet and schedule the transmission in the next slot for each device
        new_pkts = [bernoulli.rvs(alpha/device_nb) for i in range(device_nb)]
        sim_history[slot+1] = np.array([new_pkts[id] if k == 0 else k for id, k in enumerate(sim_history[slot+1])])


    statistics = np.bincount(
        sim_history[warm_t:sim_duration, ::].reshape(1, device_nb*(sim_duration-warm_t))[0]
    )[1:]
    vector_p = [sum(statistics[i:])*1.0/sum(statistics) for i in range(len(statistics))]
    # The above code is equivalent to the following code:
    # truncated_history = sim_history[WARM_UP:SIM_DURATION, ::]
    # statistics = np.bincount(truncated_history.reshape(1, NB_DEVICE*SIM_DURATION)[0])[1:]

    end_t = int(time())
    time_elapsed = float(end_t-start_t)/60.0
    print "Execution time: ", time_elapsed
    return alpha, statistics, vector_p


if __name__ == "__main__":
    DEVICE_NB, SIM_DURATION = 500, 20000
    BACKOFF = 100
    THERSHOLD = 1.0
    L = 1
    M = 1
    MAX_TRANS = 5
    ALPHA = 0.26
    WARM_T = 0

    print run_simulation(ALPHA, MAX_TRANS, DEVICE_NB, THERSHOLD, L, M, BACKOFF, SIM_DURATION, WARM_T)


# start = int(time())
# print history.shape
# for slot in range(history.shape[0]-1):
#     # First judge which transmission is successful, which one should be backlogged.
#     total_p = sum([L**(k-1)*M**(MAX_TRANS-k) for k in history[slot] if k != 0])
#     if total_p != 0:
#     # that means the current slot is not idle
#         for devic_id, x in enumerate(history[slot]):
#             if total_p > L**(x-1)*M**(MAX_TRANS-x) / THERSHOLD:
#                 if x != 0 and x != MAX_TRANS:
#                     # MAX_TRANS has not been reached. Execute backoff procedure
#                     # The new slot index is the sum of current one and backoff length
#                     new_slot = int(np.random.exponential(scale=BACKOFF)) + 1 + slot
#                     # print new_slot
#                     if new_slot <= SIM_DURATION-1:
#                     # Take care that the selected new slot should not be out of range.
#                     # Also we should note that selected new slot has not yet scheduled for another retransmission
#                     # transmission trial should be incremented by 1
#                         history[(new_slot, devic_id)] = x+1
#                     # Do not forget to 清零 for this slot.
#                     history[(slot, devic_id)] = 0
#                 elif x == MAX_TRANS:
#                     # The case where MAX_TRANS has been reached. Failure of this packet transmission
#                     history[(slot, devic_id)] = MAX_TRANS + 1
#
#     # Second, generate packet and schedule the transmission in the next slot for each device
#     new_pkts = [bernoulli.rvs(INTENSITY/NB_DEVICE) for i in range(NB_DEVICE)]
#     # print new_pkts
#     history[(slot+1, devic_id)] = [new_pkts[devic_id] if k == 0 else k for devic_id, k in enumerate(history[slot+1])]
#     for devic_id, k in enumerate(history[slot+1]):
#         if k == 0:
#             history[(slot+1, devic_id)] = new_pkts[devic_id]
#
#
# partial_history = history[WARM_UP:SIM_DURATION, ::]
#
# statistics = np.bincount(partial_history.reshape(1, NB_DEVICE*SIM_DURATION)[0])[1:]
# vector_p = [sum(statistics[i:])*1.0/sum(statistics) for i in range(len(statistics))]
# print statistics, vector_p
# # print history
# end = int(time())
#
# time_elapsed = float(end-start)/60.0

# print "Execution time: ", time_elapsed
