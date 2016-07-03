# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
from scipy.stats import bernoulli
from time import time
import multiprocessing as mp
import csv
import os
import json
from time import strftime

# Although many processes could involve in the treatment of log file, we have only one process responsible for
# writing logfile-related entry(a row in the context of CSV file) into the given CSV file. The advantage of this
# is to avoid the concurrence about the write access to target CSV file.
def listener(sim_result_f, q):
    with open(sim_result_f, 'a') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        while 1:
            csv_row = q.get()
            # 如果从 (被N多个进程共享的) result queue 中提取到了 字符串 ‘TERMINATE’，那么结束listener进程
            if csv_row == "TERMINATE":
                break
            spamwriter.writerow(csv_row)


# In our process pool, except the one in charge of writing records into CSV file, all rest processes are used
# to treat log file stored in a certain directory. Every time a process processes a log file, it store the retrieved
# information into a QUEUE data structure, which will be served by listener process.
def worker(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t, q):
    csv_row = []
    sim_result = run_simulation2(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t)
    csv_row.extend(sim_result[2])
    csv_row.append(alpha)
    q.put(csv_row)


def run_simulation2(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t):
    '''
        In this method, during the backoff slot, a new packet is allowed to be transmitted until the deliveray or failure
         of last packet. The generated packets generate will be buffered.
    '''
    start_t = int(time())
    seed = hash(start_t + os.getpid()*13)
    np.random.seed(seed)
    pkt_buffer = [[] for i in range(device_nb)]
    backlog = [0 for i in range(device_nb)]
    sim_history = np.zeros((sim_duration, device_nb), dtype=np.int)
    for slot in range(sim_duration):
        # First generate new packets
        new_pkts = [bernoulli.rvs(alpha/device_nb) for i in range(device_nb)]
        # Second determine which packets to be transmitted
        for devic_id, k in enumerate(sim_history[slot]):
            if backlog[devic_id] == 1:
                if new_pkts[devic_id] != 0:
                    # The case that next slot is schedulded for retransmission. and generate a new packet
                    # We put this newly generated packet into packet buffer
                    pkt_buffer[devic_id].append(1)
            else:
                # The case not backlogged.
                # First send buffered packets if any. Otherwise send the newly generated packet if any.
                if len(pkt_buffer[devic_id]) != 0:
                    sim_history[(slot, devic_id)] = pkt_buffer[devic_id].pop(0)
                    if new_pkts[devic_id] != 0:
                        # Since the newly generated packe has no chance to be transmitted, it has to be store into buffer.
                        pkt_buffer[devic_id].append(1)
                else:
                    # The case the packet buffer is empty.
                    if new_pkts[devic_id] != 0:
                        sim_history[(slot, devic_id)] = 1

        total_p = sum([l**(k-1)*m**(max_trans-k) for k in sim_history[slot] if k != 0])
        if total_p != 0:
            # that means the current slot is not idle
            for device_id, x in enumerate(sim_history[slot]):
                if total_p > l**(x-1)*m**(max_trans-x) / threshold:
                    # the case of transmission failure
                    backlog[device_id] = 1
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
                else:
                    # the case of successful transmission. Do not forget to flag for this device as non-backlogged.
                    backlog[device_id] = 0

    statistics = np.bincount(
        sim_history[warm_t:sim_duration, ::].reshape(1, device_nb*(sim_duration-warm_t))[0]
    )[1:]
    vector_p = [sum(statistics[i:])*1.0/sum(statistics) for i in range(len(statistics))]
    # The above code is equivalent to the following code:
    # truncated_history = sim_history[WARM_UP:SIM_DURATION, ::]
    # statistics = np.bincount(truncated_history.reshape(1, NB_DEVICE*SIM_DURATION)[0])[1:]

    end_t = int(time())
    time_elapsed = float(end_t-start_t)/60.0
    print "Execution time ", time_elapsed, "for this task ", alpha, "Result:", statistics, vector_p
    return alpha, statistics, vector_p


def run_simulation(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t):
    '''
        In this method, during the backoff slot, it is possible to generate new packets.
        Especially the number of backoff slot respects exponential distribution.
        I have ran several simulations for l=1 m=1, threshould=1. For an intensity of 0.28.
        The packet loss rate is systematically greater than that of analytical result.
    '''
    start_t = int(time())
    seed = hash(start_t + os.getpid()*13)
    np.random.seed(seed)
    sim_history = np.zeros((sim_duration, device_nb), dtype=np.int)
    for slot in range(sim_duration-1):
        # First generate new packets.
        # new_pkts = [bernoulli.rvs(alpha/device_nb) for i in range(device_nb)]
        # Which device has packets to transmit?
        # sim_history[slot] = np.array([new_pkts[id] if k == 0 else k for id, k in enumerate(sim_history[slot])])
        sim_history[slot] = np.array([bernoulli.rvs(alpha/device_nb) if k == 0 else k for k in sim_history[slot]])
        total_p = sum([l**(k-1)*m**(max_trans-k) for k in sim_history[slot] if k != 0])
        if total_p != 0:
            # that means the current slot is not idle
            for device_id, x in enumerate(sim_history[slot]):
                if total_p > l**(x-1)*m**(max_trans-x) / threshold:
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
    statistics = np.bincount(
        sim_history[warm_t:sim_duration, ::].reshape(1, device_nb*(sim_duration-warm_t))[0]
    )[1:]
    vector_p = [sum(statistics[i:])*1.0/sum(statistics) for i in range(len(statistics))]
    # The above code is equivalent to the following code:
    # truncated_history = sim_history[WARM_UP:SIM_DURATION, ::]
    # statistics = np.bincount(truncated_history.reshape(1, NB_DEVICE*SIM_DURATION)[0])[1:]

    end_t = int(time())
    time_elapsed = float(end_t-start_t)/60.0
    print "Execution time ", time_elapsed, "for this task ", alpha, "Result:", statistics, vector_p
    return alpha, statistics, vector_p

def main(config_f):
    #must use Manager queue here, or will not work
    with open(config_f) as json_file:
        json_config = json.load(json_file)

    SIM_DURATION = json_config['SIM_DURATION']
    DEVICE_NB = json_config['N']
    SIM_STEP = json_config['sim_step']
    BACKOFF = json_config["BACKOFF"]
    THERSHOLD = json_config['THRESLD']
    L = json_config['l']
    M = json_config['m']
    MAX_TRANS = json_config['MAX_TRANS']
    ALPHA_START = json_config['alpha_start']
    ALPHA_END = json_config['alpha_end']
    WARM_T = json_config['WARM_UP']

    ALPHA_INTERVAL = [ALPHA_START + i*SIM_STEP for i in range(int((ALPHA_END-ALPHA_START)/SIM_STEP)+1)]


    n = 5
    if len(ALPHA_INTERVAL) == 1:
        ALPHA_INTERVAL = [ALPHA_START for i in range(n)]

    sim_result_f = \
        "dataset_expon_backoff/sim_simd={0}_N={1}_threshold={2}_l={3}_m={4}_backoff={5}_start={6}_simstep={7}_{8}.csv"\
            .format(SIM_DURATION, DEVICE_NB, THERSHOLD, L, M, BACKOFF, ALPHA_START, SIM_STEP, strftime("%H%M%S"))

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count()+2)

    f_handler = open(sim_result_f, 'w')
    #put listener to work first
    watcher = pool.apply_async(listener, (sim_result_f, q,))

    #fire off workers
    jobs = []

    for alpha in ALPHA_INTERVAL:
        job = pool.apply_async(worker, (alpha, MAX_TRANS, DEVICE_NB, THERSHOLD, L, M, BACKOFF, SIM_DURATION, WARM_T, q))
        jobs.append(job)

    #collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #Now we are done, kill the listener
    q.put("TERMINATE")
    pool.close()
    # Do not forget to close file at the end.
    f_handler.close()

if __name__ == "__main__":

    start_t = int(time())
    config_f = os.path.join('exp_configs', 'case_l=1_m=1_threshold=0dB_N=500.json')
    print "Now do simulation with configuration file: ", config_f

    main(config_f)
    end_t = int(time())
    time_elapsed = float(end_t-start_t)/60.0

    print "Total Execution time: ", time_elapsed
    # First read


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


    # with open(config_f) as json_file:
    #     json_config = json.load(json_file)
    #
    # SIM_DURATION = json_config['SIM_DURATION']
    #
    # DEVICE_NB = json_config['N']
    # SIM_STEP = json_config['sim_step']
    # WARM_UP = json_config['WARM_UP']
    # BACKOFF = json_config["BACKOFF"]
    # SIM_DURATION = json_config['SIM_DURATION']
    # THERSHOLD = json_config['THRESLD']
    # L = json_config['l']
    # M = json_config['m']
    # MAX_TRANS = json_config['MAX_TRANS']
    # ALPHA_1 = 0.25
    # WARM_T = json_config['WARM_UP']
    #
    # print run_simulation(ALPHA_1, MAX_TRANS, DEVICE_NB, THERSHOLD, L, M, BACKOFF, SIM_DURATION, WARM_T)
