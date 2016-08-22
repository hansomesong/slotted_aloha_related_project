# -*- coding: utf-8 -*-
# Compared with script simplified_simulator.py, this script takes the definition of SINR as the received power and other
# recevied power
__author__ = 'qsong'

import numpy as np
from scipy.stats import bernoulli
from time import time
import multiprocessing as mp
import csv
import os
import json
from time import strftime
import glob
import pprint
from scipy.stats import itemfreq

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
def worker(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing, sigma_shadowing, q):
    csv_row = []
    sim_result = run_simulation(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing, sigma_shadowing)
    csv_row.extend(sim_result[2])
    csv_row.append(alpha)
    q.put(csv_row)

def run_simulation(alpha, max_trans, device_nb, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing, sigma_shadowing):
    '''
        In this method, during the backoff slot, it is possible to generate and transmit new packets.
        Only the packets generated during the slot scheduled for retransmission will be abandoned.
        In addition, when choosing the backoff slot, make sure that the choosed slot is not scheduled
        for retransmission. Otherwise continue to choose the backoff slot until the allowed one.

        I have ran several simulations for l=1 m=1, threshould=1. For an intensity of 0.28.
        The packet loss rate is systematically greater than that of analytical result.

        ipc_way: referst to imperfect power control way. If set as "DIFF", each retransmission power error may be different.
        If set as "SAME", each retransmission power error is identical.
    '''
    # LM = []
    # if ipc_way == "DIFF":
    #     LM = [1.0*l**(k)*m**(max_trans-k-1) for k in range(max_trans)]
    # elif ipc_way == "SAME":
    #     LM = [(1.0*l/m)**k for k in range(max_trans)]
    #     IPC_LEVEL = np.random.lognormal(mu_ipc, sigma_shadowing, device_nb)
    # else:
    #     print "Wrong configuration for Imperfect power control manner"
    #     exit()
    BETA = np.log(10)/10.0
    LM = [1.0*l**(k)*m**(max_trans-k-1) for k in range(max_trans)]
    start_t = int(time())
    seed = hash((start_t + os.getpid()*13)*0.0000001)
    np.random.seed(seed)
    sim_history = np.zeros((sim_duration, device_nb), dtype=np.int)
    # shadowings = np.random.lognormal(BETA*mu_shadowing, BETA*sigma_shadowing, device_nb)

    for slot in range(sim_duration-1):
        # First generate new packets.
        # Which device has packets to transmit?
        # In the following, k refers to the k th transmision instead of retransmission
        sim_history[slot] = np.array([bernoulli.rvs(alpha/device_nb) if k == 0 else k for k in sim_history[slot]])

        # With which transmit power they can sue?
        power_levels = np.array([LM[k-1] if k != 0 else k*1.0 for k in sim_history[slot]])
        shadowings = np.random.lognormal(BETA*mu_shadowing, BETA*sigma_shadowing, device_nb)
        fadings = np.random.exponential(scale=mu_fading, size=device_nb)
        power_levels *= fadings*shadowings
        # print power_levels
        # power_levels *= shadowings
        # nb_pkts = sum([1 for k in sim_history[slot] if k != 0])

        # 我们采用 SINR门限值*干扰值 和 接收功率 比较的方式，加速仿真的执行
        # (计算实际接收信噪比的思路会不可避免地考虑信道中只有一个传输的情况，导致程序的执行效率低下)
        total_p = sum(power_levels)

        curr_trans_results = [
            (10**(0.1*threshold)) * (total_p - power_levels[device_id]) > power_levels[device_id]
            if k != 0 else False for device_id, k in enumerate(sim_history[slot])
        ]

        # print curr_trans_results
        for device_id, curr_trans_result in enumerate(curr_trans_results):
            # 如果 SINR门限值*干扰值 大于 接收功率，则需要对这个包执行重传操作
            if curr_trans_result:
                # 需要知道，这是第几次传输失败了
                x = sim_history[slot][device_id]
                if x != max_trans:
                    # max_trans has not been reached. Execute backoff procedure
                    # The new slot index is the sum of current one and backoff length
                    while 1:
                        new_slot = int(np.random.exponential(scale=backoff)) + 1 + slot
                        if new_slot <= sim_duration-1:
                        # Take care that the selected new slot should not be out of range.
                            # Also we should note that selected new slot has not yet scheduled
                            # for another retransmission
                            if sim_history[(new_slot, device_id)] == 0:
                            # transmission trial should be incremented by 1
                                sim_history[(new_slot, device_id)] = x+1
                                break
                        else:
                            break
                    # Do not forget to 清零 for this slot.
                    sim_history[(slot, device_id)] = 0
                elif x == max_trans:
                    # The case where max_trans has been reached. Failure of this packet transmission
                    sim_history[(slot, device_id)] = max_trans + 1

        # if nb_pkts > 1:
        #     # that means the current slot is selected by more than one packets.
        #     # In this case, calculate the SINR to determine which one can be received, which one should be backlogged.
        #     # For other cases, slot Idle or occupied by one packet, do nothing.
        #     total_p = sum(power_levels)
        #     for device_id, x in enumerate(sim_history[slot]):
        #         if x != 0:
        #             # The case where device has no packet to transmit, we do not care...
        #             if 10**(0.1*threshold) > power_levels[device_id] / (total_p - power_levels[device_id]):
        #                 if x != max_trans:
        #                     # max_trans has not been reached. Execute backoff procedure
        #                     # The new slot index is the sum of current one and backoff length
        #                     while 1:
        #                         new_slot = int(np.random.exponential(scale=backoff)) + 1 + slot
        #                         if new_slot <= sim_duration-1:
        #                         # Take care that the selected new slot should not be out of range.
        #                             # Also we should note that selected new slot has not yet scheduled
        #                             # for another retransmission
        #                             if sim_history[(new_slot, device_id)] == 0:
        #                             # transmission trial should be incremented by 1
        #                                 sim_history[(new_slot, device_id)] = x+1
        #                                 break
        #                         else:
        #                             break
        #                     # Do not forget to 清零 for this slot.
        #                     sim_history[(slot, device_id)] = 0
        #                 elif x == max_trans:
        #                     # The case where max_trans has been reached. Failure of this packet transmission
        #                     sim_history[(slot, device_id)] = max_trans + 1

        # Second, generate packet and schedule the transmission in the next slot for each device

    # print "History:", sim_history[warm_t:sim_duration, ::].reshape(1, device_nb*(sim_duration-warm_t))
    # statistics = np.bincount(
    #     sim_history[warm_t:sim_duration, ::].reshape(1, device_nb*(sim_duration-warm_t))[0]
    # )[1:]

    # 统计第 i 次传输出现的次数，并据此计算至少需要传输 i 次的频数，如果实验足够长，频数将趋近于概率
    statistics = itemfreq(sim_history[warm_t:sim_duration, ::].reshape(1, device_nb*(sim_duration-warm_t))[0])[1:]

    vector_p = [0 for i in range(max_trans+1)]

    total_transmission = sum([element[1] for element in statistics])

    for element in statistics:
        ith_trans = element[0]
        satisfied_trans = sum([element[1] for element in statistics if element[0] >= ith_trans])
        vector_p[ith_trans - 1] = satisfied_trans*1.0/total_transmission


    # Initialize the return probability vector with 0, the length of this vector is 1+MAX_ALLOWED_RETRANSMISSION_NB
    # vector_p = [sum(statistics[i:])*1.0/sum(statistics) for i in range(len(statistics))]
    end_t = int(time())
    time_elapsed = float(end_t-start_t)/60.0

    pprint.pprint(sim_history)
    print "Execution time", int(time_elapsed), "for this task", alpha, "with seed ", seed, "Result:", statistics, vector_p
    return alpha, list(statistics), vector_p

def main(com_config_f, logs_directory):
    '''
    :param com_config_f: 公共配置文件，存储着通用的仿真参数。
    :param logs_directory:
    :return: nothing
    '''
    #must use Manager queue here, or will not work
    # 注意：一定要使用 Manager queue，否则无法工作，至于原因，我不晓得
    with open(com_config_f) as json_file:
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
    MU_FADING = json_config['MU_FADING']
    MU_SHADOWING = json_config['MU_SHADOWING']
    SIGMA_SHADOWING = json_config['SIGMA_SHADOWING']

    n = 10
    if len(ALPHA_INTERVAL) == 1:
        ALPHA_INTERVAL = [ALPHA_START for i in range(n)]

    # 将仿真结果存储在 sim_result_f 指向的文件中
    sim_result_f = os.path.join(
        logs_directory,
        "simd={0}_N={1}_threshold={2}_l={3}_m={4}_backoff={5}_start={6}_end={7}_simstep={8}_sigmas={9}_{10}.csv".format(
            SIM_DURATION, DEVICE_NB, THERSHOLD, L, M, BACKOFF, ALPHA_START, ALPHA_END, SIM_STEP, SIGMA_SHADOWING, strftime("%Y%m%d%H%M%S")
        )
    )

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count()+3)

    f_handler = open(sim_result_f, 'w')
    #put listener to work first
    watcher = pool.apply_async(listener, (sim_result_f, q,))

    #fire off workers
    jobs = []

    for alpha in ALPHA_INTERVAL:
        job = pool.apply_async(
            worker,
            (alpha, MAX_TRANS, DEVICE_NB, THERSHOLD, L, M, BACKOFF, SIM_DURATION, WARM_T, MU_FADING, MU_SHADOWING, SIGMA_SHADOWING, q)
        )
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

    all_logs = glob.glob(os.path.join('sim_configs', 'fading_shadowing', "*.json"))
    logs_directory = 'logs'
    print all_logs

    for config_f in all_logs:

        # config_f = os.path.join('sim_configs', 'case_K=5_l=1_m=1_threshold=3dB.json')
        # The simulation result will be logged into files of type CSV, in folder logs.
        # First check the existence of this folder and creat it if necessary.


        if not os.path.exists(logs_directory):
            os.makedirs(logs_directory)
        print "Now do simulation with configuration file: ", config_f

        main(config_f, logs_directory)
        end_t = int(time())
        time_elapsed = float(end_t - start_t)/60.0

        print "Total Execution time: ", time_elapsed

