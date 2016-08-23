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
        # 如果 shadowing 的方差不是0，那么生成一系列 log-normal 随机数，否则生成同等长度的1
        if sigma_shadowing > 0:
            shadowings = np.random.lognormal(BETA*mu_shadowing, BETA*sigma_shadowing, device_nb)
        else:
            shadowings = np.ones(device_nb)

        # 对于满足指数分布的 fading effect, 也是同理
        if mu_fading > 0:
            fadings = np.random.exponential(scale=mu_fading, size=device_nb)
        else:
            fadings = np.ones(device_nb)

        # 计算考虑诸多因素之后的，接受功率
        power_levels *= fadings*shadowings
        total_p = sum(power_levels)
        # 我们采用 SINR门限值*干扰值 和 接收功率 比较的方式，加速仿真的执行
        # (计算实际接收信噪比的思路会不可避免地考虑信道中只有一个传输的情况，导致程序的执行效率低下)
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
    print "Time:", int(time_elapsed), "Alpha:", alpha, "Seed:", seed, "Result:", [e[1] for e in statistics], vector_p
    return alpha, list(statistics), vector_p

def main(sim_config_dict, logs_directory):
    '''
    :param com_config_f: 公共配置字典，存储着通用的仿真参数。
    :param logs_directory:
    :return: nothing
    '''
    #must use Manager queue here, or will not work
    # 注意：一定要使用 Manager queue，否则无法工作，至于原因，我不晓得

    SIM_DURATION = sim_config_dict['SIM_DURATION']
    DEVICE_NB = sim_config_dict["DEVICE_NB"]
    BACKOFF = sim_config_dict["BACKOFF"]
    THERSHOLD = sim_config_dict['THRESLD']
    L = sim_config_dict['L']
    M = sim_config_dict['M']
    MAX_TRANS = sim_config_dict['MAX_TRANS']
    ALPHA = sim_config_dict['ALPHA']
    WARM_T = sim_config_dict['WARM_UP']
    MU_FADING = sim_config_dict['MU_FADING']
    MU_SHADOWING = sim_config_dict['MU_SHADOWING']
    SIGMA_SHADOWING = sim_config_dict['SIGMA_SHADOWING']
    # 针对每个 alpha 值，仿真重复次数
    SIM_REPEAT_NB = sim_config_dict['SIM_REPEAT_NB']

    ALPHAS = [ALPHA for i in range(SIM_REPEAT_NB)]

    # 将仿真结果存储在 sim_result_f 指向的文件中
    sim_result_f = os.path.join(
        logs_directory,
        "simd={0}_N={1}_threshold={2}dB_l={3}_m={4}_backoff={5}_alpha={6}_mufading={7}_mushadowing={8}_sigmashadowing={9}_tmp={10}.csv".format(
            SIM_DURATION, DEVICE_NB, THERSHOLD, L, M, BACKOFF, ALPHA, MU_FADING, MU_SHADOWING, SIGMA_SHADOWING, strftime("%Y%m%d%H%M%S")
        )
    )

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count()+3)

    f_handler = open(sim_result_f, 'w')

    # put listener to work first
    watcher = pool.apply_async(listener, (sim_result_f, q,))

    #fire off workers
    jobs = []

    for alpha in ALPHAS:
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
    logs_directory = 'logs'
    # The simulation result will be logged into files of type CSV, in folder logs.
    # First check the existence of this folder and creat it if necessary.
    if not os.path.exists(logs_directory):
        os.makedirs(logs_directory)
    sim_config_f = os.path.join('sim_configs', 'fading_shadowing', 'case_K=5_l=1_m=1_threshold=-3dB.json')
    print "Now do simulation with configuration file: ", sim_config_f

    sim_config_dict = {}

    with open(sim_config_f) as json_file:
        json_config = json.load(json_file)

    sim_config_dict["BACKOFF"] = json_config["BACKOFF"]
    sim_config_dict['THRESLD'] = json_config['THRESLD']
    sim_config_dict["L"] = json_config['L']
    sim_config_dict["M"] = json_config['M']
    sim_config_dict["MAX_TRANS"] = json_config['MAX_TRANS']
    sim_config_dict["WARM_UP"] = json_config['WARM_UP']
    sim_config_dict["MU_FADING"] = json_config['MU_FADING']
    sim_config_dict["MU_SHADOWING"] = json_config['MU_SHADOWING']
    sim_config_dict["SIGMA_SHADOWING"] = json_config['SIGMA_SHADOWING']
    # 针对每个 alpha 值，仿真重复次数
    sim_config_dict["SIM_REPEAT_NB"] = json_config['SIM_REPEAT_NB']
    sim_config_dict["SIM_INCRE_STEP"] = json_config['SIM_INCRE_STEP']

    ALPHA_START = json_config['ALPHA_START']
    ALPHA_END = json_config['ALPHA_END']
    SIM_INCRE_STEP = json_config['SIM_INCRE_STEP']
    SIM_DURATION = json_config['SIM_DURATION']
    DEVICE_NB = json_config['DEVICE_NB']
    ALPHA_INTERVAL = np.arange(ALPHA_START, ALPHA_END, SIM_INCRE_STEP)

    for order, ALPHA in enumerate(ALPHA_INTERVAL, 1):
        # 我们一般从 packet loss rate 大致为 0.01 的 ALPHA 值开始，所以前四次仿真，我们把 simulation duration 设为 10000
        # 仿真需要的设备数目因而也不需要很高 (因为 ALPHA 取值不高)
        sim_config_dict["ALPHA"] = ALPHA
        if order < 5:
            sim_config_dict["DEVICE_NB"] = DEVICE_NB[0]
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[0]
        elif 5 <= order < 16:
            sim_config_dict["DEVICE_NB"] = DEVICE_NB[1]
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[1]
        else:
            sim_config_dict["DEVICE_NB"] = DEVICE_NB[2]
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[2]
        print sim_config_dict
        # 将填充好的仿真参数字典传递给 main(), 开启多进程下的仿真
        main(sim_config_dict, logs_directory)

    end_t = int(time())
    time_elapsed = float(end_t - start_t)/60.0
    print "Total Execution time: ", time_elapsed

