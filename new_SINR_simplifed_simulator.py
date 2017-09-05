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
from statsmodels.distributions.empirical_distribution import ECDF
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
def worker(alpha, max_trans, binomial_p, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing, sigma_shadowing, q):
    csv_row = []
    sim_result = run_simulation(alpha, max_trans, binomial_p, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing, sigma_shadowing)
    csv_row.extend(sim_result[2])
    csv_row.append(alpha)
    q.put(csv_row)

def run_simulation(alpha, max_trans, binomial_p, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing, sigma_shadowing, output_statistics=False):
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
    start_t = int(time())
    # Generate the random number generator seed, make sure each process uses a different seed otherwise
    # the simulation result will be the same when multiple processes run the this method.
    seed = hash((start_t + os.getpid()*13)*0.0000001)
    # seed = 2160133643
    # seed = 1792749438
    # seed = 2954458010
    np.random.seed(seed)

    BETA = np.log(10)/10.0
    # BACK_OFFS = [backoff*np.power(2, i) for i in range(max_trans)]
    # BACK_OFFS = [8, 30, 60, 100]
    BACK_OFFS = [backoff for i in range(max_trans)]
    # The involved device number will be determined together by binomial probability and alpha
    # The involved device number should be not less than 400.
    # The probability p in binomial distribution should less than 0.01 to assure that binomial distribution is close to poisson distribution
    device_nb = int(alpha/binomial_p)

    # Now Add one zero in the beginning of LM list, to facilitate the generation of
    # Then each element in LM list represents the transmit power for the kth transmission (Not retransmission).
    # For example, if max_trans = 5, l=2, m=1, then
    # LM = [0, 1, 2, 4, 8, 16] The 0th transmission transmit power is of course 0
    LM = [1.0*np.power(l, k-1)*np.power(m, max_trans-k) if k != 0.0 else 0.0 for k in range(max_trans+1)]

    # sim_history now is 3-d array.
    # The third dimension is used to represent the transmission index and received power levels history
    sim_history = np.zeros((sim_duration, device_nb, max_trans+1), dtype=np.float)
    for slot in range(sim_duration):
        # First generate new packets.
        # Which device has packets to transmit?
        # In the following, k refers to the k th transmision instead of retransmission
        sim_history[slot, :, 0] = np.array(
            [bernoulli.rvs(binomial_p) if k == 0.0 else k for k in sim_history[slot, :, 0]]
        )
        # With which receive power they can sue?
        rec_power_levels = np.array([LM[int(k)] for k in sim_history[slot, :, 0]])
        # 如果 shadowing 的方差不是0，那么生成一系列 log-normal 随机数，否则生成同等长度的1
        if sigma_shadowing > 0:
            shadowings = np.random.lognormal(BETA*mu_shadowing, BETA*sigma_shadowing, device_nb)
        else:
            shadowings = np.ones(device_nb)

        if mu_fading > 0:
            fadings = np.random.exponential(scale=mu_fading, size=device_nb)
        else:
            fadings = np.ones(device_nb)

        # With impact of factors such as fading and power control error
        rec_power_levels *= fadings*shadowings

        # 将每个slot，每个设备消耗的能量，存进 energy_history. 我们假设整个slot都以恒定的功率传输数据(不知道这个假设效果如何)
        for device_id, trans_index in enumerate(sim_history[slot, :, 0]):
            sim_history[slot, device_id, int(trans_index)] = rec_power_levels[device_id]

        total_p = sum(rec_power_levels)
        # 我们采用 SINR门限值*干扰值 和 接收功率 比较的方式，加速仿真的执行
        # (计算实际接收信噪比的思路会不可避免地考虑信道中只有一个传输的情况，导致程序的执行效率低下)
        curr_trans_results = [
            (total_p - rec_power_levels[device_id])/rec_power_levels[device_id] > 1.0/(10**(0.1*threshold))
            if k != 0 else False for device_id, k in enumerate(sim_history[slot, :, 0])
        ]

        # print curr_trans_results
        for device_id, curr_trans_result in enumerate(curr_trans_results):
            # 如果 curr_trans_result = True，则需要对这个包执行重传操作
            if curr_trans_result:
                # 需要知道，这是第几次传输失败了
                # Don't forget to convert x into int type.
                x = np.int(sim_history[slot, device_id, 0])
                if x != max_trans:
                    # max_trans has not been reached. Execute backoff procedure
                    # The new slot index is the sum of current one and backoff length
                    while 1:
                        new_slot = int(np.random.exponential(scale=BACK_OFFS[x-1])) + 1 + slot
                        if new_slot < sim_duration:
                        # Take care that the selected new slot should not be out of range.
                        # Scheduled packets out of simualtion is not considered in the following statisticals calculation stage
                            # Also we should note that selected new slot has not yet scheduled for any others retransmission
                            if sim_history[new_slot, device_id, 0] == 0:
                                # copy the content of sim_history[slot, device_id, 0] to
                                # sim_history[new_slot, device_id, 0]
                                sim_history[new_slot, device_id] = sim_history[slot, device_id]
                                # transmission trial should be incremented by 1
                                sim_history[new_slot, device_id, 0] = x+1
                                # Retransmission scheduling work is done. Break from while loop...
                                break
                        # 为什么这里 我要break呢？我觉得应该注释掉啊。。。
                        # 我知道了 是因为到 simulation 的最后阶段，几乎已经不太可能在 simulation 时间范围内找到合适的 new_slot了
                        # 所以 这种情况直接视为 已经在未来找到了合适的slot做 retransmission, 选择直接 break
                        else:
                            break
                    # Do not forget to reinitialize for this (slot, device) pair .
                    sim_history[slot, device_id] = np.zeros(max_trans+1, dtype=np.float)
                elif x == max_trans:
                    # The case where max_trans has been reached. Failure of this packet transmission
                    sim_history[slot, device_id, 0] = max_trans + 1

    # To convert the type of an array, use the .astype() method (preferred) or the type itself as a function.
    # For example: z.astype(float)
    trans_index_array = sim_history[warm_t:sim_duration, :, 0].reshape(device_nb*(sim_duration-warm_t))
    trans_index_array = trans_index_array.astype(int)
    # print "trans_index_array", [e for e in trans_index_array if e != 0]
    # 统计第 i 次传输出现的次数，并据此计算至少需要传输 i 次的频数，如果实验足够长，频数将趋近于概率
    statistics = [[x, 0] for x in range(1, max_trans+2, 1)]

    for entry_1 in itemfreq(trans_index_array)[1:]:
        statistics[entry_1[0]-1][1] = entry_1[1]

    # print "statistics", statistics

    # Convert all received power level history into a one dimension array and calculate the total consumed energies
    total_energies = sum(sim_history[warm_t:sim_duration, :, 1:].reshape(max_trans*device_nb*(sim_duration-warm_t)))

    # 发现了两个bug: 一个是在处理failed packet的时候，错误地把所有的received power level 之和弄成了30
    # 第二个bug是，程序之前并没有处理过simulation的最后一个slot, 结果统计的时候却把这一行尚处于buffered状态的packets都按成功处理了。。。
    # (修改之后 packet loss rate 最后应该会提高一些。。。)
    # for element in sim_history[warm_t:sim_duration]:
    #     print [(hha[0], sum(hha[1:])) for hha in element]
    # print sim_history[warm_t:sim_duration, :, :]
    statistics_vector = [e[1] for e in statistics]

    vector_p = [0 for i in range(max_trans+1)]

    total_transmission = sum(statistics_vector)
    # The following calculation requires that the last element is always the numerb of failed packets!
    succ_transmission = sum(statistics_vector[0:-1])

    for element in statistics:
        ith_trans = element[0]
        satisfied_trans = sum([element[1] for element in statistics if element[0] >= ith_trans])
        vector_p[ith_trans - 1] = satisfied_trans*1.0/total_transmission

    # print "total_energies", total_energies
    # print "succ_transmission", succ_transmission

    energy_efficiency = total_energies/succ_transmission
    statistics_vector.append(energy_efficiency)
    # calculate the throughput: amount of packets succesfully delivered per slot
    statistics_vector.extend(vector_p)
    # Get end time and calculate elapsed time
    end_t = int(time())
    time_elapsed = float(end_t-start_t)/60.0
    print "Time:", int(time_elapsed), "Alpha:", alpha, "Seed:", seed, "Result:", statistics_vector
    return alpha, list(statistics), statistics_vector

def main(sim_config_dict, logs_directory):
    '''
    :param com_config_f: 公共配置字典，存储着通用的仿真参数。
    :param logs_directory:
    :return: nothing
    '''
    #must use Manager queue here, or will not work
    # 注意：一定要使用 Manager queue，否则无法工作，至于原因，我不晓得

    SIM_DURATION = sim_config_dict['SIM_DURATION']
    BINOMIAL_P = sim_config_dict["BINOMIAL_P"]
    # DEVICE_NB = sim_config_dict["DEVICE_NB"]
    ALPHA = sim_config_dict['ALPHA']

    BACKOFF = sim_config_dict["BACKOFF"]
    THERSHOLD = sim_config_dict['THRESLD']

    MAX_TRANS, L, M = sim_config_dict['MAX_TRANS'], sim_config_dict['L'], sim_config_dict['M']
    WARM_T = sim_config_dict['WARM_UP']
    MU_FADING = sim_config_dict['MU_FADING']
    MU_SHADOWING = sim_config_dict['MU_SHADOWING']
    SIGMA_SHADOWING = sim_config_dict['SIGMA_SHADOWING']
    # 针对每个 alpha 值，仿真重复次数
    SIM_REPEAT_NB = sim_config_dict['SIM_REPEAT_NB']
    ALPHAS = [ALPHA for i in range(SIM_REPEAT_NB)]
    DEVICE_NB = int(ALPHA/BINOMIAL_P)
    # 将仿真结果存储在 sim_result_f 指向的文件中
    sim_result_f = os.path.join(
        logs_directory,
        "simd={0}_warm={1}_maxtrans={2}_N={3}_threshold={4}dB_l={5}_m={6}_backoff={7}_alpha={8}_mufading={9}_mushadowing={10}_sigmashadowing={11}_tmp={12}.csv".
        format(
            SIM_DURATION,
            WARM_T,
            MAX_TRANS,
            DEVICE_NB,
            THERSHOLD,
            L,
            M,
            BACKOFF,
            ALPHA,
            MU_FADING,
            MU_SHADOWING,
            SIGMA_SHADOWING,
            strftime("%Y%m%d%H%M%S")
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
            (alpha, MAX_TRANS, BINOMIAL_P, THERSHOLD, L, M, BACKOFF, SIM_DURATION, WARM_T, MU_FADING, MU_SHADOWING, SIGMA_SHADOWING, q)
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
    SIM_CONFIG_DIR = 'sim_configs'
    # SIM_PART = 'fading_shadowing'
    SIM_PART = 'shadowing'
    # SIM_PART = 'perfect'

    SIM_CONFIG_FILE = 'case_K=5_l=1_m=1_threshold=-3dB.json'
    # The simulation result will be logged into files of type CSV, in folder logs.
    # First check the existence of this folder and creat it if necessary.
    if not os.path.exists(logs_directory):
        os.makedirs(logs_directory)
    sim_config_f = os.path.join(SIM_CONFIG_DIR, SIM_PART, SIM_CONFIG_FILE)
    print "Simulation for:", SIM_PART
    print "Now do simulation with configuration file: ", sim_config_f

    sim_config_dict = {}

    with open(sim_config_f) as json_file:
        json_config = json.load(json_file)

    sim_config_dict["BACKOFF"] = json_config["BACKOFF"]
    sim_config_dict['THRESLD'] = json_config['THRESLD']
    sim_config_dict["L"] = json_config['L']
    sim_config_dict["M"] = json_config['M']
    sim_config_dict["MAX_TRANS"] = json_config['MAX_TRANS']
    SIM_DURATION = json_config['SIM_DURATION']
    # WARM_UP attribute in JSON file is a percent value, for example, 10%
    sim_config_dict["MU_FADING"] = json_config['MU_FADING']
    sim_config_dict["MU_SHADOWING"] = json_config['MU_SHADOWING']
    sim_config_dict["SIGMA_SHADOWING"] = json_config['SIGMA_SHADOWING']
    # 针对每个 alpha 值，仿真重复次数
    sim_config_dict["SIM_REPEAT_NB"] = json_config['SIM_REPEAT_NB']
    sim_config_dict["SIM_INCRE_STEP"] = json_config['SIM_INCRE_STEP']
    sim_config_dict["BINOMIAL_P"] = json_config["BINOMIAL_P"]

    ALPHA_START = json_config['ALPHA_START']
    ALPHA_END = json_config['ALPHA_END']
    SIM_INCRE_STEP = json_config['SIM_INCRE_STEP']
    # DEVICE_NB = json_config['DEVICE_NB']
    ALPHA_INTERVAL = np.arange(ALPHA_START, ALPHA_END+0.00001, SIM_INCRE_STEP) # Do not forget to include the ALPHA_END

    for order, ALPHA in enumerate(ALPHA_INTERVAL, 1):
        # 我们一般从 packet loss rate 大致为 0.01 的 ALPHA 值开始，所以前四次仿真，我们把 simulation duration 设为 10000
        # 仿真需要的设备数目因而也不需要很高 (因为 ALPHA 取值不高)
        sim_config_dict["ALPHA"] = ALPHA
        if order < 5:
            # sim_config_dict["DEVICE_NB"] = DEVICE_NB[0]
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[0]
            sim_config_dict["WARM_UP"] = int(json_config['WARM_UP']*0.01*SIM_DURATION[0])

        elif 5 <= order < 16:

            # sim_config_dict["DEVICE_NB"] = DEVICE_NB[1]
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[1]
            sim_config_dict["WARM_UP"] = int(json_config['WARM_UP']*0.01*SIM_DURATION[1])

        else:
            # sim_config_dict["DEVICE_NB"] = DEVICE_NB[-1]
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[-1]
            sim_config_dict["WARM_UP"] = int(json_config['WARM_UP']*0.01*SIM_DURATION[-1])

        print sim_config_dict
        # 将填充好的仿真参数字典传递给 main(), 开启多进程下的仿真
        main(sim_config_dict, logs_directory)

    end_t = int(time())
    time_elapsed = float(end_t - start_t)/60.0
    print "Total Execution time: ", time_elapsed

