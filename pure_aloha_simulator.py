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
import pandas as pd
from scipy.stats import itemfreq

# One process for writing an entry(a row in the context of CSV file) into the given CSV file.
# The advantage is to avoid the concurrence about the write access to target CSV file.
def listener(sim_result_f, q):
    # 经过一晚上的奋战，终于找到BUG所在了！！！
    while 1:
            df_row = q.get()
            # 如果从 (被N多个进程共享的) result queue 中提取到了 字符串 ‘TERMINATE’，那么结束listener进程
            # 注意：无法直接比较对象df_row和"TERMINATE"，会抛出异常退出的！而没有try-catch结构，主进程(本脚本)无法知道异常出现
            if isinstance(df_row, str):
                # if the obtained object from Queue is of type str. prepare to terminate this process.
                # stats_df = pd.read_csv(sim_result_f).sort_index()
                # # overwrite initial csv content with sorted statistics
                # stats_df.to_csv(sim_result_f, mode='w', header=False)
                break
            # I think, within to_csv method, there exists some file open/close operation!!
            # No need to explicitly do open/close operations.
            df_row.to_csv(sim_result_f, mode='a', index=False, header=False)

# All processes other than listener are used to run a simulation
def worker(sim_settings, q):
    sim_result = run_simulation(sim_settings)
    q.put(sim_result)

# This function will be vectorized by Numpy in the corp of function run_simluation(**)
# and run on a matrix where the value in each cell is the square of distance
# That's why, we multiply 0.5 with input path loss exponent.
def path_loss_law(dist_square, path_loss_exponent):
    # $r^{-\gamma}$, where r is distance, \gamma is path-loss component
    return np.power(dist_square, -0.5*path_loss_exponent)


def calculate_sinr(rec_power_vector, threshold, curr_slot_trans_index):
    # curr_slot_trans_index = sim_history[slot, :, 0]
    # 我们采用 SINR门限值*干扰值 和 接收功率 比较的方式，加速仿真的执行
    total_p = np.sum(rec_power_vector)
    trans_state_result = [
        (total_p - rec_power_vector[d_id]) / rec_power_vector[d_id] > 1.0 / (10 ** (0.1 * threshold))
        if k != 0 else False for d_id, k in enumerate(curr_slot_trans_index)]
    return trans_state_result

def process_simultenaous_trans(slot, curr_trans_results, sim_history, max_trans, sim_duration, BACK_OFFS):
    """
    :param slot: the current slot index, for example, 1,2,3,4,5, ...
    :param curr_trans_results: a vector of boolean variable, indicating the current transmission result for each device
    :param sim_history: a matrix recording the simulation results. 0=> no transmission in this slot
    :return: None.
    """
    for device_id, curr_trans_result in enumerate(curr_trans_results):
        # 如果 curr_trans_result = True，则需要对这个包执行重传操作
        if curr_trans_result:
            # 需要知道，这是第几次传输失败了
            # Don't forget to convert x into int type.
            x = np.int(sim_history[slot, device_id, 0])
            if x != max_trans:
                # max_trans has not been reached. Execute backoff procedure
                while 1:
                    # The new slot index is the sum of current one and backoff length
                    new_slot = int(np.random.exponential(scale=BACK_OFFS[x-1])) + 1 + slot
                    if new_slot < sim_duration:
                    # Take care that the selected new slot should not be out of range.
                    # Scheduled packets out of simualtion is not considered in statisticals calculation stage
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

def retransmission():
    pass

def crs_prt(pkt_1, pkt_2):
    """
        Calculate the cross part with other packet
    """
    return max(min(pkt_1[2], pkt_2[2]) - max(pkt_1[1], pkt_2[1]), 0)

def deploy_nodes(width, alpha, intensity_bs):
    # The involved device number will be one sample from a spatial PPP over a finit disk area
    # Generate the needed device nb and base station in this simulation
    AREA_SURFACE = np.pi*np.power(width, 2)
    device_nb = int(np.random.poisson(alpha*AREA_SURFACE, 1))
    # We should have at least one BS
    bs_nb = max(int(np.random.poisson(intensity_bs*AREA_SURFACE, 1)), 1)
    # Uniformelly distribute devices and base stations using polar coordinates.
    # We assume that always one device at the origin
    device_rho = np.concatenate(([0.0], width*np.sqrt(np.random.uniform(0, 1, device_nb-1))))
    device_arguments = np.random.uniform(-np.pi, np.pi, device_nb)
    coordinates_devices_array = zip(device_rho, device_arguments)
    bs_rho = width*np.sqrt(np.random.uniform(0, 1, bs_nb))
    bs_arguments = np.random.uniform(-np.pi, np.pi, bs_nb)
    coordinates_bs_array = zip(bs_rho, bs_arguments)
    return device_nb, bs_nb, coordinates_devices_array, coordinates_bs_array

def nodes_location_process(device_nb, bs_nb, coordinates_devices_array, coordinates_bs_array, path_loss):
    # An intermediare data structure for the device-to-base_station distance.
    d2bs_dist_matrix = np.zeros((bs_nb, device_nb))
    for index, line in enumerate(d2bs_dist_matrix):
        # index => index of BS
        curr_bs_rho, curr_bs_arg = coordinates_bs_array[index][0], coordinates_bs_array[index][1]
        # The value is the square of actual distance, whatever, just for selection of nearest BS
        d2bs_dist_matrix[index] = np.array(
            [np.power(coordinate[0], 2)+np.power(curr_bs_rho, 2)-2*coordinate[0]*curr_bs_rho*np.cos(coordinate[1]-curr_bs_arg)
                 for coordinate in coordinates_devices_array])

    # Allocate the nearest base station, according to the min of distance of each row
    device_base_table = d2bs_dist_matrix.argmin(axis=0)
    # Numpy provides np.vectorize to turn Pythons which operate on numbers into functions that operate on numpy arrays
    vectorized_path_loss_f = np.vectorize(path_loss_law)
    # A maxtrix of size: DEVICE_NB * BS_NB
    path_loss_matrix = vectorized_path_loss_f(d2bs_dist_matrix, path_loss)
    path_loss_matrix = np.transpose(path_loss_matrix)
    return device_base_table, path_loss_matrix


def run_simulation(sim_config_dict):
    """
        In this method, during the backoff slot, it is possible to generate and transmit new packets.
        Only the packets generated during the slot scheduled for retransmission will be abandoned.
        In addition, when choosing the backoff slot, make sure that the choosed slot is not scheduled
        for retransmission. Otherwise continue to choose the backoff slot until the allowed one.
    """
    # Create a minimum float number, which is extremely to zero, cause np.lognormal does not accept 0
    start_t = int(time())
    seed = hash(hash(os.urandom(os.getpid()))) % 4294967295 # np.seed [0, 4294967295]
    np.random.seed(seed)

    F_EPS = np.finfo(float).eps
    METHOD = sim_config_dict['METHOD']
    binomial_p, alpha, intensity_bs = sim_config_dict["BINOMIAL_P"], sim_config_dict['ALPHA'], sim_config_dict["INTENSITY_BS"]
    sim_duration, warm_t, backoff = sim_config_dict['SIM_DURATION'], sim_config_dict['WARM_UP'], sim_config_dict["BACKOFF"]
    threshold = sim_config_dict['THRESLD']
    max_trans, l, m = sim_config_dict['MAX_TRANS'], sim_config_dict['L'], sim_config_dict['M']
    mu_fading,  mu_shadowing, sigma_dB = sim_config_dict['MU_FADING'], sim_config_dict['MU_SHADOWING'], sim_config_dict['SIGMA_SHADOWING']
    sigma_dB += F_EPS
    width, path_loss = sim_config_dict["WIDTH"], sim_config_dict['PATH_LOSS']
    BETA = np.log(10)/10.0
    # The vector format of back off values allows to implement different backoff for each retransmission
    BACK_OFFS = [backoff for i in range(max_trans)]

    device_nb, bs_nb, coordinates_devices_array, coordinates_bs_array = deploy_nodes(width, alpha, intensity_bs)
    device_base_table, path_loss_matrix = nodes_location_process(
        device_nb, bs_nb, coordinates_devices_array, coordinates_bs_array, path_loss
    )

    # Now Add one zero in the beginning of LM list, to facilitate the generation of
    # Then each element in LM list represents the transmit power for the kth transmission (Not retransmission).
    # For example: max_trans = 5, l=2, m=1 => LM = [0, 1, 2, 4, 8, 16] The 0th transmission transmit power is 0
    LM = [1.0*np.power(l, k-1)*np.power(m, max_trans-k) if k != 0.0 else 0.0 for k in range(max_trans+1)]

    # sim_history: a matrix of size: sim_duration * device_nb
    sim_history = np.zeros((sim_duration, device_nb), dtype=np.float)
    # rec_power_history: a matrix of size: sim_duration * device_nb * bs_nb
    # Store the received power level from each device to each BS
    rec_power_history = np.zeros((sim_duration, device_nb, bs_nb), dtype=np.float)
    # trans_power_history: a matrix of size: sim_duration * device_nb * bs_nb,  useful for energy-efficiency analysis
    trans_power_history = np.zeros((sim_duration, device_nb, bs_nb), dtype=np.float)

    # Packet is represented by a 3-elements list: trans_index, start_t, end_t
    last_channel_effect = np.zeros((device_nb, bs_nb), dtype=np.float)

    # Packet is represented by a 3-elements list: trans_index, start_t, end_t. Time-related info
    pkts_history = np.zeros((sim_duration, device_nb, 3), dtype=np.float)

    for slot_index in range(sim_duration-1):
        # At the beginning of each loop, determine which devices will generate new packets
        # Register the device ids generating new packets in current slot
        new_pkt_d_ids = [
            device_id for device_id in range(device_nb)
            if trans_power_history[slot_index, device_id, 0] == 0.0 and bernoulli.rvs(binomial_p) == 1
        ]
        # print "new_pkt_d_ids", new_pkt_d_ids
        for device_id in new_pkt_d_ids:
            trans_power_history[slot_index, device_id, :] = np.array([1]*bs_nb)
            trans_power_history[slot_index+1, device_id, :] = np.array([1]*bs_nb)

            start_t = np.random.uniform(slot_index, slot_index+1)
            pkts_history[slot_index][device_id] = [1, start_t, slot_index+1]
            pkts_history[slot_index+1][device_id] = [1, slot_index+1, start_t+1]

        #update rec_power_history
        channel_variations = np.array([
            np.random.lognormal(BETA*mu_shadowing, BETA*sigma_dB, bs_nb)*np.random.exponential(mu_fading, bs_nb)
            if sum(rand_effect) == 0 else rand_effect for device_id, rand_effect in enumerate(last_channel_effect)
        ])
        # At the end of each iteration, do not to save fading and shadowing values in this iteration
        last_channel_effect = np.array(
            [rand_effect if device_id in new_pkt_d_ids else [0.0]*bs_nb
             for device_id, rand_effect in enumerate(channel_variations)]
        )

        # fadings = np.random.exponential(scale=mu_fading, size=(device_nb, bs_nb))
        rec_power_history[slot_index] = trans_power_history[slot_index]*channel_variations*path_loss_matrix

        # If there exit finished transmission. If yes. proceed them... Else go to next iteration.
        emitted_device_ids=[device_id for device_id, pkt in enumerate(pkts_history[slot_index]) if pkt[0]!= 0 and pkt[2] < slot_index+1]
        if len(emitted_device_ids) > 0:
            for device_id in emitted_device_ids:
                # First calculate cross part. This is very important for mean cumulative interference calcul.
                ref_pkt_1 = pkts_history[slot_index-1][device_id]
                ref_pkt_2 = pkts_history[slot_index][device_id]
                crs_prt_1 = np.array([crs_prt(ref_pkt_1, pkt) for pkt in pkts_history[slot_index-1]])
                crs_prt_2 = np.array([crs_prt(ref_pkt_2, pkt) for pkt in pkts_history[slot_index]])
                crs_prt_1_matrix = np.tile(crs_prt_1, (bs_nb, 1)).transpose()
                crs_prt_2_matrix = np.tile(crs_prt_2, (bs_nb, 1)).transpose()

                # The cumulative interference depends on "mean mode" or "max mode"
                # if ITF_MODE == "MEAN":
                cumu_itf_matrix = rec_power_history[slot_index-1]*crs_prt_1_matrix \
                                    +rec_power_history[slot_index]*crs_prt_2_matrix                # else:

                ref_rec_power_vector = rec_power_history[slot_index, device_id]
                # Sum on the basis of column, i.e., cumulative power for each BS. Thus len(cumu_itf_vector) == bs_nb
                # To avoid 0 in the denominator, add an extremly small value
                cumu_itf_vector = cumu_itf_matrix.sum(axis=0) + F_EPS - ref_rec_power_vector
                # The failure state of a device at all BS.
                # For BS_RX_DIVERS, if any of this list is false (i.e., successful transmission)
                # For BS_NST_ATT, we just care about the state of attached BS.
                ref_rec_sinr_vector = [sinr < np.power(10, threshold/10) for sinr in ref_rec_power_vector/cumu_itf_vector]

                # Process with received SINR at the BS side.
                if METHOD == "BS_NST_ATT":
                    if ref_rec_sinr_vector[device_base_table[device_id]]:
                        # True => sinr < thresold => failure => retransmission procedure
                        #Todo: to be implemented. This method. For one-shot acess, we can leave it empty...
                        retransmission()
                        sim_history[slot_index][device_id] = 2.0

                    else:
                        # successful transmission. Do some clean work.
                        sim_history[slot_index][device_id] = 1.0

                elif METHOD == "BS_RX_DIVERS":
                    if False in ref_rec_sinr_vector:
                    # If there exist False in ref_rec_sinr_vector => a certain BS has received the packet!!
                        sim_history[slot_index][device_id] = 1.0
                    else:
                        retransmission()
                        sim_history[slot_index][device_id] = 2.0

                else:
                    print "No simulated method is specified (BS_NST_ATT or BS_RX_DIVERS), program exit. Please check you simulation configuration JSON file."
                    exit()

                # For each finished packet, do not forget to reset transmit power in last slot as 0
                trans_power_history[slot_index-1][device_id] = 0.0


    statistics_vector = sim_result_statistics(
        sim_config_dict, device_nb, coordinates_devices_array, sim_history, trans_power_history
    )
    end_t = int(time())
    time_elapsed = float(end_t-start_t)/60.0
    print "Time:", int(time_elapsed), "Alpha:", alpha, "Seed:", seed, "Result:"
    print statistics_vector.to_string(index=False, header=False)
    return statistics_vector

def sim_result_statistics(sim_config_dict, device_nb, coordinates_devices_array, sim_history, trans_power_history):
    # Remove the statistics in the border area. Particularly important for BS reception diversity.
    # print "In statistics.s..."
    # print sim_history[:,0]
    # print trans_power_history[:, 0, 0]
    width = sim_config_dict["WIDTH"]
    sim_duration, warm_t = sim_config_dict['SIM_DURATION'], sim_config_dict['WARM_UP']
    max_trans = sim_config_dict['MAX_TRANS']
    global_statistics = []
    stat_percents = sim_config_dict["STAT_PERCENTS"]
    for stat_percent in stat_percents:
        accepted_device_index = np.array([1.0 if coordinate[0] <= width*stat_percent else 0.0 for coordinate in coordinates_devices_array])
        trans_index_array = sim_history[warm_t:sim_duration, :]*np.tile(accepted_device_index, (sim_duration-warm_t, 1))
        # To convert the type of an array, use the .astype() method (preferred) or the type itself as a function.
        # For example: z.astype(float)
        trans_index_array = trans_index_array.reshape(device_nb*(sim_duration-warm_t))
        trans_index_array = trans_index_array.astype(int)
        # Count the occurrence of transmission trial with index x.
        # e..g, x=1 refers to the first transmission trial, x=max_trans is the last transmission trial.
        # All failed transmssion will be with index max_trans+1
        statistics = [[x, 0] for x in range(1, max_trans+2, 1)]
        statistics_vector = [int(stat_percent*100), sim_config_dict['ALPHA']]
        item_pairs = itemfreq(trans_index_array)[1:]  # count for items >= 1
        # Attention! Must iterate for item_pairs instead of statistics: item_pairs may just contain counts for 0 and 1
        for trans_index, element in enumerate(item_pairs):
            statistics[trans_index][1] = element[1]

        item_counts = [e[1] for e in statistics]
        vector_p = [0 for i in range(max_trans+1)]
        statistics_vector.extend(item_counts)

        total_transmission = sum(item_counts)
        # The following calculation requires that the last element is always the numeber of failed packets!
        succ_transmission = sum(item_counts[0:-1])
        for element in statistics:
            ith_trans = element[0]
            satisfied_trans = sum([element[1] for element in statistics if element[0] >= ith_trans])
            vector_p[ith_trans - 1] = satisfied_trans*1.0/total_transmission
        statistics_vector.extend(vector_p)

        # Convert all received power level history into a one dimension array and calculate the total consumed energies
        # We iterate max_trans times (i.e., max_trans+1-1) to summarize the consumed energies by devices within
        # statistical region.
        # TODO: There must be some way to avoid the following loop...

        energy_index_array = trans_power_history[warm_t:sim_duration, :, 0]*np.tile(accepted_device_index, (sim_duration-warm_t, 1))
        energy_index_array = energy_index_array.reshape(device_nb*(sim_duration-warm_t))
        total_energies = energy_index_array.sum()

        energy_efficiency = total_energies/succ_transmission
        statistics_vector.append(energy_efficiency)
        # calculate the throughput: amount of packets succesfully delivered per slot
        global_statistics.append(statistics_vector)

    statistics_df = pd.DataFrame(global_statistics)
    return statistics_df

def main(sim_config_dict, logs_directory):
    """
    :param sim_config_dict:     dictionary, simulation settings.
    :param logs_directory:      target directory to save the simualtion-generated csv logs
    :return: nothing
    """
    # must use Manager queue here, or will not work
    # 注意：一定要使用 Manager queue，否则无法工作，至于原因，我不晓得
    max_trans = sim_config_dict['MAX_TRANS']
    if max_trans == 1:
        # If one-shot access
        # Some parameters, such as warm time, l, m, back-off time, no need to know...
        output_csv_f_name = \
            "METHOD={METHOD}_TRLD={THRESLD}_ALPHA={ALPHA}_BSDENSITY={INTENSITY_BS}_FADING={MU_FADING}_SHADOW={SIGMA_SHADOWING}_R={WIDTH}_".format(**sim_config_dict)
        output_csv_f_name +="TMP={}.csv".format(strftime("%Y%m%d%H%M%S"))

    else:
        output_csv_f_name ="METHOD={METHOD}_SIMD={SIM_DURATION}_WARM={WARM_UP}_MAXTR={MAX_TRANS}_BSITSY={INTENSITY_BS}_TRLD={THRESLD}dB_L={L}_M={M}_BACKOFF={BACKOFF}_".format(**sim_config_dict)
        output_csv_f_name += "ALPHA={ALPHA}_FADING={MU_FADING}_SHADOW={SIGMA_SHADOWING}_R={WIDTH}_".format(**sim_config_dict)
        output_csv_f_name +="TMP={}.csv".format(strftime("%Y%m%d%H%M%S"))

    sim_result_f = os.path.join(logs_directory, output_csv_f_name)

    manager = mp.Manager()
    q = manager.Queue()
    # sim_config_dict['"NB_PROCESS"'] at least 2, one producer, one consumer.
    pool = mp.Pool(sim_config_dict["NB_PROCESS"])

    # Create a csv file from scratch, use csv module to write the label of each column into this file.
    # This should be finished in main processing. Other worker processes just writing data, don't care column label
    # Thus, make sure that in method sim_result_statistics(), the list has the same order with column label.
    df_column = ["STAT_PERCENT", "ALPHA"]
    df_column.extend(["{0}nb".format(element+1) for element in range(max_trans+1)])
    df_column.extend(["{0}pr".format(element+1) for element in range(max_trans+1)])
    # Be careful!!! If the element order of statistics_vector changed, do not forget to change the column label order!!!
    df_column.append("EE")
    # Note that, if MAX_TRANS=2, "2pr" is the column of outage probability.
    with open(sim_result_f, 'w') as f:
        csv_w = csv.writer(f, delimiter=',')
        csv_w.writerow(df_column)

    # put listener to work first
    watcher = pool.apply_async(listener, (sim_result_f, q,))

    #fire off workers
    jobs = []

    for i in range(sim_config_dict['SIM_REPEAT_NB']):
        job = pool.apply_async(
            worker,
            (sim_config_dict, q)
        )
        jobs.append(job)

    #collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #Now we are done, kill the listener
    q.put("TERMINATE")
    pool.close()
    # Do not forget to close file at the end.
    # f_handler.close()

if __name__ == "__main__":
    START_T= int(time())
    logs_directory = 'logs'
    SIM_CONFIG_DIR = 'sim_configs'
    # SIM_PART = 'fading_shadowing'
    # SIM_PART = 'perfect'
    SIM_PART = 'bs_rx_divers'
    SIM_CONFIG_FILE = 'case_AVGITF_K=1threshold=3dB.json'
    # Check the existence of "logs" folder and create it if necessary.
    if not os.path.exists(logs_directory):
        os.makedirs(logs_directory)

    sim_config_f = os.path.join(SIM_CONFIG_DIR, SIM_PART, SIM_CONFIG_FILE)
    print "Simulation for:", SIM_PART
    print "Now do simulation with configuration file: ", sim_config_f

    sim_config_dict = {}

    with open(sim_config_f) as json_file:
        sim_config_dict = json.load(json_file)

    SIM_DURATION = sim_config_dict['SIM_DURATION']
    ALPHA_START = sim_config_dict['ALPHA_START']
    ALPHA_END = sim_config_dict['ALPHA_END']
    SIM_INCRE_STEP = sim_config_dict['SIM_INCRE_STEP']
    ALPHA_INTERVAL = np.arange(ALPHA_START, ALPHA_END+0.00001, SIM_INCRE_STEP) # Do not forget to include the ALPHA_END

    #在脚本区，遍历从ALPHA_START到ALPHA_END中的每一个ALPHA值，进行仿真，因而字典sim_config_dict在每次循环中，有几个参数都会被修改。
    for order, ALPHA in enumerate(ALPHA_INTERVAL, 1):
        # 我们一般从 packet loss rate 大致为 0.01 的 ALPHA 值开始，所以前四次仿真，我们把 simulation duration 设为 10000
        # 仿真需要的设备数目因而也不需要很高 (因为 ALPHA 取值不高)
        sim_config_dict["ALPHA"] = ALPHA
        if order < 5:
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[0]
            sim_config_dict["WARM_UP"] = int(sim_config_dict['WARM_UP']*0.01*SIM_DURATION[0])

        elif 5 <= order < 16:
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[1]
            sim_config_dict["WARM_UP"] = int(sim_config_dict['WARM_UP']*0.01*SIM_DURATION[1])

        else:
            sim_config_dict["SIM_DURATION"] = SIM_DURATION[-1]
            sim_config_dict["WARM_UP"] = int(sim_config_dict['WARM_UP']*0.01*SIM_DURATION[-1])

        # Just print the simulation settings.
        if sim_config_dict['MAX_TRANS'] == 1:
            # If one-shot access
            # Some parameters, such as warm time, l, m, back-off time, no need to know...
            output_csv_f_name = \
                "METHOD={METHOD}_TRLD={THRESLD}_ALPHA={ALPHA}_BSDENSITY={INTENSITY_BS}_FADING={MU_FADING}_SHADOW={SIGMA_SHADOWING}_R={WIDTH}".format(**sim_config_dict)

        else:
            output_csv_f_name ="METHOD={METHOD}_SIMD={SIM_DURATION}_WARM={WARM_UP}_MAXTR={MAX_TRANS}_BSITSY={INTENSITY_BS}_TRLD={THRESLD}dB_l={L}_m={M}_backoff={BACKOFF}_".format(**sim_config_dict)
            output_csv_f_name += "alpha={ALPHA}_mufading={MU_FADING}_shadowing={SIGMA_SHADOWING}_R={WIDTH}_".format(**sim_config_dict)

        print output_csv_f_name
        main(sim_config_dict, logs_directory)

    END_T = int(time())
    time_elapsed = float(END_T - START_T)/60.0
    print "Total Execution time: ", time_elapsed

