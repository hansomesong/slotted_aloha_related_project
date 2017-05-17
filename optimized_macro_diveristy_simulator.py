# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
from scipy.stats import bernoulli
from time import time
import multiprocessing as mp
import csv
import os
import sys
import json
from time import strftime
import pandas as pd
from scipy.stats import itemfreq

# One process for writing an entry(a row in the context of CSV file) into the given CSV file.
# The advantage is to avoid the concurrence about the write access to target CSV file.
def listener(sim_result_f, q):
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
# and run on a matrix where the value in each cell is the distance between a BS and device
def path_loss_law(dist, path_loss_exponent):
    # $r^{-\gamma}$, where r is distance, \gamma is path-loss component
    return np.power(dist, -1.0*path_loss_exponent)


def calculate_sinr(rec_power_vector, threshold, curr_slot_trans_index):
    """

    :param rec_power_vector:                a list of received power at each BS;
    :param threshold:                       capture ration, unit dB, typically 3dB
    :param curr_slot_trans_index:           a list of transmission index, 0=> no transmission, 1=>first trial
    :return:
            trans_state_result:             a list of bool type. True=> transmission False=> Success or no transmission
    """
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
        # curr_trans_result = True <=> SINR < threshold <=> Transmission Failure <=> Retransmission
        if curr_trans_result:
            # Need to know which transmission trial failed? the first one? second one?...
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


def is_retransmited(a):
    """
        This method is uniquely useful for BS reception diversity approach, where just one of BS in the simulation
        says no need to do retransmission. This method returns a False, which means this transmission is successful.
    :param a:
    :return:
    """
    """
    :param a:
    :return:
    """
    result = True
    for state in a:
        result = result and state
    return result

def output_simulation_csv():
    """
        This method is used to write the simulation details into a TxT file
        The format of output file is as follows:
        #BS #DEVICE #AREA
        BS_COORDINATES
        DEVICES_COORAN

        To be finised...
    """
    return 0

def get_involved_nodes_nb(width, alpha, intensity_bs, coord_type="Polar"):
    """
    :param width:                       the radius of simulation area
    :param alpha:                       the active device intensity (Note: not transmitting device intensity)
    :param intensity_bs:                the Base Station spatial intensity
    :return:
    """
    if coord_type == "Polar":
        AREA_SURFACE = np.pi*np.power(width, 2)
        # device_nb = int(np.random.poisson(alpha*AREA_SURFACE, 1))                   # At least one device
        # bs_nb = max(int(np.random.poisson(intensity_bs*AREA_SURFACE, 1)), 1)        # At least one BS
        device_nb = int(alpha*AREA_SURFACE)                   # At least one device
        bs_nb = int(intensity_bs*AREA_SURFACE)   # At least one BS
    elif coord_type == "Cartesian":
        AREA_SURFACE = np.power(width, 2)
        device_nb = int(np.random.poisson(alpha*AREA_SURFACE, 1))                   # At least one device
        bs_nb = max(int(np.random.poisson(intensity_bs*AREA_SURFACE, 1)), 1)        # At least one BS
    else:
        sys.exit("Error! Unknown coordinates type. Either Polar or Cartesian")
    return device_nb, bs_nb

def deploy_nodes(width, device_nb, bs_nb, coord_type="Polar"):
    """
    :param width:                       the radius of simulation area
    :param device_nb:                   the active device number (Note: not transmitting device number)
    :param bs_nb:                       the Base Station number
    :param coord_type:                  the coordination type of simulation area. either Polar Coordinates  or Cartesian Coordinates
    :return:
            coordinates_devices_array: the list of tuple, each tuple represents for the coordination pair for a device
            coordinates_bs_array:      the list of tuple, each tuple represents for the coordination pair for a BS
    """
    # The involved device number will be one sample from a spatial PPP over a finite disk area
    # Generate the needed device nb and base station in this simulation
    if coord_type == "Polar":
        # Uniformly distribute devices and base stations using polar coordinates.
        # We assume that always one device at the origin
        device_rho = np.concatenate(([0.0], width*np.sqrt(np.random.uniform(0, 1, device_nb-1))))
        device_arguments = np.random.uniform(-np.pi, np.pi, device_nb)
        coordinates_devices_array = zip(device_rho, device_arguments)
        bs_rho = width*np.sqrt(np.random.uniform(0, 1, bs_nb))
        bs_arguments = np.random.uniform(-np.pi, np.pi, bs_nb)
        coordinates_bs_array = zip(bs_rho, bs_arguments)
    elif coord_type == "Cartesian":
        device_x = np.random.uniform(0, width, device_nb)
        device_y = np.random.uniform(0, width, device_nb)
        coordinates_devices_array = zip(device_x, device_y)
        bs_x = np.random.uniform(0, width, bs_nb)
        bs_y = np.random.uniform(0, width, bs_nb)
        coordinates_bs_array = zip(bs_x, bs_y)
    else:
        sys.exit("Error! Unknown coordinates type. Either Polar or Cartesian")

    return coordinates_devices_array, coordinates_bs_array

def nodes_location_process(device_nb, bs_nb, coordinates_devices_array, coordinates_bs_array, path_loss, coord_type="Polar"):
    """

    :param device_nb:                       the number of devices involved in simulation
    :param bs_nb:                           the number of BS involved in simulation
    :param coordinates_devices_array:       the list of tuple for device coordinates
    :param coordinates_bs_array:            the list of tuple for Base Station coordinates
    :param path_loss:                       the path-loss function, the most simple form is r^{-gamma}
    :param coord_type:                      the coordinates type, either polar or Cartesian
    :return: path_loss_matrix:              a matrix of dimension device_nb * bs_nb,
                                            An intermediare data structure for the device-to-base_station distance.
    """
    d2bs_dist_matrix = \
        calculate_2bs_dist_matrix(bs_nb, device_nb, coordinates_bs_array, coordinates_devices_array, coord_type)

    # Allocate the geographically nearest base station, according to the min of distance of each column
    # device_base_table save the index of BS to which the device attaches
    device_base_table = d2bs_dist_matrix.argmin(axis=0)
    # Numpy provides np.vectorize to turn Pythons which operate on numbers into functions that operate on numpy arrays
    vectorized_path_loss_f = np.vectorize(path_loss_law)
    # A maxtrix of size: DEVICE_NB * BS_NB
    path_loss_matrix = vectorized_path_loss_f(d2bs_dist_matrix, path_loss)
    #TODO: why we need to transpose matrix "path_loss_matrix"? I think not necessary. comment it!
    # path_loss_matrix = np.transpose(path_loss_matrix)
    return d2bs_dist_matrix, device_base_table, path_loss_matrix

def calculate_2bs_dist_matrix(bs_nb, device_nb, coordinates_bs_array, coordinates_devices_array, coord_type="Polar"):
    """
        This method is used to calculate device-to-BS distance matrix.
        Each row is a distance vector for a BS to each device
        Each column is a distance vector for a device to each BS
    :param bs_nb:
    :param device_nb:
    :param coordinates_bs_array:
    :param coordinates_devices_array:
    :param coord_type:
    :return: d2bs_dist_matrix:  matrix of distance, each cell is a distance between a device and BS pair.
    """
    # An intermediare data structure for the device-to-base_station distance.
    d2bs_dist_matrix = np.zeros((bs_nb, device_nb))
    if coord_type == "Polar":
        for index, line in enumerate(d2bs_dist_matrix):
            # index => index of BS
            curr_bs_rho, curr_bs_arg = coordinates_bs_array[index][0], coordinates_bs_array[index][1]
            # The value is the square of actual distance, whatever, just for selection of nearest BS
            tmp_array = np.array(
                [np.power(coordinate[0], 2)+np.power(curr_bs_rho, 2)-2*coordinate[0]*curr_bs_rho*np.cos(coordinate[1]-curr_bs_arg)
                     for coordinate in coordinates_devices_array])
            d2bs_dist_matrix[index] = np.sqrt(tmp_array)

    elif coord_type == "Cartesian":
        for index, line in enumerate(d2bs_dist_matrix):
            # index => index of BS
            curr_bs_x, curr_bs_y = coordinates_bs_array[index][0], coordinates_bs_array[index][1]
            # The value is the square of actual distance, whatever, just for selection of nearest BS
            tmp_array = np.array(
                [np.power(coordinate[0] - curr_bs_x, 2)+np.power(coordinate[1]-curr_bs_y, 2) for coordinate in coordinates_devices_array])
            d2bs_dist_matrix[index] = np.sqrt(tmp_array)

    return d2bs_dist_matrix



def f(recv_p, cumu_itf, threshold):
    if recv_p > 0:
        return cumu_itf / recv_p > 1.0 / (10 ** (0.1 * threshold))
    else:
        return False

def run_simulation(sim_config_dict):
    """
        In this method, during the backoff slot, it is possible to generate and transmit new packets.
        Only the packets generated during the slot scheduled for retransmission will be abandoned.
        In addition, when choosing the backoff slot, make sure that the selected slot is not scheduled
        for retransmission. Otherwise continue to choose the backoff slot until the allowed one.
    :param sim_config_dict:                 the dict which contains all simulation settings.
    :return:
    """
    start_t = int(time())
    seed = hash(hash(os.urandom(os.getpid()))) % 4294967295 # np.seed [0, 4294967295]
    np.random.seed(seed)

    # Create an extremely small float number and add it to shadowing standard error => np.lognormal does not accept 0
    F_EPS = np.finfo(float).eps
    METHOD = sim_config_dict['METHOD']
    binomial_p, alpha, intensity_bs = sim_config_dict["BINOMIAL_P"], sim_config_dict['ALPHA'], sim_config_dict["INTENSITY_BS"]
    sim_duration, warm_t, backoff = sim_config_dict['SIM_DURATION'], sim_config_dict['WARM_UP'], sim_config_dict["BACKOFF"]
    threshold = sim_config_dict['THRESLD']
    max_trans, l, m = sim_config_dict['MAX_TRANS'], sim_config_dict['L'], sim_config_dict['M']
    mu_fading,  mu_shadowing, sigma_dB = sim_config_dict['MU_FADING'], sim_config_dict['MU_SHADOWING'], sim_config_dict['SIGMA_SHADOWING']
    width, path_loss = sim_config_dict["WIDTH"], sim_config_dict['PATH_LOSS']
    coord_type = sim_config_dict['COORD_TYPE']

    sigma_dB += F_EPS
    BETA = np.log(10)/10.0
    sigma = BETA*sigma_dB

    device_nb, bs_nb = get_involved_nodes_nb(width, alpha, intensity_bs, coord_type)
    coordinates_devices_array, coordinates_bs_array = deploy_nodes(width, device_nb, bs_nb, coord_type)

    d2bs_dist_matrix, device_base_table, path_loss_matrix = \
        nodes_location_process(device_nb, bs_nb, coordinates_devices_array, coordinates_bs_array, path_loss, coord_type)
    #TODO: we can consider to introduce a factor as power control component...

    sim_history = np.ones(sim_duration, dtype=int)

    for slot in range(sim_duration):
        # first determine the location of device and BS
        coordinates_devices_array, coordinates_bs_array = deploy_nodes(width, device_nb, bs_nb, coord_type)

        d2bs_dist_matrix, device_base_table, path_loss_matrix = \
        nodes_location_process(device_nb, bs_nb, coordinates_devices_array, coordinates_bs_array, path_loss, coord_type)

        # Generate and calculate the cumulative interference from all devices
        trans_power_levels = np.random.binomial(1, binomial_p, device_nb)
        trans_power_matrix = np.tile(trans_power_levels, (bs_nb, 1))                          # matrix bs_ns * device_nb
        shadowings = np.random.lognormal(BETA*mu_shadowing, sigma, size=(bs_nb, device_nb))
        fadings = np.random.exponential(scale=mu_fading, size=(bs_nb, device_nb))
        rec_power_matrix = trans_power_matrix*fadings*shadowings*path_loss_matrix
        cumu_itf_levels = np.sum(rec_power_matrix, axis=1)  # A vector of bs_nb elements.Interference received at each BS

        shadowings0 = np.random.lognormal(BETA*mu_shadowing, sigma, bs_nb)
        fadings0 = np.random.exponential(mu_fading, bs_nb)
        #TODO: it just need one vector instead of matrix here.
        # path_loss_matrix[:, 0] is for the device at the origin to each BS
        rec_power_levels0 = shadowings0 * fadings0 * path_loss_matrix[:, 0]


        f2 = np.vectorize(f, otypes=[np.bool])
        curr_trans_levels = f2(rec_power_levels0, cumu_itf_levels, threshold)


        if METHOD in ["BS_NST_ATT", "BS_BEST_ATT"]:
        # The nearest-base-station approach
        # The BS selected is the one to which the physical distance between device and BS is smallest
        # Allocate the nearest base station, according to the min of distance of each row
            # device_base_table[i] => index of attached BS, index of row
            if curr_trans_levels[device_base_table[0]]:
                sim_history[slot] = 2

        elif METHOD == "BS_RX_DIVERS":
        # The fire-and-forget approach
        # curr_trans_results will be a scalar, true of false
            curr_trans_results = is_retransmited(curr_trans_levels)
            if curr_trans_results:
                # true => transmission failure, need retransmission
               sim_history[slot] = 2

        else:
            print "No simulated method is specified (BS_NST_ATT or BS_RX_DIVERS), program exit. " \
                  "Please check you simulation configuration JSON file."
            exit()

    statistics_vector = sim_result_statistics(sim_config_dict, device_nb, coordinates_devices_array, sim_history)
    end_t = int(time())
    time_elapsed = float(end_t-start_t)/60.0
    print "Time:", int(time_elapsed), "Alpha:", alpha, "Seed:", seed, "Result:"
    print statistics_vector.to_string(index=False, header=False)
    return statistics_vector

def sim_result_statistics(sim_config_dict, device_nb, coordinates_devices_array, sim_history):
    # Remove the statistics in the border area. Particularly important for BS reception diversity.
    width = sim_config_dict["WIDTH"]
    sim_duration, warm_t = sim_config_dict['SIM_DURATION'], sim_config_dict['WARM_UP']
    max_trans = sim_config_dict['MAX_TRANS']
    sim_history = sim_history.astype(int)
    statistics = [[x, 0] for x in range(1, max_trans+2, 1)]
    item_pairs = itemfreq(sim_history)  # count for items >= 1. The element in sim_history is atleat 1
    statistics_vector = [0, sim_config_dict['ALPHA']]

    for element in item_pairs:
        trans_index = element[0] - 1
        statistics[trans_index][1] = element[1]

    item_counts = [e[1] for e in statistics]
    vector_p = [0 for i in range(max_trans+1)]
    statistics_vector.extend(item_counts)

    total_transmission = sum(item_counts)
    for element in statistics:
        ith_trans = element[0]
        satisfied_trans = sum([element[1] for element in statistics if element[0] >= ith_trans])
        vector_p[ith_trans - 1] = satisfied_trans*1.0/total_transmission

    statistics_vector.extend(vector_p)
    # if we use syntax such as: pd.DataFrame(statistics_vector).
    # the content of statistic_vector will be written into a column
    statistics_df = pd.DataFrame([statistics_vector])
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
            "SLOT_METHOD={METHOD}_TRLD={THRESLD}_P={BINOMIAL_P}_ALPHA={ALPHA}_BSDENSITY={INTENSITY_BS}_FADING={MU_FADING}_SHADOW={SIGMA_SHADOWING}_R={WIDTH}_".format(**sim_config_dict)
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
    SIM_CONFIG_FILE = 'case_K=1_l=1_m=1_threshold=3dB.json'
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
                "Slotted_METHOD={METHOD}_TRLD={THRESLD}_P={BINOMIAL_P}_ALPHA={ALPHA}_BSDENSITY={INTENSITY_BS}_" \
                "FADING={MU_FADING}_SHADOW={SIGMA_SHADOWING}_R={WIDTH}".format(**sim_config_dict)

        else:
            output_csv_f_name ="METHOD={METHOD}_SIMD={SIM_DURATION}_WARM={WARM_UP}_MAXTR={MAX_TRANS}_BSITSY={INTENSITY_BS}_TRLD={THRESLD}dB_l={L}_m={M}_backoff={BACKOFF}_".format(**sim_config_dict)
            output_csv_f_name += "alpha={ALPHA}_mufading={MU_FADING}_shadowing={SIGMA_SHADOWING}_R={WIDTH}_".format(**sim_config_dict)

        print output_csv_f_name
        main(sim_config_dict, logs_directory)

    END_T = int(time())
    time_elapsed = float(END_T - START_T)/60.0
    print "Total Execution time: ", time_elapsed

