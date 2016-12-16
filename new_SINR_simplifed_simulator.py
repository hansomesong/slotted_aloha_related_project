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
import glob
import pprint
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import itemfreq

# Only one process responsible for writing an entry(a row in the context of CSV file) into the given CSV file.
# The advantage is to avoid the concurrence about the write access to target CSV file.
def listener(sim_result_f, q):
    with open(sim_result_f, 'a') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        while 1:
            csv_row = q.get()
            # 如果从 (被N多个进程共享的) result queue 中提取到了 字符串 ‘TERMINATE’，那么结束listener进程
            if csv_row == "TERMINATE":
                break
            spamwriter.writerow(csv_row)

# Except the one in charge of writing records into CSV file, all rest processes are used
# to run a simulation
def worker(alpha, max_trans, binomial_p, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing, sigma_shadowing, width, intensity_bs, path_loss, q):
    csv_row = []
    sim_result = run_simulation(
        alpha, max_trans, binomial_p, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing, sigma_shadowing, width, intensity_bs, path_loss
    )
    csv_row.extend(sim_result[2])
    csv_row.append(alpha)
    q.put(csv_row)

# This function will be vectorized by Numpy in the corp of function run_simluation(**)
# and run on a matrix where the value in each cell is the square of distance
# That's why, we multiply 0.5 with input path loss exponent.
def path_loss_law(dist_square, path_loss_exponent):
    # The reason to multiply 1000 before is to avoid that channel gain is much too small.
    return 1e3*np.power(dist_square, -0.5*path_loss_exponent)

def calculate_sinr(rec_power_vector, threshold, curr_slot_trans_index):
    # curr_slot_trans_index = sim_history[slot, :, 0]
    # 我们采用 SINR门限值*干扰值 和 接收功率 比较的方式，加速仿真的执行
    # (计算实际接收信噪比的思路会不可避免地考虑信道中只有一个传输的情况，导致程序的执行效率低下)
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


def is_retransmited(a):
    result = True
    for state in a:
        result = result and state
    return result

def run_simulation(alpha, max_trans, binomial_p, threshold, l, m, backoff, sim_duration, warm_t, mu_fading, mu_shadowing,
                   sigma_shadowing, width, intensity_bs, path_loss, output_statistics=False
    ):
    """
        In this method, during the backoff slot, it is possible to generate and transmit new packets.
        Only the packets generated during the slot scheduled for retransmission will be abandoned.
        In addition, when choosing the backoff slot, make sure that the choosed slot is not scheduled
        for retransmission. Otherwise continue to choose the backoff slot until the allowed one.
    """
    start_t = int(time())
    seed = hash(hash(os.urandom(os.getpid()))) % 4294967295 # np.seed [0, 4294967295]
    np.random.seed(seed)

    BETA = np.log(10)/10.0
    # The vector format of back off values allows to implement different backoff for each retransmission
    BACK_OFFS = [backoff for i in range(max_trans)]
    # The involved device number will be one sample from a spatial PPP over a finit disk area
    # Generate the needed device nb and base station in this simulation
    AREA_SURFACE = np.pi*np.power(width, 2)
    device_nb = int(np.random.poisson(alpha*AREA_SURFACE, 1))
    bs_nb = max(int(np.random.poisson(intensity_bs*AREA_SURFACE, 1)), 1)

    # Uniformelly distribute devices and base stations using polar coordinates
    device_rho = width*np.sqrt(np.random.uniform(0, 1, device_nb))
    device_arguments = np.random.uniform(-np.pi-10e-4, np.pi+10e-4, device_nb)
    coordinates_devices_array = zip(device_rho, device_arguments)
    bs_rho = width*np.sqrt(np.random.uniform(0, 1, bs_nb))
    bs_arguments = np.random.uniform(-np.pi-10e-4, np.pi+10e-4, bs_nb)
    coordinates_bs_array = zip(bs_rho, bs_arguments)

    # An intermediare data structure for the device-to-base_station distance.
    d2bs_dist_matrix = np.zeros((bs_nb, device_nb))
    for index, line in enumerate(d2bs_dist_matrix):
        curr_bs_rho, curr_bs_arg = coordinates_bs_array[index][0], coordinates_bs_array[index][1]
        # The value is the square of actual distance, whatever, just for selection of nearest BS
        d2bs_dist_matrix[index] = np.array(
            [np.power(coordinate[0], 2)+np.power(curr_bs_rho, 2)-2*coordinate[0]*curr_bs_rho*np.cos(coordinate[1]-curr_bs_arg)
                 for coordinate in coordinates_devices_array])

    # Allocate the nearest base station, according to the min of distance of each row
    device_base_table = d2bs_dist_matrix.argmin(axis=0)
    # Numpy provides np.vectorize to turn Pythons which operate on numbers into functions that operate on numpy arrays
    vectorized_path_loss_f = np.vectorize(path_loss_law)
    path_loss_matrix = vectorized_path_loss_f(d2bs_dist_matrix, path_loss)

    # print "device_rho:", coordinates_devices_array
    # print coordinates_bs_array
    # print bs_nb, device_nb, "over", AREA_SURFACE
    # print d2bs_dist_matrix
    # print path_loss_matrix

    # Now Add one zero in the beginning of LM list, to facilitate the generation of
    # Then each element in LM list represents the transmit power for the kth transmission (Not retransmission).
    # For example: max_trans = 5, l=2, m=1 => LM = [0, 1, 2, 4, 8, 16] The 0th transmission transmit power is 0
    LM = [1.0*np.power(l, k-1)*np.power(m, max_trans-k) if k != 0.0 else 0.0 for k in range(max_trans+1)]

    # sim_history now is 3-d array.
    # The third dimension is used to represent the transmission index and received power levels history
    sim_history = np.zeros((sim_duration, device_nb, max_trans+1), dtype=np.float)

    # 根据是否考虑fading, shadowing，分类讨论，避免在遍历sim_history时候，每次都要执行if-else语句
    if mu_fading > 0 and sigma_shadowing == 0:
        for slot in range(sim_duration):
            # First generate new packets. k refers to the k th transmision instead of retransmission
            sim_history[slot, :, 0] = np.array(
                [bernoulli.rvs(binomial_p) if k == 0.0 else k for k in sim_history[slot, :, 0]]
            )
            # print "In slot ", slot, sim_history[slot, :, 0]
            rec_power_levels = np.tile(np.array([LM[int(k)] for k in sim_history[slot, :, 0]]), (bs_nb, 1))

            fadings = np.random.exponential(scale=mu_fading, size=(bs_nb, device_nb))
            # 将每个slot，每个设备消耗的能量，存进 energy_history. 我们假设整个slot都以恒定的功率传输数据(不知道这个假设效果如何)
            for device_id, trans_index in enumerate(sim_history[slot, :, 0]):
                sim_history[slot, device_id, int(trans_index)] = rec_power_levels[0, device_id]

            # With impact of factors such as fading and power control error
            # Also take into account the path loss. Size of rec_power_levels:(bs_nb, device_nb)
            rec_power_levels *= fadings*path_loss_matrix

            # If the element is True, means the transmission of corresponding device is failed. Need retransmission or drop
            # If the element is False, means no message to transmit or successful transmission
            # This vector will be executed "and" operation with each base station associated transmission state vector
            # If one transmission state is "False", even others are all True,
            # then this message is received by some BS with success. Not need retransmission
            curr_trans_matrix = np.apply_along_axis(calculate_sinr, 1, rec_power_levels, threshold, sim_history[slot, :, 0])
            # The nearest-base-station approache
            # curr_trans_results = np.array([curr_trans_matrix[device_base_table[i], i] for i in range(device_nb)])

            # The fire-and-forget approach
            curr_trans_results = np.apply_along_axis(is_retransmited, 0, curr_trans_matrix)

            # Iterate curr_trans_results to proceed retransmission...
            process_simultenaous_trans(slot, curr_trans_results, sim_history, max_trans, sim_duration, BACK_OFFS)
            # print rec_power_levels
            # print curr_trans_matrix
            # print "In slot ", slot, "transmission result", curr_trans_results
            # print sim_history[slot, :, 0]

    elif mu_fading > 0 and sigma_shadowing > 0:
        print "in fading and shadowing case..."
        for slot in range(sim_duration):
            # First generate new packets. k refers to the k th transmision instead of retransmission
            sim_history[slot, :, 0] = np.array(
                [bernoulli.rvs(binomial_p) if k == 0.0 else k for k in sim_history[slot, :, 0]]
            )
            rec_power_levels = np.tile(np.array([LM[int(k)] for k in sim_history[slot, :, 0]]), (bs_nb, 1))
            shadowings = np.random.lognormal(BETA*mu_shadowing, BETA*sigma_shadowing, size=(bs_nb, device_nb))
            fadings = np.random.exponential(scale=mu_fading, size=(bs_nb, device_nb))
            # 将每个slot，每个设备消耗的能量，存进 energy_history. 我们假设整个slot都以恒定的功率传输数据(不知道这个假设效果如何)
            for device_id, trans_index in enumerate(sim_history[slot, :, 0]):
                sim_history[slot, device_id, int(trans_index)] = rec_power_levels[0, device_id]

            # With impact of factors such as fading and power control error
            # Also take into account the path loss
            # Actually rec_power_levels is a bs_nb X device_nb matrix
            # print rec_power_levels.shape, fadings.shape, path_loss_matrix.shape, shadowings.shape
            rec_power_levels *= fadings*shadowings*path_loss_matrix

            curr_trans_matrix = np.apply_along_axis(calculate_sinr, 1, rec_power_levels, threshold, sim_history[slot, :, 0])
            # The nearest-base-station approache
            # curr_trans_results = np.array([curr_trans_matrix[device_base_table[i], i] for i in range(device_nb)])

            # The fire-and-forget approach
            curr_trans_results = np.apply_along_axis(is_retransmited, 0, curr_trans_matrix)
            # Iterate curr_trans_results to proceed retransmission...
            process_simultenaous_trans(slot, curr_trans_results, sim_history, max_trans, sim_duration, BACK_OFFS)

    # Simulation finished here. Now do statistics work.
    # Before the statistics, remove the statistics in the border area.
    # The following can be grouped into a separate function
    accepted_device_index = np.array([1.0 if coordinate[0] <= width*1.0 else 0.0 for coordinate in coordinates_devices_array])
    trans_index_array = sim_history[warm_t:sim_duration, :, 0]*np.tile(accepted_device_index, (sim_duration-warm_t, 1))
    # print sim_history[warm_t:sim_duration, :, 0]
    # print coordinates_devices_array
    # print coordinates_bs_array, "BS"
    # print "seleted: ", np.tile(accepted_device_index, (sim_duration-warm_t, 1))
    # To convert the type of an array, use the .astype() method (preferred) or the type itself as a function.
    # For example: z.astype(float)
    trans_index_array = trans_index_array.reshape(device_nb*(sim_duration-warm_t))
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

    WIDTH = sim_config_dict["WIDTH"]
    INTENSITY_BS = sim_config_dict["INTENSITY_BS"]
    PATH_LOSS = sim_config_dict['PATH_LOSS']

    # 针对每个 alpha 值，仿真重复次数
    SIM_REPEAT_NB = sim_config_dict['SIM_REPEAT_NB']
    ALPHAS = [ALPHA for i in range(SIM_REPEAT_NB)]
    # DEVICE_NB = int(ALPHA/BINOMIAL_P)
    # 将仿真结果存储在 sim_result_f 指向的文件中
    sim_result_f = os.path.join(
        logs_directory,
        "simd={0}_warm={1}_maxtrans={2}_bsintensity={3}_threshold={4}dB_l={5}_m={6}_backoff={7}_alpha={8}_mufading={9}_mushadowing={10}_sigmashadowing={11}_R={12}_tmp={13}.csv".
        format(
            SIM_DURATION,
            WARM_T,
            MAX_TRANS,
            INTENSITY_BS,
            THERSHOLD,
            L,
            M,
            BACKOFF,
            ALPHA,
            MU_FADING,
            MU_SHADOWING,
            SIGMA_SHADOWING,
            WIDTH,
            strftime("%Y%m%d%H%M%S")
        )
    )

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(6)

    f_handler = open(sim_result_f, 'w')

    # put listener to work first
    watcher = pool.apply_async(listener, (sim_result_f, q,))

    #fire off workers
    jobs = []

    for alpha in ALPHAS:
        job = pool.apply_async(
            worker,
            (alpha, MAX_TRANS, BINOMIAL_P, THERSHOLD, L, M, BACKOFF, SIM_DURATION, WARM_T, MU_FADING, MU_SHADOWING, SIGMA_SHADOWING, WIDTH, INTENSITY_BS, PATH_LOSS, q)
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
    SIM_PART = 'fading'
    # SIM_PART = 'perfect'

    SIM_CONFIG_FILE = 'case_K=1_l=1_m=1_threshold=3dB.json'
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
    sim_config_dict["INTENSITY_BS"] = json_config["INTENSITY_BS"]
    sim_config_dict["WIDTH"] = json_config["WIDTH"]
    sim_config_dict["PATH_LOSS"] = json_config["PATH_LOSS"]

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

