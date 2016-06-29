# -*- coding: utf-8 -*-
__author__ = 'qsong'
import json
import os


from single_channel_aloha_analytical_method import *
from single_channel_aloha_simulator import *

def run_analytical_simulation(config_f):
    # The case of slotted aloha with threshold -3dB, without power increment
    with open(config_f) as json_file:
        json_config = json.load(json_file)

    ## The section for running analytical method
    np.set_printoptions(precision=4)
    MAX_TRANS = json_config['MAX_TRANS']
    alpha_start = json_config['alpha_start']
    l = json_config['l']
    m = json_config['m']
    DELTA = json_config['DELTA']
    THRESLD = json_config['THRESLD']
    alpha_end =json_config['alpha_end']
    ana_step = json_config['ana_step']

    P = np.zeros(MAX_TRANS)
    P[0] = 1

    ana_result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, THRESLD, ana_step)
    ana_result_f = os.path.join("data_files", "analytical_result_threshold_{0}_l={1}_m={2}.csv".format(THRESLD, l, m))
    with open(ana_result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')
        for n, vector_p in enumerate(ana_result, 1):
            print n, vector_p
            spamwriter.writerow(vector_p)


    ## The section for running simulation
    SIM_DURATION = json_config['SIM_DURATION']
    # The POWER_LEVELS: the possible transmit power (received power for packet)
    POWER_LEVELS = [l**k*m**(MAX_TRANS-1-k) for k in range(MAX_TRANS)]
    N = json_config['N']
    sim_step = json_config['sim_step']
    WARM_UP = json_config['WARM_UP']
    SIM_NB = json_config['SIM_NB']


    sim_result_f = os.path.join("data_files", "sim_result_sim_{0}_N_{1}_threshold_{2}_l={3}_m={4}.csv".format(SIM_NB, N, THRESLD, l, m))
    sim_result = iteration(alpha_start, alpha_end, sim_step, THRESLD, N, POWER_LEVELS, MAX_TRANS, SIM_NB, WARM_UP, SIM_DURATION)
    with open(sim_result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')
        for n, row in enumerate(sim_result, 1):
            print n, row
            spamwriter.writerow(row)



if __name__ == "__main__":
    # run_analytical_simulation('case_l=1_m=1_threshold=-3dB_N=500.json')
    # # The case of slotted aloha with threshold -3dB, without power increment
    config_f = os.path.join('exp_configs', 'case_l=1_m=1_threshold=0dB_N=500.json')
    print "Now do simulation with configuration file: ", config_f
    with open(config_f) as json_file:
        json_config = json.load(json_file)

    ## The section for running analytical method
    np.set_printoptions(precision=4)
    MAX_TRANS = json_config['MAX_TRANS']
    alpha_start = json_config['alpha_start']
    l = json_config['l']
    m = json_config['m']
    DELTA = json_config['DELTA']
    THRESLD = json_config['THRESLD']
    alpha_end =json_config['alpha_end']
    ana_step = json_config['ana_step']
    print ana_step

    P = np.zeros(MAX_TRANS)
    P[0] = 1

    ana_result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, THRESLD, ana_step)
    ana_result_f = "analytical_result_threshold={0}_l={1}_m={2}.csv".format(THRESLD, l, m)
    with open(ana_result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')
        for n, vector_p in enumerate(ana_result, 1):
            print n, vector_p
            spamwriter.writerow(vector_p)


    ## The section for running simulation
    SIM_DURATION = json_config['SIM_DURATION']
    # The POWER_LEVELS: the possible transmit power (received power for packet)
    POWER_LEVELS = [l**k*m**(MAX_TRANS-1-k) for k in range(MAX_TRANS)]
    N = json_config['N']
    sim_step = json_config['sim_step']
    WARM_UP = json_config['WARM_UP']
    SIM_NB = json_config['SIM_NB']


    # sim_result_f = "sim_result_sim_{0}_N_{1}_threshold_{2}_l={3}_m={4}.csv".format(SIM_NB, N, THRESLD, l, m)
    # sim_result = iteration(alpha_start, alpha_end, sim_step, THRESLD, N, POWER_LEVELS, MAX_TRANS, SIM_NB, WARM_UP, SIM_DURATION)
    # with open(sim_result_f, 'w') as f_handler:
    #     spamwriter = csv.writer(f_handler, delimiter=',')
    #     for n, row in enumerate(sim_result, 1):
    #         print n, row
    #         spamwriter.writerow(row)


    # 真他妈的 Bizart啊。。。我直接设定 intensity = 0.8 算出来的 丢包率是 0.9 从0.1开始，现在就接近于0了。。。什么世道

    sim_result_f = "sim_result_simd={0}_N={1}_threshold={2}_l={3}_m={4}.csv".format(SIM_DURATION, N, THRESLD, l, m)

    with open(sim_result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')

        while alpha_start <= alpha_end:
            result = []

            pool = multiprocessing.Pool(SIM_NB)

            # Populate the task list
            tasks =[]
            for n in range(SIM_NB):
                devices = [Device(i, alpha_start/N, POWER_LEVELS, MAX_TRANS) for i in range(N)]
                channel = Channel(devices, MAX_TRANS)
                tasks.append((channel, THRESLD, int(time()+n*100), WARM_UP, SIM_DURATION, ))

            # Start all tasks
            result = [pool.apply_async(run_simulation, t) for t in tasks]

            # Iterate the final result
            tmp_array = np.array([prob_vector.get() for prob_vector in result])
            # for prob_vector in result:
            #     print prob_vector.get()

            # Calculate the average of all simulations
            row = [float("{0:.4g}".format(p)) for p in np.mean(np.array(tmp_array), axis=0)]
            print "result", alpha_start, row

            row.append(alpha_start)

            spamwriter.writerow(row)

            alpha_start += sim_step





