# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
from scipy.stats import bernoulli
import csv
from time import time
import multiprocessing

class Packet(object):
    '''
        The class Packet modelize the packet sent by device to the eNB. It should contain
        three attributes:
        @pkt_id         :the unique identifier of this packet;
        @device_id      :the unique identifier of the device sending this packet;
        @tr_id          :transmission index of this packet, meaning that this packet has been transmitted how
                            much times.
        @power          :the transmit power associated with this transmission
        @timer          :Attribute timer refer to the retransmission timer. -1 means the timer is not still armed.
                        0 refers to the expiration of timer
    '''
    def __init__(self, pkt_id, device_id, tr_id, power, timer):
        self.pkt_id = pkt_id
        self.device_id = device_id
        self.tr_id = tr_id
        self.power = power
        self.timer = timer

    def __str__(self):
        return "Packet({0}, {1}, {2}, {3}, {4})".format(self.pkt_id, self.device_id, self.tr_id, self.power, self.timer)

    def __eq__(self, other):
        return (isinstance(other, self.__class__)) \
               and (self.__dict__ == other.__dict__) \
               and self.pkt_id == other.pkt_id \
               and self.device_id == other.device_id \
               and self.tr_id == other.tr_id \
               and self.power == other.power \
               and self.timer == other.timer

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "Packet({0}, {1}, {2}, {3}, {4})".format(self.pkt_id, self.device_id, self.tr_id, self.power, self.timer)


class Device(object):
    """
        The device should keep two separate buffer, one is for fresh transmission packet,
        the other is for collided packets.

        We should choose one transmit strategy: always process one packet until its succesful transmission
        or process other new packet transmission and

        一些不得不重视的陷阱：仔细想来，failed packet 应该是 Device 来处理比较好，确保合法的 packet 才可以发送到 Channel 上去
    """

    def __init__(self, index, proba, POWER_LEVELS, MAX_TRANS, BACKOFF):
        self.index = index
        self.proba = proba
        self.buffer = []
        self.packets = []
        self.fail_pkts = []
        self.succ_pkts = []
        self.sent_pkts = []
        self.POWER_LEVELS = POWER_LEVELS
        self.MAX_TRANS = MAX_TRANS
        self.BACKOFF = BACKOFF


    def __str__(self):
        return "Device: {0} Pakcet buffer length: {1} Buffer length: {2}".format(
            self.index, len(self.packets), len(self.buffer)
        )

    def generate_pkt(self, round_index):
        # No matter whether the current round is scheduled for retransmission,
        # we should check if new packet is arrived and save it into the buffer
        # Each device is modelized one packet source respecting to Bernoulli distribution.
        # The probability is defined as $\lambda/n$, where $\lambda$ is the fresh random access
        # poisson distribution mean, n is the total device number.
        if bernoulli.rvs(self.proba) == 1:
            pkt = Packet(round_index, self.index, 1, self.POWER_LEVELS[0], -1)
            self.packets.append(pkt)

        # Then we should decide transmit a new packet or a collided packet.
        # First judge whether the current round is a scheduled round for retransmission

    def get_pkt(self, round_index):
        '''
            If current round is scheduled for retransmission, return the collided packet
            buffered in corresponding buffer;
            If current round is unfortunately scheduled as retransmission occasion for multiple collided packets,
            return always the frist packet; the remainings will be processed in next round;
            If current round is not scheduled for retransmission, check the buffer for fresh
            packet, if not empty, try to transmit the first packet in this buffer.

            注意：发送一个 packet 之后，这个被发送的 packet 应该存储在 channel 的 sent_packet list 中来
            不过 要不要在 device里面也存一下呢？同理发送成功的 packet 要不要也存在 device里呢？存一下好了。。。这样
            可以统计每一个设备的 发送成功率 以及 packet loss rate...
        '''
        # First judge whether the current round is a scheduled round for retransmission
        # If the retransmission timer is not zero (timer is armed, cause each round this timer decrements by one),
        # in this round, we can transmit fresh packe
        # if the buffer for fresh packets is not empty.

        # If the collided packets buffer is empty and fresh packets buffer is not empty. It is time to transmit the
        # first packet queued in fresh packets buffer

        # We need to know which packet to be transmitted.

        # 4 combination probabilities of buffer length and fresh packets length
        pkt = None

        if len(self.buffer) == 0 and len(self.packets) != 0:
                # case only fresh packets are available..
                # first, pop one packet from fresh packets buffer
                pkt = self.packets.pop()

        elif len(self.buffer) != 0:
            # case collided packets are available...
            # we first check if there are some packets in collided packets buffer. If yes we send
            pkts = [element for element in self.buffer if element.timer == 0]
            if len(pkts) != 0:
                # print "case where collided packet is available...", [str(element) for element in self.buffer]

                # there exist some packets in collided packets buffer whose timer expired
                # pop the first element of them
                pkt = pkts[0]
                self.buffer.pop(self.buffer.index(pkt))
                for element in self.buffer:
                    if element not in pkts:
                        element.timer -= 1

            else:
                # print "case where no collided packets available...", [str(element) for element in self.buffer]
                # Decrement one for all collided packets
                for element in self.buffer:
                    element.timer -= 1
                # if len(self.packets) != 0:
                #     pkt = self.packets.pop()
                # else:
                    # case neither collided nor fresh packets are available
                    # print "\t\tDo nothing in transmit method of device ", self.index

        if pkt != None:
            # print "SLOT {0}:\t\t{1} is transmited".format(round_index, pkt)
            self.sent_pkts.append(pkt)
            # self.channel.rec_pkt(pkt, round_index)

        # 没有任何 packet 可以发送的时候，我们什么也不做
        # 最后别忘了产生新的 packets, 留待下一个 round 去发送
        self.generate_pkt(round_index)
        return pkt



    def backoff(self, packet, round_index):
        '''
            注意：Backoff 算法应该同时包含了对于达到 最大传输次数的 包的处理工作。保证 Device 各种 Buffer 中的 packet 都
            是需要发的包

            其实本质上来讲，backoff 应该算是一个 packet generator 方法。即：根据 input 给定的 packet 的生成一个新的packet，
            只放在重传的 Buffer 里面去
        '''
        # Note that if packet transmission index reaches to 5, abandon this packet
        # backoff does not taken into accout the packet with transmission index with MAX_TRANS
        if packet.tr_id != self.MAX_TRANS:
            # 如果需要被 backoff 的 packet 还没有达到最大传输数目，则创建一个新的packet放到buffer中去
            # In fact, the backoff window should be a CONSTANT system parameter
            p = Packet(
                    packet.pkt_id,
                    packet.device_id,
                    packet.tr_id+1,
                    self.POWER_LEVELS[packet.tr_id],
                    # np.random.randint(1, self.BACKOFF))
                    int(np.random.exponential(scale=self.BACKOFF))
            )
            self.buffer.append(p)

            # print "SLOT {0}:\t\t{1} is transmitted but backlogged, scheduled to make the {2} transmission in SLOT {3}.".format(round_index, packet, p.tr_id, round_index+p.timer+1)

        else:
            # 需要被 backoff 的 packet 已经达到了最大重传次数，直接放弃
            # 放入相应的 list 之中去, 无需其他操作
            self.fail_pkts.append((round_index, packet))
            # print "SLOT {0}:\t\t{1} is transmitted but failed and dropped.".format(round_index, packet)

    def ack(self, pkt, round_index):
        '''
            这个方法主要就是把 packet 放入到被成功发送的 list 中去
        '''
        self.succ_pkts.append((round_index, pkt))
        # print "SLOT {0}:\t\t{1} is transmitted with success.".format(round_index, pkt)

class Channel(object):
    '''
        The class Channel simulates the physical channel which conveys the packet from devcies
        to eNodeB.
        经过观察分析发现，Xavier建议的 warm-up period非常有必要。比如前 100 个实验 rounds里确实出现了非常多的成功的 transmission
        而我们的 analytical model 则考虑的是稳态的情况。
    '''
    def __init__(self, devices, MAX_TRANS):
        self.devices = devices
        self.history = []
        # A variable to store all succesfully sent packets.
        self.succ_pkts =[]
        # Store all the dropped packets
        self.fail_pkts = []
        self.MAX_TRANS = MAX_TRANS


    def rec_pkt(self, round_index, threld):
        rec_pkts = [device.get_pkt(round_index) for device in self.devices]
        # Remove the None of above list
        rec_pkts = [pkt for pkt in rec_pkts if pkt != None]
        # 将实验某 round 期间受到的 pkts 存储在一个list中
        self.history[round_index] = rec_pkts

        if len(rec_pkts) != 0:  # namely  len(rec_pkts) >= 1
        # 如果 cumulative interference 以后不包含本身的话，要把1 还有 >1的情况分开讨论
            total_p = sum([e.power for e in rec_pkts])
            for pkt in rec_pkts:
                if total_p <= pkt.power/threld:
                    self.devices[pkt.device_id].ack(pkt, round_index)
                else:
                    self.devices[pkt.device_id].backoff(pkt, round_index)

    def get_sent_pkts(self, warm_up_t):
        result = []
        for i, rec_pkts in enumerate(self.history):
            if i > warm_up_t:
                result.extend(rec_pkts)
        return result

    def get_buffered_pkts(self):
        # Do not need to consider the warm up time issue, since at the end of the simulation, the packets in buffer
        # are surely sent after the warm up period
        result = []
        for device in self.devices:
            result.extend(device.buffer)
        return result


    def get_fail_pkts(self, warm_up_t):
        result = []
        for device in self.devices:
            result.extend([element[1] for element in device.fail_pkts if element[0] > warm_up_t])
        return result

    def get_succ_pkts(self, warm_up_t):
        result = []
        for device in self.devices:
            result.extend([element[1] for element in device.succ_pkts if element[0] > warm_up_t])
        return result


    def create_one_slot(self):
        self.history.append([])



    def statistics(self, warm_up_t):
        #Statistics 2: 在成功还有失败的包里统计
        succ = self.get_succ_pkts(warm_up_t)
        fail = self.get_fail_pkts(warm_up_t)
        buffered = self.get_buffered_pkts()
        t_pkt_nb = len(succ) + len(fail) + len(buffered)
        all = self.get_sent_pkts(warm_up_t)

        result1 = [
            (sum([1 for pkt in succ if pkt.tr_id >= i+1])+len(fail)
             +sum([1 for pkt in buffered if pkt.tr_id >= i+1]))*1.0/t_pkt_nb for
                i in range(self.MAX_TRANS)
        ]

        # result2 = [len([pkt for pkt in all if pkt.tr_id >= i+1])*1.0/len(all) for i in range(MAX_TRANS)]

        p_loss = len(fail)*1.0/(len(fail)+len(succ))

        result1.append(p_loss)

        return result1


def run_simulation(channel, threld, seed, warm_up_t, sim_duration):
    '''
        @sim_para:
        @seed: a random number generator seed to guarantee that each simulation will generate a different r.v.s
    '''
    np.random.seed(seed)

    for i in range(sim_duration):
        channel.create_one_slot()
        channel.rec_pkt(i, threld)
    # print len(channel.history), [len(e) for e in channel.history]
    # print len(channel.history[warm_up_t:]), [len(e) for e in channel.history[warm_up_t:]]
    # print channel.history
    # print "SUCC:", len(channel.get_succ_pkts(warm_up_t)), [str(pkt) for pkt in channel.get_succ_pkts(warm_up_t)]
    # print "FAIL:", len(channel.get_fail_pkts(warm_up_t)), [str(pkt) for pkt in channel.get_fail_pkts(warm_up_t)]
    return channel.statistics(warm_up_t)


def iteration(alpha_start, alpha_end, sim_step, THRESLD, N, POWER_LEVELS, MAX_TRANS, SIM_NB, WARM_UP, SIM_DURATION):
    results = []
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

        # Calculate the average of all simulations
        row = [float("{0:.4g}".format(p)) for p in np.mean(np.array(tmp_array), axis=0)]
        row.append(alpha_start)
        result.append(row)
        results.append(result)
        alpha_start += sim_step
        print alpha_start

    return results


if __name__ == "__main__":
    start = time()
    # The basic system parameters for simulation.
    # MAX_TRANS: the maximum allowed transmission number of a packet
    MAX_TRANS = 5
    l = 1.0
    m = 1.0
    SIM_DURATION = 5000
    # The POWER_LEVELS: the possible transmit power (received power for packet)
    POWER_LEVELS = [l**k*m**(MAX_TRANS-1-k) for k in range(MAX_TRANS)]
    THLD = -3.0
    THLD = 10**(THLD/10.0)
    THLD = 0.5
    N = 500
    intensity = 0.1
    end_intensity = 2.05
    step = 0.05
    WARM_UP = 1000
    SIM_NB = 1

    # 真他妈的 Bizart啊。。。我直接设定 intensity = 0.8 算出来的 丢包率是 0.9 从0.1开始，现在就接近于0了。。。什么世道

    result_f = "sim_{0}_{1}_result_with_{2}_test.csv".format(N, SIM_NB, THLD)

    with open(result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')


        while intensity <= end_intensity:
            proba = intensity/N
            result = []

            pool = multiprocessing.Pool(SIM_NB)

            # Populate the task list
            tasks =[]
            for n in range(SIM_NB):
                devices = [Device(i, proba, POWER_LEVELS, MAX_TRANS) for i in range(N)]
                channel = Channel(devices, MAX_TRANS)
                tasks.append((channel, THLD, int(time()+n*100), WARM_UP, SIM_DURATION, ))

            # Start all tasks
            result = [pool.apply_async(run_simulation, t) for t in tasks]


            # Iterate the final result
            tmp_array = np.array([prob_vector.get() for prob_vector in result])
            # for prob_vector in result:
            #     print prob_vector.get()

            # Calculate the average of all simulations
            row = [float("{0:.4g}".format(p)) for p in np.mean(np.array(tmp_array), axis=0)]
            print "result", intensity, row

            row.append(intensity)

            spamwriter.writerow(row)

            intensity += step


    end = time()

    print "Execution time in minutes: ", end-start, " namely: ", float(end-start)/60.0, "minutes."
