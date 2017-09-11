__author__ = 'qsong'

from scipy.special import gamma as gamma_f
import numpy as np
from scipy import optimize
from scipy.misc import comb as comb
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
params = {
    'legend.fontsize': 15,
    "lines.markersize" : 4,
    'lines.linewidth' : 1,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 20,
    'legend.numpoints': 1
}
plt.rcParams.update(params)

from analytical_model import sgam

# Some common constant declaration
EPS = np.finfo(np.float64).eps
BETA = np.log(10.0)/10.0

FIGSIZE = (10, 6)

def outage_distance(p_outage, lambda_b, sigma_dB, thetha_dB, gamma):

    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)

    return np.sqrt(np.log(1.0/(p_outage)) / (lambda_b*np.pi*fm_shadowing))

def func(lambda_m, p, thetha_dB, gamma, r, lambda_b, p_outage, pure=False, itf_mean=True):
    """

        Fixed function to be solved.
    """
    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    L = p*lambda_m/lambda_b
    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = factor*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    term_1 = np.power(p*np.pi*A*np.power(THETA, 2.0/gamma)*fm_shadowing*np.power(r, 2), -1)
    term_2 = np.log(1-np.exp(-p*lambda_m*np.pi*A*np.power(THETA, 2.0/gamma)*fm_shadowing*np.power(r, 2)*np.power(r, 2))) - np.log(p_outage)
    term_3 = np.log(A*np.power(THETA, 2.0/gamma)*L*term_2)
    result = -1.0*term_1*term_3
    return result


def func_plr(plr, p, thetha_dB, gamma, lambda_b, lambda_m, N, pure=False, itf_mean=True):

    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    thetha_term = np.power(THETA, 2.0/gamma)
    L = p*lambda_m/lambda_b
    G = (1-plr)*L/(1-np.power(plr, 1.0/N))
    # G = *L
    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = factor*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    return np.exp(np.sum([comb(N, i)*np.power(-1, i)*np.power(i*A*thetha_term*G, -1) for i in range(1, N+1)]))


def func_best_att_plr(plr, p, thetha_dB, gamma, lambda_b, lambda_m, N, pure=False, itf_mean=True):

    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    thetha_term = np.power(THETA, 2.0/gamma)
    L = p*lambda_m/lambda_b
    G = (1-plr)*L/(1-np.power(plr, 1.0/N))
    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = factor*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    return np.sum([comb(N, i)*np.power(-1, i)*np.power(i*A*thetha_term*G + 1, -1) for i in range(0, N+1)])


def single_attach(p, thetha_dB, gamma, r, p_outage, pure=False, itf_mean=True):
    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = factor*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    term_2 = -1.0*np.log(1-p_outage)
    term_1 = np.power(p*np.pi*A*np.power(THETA, 2.0/gamma)*fm_shadowing*np.power(r, 2), -1)
    return term_1 * term_2


def macro_diversity_fixed_retransmission(p, thetha_dB, gamma, lambda_b, lambda_m, N, pure=False, itf_mean=True):
    THETA = np.power(10, thetha_dB/10)
    thetha_term = np.power(THETA, 2.0/gamma)
    L = p*lambda_m/lambda_b
    G = N*L
    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = factor*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    return np.exp(np.sum([comb(N, i)*np.power(-1, i)*np.power(i*A*thetha_term*G, -1) for i in range(1, N+1)]))

if __name__ == "__main__":

    MARKER_EVERY = 5
    lambda_b, sigma_dB, thetha_dB, gamma = 0.08, 8, 3.0, 4
    THETA = np.power(10, thetha_dB/10)
    p = 0.008


    fig2, axes2 = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    fig1, axes1 = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)
    fig0, axes0 = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)



    lambda_m_list  = np.linspace(0.5, 5, 100)
    # lambda_m_list  = [1.0]
    N_vector = [1, 2, 3, 4]

    plr_list_div_array = []
    plr_list_best_att_array =[]
    plr_list_div_fix_array = []

    for N in N_vector:
        print "N %s" % N
        plr_list_div = []
        plr_list_best_att =[]
        plr_list_div_fix = []
        for lambda_m in lambda_m_list:

            tmp0 = optimize.fixed_point(
                func_plr,
                [0.1],
                args=(p, thetha_dB, gamma, lambda_b, lambda_m, N, True, False)
            )
            plr_list_div.append(tmp0[0])

            tmp1 = optimize.fixed_point(
                func_best_att_plr,
                [0.1],
                args=(p, thetha_dB, gamma, lambda_b, lambda_m, N, False, True)
            )
            plr_list_best_att.append(tmp1[0])
            plr_list_div_fix.append(macro_diversity_fixed_retransmission(p, thetha_dB, gamma, lambda_b, lambda_m, N, True, False))

        print "macro diveristy", plr_list_div
        print "best attacg", plr_list_best_att
        print "fixed macro diversity", plr_list_div_fix
        plr_list_div_array.append(plr_list_div)
        plr_list_div_fix_array.append(plr_list_div_fix)
        plr_list_best_att_array.append(plr_list_best_att)


    index_N = 0


    axes0.plot(
        p*lambda_m_list/lambda_b,
        plr_list_best_att_array[index_N],
        color='b',  marker='o', linestyle='-', markevery=MARKER_EVERY, label=r"Best BS attach, slotted ALOHA, $N_{max}$ = %s " % N_vector[index_N]
    )
    axes1.plot(
        p*lambda_m_list/lambda_b,
        plr_list_div_array[index_N],
        color='b',  marker='o', linestyle='-', markevery=MARKER_EVERY, label=r"SC Macro diversity, pure ALOHA"+"\n"+"maximum interference, $N_{max}$ = %s " % N_vector[index_N]
    )
    axes2.plot(
        p*lambda_m_list/lambda_b,
        plr_list_div_fix_array[index_N],
        color='b',  marker='o', linestyle='-', markevery=MARKER_EVERY, label=r"SC Macro diversity, pure ALOHA"+"\n"+"maximum interference, $N_{max}$ = %s " % N_vector[index_N]
    )

    index_N = 1
    axes0.plot(
        p*lambda_m_list/lambda_b,
        plr_list_best_att_array[index_N],
        color='g',  marker='s', linestyle='--', markevery=MARKER_EVERY, label=r"Best BS attach, slotted ALOHA, $N_{max}$ = %s " % N_vector[index_N]
    )
    axes1.plot(
        p*lambda_m_list/lambda_b,
        plr_list_div_array[index_N],
        color='g',  marker='s', linestyle='--', markevery=MARKER_EVERY, label=r"SC Macro diversity, pure ALOHA"+"\n"+"maximum interference, $N_{max}$ = %s " % N_vector[index_N]
    )
    axes2.plot(
        p*lambda_m_list/lambda_b,
        plr_list_div_fix_array[index_N],
        color='g',  marker='s', linestyle='--', markevery=MARKER_EVERY, label=r"SC Macro diversity, pure ALOHA"+"\n"+"maximum interference, $N_{max}$ = %s " % N_vector[index_N]
    )



    index_N = 2
    axes0.plot(
        p*lambda_m_list/lambda_b,
        plr_list_best_att_array[index_N],
        color='c',  marker='^', linestyle='-.', markevery=MARKER_EVERY-1, label=r"Best BS attach, slotted ALOHA, $N_{max}$ = %s " % N_vector[index_N]
    )
    axes1.plot(
        p*lambda_m_list/lambda_b,
        plr_list_div_array[index_N],
        color='c',  marker='^', linestyle='-.', markevery=MARKER_EVERY-1, label=r"SC Macro diversity, pure ALOHA"+"\n"+"maximum interference, $N_{max}$ = %s " % N_vector[index_N]
    )
    axes2.plot(
        p*lambda_m_list/lambda_b,
        plr_list_div_fix_array[index_N],
        color='c',  marker='^', linestyle='-.', markevery=MARKER_EVERY-1, label=r"SC Macro diversity, pure ALOHA"+"\n"+"maximum interference, $N_{max}$ = %s " % N_vector[index_N]
    )
    #



    index_N = 3
    axes0.plot(
        p*lambda_m_list/lambda_b,
        plr_list_best_att_array[index_N],
        color='k',  marker='v', linestyle=':', markevery=MARKER_EVERY+3, label=r"Best BS attach, slotted ALOHA, $N_{max}$ = %s " % N_vector[index_N]
    )
    axes1.plot(
        p*lambda_m_list/lambda_b,
        plr_list_div_array[index_N],
        color='k',  marker='v', linestyle=':',  markevery=MARKER_EVERY+3, label=r"SC Macro diversity, pure ALOHA"+"\n"+"maximum interference, $N_{max}$ = %s " % N_vector[index_N]
    )
    axes2.plot(
        p*lambda_m_list/lambda_b,
        plr_list_div_fix_array[index_N],
        color='k',  marker='v', linestyle=':',  markevery=MARKER_EVERY+3, label=r"SC Macro diversity, pure ALOHA"+"\n"+"maximum interference, $N_{max}$ = %s " % N_vector[index_N]
    )
    #


    # axes.plot(
    #     p*lambda_m_list/lambda_b,
    #     sgam.bs_rx_div_op(lambda_m_list, lambda_b, gamma, p, thetha_dB, 8, True, False),
    #     color='k',  marker='', linestyle='--', linewidth=2, label="Macro diversity, analytical"
    # )

    print sgam.bs_rx_div_op(1.0, 0.08, gamma, 0.008, thetha_dB, 8, True, False),
    axes0.set_yscale("log")
    axes0.set_xscale('linear')
    axes0.set_xlabel(r"Fresh Normalized Load")
    axes0.set_ylabel(r"Packet Loss Rate")
    axes0.set_ylim(1e-3, 0.8)
    axes0.set_xlim(0.02, 0.5)
    axes0.grid()

    axes1.set_yscale("log")
    axes1.set_xscale('linear')
    axes1.set_xlabel(r"Fresh Normalized Load")
    axes1.set_ylabel(r"Packet Loss Rate")
    axes1.set_ylim(1e-3, 0.8)
    axes1.set_xlim(0.02, 0.5)
    axes1.grid()

    axes2.set_yscale("log")
    axes2.set_xscale('linear')
    axes2.set_xlabel(r"Fresh Normalized Load")
    axes2.set_ylabel(r"Packet Loss Rate")
    axes2.set_ylim(1e-3, 0.8)
    axes2.set_xlim(0.02, 0.5)
    axes2.grid()


    axes0.legend(loc='best')
    axes1.legend(loc='best')
    axes2.legend(loc='best')
    #
    # shift = max([t.get_window_extent().width for t in axes2.legend.get_texts()])
    #
    # for t in axes2.legend.get_texts():
    #     t.set_ha('left') # ha is alias for horizontalalignment
    #     t.set_position((shift, 0))

    fig2.savefig('fixed_normalized_load_vs_packet_loss_rate_macro_diversity_N=%s.eps' % N, bbox_inches='tight', format='eps', dpi=300)
    fig1.savefig('normalized_load_vs_packet_loss_rate_macro_diversity_N=%s.eps' % N, bbox_inches='tight', format='eps', dpi=300)
    fig0.savefig('normalized_load_vs_packet_loss_rate_best_bs_N=%s.eps' % N, bbox_inches='tight', format='eps', dpi=300)

    plt.show()






