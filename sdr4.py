import numpy as np
from matplotlib import pyplot as plt
import sdr1
import sdr2
import sdr3


def autocorr(sig1):
    tsig = sig1
    result = np.zeros(np.size(sig1))+1j*np.zeros(np.size(sig1))
    for k in range(np.size(sig1)):
        tsig = sig1
        tsig = np.roll(tsig,k)
        for l in range(k):
            tsig[l] = 0 + 1j*0
        result[k] = np.sum(np.multiply(sig1,np.conj(tsig)))
    return result         
    

def gen_msg_syms():
    tx_msg_base = [0,1,0,0,0,1,1,0,1,1,0,1,0,1,1,1,1,0,0,0,0,1,0,1,1,0,1,1,0,0,0,1]
    #tx_msg_base = [0,0,1,1]
    tx_msg = np.array([])
    # Make tx_msg by appending rept times
    rept = 160
    for i in range(rept):
        tx_msg = np.append(tx_msg,tx_msg_base)
    #print(np.size(tx_msg))
    tx_msg_I = np.array([])
    tx_msg_Q = np.array([])
    for i in range(int(np.size(tx_msg)/2)):
        tx_msg_It,tx_msg_Qt = sdr2.mapBits2QAMSyms(4,tx_msg[2*i:2*(i+1)])
        tx_msg_I = np.append(tx_msg_I,tx_msg_It)
        tx_msg_Q = np.append(tx_msg_Q,tx_msg_Qt)
    return tx_msg_I,tx_msg_Q


#Feed forward carrier phase recovery
def ff_carr_phs_recovery(ref_cplx,rx_cplx):
    ref_cplx = np.conj(ref_cplx)
    #print(ref_cplx,rx_cplx,np.multiply(ref_cplx,rx_cplx))
    est_phs = np.angle(np.average(np.multiply(ref_cplx,rx_cplx)))    
    #print(est_phs*180/np.pi)
    return est_phs
    
#Delay and multiply  coarse freqency offset estimation
def ff_coarse_cfo_estimation(rx_cplx,L):
    base_len = 16
    cntr = np.size(rx_cplx) - L*base_len
    phase_acc = 0 +1j*0
    for i in range(cntr):
        phase_acc = phase_acc + np.multiply(rx_cplx[i+(L*base_len)],np.conj(rx_cplx[i]))
    
    phase_acc_avg = phase_acc/cntr
    est_phs_cfo = np.angle(phase_acc_avg)
    return est_phs_cfo/(base_len*L)
    
def ff_fine_cfo_estimation(ref_cplx,rx_cplx):
    ref_cplx = np.conj(ref_cplx)
    conprdct = np.multiply(ref_cplx,rx_cplx)
    ac_result = autocorr(conprdct)
    res = 0.0    
    for i in range(np.size(ac_result)-1):
        res = res + (np.angle(ac_result[i+1])/(i+1))   
        #print(res)
    return (res/np.size(ac_result))
        
    
    




