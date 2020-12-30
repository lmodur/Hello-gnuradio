import numpy as np
import sdr4
from matplotlib import pyplot as plt

def createThermalNoise(nsamps,bw,NF):
    noise = np.random.normal(0,1,nsamps)+1j*np.random.normal(0,1,nsamps)
    noise = np.sqrt(0.5)*noise
    noise_psd = -174 #dBm/Hz
    noise_pwr = noise_psd + 10*np.log10(bw) + NF
    noise_pwr_scal = 10**((noise_pwr-30)/10.0)
    noise = np.sqrt(noise_pwr_scal)*noise
    return noise
    
    
def measureSigPwr(sig):
    pwrdBm = 10*np.log10(np.average(np.abs(sig)**2)*10**3)
    print(pwrdBm)
    return pwrdBm

def scaleSigPwr(sig,tgt_pwr):
    pwrdBm = measureSigPwr(sig)
    scale = np.sqrt(10**((tgt_pwr - pwrdBm)/10.0))
    sig = sig*scale
    measureSigPwr(sig)
    return sig
    
def signalDetect(sig,noise_pwr):
    SD = 4 + noise_pwr
    detect = False
    if measureSigPwr(sig) >= SD:
        detect = True
    return detect


def packetDetect(r,L):
    base_len = 13
    D = base_len*L 
    corrn = np.zeros(np.size(r)) + 1j*np.zeros(np.size(r))
    corrd = np.zeros(np.size(r)) + 1j*np.zeros(np.size(r))
    res = np.array([])    
    rd = np.roll(r,D)
    rd[0:D] = np.zeros(D) + 1j*np.zeros(D)
    corrn = np.multiply(r,np.conj(rd))
    corrd = np.multiply(r,np.conj(r))
    
    sum_n = 0 + 1j*0
    sum_d = 0 + 1j*0
    for k in range(np.size(r)-D):
        for j in range(D):
            sum_n = sum_n + corrn[k+j]
            sum_d = sum_d + corrd[k+j]
        if sum_d != 0:
            res = np.append(res,np.abs(sum_n)/np.abs(sum_d))
        sum_n = 0 + 1j*0
        sum_d = 0 + 1j*0
        
    thresh = 0.7
    fstart = 0
    for z in range(np.size(res)):
        if np.average(res[z:z+(D*3)]) > thresh:
            fstart = z
            break
    return res,fstart

def frameAlignment(rx,ref):
    corr = np.zeros(np.size(rx)) + 1j*np.zeros(np.size(rx))
    ref = np.insert(ref,np.size(ref),np.zeros(np.size(rx)-np.size(ref)))
    #plt.plot(np.abs(rx))
    #plt.plot(np.abs(ref))
    for k in range(np.size(rx)):
        corr[k] = np.sum(np.multiply(rx,np.conj(ref)))
        ref = np.roll(ref,1)
    return corr
     
    
    
    
    
    
