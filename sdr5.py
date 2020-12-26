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


def packetDetect(r,W,L):
    D = 16*L
    Lp = W*L 
    corrn = np.zeros(np.size(r)) + 1j*np.zeros(np.size(r))
    corrd = np.zeros(np.size(r)) + 1j*np.zeros(np.size(r))
    res = np.zeros(np.size(r))
    resd = np.zeros(np.size(r))
    for n in range(np.size(r)):
        for k in range(Lp):
            if n+k+D < np.size(r):
                corrn[n] = corrn[n] + r[n+k]*np.conj(r[n+k+D])
                corrd[n] = corrd[n] + r[n+k+D]*np.conj(r[n+k+D])
        if corrd[n] != 0:
            res[n] = np.abs(corrn[n])/np.abs(corrd[n])
            resd[n] = res[n] - res[n-1]   
    return res

def frameAlignment(rx,ref):
    corr = np.zeros(np.size(rx)) + 1j*np.zeros(np.size(rx))
    ref = np.insert(ref,np.size(ref),np.zeros(np.size(rx)-np.size(ref)))
    #plt.plot(np.abs(rx))
    #plt.plot(np.abs(ref))
    for k in range(np.size(rx)):
        corr[k] = np.sum(np.multiply(rx,np.conj(ref)))
        ref = np.roll(ref,1)
    return corr
     
    
    
    
    
    
