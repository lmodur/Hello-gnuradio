import numpy as np
from scipy import signal as signal
from scipy import fftpack as fft
from matplotlib import pyplot as plt
import string
import random

def str2Bytes(str):
    arr = bytes(str,'ascii')
    return arr
    
def Bytes2str(val):
    return chr(int(val))

def unpackByte2Bits(num):
    out = [1 if num & (1 << (7-n)) else 0 for n in range(8)]
    return np.array(out)
    
def symbolmapQAM4(bitstream):
    val1 = 1/np.sqrt(2)
    valm1 = -1/np.sqrt(2)
    qam4LUT_I = {0:val1,1:valm1,3:val1,2:valm1}
    qam4LUT_Q = {0:val1,1:val1,3:valm1,2:valm1}
    j = 0
    k = 0
    Isymstream = np.zeros(int(np.size(bitstream)/2))
    Qsymstream = np.zeros(int(np.size(bitstream)/2))
    while j < np.size(bitstream):
        idx = 2*bitstream[j] + bitstream[j+1]
        j = j + 2
        Isymstream[k] = qam4LUT_I[int(idx)]
        Qsymstream[k] = qam4LUT_Q[int(idx)]
        k = k + 1
    return Isymstream,Qsymstream
    
def symbolunmapQAM4(Isymstream,Qsymstream):
    val1 = 1/np.sqrt(2)
    valm1 = -1/np.sqrt(2)
    k = 0
    bitstream = np.zeros(np.size(Isymstream)*2)
    for i in range(np.size(Isymstream)):
        if Isymstream[i] == val1 and Qsymstream[i] == val1:
            bitstream[k] = 0
            bitstream[k+1] = 0
        elif Isymstream[i] == valm1 and Qsymstream[i] == val1:
            bitstream[k] = 0
            bitstream[k+1] = 1
        elif Isymstream[i] == val1 and Qsymstream[i] == valm1:
            bitstream[k] = 1
            bitstream[k+1] = 1
        elif Isymstream[i] == valm1 and Qsymstream[i] ==valm1:
            bitstream[k] = 1
            bitstream[k+1] = 0
        k = k + 2
    return bitstream
                        
def symbolmapPAM(bitstream,Mary):
    pam2LUT = {0:-1,1:1}
    pam4LUT = {2:-3,0:-1,1:1,3:3}
    cnt = np.size(bitstream)
    if Mary == 1:
        LUT = pam2LUT
    elif Mary == 2:
        LUT = pam4LUT
        cnt = cnt/2
    symstream = np.zeros(int(cnt))
    if Mary == 1:
        for i in range(np.size(bitstream)):
            symstream[i] = LUT[bitstream[i]]
    if Mary == 2:
        j = 0
        k = 0
        while j < np.size(bitstream):
            idx = 2*bitstream[j] + bitstream[j+1]
            j = j + 2
            symstream[k] = LUT[int(idx)]
            k = k + 1
    return symstream
    
def symbolDetectQAM4(Isymstream,Qsymstream):
    val1 = 1/np.sqrt(2)
    valm1 = -1/np.sqrt(2)
    syms =np.array([val1,valm1])
    Idetstream = np.zeros(np.size(Isymstream))
    Qdetstream = np.zeros(np.size(Qsymstream))
    
    for i in range(np.size(Isymstream)):
        idxI = np.argmin(np.abs(syms - Isymstream[i]))
        idxQ = np.argmin(np.abs(syms - Qsymstream[i]))
        Idetstream[i] = syms[idxI]
        Qdetstream[i] = syms[idxQ]
    return Idetstream,Qdetstream
        
def symbolDetect(symstream,Mary):
    pam2syms = np.array([-1,1])
    pam4syms = np.array([-3,-1,1,3])
    if Mary == 1:
        syms = pam2syms
    else:
        syms = pam4syms
        
    detstream = np.zeros(np.size(symstream))
    for i in range(np.size(symstream)):
        idx = np.argmin(np.abs(syms - symstream[i]))
        detstream[i] = syms[idx]
    return detstream
    
def symbolunmapPAM(symstream,Mary):
    pam2LUT = {-1:0,1:1}
    pam4LUT = {-3:2,-1:0,1:1,3:3}
    cnt = np.size(symstream)
    if Mary == 1:
        LUT = pam2LUT
    elif Mary == 2:
        LUT = pam4LUT
        cnt = cnt*2
    bitstream = np.zeros(int(cnt))
    if Mary == 1:
        for i in range(np.size(symstream)):
            bitstream[i] = LUT[symstream[i]]
    if Mary == 2:
        j = 0
        k = 0
        while k < np.size(symstream): 
            idx = LUT[symstream[k]]
            k = k + 1
            bitstream[j] = idx//2
            bitstream[j+1] =  idx - 2*(idx//2)
            j = j + 2
    return bitstream
    
    
    
def pulseshaper(taps,beta=0,type='boxcar'):
    pulse_coeffs = np.ones(taps)
    if type == 'rrc':
        n = 0
        while n < taps:
            p = n - int(taps/2)
            numr = (np.sin(np.pi*p*(1-beta)/taps)) + (4*beta*p/taps)*(np.cos(np.pi*p*(1+beta)/taps))
            denr = (np.pi*p/taps)*(1-((4*beta*p/taps)**2))
            if numr == 0 and denr == 0:
                pulse_coeffs[n] = 1-beta+(4*beta/np.pi)
            else:
                pulse_coeffs[n] = numr/denr
            n = n + 1   
    return pulse_coeffs/(np.sqrt(taps))
    
def packBits2Bytes(num):
    out = 0
    for i in range(np.size(num)):
        out = out + (2**(7-i))*num[i]
    return int(out)
    
def addNoise(sig,SNRdB):
    sig_var = 10.0*np.log10(sig@sig.T/(np.size(sig)))
    #print(sig_var)
    noise_var = 10**((sig_var - SNRdB)/10.0)
    noise = np.sqrt(noise_var)*np.random.normal(0,1,np.size(sig))
    noise_pwr = 10*np.log10(noise@noise.T/np.size(noise))
    #print(noise_pwr)
    sig_n = sig + noise
    return sig_n
    
def PCMquantizer(sig):
    #3 bit quantizer
    sig = sig/np.max(sig)
    n_lvls = 2**3
    quant = 2/n_lvls 
    lvls = np.arange(-1+quant/2,1,quant)
    print(lvls)
    idx = np.argmin(np.abs(sig -lvls.reshape(-1,1)),axis=0)
    
    qsig = np.zeros(np.size(sig))
    for i in range(np.size(sig)):
        qsig[i] = lvls[idx[i]]
    return qsig
    
def calcBER(x_bits,y_bits):
    ber = 0
    err_bits = x_bits - y_bits
    for val in err_bits:
        if val:
            ber = ber + 1
    return (ber/np.size(x_bits))
                
    
symsize = 11
Ts = 1e3
Tb = Ts/symsize

str = "The quick brown fox jumps over a lazy dog"

pkt = str2Bytes(str)

x = np.array([])
for pkts in pkt:
    x=np.append(x,unpackByte2Bits(pkts))

# Map to symbols
M = 2
symstream = symbolmapPAM(x,M)    



#Upsample using pulse shaper
beta = 1
x_up = np.zeros(np.size(symstream)*symsize)
pulse_coeffs = pulseshaper(symsize,beta,'rrc')
#plt.plot(pulse_coeffs)

#convolve with pulse coeffs
for i in range(np.size(symstream)):
    x_up[i*symsize] = symstream[i] 
x_up = signal.convolve(x_up,pulse_coeffs) 

X_UP = fft.fft(x_up)
freqX = fft.fftfreq(np.size(x_up),1/Ts)

#plt.plot(freqX,20*np.log10(np.abs(X_UP)))


#plt.subplot(2,1,1)
t = np.arange(0,np.size(x_up)*Tb,Tb)
#plt.plot(t,x_up,label='without noise')

y = addNoise(x_up,15)
#plt.plot(t,y,label='with noise')

#plt.legend(loc='upper right')


#Decode the PAM signal --Matched filter and downsample
y_recvd = np.convolve(y,pulse_coeffs)
Y_recvd = fft.fft(y_recvd)
freqY =fft.fftfreq(np.size(y_recvd),1/Ts)
#plt.plot(freqY,20*np.log10(np.abs(Y_recvd)))
y_syms = symbolDetect(y_recvd[symsize-1::symsize],M)



y_bits = symbolunmapPAM(y_syms,M)


#print(y_bits)
packedBytes = np.zeros(int(np.size(y_bits)/8))
for i in range(int(np.size(y_bits)/8)):
    packedBytes[i] = packBits2Bytes(y_bits[i*8:i*8+8])

rx_str = ''
for val in packedBytes:
    #print(val)
    rx_str =rx_str + Bytes2str(val)
print("Tx text is {} \n and Rx txt is {}".format(str,rx_str))

## Test PCM quantizer
fs = 1e3
f0 = 100
t = np.arange(0,0.5,1.0/fs)
sigx = np.sin(2*np.pi*f0*t)
#sigx = sigx.astype(np.float32)
#sigy = PCMquantizer(sigx)

imp = np.array([1])
z = signal.convolve(imp,pulse_coeffs)
zz = signal.convolve(z,pulse_coeffs)
#plt.plot(zz)


    ## Test QAM Tx and Rx loopback
Ts = 1e3
Tb = 2*Ts
L = 11
nchars = 1e4
SNR = np.arange(-3,10,1)
ber = np.zeros(np.size(SNR))
cntr = 0
for k in SNR:
    str = ''.join(random.choice(string.ascii_lowercase) for _ in range(int(nchars)))
    pkts = str2Bytes(str)
    x_bits = np.array([])
    for pkt in pkts:
        x_bits = np.append(x_bits,unpackByte2Bits(pkt))
        
    #print(x)
    Isymstream,Qsymstream = symbolmapQAM4(x_bits)
    #print(Isymstream)
    #print(Qsymstream)
    #Upsample using pulse shaper
    beta = 1
    xI_up = np.zeros(np.size(Isymstream)*L)
    xQ_up = np.zeros(np.size(Qsymstream)*L)
    pulse_coeffs = pulseshaper(symsize,beta,'rrc')
    #plt.plot(pulse_coeffs)
    #convolve with pulse coeffs
    for i in range(np.size(Isymstream)):
        xI_up[i*L] = Isymstream[i]
        xQ_up[i*L] = Qsymstream[i]
    xI_up = signal.convolve(xI_up,pulse_coeffs)
    xQ_up = signal.convolve(xQ_up,pulse_coeffs)
    
    #yI = addNoise(xI_up,SNR)
    #yQ = addNoise(xQ_up,SNR)
    
    #plt.plot(xI_up)
    #plt.plot(yI)
    #plt.plot(xQ_up)
    #plt.plot(yQ)
    
    yI = signal.convolve(xI_up,pulse_coeffs)
    yQ = signal.convolve(xQ_up,pulse_coeffs)
    #plt.plot(yI_recvd)
    #plt.plot(yQ_recvd)
    yI_recvd = addNoise(yI,k)
    yQ_recvd = addNoise(yQ,k)
    
    rxIsymstream,rxQsymstream = symbolDetectQAM4(yI_recvd[L-1:np.size(xI_up):L],yQ_recvd[L-1:np.size(xQ_up):L])
    
    
    #plt.scatter(yI_recvd[L-1::L],yQ_recvd[L-1::L])
    #plt.show()
    
    y_bits = symbolunmapQAM4(rxIsymstream,rxQsymstream)
    
    ber[cntr] = calcBER(x_bits,y_bits)
    cntr = cntr + 1

plt.semilogy(SNR,ber)
plt.show()

packedBytes = np.zeros(int(np.size(y_bits)/8))
for i in range(int(np.size(y_bits)/8)):
    packedBytes[i] = packBits2Bytes(y_bits[i*8:i*8+8])

rx_str = ''
for val in packedBytes:
    #print(val)
    rx_str =rx_str + Bytes2str(val)
#print("Tx text is {} \n and Rx txt is {}".format(str,rx_str))




