import string
import numpy as np
from matplotlib import pyplot as plt 
import sdr1

def bin2gray(num):
    return num^(num>>1)

def gray2bin(num):
    mask = num
    while mask !=0:
        mask = mask >> 1
        num = num ^ mask
    return num
    
def genQuantLevels(M):
    if M == 4:
        rI = int(np.sqrt(M))
        rQ = rI
    elif M == 16:
        rI = int(np.sqrt(M))
        rQ = rI
    elif M == 64:
        rI = int(np.sqrt(M))
        rQ = rI
    elif M == 128:
        rI = 16
        rQ = 8
    elif M == 256:
        rI = int(np.sqrt(M))
        rQ = rI
        
    m_bits = int(np.log10(M)/np.log10(2))
    dI = np.arange(1,rI+1,1)
    dQ = np.arange(1,rQ+1,1)
    ref_syms_I = 2*dI - 1 - rI
    ref_syms_Q = 2*dQ - 1 - rQ
    ref_syms_f = np.flip(ref_syms_Q)
    #Constellation normalization
    normf = 0
    for j in range(np.size(ref_syms_I)):
        for k in range(np.size(ref_syms_Q)):
            normf = normf + ref_syms_I[j]**2 + ref_syms_Q[k]**2 
    normf = np.sqrt(normf/M)
    ref_syms_I = ref_syms_I/normf
    ref_syms_Q = ref_syms_Q/normf
    ref_syms_f = ref_syms_f/normf
    return ref_syms_I,ref_syms_Q,ref_syms_f                           
                                            
    
def mapBits2QAMSyms(M,in_bits):
    ref_syms_I,ref_syms_Q,ref_syms_f = genQuantLevels(M)
    m_bits = int(np.log10(M)/np.log10(2))
    temp = np.insert(in_bits,0,np.zeros(8-np.size(in_bits)))
    #print(temp)
    temp1 = gray2bin(sdr1.packBits2Bytes(temp))
    #print(temp1)
    rI = np.size(ref_syms_I)
    rQ = np.size(ref_syms_Q)
    isel = int(temp1/rQ)
    qsel = temp1 - isel*rQ
    I_val = ref_syms_I[isel]
    if (isel - 2*(isel//2) == 0):
        Q_val = ref_syms_Q[qsel]
    else:
        Q_val = ref_syms_f[qsel]
    return I_val,Q_val
    
    
def quantizeQAMsym(M,Ival,Qval):
    ref_syms_I,ref_syms_Q,ref_syms_f = genQuantLevels(M)           
    cx_rx_sym = Ival + 1j*Qval
    cx_ref_sym = np.array([])
    for i in range(np.size(ref_syms_I)):
        for j in range(np.size(ref_syms_Q)):
            cx_ref_sym = np.append(cx_ref_sym,ref_syms_I[i]+1j*ref_syms_Q[j])
    #print(cx_ref_sym)
    minidx = np.argmin(np.abs(cx_ref_sym - cx_rx_sym))
    minI = int(minidx/np.size(ref_syms_Q))
    minQ = int(minidx  - minI*np.size(ref_syms_Q))
    #print(minidx,minI,minQ)
    return ref_syms_I[minI],ref_syms_Q[minQ]
    
def unmapQAMSym2Bits(M,Ival,Qval):
    ref_syms_I,ref_syms_Q,ref_syms_f = genQuantLevels(M)
    whichI = np.where(ref_syms_I == Ival)
    if(whichI[0] - 2*(whichI[0]//2)) == 0: 
        whichQ = np.where(ref_syms_Q == Qval)
    else:
        whichQ = np.where(ref_syms_f == Qval)
    #print(whichI,whichQ)
    binval = whichI[0]*np.size(ref_syms_Q) + whichQ[0]
    return bin2gray(binval)
    
def mapBits2BPSK(in_bit):
    I_val = np.zeros(1)
    if in_bit == 0:
        I_val = -1.0/np.sqrt(2)
    else:
        I_val = 1.0/np.sqrt(2)
    return I_val
    
def quantizeBPSK(Ival):
    ref1 = 1.0/np.sqrt(2)
    ref2 = -1.0/np.sqrt(2)
    if (Ival-ref1)**2 < (Ival-ref2)**2:
        return ref1
    else:
        return ref2

def unmapBPSK2Bits(Ival):
    ref1 = 1.0/np.sqrt(2)
    ref2 = -1.0/np.sqrt(2)
    if Ival == ref1:
        return 1
    else:
        return 0
    












