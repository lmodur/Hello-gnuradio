import string
import numpy as np

def str2Bytes(str):
    arr = bytes(str,'ascii')
    return arr
    
def Bytes2str(val):
    return chr(int(val))

def unpackByte2Bits(num):
    out = [1 if num & (1 << (7-n)) else 0 for n in range(8)]
    return np.array(out)

def packBits2Bytes(num):
    out = 0
    for i in range(np.size(num)):
        out = out + (2**(7-i))*num[i]
    return int(out)
    
def calcBER(x_bits,y_bits):
    ber = 0
    err_bits = x_bits - y_bits
    for val in err_bits:
        if val:
            ber = ber + 1
    return (100*ber/np.size(x_bits))

    
def scrambler(tx_bits):
    lsr = np.array([1,0,1,1,1,0,1])
    for i in range(np.size(tx_bits)):
        out = 1
        if lsr[3] == lsr[6]:
            out = 0
        lsr = np.roll(lsr,1)
        lsr[0] = out
        if tx_bits[i] == out:
            tx_bits[i] = 0
        else:
            tx_bits[i] = 1
    return tx_bits
    