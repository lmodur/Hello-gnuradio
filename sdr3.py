import numpy as np
from matplotlib import pyplot as plt
from scipy import signal as signal

h = np.array([1,2,3,4,5,6])
x = np.array([1,2,3,4,5,6,7,8])

def tapped_delay_line(coeffs,inp):
    d = np.zeros(np.size(coeffs))
    outp = np.zeros(np.size(inp))
    for i in range(np.size(inp)):
        d = np.roll(d,1)
        d[0] = inp[i]
        outp[i] = np.sum(np.multiply(coeffs,d))
    return outp

outp = tapped_delay_line(h,x)
#print(outp)


def polyphase_interpolate(L,coeffs,inp):
    n_r = np.size(coeffs)/L
    c = coeffs.reshape(int(n_r),L)
    c = c.T
    
    temp_out = np.zeros((L,np.size(inp)))
    for i in range(L):
        temp_out[i,:] = tapped_delay_line(c[i,:],inp)
    temp_out = temp_out.T
    outp = temp_out.flatten()
    return outp    


def polyphase_decimate(D,coeffs,inp):
    n_r =  np.size(coeffs)/D
    c = coeffs.reshape(int(n_r),D)
    c = c.T
    
    n_z = np.ceil(np.size(inp)/D) - int(np.size(inp)/D)
    inp = np.append(inp,np.zeros(int(n_z)))
    inp_t = inp.reshape(int(np.size(inp)/D),D)
    inp_t = inp_t.T

    outp = np.zeros(int(np.size(inp)/D))
    for i in range(D):
        t1 = np.roll(tapped_delay_line(c[i,:],inp_t[i,:]),i)
        #print(t1)
        if i > 0 and i < np.size(t1):
            t1[0:i] = np.zeros(i)
        outp = t1 + outp
    return outp
        
def pulseshaper(taps,beta=0,type='boxcar'):
    pulse_coeffs = np.ones(taps)
    if type == 'rrc':
        n = 0
        while n < taps-1:
            p = n - int(taps/2)
            numr = (np.sin(np.pi*p*(1-beta)/taps)) + (4*beta*p/taps)*(np.cos(np.pi*p*(1+beta)/taps))
            denr = (np.pi*p/taps)*(1-((4*beta*p/taps)**2))
            if numr == 0 and denr == 0:
                pulse_coeffs[n] = 1-beta+(4*beta/np.pi)
            else:
                pulse_coeffs[n] = numr/denr
            n = n + 1  
        #Normalize energy to 1
        energy = 0
        for i in range(np.size(pulse_coeffs)):
            energy = energy + pulse_coeffs[i]**2
    return pulse_coeffs/(np.sqrt(energy))


    

