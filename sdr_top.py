import math
import string
import random
import numpy as np
from matplotlib import pyplot as plt
import sdr1
import sdr2
import sdr3
import sdr4
import sdr5

#Convert a text to bytes and bits
tx_string = 'A quick brown fox jumps over the lazy dog'
#tx_string = 'G'
#tx_string = ''.join(random.choices(string.ascii_uppercase,k=1000))

tx_bytes = sdr1.str2Bytes(tx_string)
tx_bits = np.array([])

# Map bits to QAM symbols- dimension M
M = 4
m_bits = int(np.log10(M)/np.log10(2)) 
tx_syms_I = np.array([])
tx_syms_Q = np.array([])

tx_bits = np.array([])
for tx_byte in tx_bytes:
    tx_bits = np.append(tx_bits,sdr1.unpackByte2Bits(tx_byte))

#Scramble data bits
#print(tx_bits)
tx_bits = sdr1.scrambler(tx_bits)

for i in range(math.ceil(np.size(tx_bits)/m_bits)):
        temp_I,temp_Q = sdr2.mapBits2QAMSyms(M,tx_bits[i*m_bits:(i+1)*m_bits])    
        tx_syms_I = np.append(tx_syms_I,temp_I)
        tx_syms_Q = np.append(tx_syms_Q,temp_Q)

        
#print(tx_syms_I,tx_syms_Q)
#plt.scatter(tx_syms_I,tx_syms_Q)
#plt.show()

#Generate message (preamble)  symbols and prepend to data symbols
tx_msg_I,tx_msg_Q = sdr4.gen_msg_syms()
tx_syms_I = np.insert(tx_syms_I,0,tx_msg_I)
tx_syms_Q = np.insert(tx_syms_Q,0,tx_msg_Q)


#Generate rrc match filter, upsample each sym,pass through matched filter
L = 2
h = sdr3.pulseshaper(L,0.5,'rrc')
utx_syms_I = np.array([])
utx_syms_Q = np.array([])

for i in range(np.size(tx_syms_I)):
    utx_syms_I = np.append(utx_syms_I, sdr3.polyphase_interpolate(L,h,np.array([tx_syms_I[i]])))                
    utx_syms_Q = np.append(utx_syms_Q,sdr3.polyphase_interpolate(L,h,np.array([tx_syms_Q[i]])))

#plt.show()

#Add AWGN channel
SNRdB = 5
nsamps = 41
bw = 20e6
NF = 5
wn = sdr5.createThermalNoise(nsamps,bw,NF)
utx_wn = utx_syms_I + 1j*utx_syms_Q 
sig_pwr = sdr5.measureSigPwr(utx_wn)
tgt_pwr = sig_pwr - SNRdB
wn = sdr5.scaleSigPwr(wn,tgt_pwr)

#add noise
wna = sdr5.createThermalNoise(np.size(utx_wn),bw,NF)
wna = sdr5.scaleSigPwr(wna,tgt_pwr)
noise_pwr = sdr5.measureSigPwr(wna)
utx_wn = utx_wn + wna

#insert noise samples before signal
utx_wn = np.insert(utx_wn,0,wn)


#Average every W samples and detect if pwr > SD
Wd = 40 #window size for signal detect
sdcnt = 0
while sdcnt < np.size(utx_wn):
    detect = sdr5.signalDetect(utx_wn[sdcnt:sdcnt+Wd],noise_pwr)
    if detect:
        break
    sdcnt = sdcnt+Wd
if detect:
    print("Signal detected after samps!!!",sdcnt-Wd)
sdcnt = sdcnt - Wd
#Pass W symbols through thru packet detector
base_len = 16
rept = 4
W = base_len*rept
threshold = 0.6
res = sdr5.packetDetect(utx_wn[sdcnt:sdcnt+Wd+(W*L)],W,L)
#plt.plot(res)
res2 = np.average(res[0:int(W*L/2)])
if res2 > threshold:
    print("Packet detected")
    
#Find frame start
ref_msg_corr = np.array([])
for z in range(base_len):              ref_msg_corr = np.append(ref_msg_corr,sdr3.polyphase_interpolate(L,h,np.array([tx_msg_I[z]]))    
+1j*sdr3.polyphase_interpolate(L,h,np.array([tx_msg_Q[z]])))
corr = sdr5.frameAlignment(utx_wn[sdcnt:sdcnt+Wd+(base_len*L)],ref_msg_corr)
#plt.plot(np.abs(corr))
framestart = np.argmax(np.abs(corr))
print('Frame start detected at',framestart)

#strip off noise samples to retain only
#samples from packet start
rx_cplx = utx_wn[sdcnt+framestart:]


#Carrier phase diff between tx and rx is theta 
#Rotate the rx syms with phase theta
carr_phs_off_theta = 20.0
carr_phs_off = carr_phs_off_theta*np.pi/180
carr_phs_rot = np.exp(1j*carr_phs_off)
rx_cplx = carr_phs_rot * rx_cplx


# Add Carrier freq offset
F0 = 0.00005 # foffset*symbol_time of symbol rate
n = np.arange(0,np.size(rx_cplx),1)
cmplx_freq = np.exp(2*1j*np.pi*F0*n/L)
rx2_cplx = np.multiply(rx_cplx,cmplx_freq)

#Estimate and correct carrier freq offset
est_phase_cfo  = sdr4.ff_coarse_cfo_estimation(rx2_cplx[0:L*np.size(tx_msg_I)],L)
print(est_phase_cfo)
#est_phase_cfo = 0



urx_syms_I = np.real(rx2_cplx)
urx_syms_Q = np.imag(rx2_cplx)

# pass thru matched filter
D = 1
t_rx_syms_I = np.array([])
rx_syms_I = np.array([])
t_rx_syms_Q = np.array([])
rx_syms_Q = np.array([])

for i in range(int(np.size(urx_syms_I)/L)):
    t_rx_syms_I = np.append(t_rx_syms_I, sdr3.polyphase_decimate(D,h,urx_syms_I[i*L:(i+1)*L]))                
    t_rx_syms_Q = np.append(t_rx_syms_Q,sdr3.polyphase_decimate(D,h,urx_syms_Q[i*L:(i+1)*L]))

rx3_cplx = t_rx_syms_I + 1j*t_rx_syms_Q




#Down sample and freq correct
D = L
k = np.arange(0,np.size(rx3_cplx),1) 
t2_rx_cplx = np.multiply(rx3_cplx,np.exp(-1j*est_phase_cfo*k))

t2_rx_syms_cplx = t2_rx_cplx[D-1::D]

t2_rx_syms_I = np.real(t2_rx_syms_cplx)
t2_rx_syms_Q = np.imag(t2_rx_syms_cplx)


#Seperate message and data symbols
rx_msg_syms_I = t2_rx_syms_I[0:np.size(tx_msg_I)]
rx_msg_syms_Q = t2_rx_syms_Q[0:np.size(tx_msg_Q)]
rx_syms_I = t2_rx_syms_I[np.size(tx_msg_I):]
rx_syms_Q = t2_rx_syms_Q[np.size(tx_msg_Q):]

#Fine freq error estimation and correction
ref_msg_cplx = tx_msg_I + 1j*tx_msg_Q
rx_msg_cplx = rx_msg_syms_I + 1j*rx_msg_syms_Q

#print(ref_msg_cplx,rx_msg_cplx)

est_freq_rot = sdr4.ff_fine_cfo_estimation(ref_msg_cplx,rx_msg_cplx)
print(est_freq_rot/L)

rxff_cplx = np.array([])
rxff_cplx = np.append(rxff_cplx,rx_msg_cplx)
rxff_cplx = np.append(rxff_cplx,(rx_syms_I+1j*rx_syms_Q))

#est_freq_rot = 0
#Correct fine freq offset
m = np.arange(0,np.size(rxff_cplx),1)
rxff_cplx = np.multiply(rxff_cplx,np.exp(-1j*est_freq_rot*m))

rx_msg_cplx = rxff_cplx[0:np.size(tx_msg_I)]
rxf_cplx = rxff_cplx[np.size(tx_msg_I):]

# Use rx_msg_syms to extract carrier phs rotation
#print(ref_msg_cplx,rx_msg_cplx)
est_phs_rot = sdr4.ff_carr_phs_recovery(ref_msg_cplx,rx_msg_cplx)
print(est_phs_rot*180/np.pi)

#est_phs_rot = 0
#Phase error correction
rx_msg_syms_I = np.real(rx_msg_cplx*np.exp(-1j*est_phs_rot))
rx_msg_syms_Q = np.imag(rx_msg_cplx*np.exp(-1j*est_phs_rot))

#rxf_cplx = rx_syms_I + 1j*rx_syms_Q
rx_syms_I = np.real(rxf_cplx*np.exp(-1j*est_phs_rot))
rx_syms_Q = np.imag(rxf_cplx*np.exp(-1j*est_phs_rot))

#print(rx_msg_syms_I,rx_msg_syms_Q)
#scatter plot
#print(rx_syms_I,rx_syms_Q)
plt.scatter(rx_syms_I,rx_syms_Q)
plt.scatter(rx_msg_syms_I,rx_msg_syms_Q)
plt.show()


rec_syms_I = np.zeros(np.size(rx_syms_I))
rec_syms_Q = np.zeros(np.size(rx_syms_Q))
rec_bytes = np.zeros(np.size(rx_syms_Q))
rec_bits = np.array([])

for i in range(np.size(rx_syms_I)):
    rec_syms_I[i],rec_syms_Q[i] = sdr2.quantizeQAMsym(M,rx_syms_I[i],rx_syms_Q[i])    
    #print(rec_syms_I[i],rec_syms_Q[i])
    rec_bytes[i] = sdr2.unmapQAMSym2Bits(M,rec_syms_I[i],rec_syms_Q[i])
    #print(rec_bytes[i])
    temp_r = sdr1.unpackByte2Bits(int(rec_bytes[i]))
    temp_rf = np.flip(temp_r)
    rec_bits = np.append(rec_bits,temp_rf[:m_bits])
    #print(temp_r)

rec_bits_size = int(np.size(rec_bits)/8)
recvd_bits = rec_bits[0:8*rec_bits_size]

recvd_bits_f = np.array([])
temp_bits_f = np.array([])
j = 0
for i in range(np.size(recvd_bits)):
    temp_bits_f = np.append(temp_bits_f,recvd_bits[i])     
    j = j + 1
    if j == m_bits or i == np.size(recvd_bits) -1:
        recvd_bits_f = np.append(recvd_bits_f,np.flip(temp_bits_f))
        j = 0
        temp_bits_f = np.array([])
    

#Unscramble bits
recvd_bits_f = sdr1.scrambler(recvd_bits_f)
#print(recvd_bits_f)    

packedBytes = np.zeros(rec_bits_size)
rx_string= ''
for i in range(np.size(packedBytes)):
    packedBytes[i] = sdr1.packBits2Bytes(recvd_bits_f[i*8:(i+1)*8])
for pkt in packedBytes:
    rx_string = rx_string+(sdr1.Bytes2str(pkt))

print(rx_string)
#print(sdr1.calcBER(tx_bits,recvd_bits_f))






 


    





