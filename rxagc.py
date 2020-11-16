import numpy as np
from matplotlib import pyplot as plt


def dBm2mV(pwr):
    return 1000*np.sqrt(0.1*(10**(pwr/10.0)))
    


lna = {0:3,1:17,2:21}
mxr = {0:3,1:11,2:12,3:13,4:14,5:15,6:16,7:17,8:18,9:19,10:20,11:21,12:22,13:23,14:24,15:25,16:26}
tia = {0:-6,1:0}
lpf = np.arange(0,25,1)

print('Gain table')
print('GT idx\t','LNA gain\t','Mxr gain\t','TIA  gain\t','LPF gain\t')

appg_lna = np.zeros(77)
appg_mxr = np.zeros(77)
appg_tia = np.zeros(77)
appg_lpf = np.zeros(77)

#lna bypass region
appg_lna[5:19] = 0
appg_mxr[5:19] = 0
appg_mxr[5:19] = 0
appg_lpf[5:19] = np.arange(0,14,1)

#High rssi region
appg_lna[34:40] = 1
appg_mxr[34:40] = 5
appg_tia[34:40] = 0
appg_lpf[34:40] = np.arange(3,9,1)

appg_lna[30:34] = 1
appg_mxr[30:34] = np.arange(1,5,1)
appg_tia[30:34] = 0
appg_lpf[30:34] = 3

appg_lna[19:30] = 1
appg_mxr[19:30] = 0
appg_tia[19:30] = 0
appg_lpf[19:30] = np.arange(0,11,1)

# Mid rssi region
appg_lna[40:62] = 1
appg_mxr[40:62] = 5
appg_tia[40:62] = 1
appg_lpf[40:62] = np.arange(3,25,1)

# Low rssi region
appg_lna[62:73] = 1
appg_mxr[62:73] = np.arange(6,17,1)
appg_tia[62:73] = 1
appg_lpf[62:73] = 24

appg_lna[73:78] = 2
appg_mxr[73:78] = 16
appg_tia[73:78] = 1
appg_lpf[73:78] = np.arange(21,25,1)



for i in range(77):
    total_gain = lna[int(appg_lna[i])] +mxr[int(appg_mxr[i])] + tia[int(appg_tia[i])] +lpf[int(appg_lpf[i])]
    print("{}\t{}\t{}\t{}\t{}\t{}".format(i,lna[appg_lna[i]],mxr[int(appg_mxr[i])],tia[appg_tia[i]],lpf[int(appg_lpf[i])],total_gain))

pwr_max = -5    
LMT_thresh = 126
ADC_thresh = 178
clip1 = [2,16,1,24,76]
clip2 = [1,16,1,24,72]
clip3 = [1,5,1,24,61]
clip4 = [1,5,0,8,39]
clip5 = [1,0,0,10,29]
clip6 = [0,0,0,13,18]

def agcCalcs(pwr_in,clip):
    wb_det = 0
    nb_det = 0
    pwr1 = pwr_in + lna[int(clip[0])] + mxr[int(clip[1])] + tia[int(clip[2])]
    pwr2 = pwr1 + lpf[int(clip[3])]
    if dBm2mV(pwr1) > LMT_thresh:
        wb_det = 1
    if dBm2mV(pwr2) > ADC_thresh:
        nb_det = 1
    return wb_det,nb_det
            
def agcClips(pwr_in):
    clip_whr = clip1
    clip_bucket = 1
    wb_det,nb_det = agcCalcs(pwr_in,clip_whr)
    #print(wb_det,nb_det,'clip1')
    if wb_det == 1 or nb_det == 1:
        clip_whr = clip2
        clip_bucket = 2
        wb_det,nb_det = agcCalcs(pwr_in,clip_whr)
        #print(wb_det,nb_det,'clip2')
        if wb_det == 1 or nb_det == 1:
            clip_whr = clip3
            clip_bucket = 3
            wb_det,nb_det = agcCalcs(pwr_in,clip_whr)
            #print(wb_det,nb_det,'clip3')
            if wb_det == 1 or nb_det == 1:
                clip_whr = clip4
                clip_bucket = 4
                wb_det,nb_det = agcCalcs(pwr_in,clip_whr)
                #print(wb_det,nb_det,'clip4')
                if wb_det == 1 or nb_det == 1:
                    clip_whr = clip5
                    clip_bucket = 5
                    wb_det,nb_det = agcCalcs(pwr_in,clip_whr)
                    #print(wb_det,nb_det,'clip5')
                    if wb_det == 1 or nb_det == 1:
                        clip_whr = clip6
                        clip_bucket = 6
                        wb_det,nb_det=agcCalcs(pwr_in,clip_whr)
                        #print(wb_det,nb_det,'clip6')
    return clip_whr,clip_bucket
    
def agcSelGain(pwr_in,clip_whr,clip_bucket):
    sel_gain = np.zeros(5)
    pwr_est = pwr_in + lna[clip_whr[0]] + mxr[clip_whr[1]] + tia[clip_whr[2]] + lpf[clip_whr[3]]
    #print(clip_whr,clip_bucket,pwr_est)
    idx = clip_whr[4] + pwr_max - pwr_est 
    sel_gain[0] = appg_lna[idx]
    sel_gain[1] = appg_mxr[idx]
    sel_gain[2] = appg_tia[idx]
    sel_gain[3] = appg_lpf[idx]
    return sel_gain
    
    
lna_gain = np.zeros(np.size(appg_lna))
mxr_gain = np.zeros(np.size(appg_mxr))
tia_gain = np.zeros(np.size(appg_tia))
lpf_gain = np.zeros(np.size(appg_lpf))
pwr_vals = np.zeros(np.size(appg_lna))



pwr_in = -76
clip_whr,clip_bucket = agcClips(pwr_in)
sel_gain =agcSelGain(pwr_in,clip_whr,clip_bucket)
#print(lpf[int(sel_gain[3])] 

cnt = 0
while pwr_in <= pwr_max:
    clip_whr,clip_bucket = agcClips(pwr_in)
    sel_gain =agcSelGain(pwr_in,clip_whr,clip_bucket)
    pwr_vals[cnt] = pwr_in
    lna_gain[cnt] = lna[int(sel_gain[0])]
    mxr_gain[cnt] = mxr[int(sel_gain[1])]
    tia_gain[cnt] = tia[int(sel_gain[2])]
    lpf_gain[cnt] = lpf[int(sel_gain[3])]
    pwr_in = pwr_in + 1
    cnt = cnt + 1
    
plt.subplot(4,1,1)
plt.plot(pwr_vals,lna_gain)

plt.subplot(4,1,2)
plt.plot(pwr_vals,mxr_gain)

plt.subplot(4,1,3)
plt.plot(pwr_vals,tia_gain)

plt.subplot(4,1,4)
plt.plot(pwr_vals,lpf_gain)

plt.show()

