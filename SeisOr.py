from obspy.core.inventory.inventory import read_inventory
from obspy import read
#import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
#from obspy.geodetics.base import gps2dist_azimuth
from scipy.stats import circmean,circstd
#import random

datadir= "C:/Users/josia/Dropbox/AmbientSeismicNoiseProject/MinDayReq/2/"

def readandfold(dir, station, orienstation):
    filename = dir+orienstation+'_'+station
    try:
        st = read(filename+'_RZ.SAC')
        st += read(filename+'_TZ.SAC')
        st += read(filename+'_ZZ.SAC')
        flag = True
    except: 
        filename =dir+station+'_'+orienstation
        st = read(filename+'_ZR.SAC')
        st += read(filename+'_ZT.SAC')
        st += read(filename+'_ZZ.SAC')
        flag = False
        
    # From one trace, pick off some constants we'll need:
    starttime = st[0].stats.starttime
    delta  = st[0].stats.delta
    npts2 = st[0].stats.npts
    npts = int((npts2-1)/2)

    # Fold the CCFs to improve SNR:
    st_folded= st.copy()
    for i_,tr in enumerate(st):
        causal =  st[i_].data[npts:-1]
        acausal = st[i_].data[npts:0:-1]
        st_folded[i_].data[0:npts] =  (causal+acausal)/2

    # Keep only the first half of the stream, as that is where
    # all the action is:
    st_folded.trim(starttime=starttime, endtime=starttime+npts*delta)
#    st_folded.filter("bandpass",freqmin=0.1,freqmax=1)
    # make sure the order is R,T,Z (or E, N, Z). Would be better to
    # use st.select() with a custom key header?
    st_folded.sort() 
#    st_folded.plot()
    ZZ = st_folded[2].data
    RZ = st_folded[0].data
    TZ = st_folded[1].data
    
    return ZZ,RZ,TZ,flag

# get station list from AVF inventory, where we included RBAZ: 
inv = read_inventory('../../Orientations/invAVF_with_RBAZ.xml')
net = inv[0]

for orienstation in net:
    if orienstation.code != 'WTAZ' and orienstation.code != 'RBAZ':
        # for each station pair, read the CCFs, rotate, and find
        # optimal rotation, based on Zha's method:
        angles=[]
        for station in net:
            if station.code != orienstation.code and station.code != 'RBAZ':
                # get the data from fn:
                [ZZ, RZ, TZ, flag] = readandfold(datadir, station.code,\
                                                   orienstation.code)

                # we are going to compare the RT with a 90 degree
                # shifted ZZ:
                ZZ90 = np.imag(hilbert(ZZ))
                Szz = np.correlate(ZZ90,ZZ90)

                # rotating CLOCKWISE, correlating with ZZ90:
                maxSrz= []
                thetas = np.linspace(0,2*np.pi,360)
                for i_,theta in enumerate(thetas):
                    RZ_rot =  np.cos(theta)*RZ - np.sin(theta)*TZ
                    TZ_rot = np.sin(theta)*RZ + np.cos(theta)*TZ
                    Srz = np.correlate(RZ_rot, ZZ90)
                    #Stz = np.correlate(TZ_rot, ZZ90)
                    Szz = np.correlate(ZZ90, ZZ90)
                    maxSrz.append(max(Srz)/max(Szz))

                # find the angle with the maximum correlation:
                maxmaxSrz= max(maxSrz)
                rotangle_index = maxSrz.index(maxmaxSrz)
                rotangle = thetas[rotangle_index]

                # MSNOISE pre-rotated, based on order of the station pair:
                if not flag:
                    rotangle = rotangle - np.pi
                
                # only use estimates with a max correlation larger than 0.3
                if  maxmaxSrz > 0.3:
                    angles.append(rotangle)
                    
        # calculate and print the (rounded) circular mean/std:
        circmeanangle = int(np.round(circmean(np.unwrap(angles))*180/np.pi))
        circmeanstd = int(np.round(circstd(np.unwrap(angles))*180/np.pi))
        print(orienstation.code, circmeanangle,circmeanstd, len(angles))
    
