from obspy.core.inventory.inventory import read_inventory
from obspy import read
import os
import numpy as np
from scipy.signal import hilbert
from scipy.stats import circmean,circstd

datadir= os.getcwd()+"/ImpulseResponses/"

def readandfold(dir, station, orienstation):
    filename = dir+orienstation+'_'+station
    try:
        st = read(filename+'_RZ.MSEED')
        st += read(filename+'_TZ.MSEED')
        st += read(filename+'_ZZ.MSEED')
        flag = True
    except: 
        filename =dir+station+'_'+orienstation
        st = read(filename+'_ZR.MSEED')
        st += read(filename+'_ZT.MSEED')
        st += read(filename+'_ZZ.MSEED')
        flag = False
        
    # From one trace, pick off some constants we'll need:
    starttime = st[0].stats.starttime
    delta  = st[0].stats.delta
    npts2 = st[0].stats.npts
    npts = int((npts2-1)/2)

    # Fold the CCFs to improve SNR (make optional?):
    st_folded= st.copy()
    for i_,tr in enumerate(st):
        causal =  st[i_].data[npts:-1]
        acausal = st[i_].data[npts:0:-1]
        st_folded[i_].data[0:npts] =  (causal+acausal)/2

    #Keep only the first half of the stream, as that is where all the action is.
    st_folded.trim(starttime=starttime, endtime=starttime+npts*delta)

    # make sure the order is R,T,Z (or E, N, Z). Maybe better to use st.select() with a custom key header?
    st_folded.sort() 
    ZZ = st_folded[2].data
    RZ = st_folded[0].data
    TZ = st_folded[1].data
    
    return ZZ,RZ,TZ,flag

# get station list from inventory:
inv = read_inventory(os.getcwd()+"\inv.xml")
net = inv[0]

for orienstation in net:
    if orienstation.code != 'WTAZ': #In our network, WTAZ only has a vertical component.
        print(orienstation.code)
        # for each station pair, read the CCFs, rotate, and find
        # optimal rotation, based on Zha's method:
        angles=[]
        for station in net:
            if station.code != orienstation.code:  
                # get the data from fn:
                [ZZ, RZ, TZ, flag] = readandfold(datadir, station.code,\
                                                   orienstation.code)

                # we are going to compare the RT with a 90 degree
                # shifted ZZ:
                ZZ90 = np.imag(hilbert(ZZ))
                Szz = np.correlate(ZZ90,ZZ90)

                # rotating CLOCKWISE, correlating with ZZ90:
                maxSrz= []
                #coherences=[]
                thetas = np.linspace(0,2*np.pi,360)
                for i_,theta in enumerate(thetas):
                    RZ_rot =  np.cos(theta)*RZ - np.sin(theta)*TZ
                    TZ_rot = np.sin(theta)*RZ + np.cos(theta)*TZ
                    Srz = np.correlate(RZ_rot, ZZ90)
                    Szz = np.correlate(ZZ90, ZZ90)
                    maxSrz.append(max(Srz)/max(Szz))

                # find the angle with the maximum correlation:
                maxmaxSrz= max(maxSrz)
                rotangle_index = maxSrz.index(maxmaxSrz)
                rotangle = thetas[rotangle_index]

                # MSNOISE pre-rotated, based on order of the station pair:
                if not flag:
                    rotangle = rotangle - np.pi
                    
                # Data Culling: discard if max corr < 0.3
                if  maxmaxSrz > 0.3:
                    angles.append(rotangle)
                    print(int(rotangle*180/np.pi))
                    
        # calculate and print the (rounded) circular mean/std:
        circmeanangle = int(np.round(circmean(np.unwrap(angles))*180/np.pi))
        circmeanstd = int(np.round(circstd(np.unwrap(angles))*180/np.pi))
        print(orienstation.code, circmeanangle,circmeanstd, len(angles))

        
