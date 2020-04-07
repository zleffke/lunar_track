#!/usr/bin/env python
#############################################
#   Title: Satellite utilities              #
# Project: TLE Match                        #
#    Date: Jan 2018                         #
#  Author: Zach Leffke, KJ4QLP              #
#############################################
import sys
import os
import math
import string
import struct
import scipy
import datetime
import pandas as pd
import numpy as np
from skyfield import api as sf
from skyfield import almanac

deg2rad = math.pi / 180
rad2deg = 180 / math.pi
c       = float(299792458)    #[m/s], speed of light
au2km   = 149597870.700 #km/au

def Doppler_Shift(center_freq, range_rate):
    #Returns Doppler Shift at receiver given fixed emitter
    #center_freq [Hz]:  center frequency of emitter
    #range_rate [m/s]:  range rate between emitter and receiver, negative means approaching
    #         c [m/s]:  speed of light
    #c = float(299792458) #[m/s], speed of light
    f_obs = (1.0 - range_rate / c) * center_freq
    f_delta = -1.0 * range_rate / c * center_freq
    doppler = {}
    doppler['center'] = f_obs
    doppler['offset'] = f_delta
    #return f_obs, f_delta
    return doppler

def Doppler_Shift_Invert(f_obs, range_rate):
    #Returns frequency of emitted signal given
    #      f_obs [Hz]:  measured or observed frequency
    #center_freq [Hz]:  center frequency of emitter
    #range_rate [m/s]:  range rate between emitter and receiver
    center_freq = f_obs / (1 - range_rate/c)
    return center_freq

def Freq_2_Lambda(frequency):
    #Frequency passed in Hertz
    c = 299792458 #[m/s]
    lam = c / frequency
    return lam

def Path_Loss(link_range, lam, n = 2):
    #link range - length of path in meters
    #lam         - operating wavelength in meters
    #n            - path loss exponent, exp = 2 for free space
    #loss = n*10*log10(4*pi*link_range / lam)
    #loss = n*10*math.log10(4*math.pi*link_range / lam)
    loss = n*10*np.log10(4*math.pi*link_range / lam)
    return loss

def lin_inv_xpndr_map(up, up_min = 145.900e6, up_max = 146.000e6, dn_min = 435.800e6, dn_max = 435.900e6):
    #Linear Inverting Transponder Uplink to Downlink Mapping
    #Default Values for FO-29 in Hz
    dn = dn_max - (up - up_min) * (dn_max - dn_min) / (up_max - up_min)
    return dn

def lin_inv_xpndr_map_reverse(dn, up_min = 145.900e6, up_max = 146.000e6, dn_min = 435.800e6, dn_max = 435.900e6):
    #Linear Inverting Transponder Downlink to Uplink Mapping
    #Default Values for FO-29 in Hz
    up = (dn_max - dn) * (up_max - up_min) / (dn_max - dn_min) + up_min
    return up

def Lunar_Rise_Set(e, gs, t0=None, t1=None):
    if t0== None:
        t0 = ts.utc(datetime.datetime.now(datetime.timezone.utc)) #Now
    if t1 == None:
        t1 = ts.utc(datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(hours=24)) #24 hours from now
    #print(t0)
    #print(t1)
    f = almanac.risings_and_settings(e, e['Moon'], gs)
    t, y = almanac.find_discrete(t0, t1, f)
    #print (t, y)
    for ti, yi in zip(t, y):
        #print (ti,yi)
        print('Lunar Rise:' if yi else ' Lunar Set:', ti.utc_datetime())
        if yi: #Rise
            rise_time = ti.utc_datetime()
        else:
            set_time = ti.utc_datetime()

    return {'rise':rise_time,
             'set':set_time,
             'vis':set_time < rise_time}

def Lunar_Illumination(e,t):
    fi = almanac.fraction_illuminated(e, 'Moon', t)
    return fi


# class satellite(object):
#     def __init__(self, ephem_sat, sat_name, norad_id):
#         self.ephem_sat  = ephem_sat     #PyEphem Satellite object for use in computations
#         self.sat_name   = sat_name #Common Name of spacecraft
#         self.norad_id   = norad_id #NORAD ID of spacecraft
#
#     def gen_doppler(self, gs, timestamp, rx_freq):
#             #input:  pyephem GS object
#         #--convert timestamp to useable datetime format
#         timestamp=[dt.datetime.utcfromtimestamp(element*1e-9) for element in timestamp]
#         #print timestamp
#         measured_freq = []
#         doppler_offset = []
#
#         for ts in timestamp:
#             gs.date = ephem.Date(ts) #convert datetime to ephem date
#             self.ephem_sat.compute(gs)
#             range_rate = self.ephem_sat.range_velocity
#             doppler = Doppler_Shift(rx_freq, range_rate)
#             #print ts, range_rate, doppler['center'], doppler['offset']
#             measured_freq.append(doppler['center'])
#             doppler_offset.append(doppler['offset'])
#
#         df = pd.DataFrame({ 'timestamp':timestamp,
#                             'doppler_offset':doppler_offset,
#                             'measured_freq':measured_freq})
#         df['doppler_offset'] = df['doppler_offset'].astype(float)
#         df['measured_freq'] = df['measured_freq'].astype(float)
#
#         df.name = self.sat_name +'('+ self.norad_id + ')'
#         return df
