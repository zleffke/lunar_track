#! /usr/bin/python3

import sys
import os
import math
import string
import argparse
import json
import subprocess
import datetime
import time
import skyfield
import numpy as np
import pandas as pd
from skyfield import api as sf
from skyfield import almanac
from skyfield.api import utc
from astroquery.jplhorizons import Horizons
from astropy import units as u

import utilities.general as gu
import utilities.plotting as pu

if __name__ == "__main__":
    """ Main entry point to start the service. """
    #--------START Command Line argument parser------------------------------------------------------
    parser = argparse.ArgumentParser(description="Compute Lunar Occultation Information for LRO using JPL Horizons Data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cwd = os.getcwd()
    cfg_fp_default = '/'.join([cwd, 'config'])
    cfg = parser.add_argument_group('Configuration File')
    cfg.add_argument('--cfg_path',
                       dest='cfg_path',
                       type=str,
                       default='/'.join([os.getcwd(), 'config']),
                       help="Camera Configuration File Path",
                       action="store")
    cfg.add_argument('--cfg_file',
                       dest='cfg_file',
                       type=str,
                       default="lro_occultation.json",
                       help="Configuration File",
                       action="store")
    cfg.add_argument('--offset',
                       dest='offset',
                       type=float,
                       default=0,
                       help="Apply Light Time Offset to UTC",
                       action="store")
    args = parser.parse_args()
    #--------END Command Line argument parser------------------------------------------------------
    subprocess.run(["reset"])
    np.set_printoptions(precision=10, linewidth=125)
    #setup configuration file info
    cfg = gu.import_configs(args)
    #Import Horizons Vector and Pointing Data
    df_vec_lro = gu.HorizonsCSVToDataFrame(cfg['data']['path'], cfg['data']['lro']['vec_file'])
    df_eph_lro = gu.HorizonsCSVToDataFrame(cfg['data']['path'], cfg['data']['lro']['eph_file'])
    df_vec_moon = gu.HorizonsCSVToDataFrame(cfg['data']['path'], cfg['data']['moon']['vec_file'])
    df_eph_moon = gu.HorizonsCSVToDataFrame(cfg['data']['path'], cfg['data']['moon']['eph_file'])
    #df = df_vec.merge(df_eph, how='inner', on='datetime_jd [TDB]')
    dfs = [df_vec_lro, df_eph_lro, df_vec_moon, df_eph_moon]
    #df = pd.concat(dfs, join='outer')
    df = df_vec_lro.merge(df_eph_lro, how='inner', on='datetime_utc')
    df = df.merge(df_vec_moon, how='inner', on='datetime_utc')
    df = df.merge(df_eph_moon, how='inner', on='datetime_utc')
    print(df.info())
    #print(df)
    #sys.exit()

    #---Earth to GS ----
    #set up Skyfield data Loader path
    load = gu.setup_skyfield_loader(cfg)
    ts = load.timescale()
    t = ts.utc(np.array(df['datetime_utc']))
    if args.offset:
        a = 'datetime_utc'
        b = 'MOON 1-Way Prop Delay [sec]'
        #print(df[[a,b]])
        df['offset_utc'] = df.apply(lambda x: x[a]+pd.Timedelta(seconds=x[b]/args.offset), axis=1)
        t = ts.utc(np.array(df['offset_utc']))
    else:
        t = ts.utc(np.array(df['datetime_utc']))
    df['utc_iso'] = t.utc_iso(places=6)

    #load almanac
    e = load('de421.bsp')
    earth = e['earth']
    moon = e['moon']
    #setup ground station
    gs = sf.Topos(latitude_degrees=cfg['gs']['latitude'],
                  longitude_degrees=cfg['gs']['longitude'],
                  elevation_m=cfg['gs']['altitude'])
    r_site = gs.at(t)
    r_site = np.array(r_site.position.km).transpose()
    mag_site = np.linalg.norm(r_site, axis=1)

    print("\n#---EARTH TO GS, r_site, Skyfield Computed----")
    print("              GS X [km]:", r_site.transpose()[0])
    print("              GS y [km]:", r_site.transpose()[1])
    print("              GS Z [km]:", r_site.transpose()[2])
    print("          GS Range [km]:", mag_site)

    #---EARTH TO MOON ----
    pos_keys = ['MOON x [km]','MOON y [km]','MOON z [km]']
    r_moon = np.array(df[pos_keys]) #geocentric position
    mag_moon = np.linalg.norm(np.array(r_moon), axis=1) #earth center to moon center
    print("\n#---EARTH TO MOON, r_moon, GEOCENTRIC HORIZONS Computed----")
    print("HORIZONS Moon X [km]:", np.array(df['MOON x [km]']).transpose())
    print("HORIZONS Moon Y [km]:", np.array(df['MOON y [km]']).transpose())
    print("HORIZONS Moon Z [km]:", np.array(df['MOON z [km]']).transpose())
    print("HORIZONS MOON R [km]:", mag_moon)

    #---EARTH TO LRO ----
    pos_keys = ['LRO x [km]','LRO y [km]','LRO z [km]']
    r_lro = np.array(df[pos_keys]) #geocentric position
    mag_lro = np.linalg.norm(np.array(r_lro), axis=1)#earth center to lro
    print("\n#---EARTH TO LRO, r_lro, HORIZONS Computed----")
    print("HORIZONS LRO X [km]:", np.array(df['LRO x [km]']).transpose())
    print("HORIZONS LRO y [km]:", np.array(df['LRO y [km]']).transpose())
    print("HORIZONS LRO Z [km]:", np.array(df['LRO z [km]']).transpose())
    print("HORIZONS LRO R [km]:", mag_lro)

    #---MOON TO LRO ----
    print("\n#---MOON TO LRO, rho_ms, HORIZONS Computed----")
    rho_ms = r_lro - r_moon #selenocentric
    mag_ms = np.linalg.norm(np.array(rho_ms), axis=1)
    lro_alt = mag_ms - gu.R_moon #altitude of LRO
    print("      mag_ms |r_ms| [km]:", mag_ms)
    print("       LRO Altitude [km]:", lro_alt)
    print("   Maximum Altitude [km]:", max(lro_alt))
    print("   Minimum Altitude [km]:", min(lro_alt))

    #---GS to Moon ----
    print("\n----RHO, Referenced to Ground Station ")
    print("#---GS TO MOON, rho_moon, HORIZONS(`) & Skyfield(*) Computed----")
    rho_moon = r_moon - r_site
    rho_moon_mag = np.linalg.norm(np.array(rho_moon), axis=1)

    print("  RHO HORIZONS Moon Azimuth [deg]`:", np.array(df['MOON Azimuth [deg]']).transpose())
    print("RHO HORIZONS Moon Elevation [deg]`:", np.array(df['MOON Elevation [deg]']).transpose())
    print("     RHO HORIZONS Moon Range [km]`:", np.array(df['MOON Range [km]']).transpose())
    print("         RHO HORIZONS Moon X [km]':", np.array(rho_moon).transpose()[0])
    print("         RHO HORIZONS Moon Y [km]':", np.array(rho_moon).transpose()[1])
    print("         RHO HORIZONS Moon Z [km]':", np.array(rho_moon).transpose()[2])
    print("     RHO HORIZONS MOON Range [km]':", rho_moon_mag)

    #--- Do SIGHT ALGORITHM -----
    print("\n#---SIGHT ALGORITHM----")
    print('Moon Radius [km]:', gu.R_moon)
    theta_1 = np.arccos(np.divide(gu.R_moon, rho_moon_mag)) * gu.rad2deg
    theta_2 = np.arccos(np.divide(gu.R_moon, mag_ms)) * gu.rad2deg
    #theta_1 = np.arccos(np.divide(gu.R_moon, max(rho_moon2_mag))) *gu.rad2deg
    #theta_2 = np.arccos(np.divide(gu.R_moon, max(mag_ms))) * gu.rad2deg
    theta_12 = theta_1 + theta_2
    print('   theta_1 [deg]:', theta_1)
    print('   theta_2 [deg]:', theta_2)
    print('  theta_12 [deg]:', theta_12)

    theta = np.arccos(np.divide(np.diag(np.dot(-1*rho_moon, rho_ms.transpose())),rho_moon_mag*mag_ms))* gu.rad2deg
    print('theta [deg]:', theta)
    LOS = (theta <= theta_12)
    df['LOS'] = LOS


    df['dAz [deg]'] = df['LRO Azimuth [deg]'] - df['MOON Azimuth [deg]']
    df['dEl [deg]'] = df['LRO Elevation [deg]'] - df['MOON Elevation [deg]']

    R_moon_deg = np.average(
                    np.arctan(
                    np.divide(gu.R_moon, df['MOON Range [km]'])))*gu.rad2deg
    print("Angular Radius of Moon, avg [deg]:", R_moon_deg)
    df1 = df[['datetime_utc','dAz [deg]', 'dEl [deg]', 'LOS']]
    df1.name = "Moon View"

    for i in np.array(df['MOON Elevation [deg]']):
        print(i)


    pu.plot_moon_view2(0, df1, o_path=None, save=0, R_moon_deg=R_moon_deg)

    sys.exit()
    occult = list(np.where(LOS[:-1] != LOS[1:])[0]) #True - setting, False - Rising

    for i in occult:
        print(i-1, df['datetime_utc'][i-1], LOS[i-1])
        print(i, df['datetime_utc'][i], LOS[i])
        print(i+1, df['datetime_utc'][i+1], LOS[i+1], '\n')
    print(occult)
    events = []
    if LOS[occult[0]]: #first event is setting, init with t[0]
        # rise_time = df['datetime_utc'][0]
        # set_time  = df['datetime_utc'][occult[0]]
        # duration  = (set_time - rise_time).total_seconds()
        # events.append({
        #     "rise_time":rise_time,
        #      "set_time":set_time,
        #      "duration":duration
        # })
        occult.pop(0)
        print(occult)

    #print(list(zip(occult[::2], occult[1::2])))
    for tup in list(zip(occult[::2], occult[1::2])):
        rise_time = df['datetime_utc'][tup[0]]
        set_time  = df['datetime_utc'][tup[1]]
        duration  = (set_time - rise_time).total_seconds() / 60.0
        events.append({
            "rise_time":rise_time,
             "set_time":set_time,
             "duration":duration
        })

    for idx, event in enumerate(events):
        print(idx, event)
        #print(idx, 'rise time', event['rise_time'])
        #print(idx, ' set_time', event['set_time'])





    # print('\ni, Time UTC, Alt, theta, theta_1, theta_2, theta_12, LOS')
    # for i, th in enumerate(theta):
    #     print (i, t[i].utc_iso(), lro_alt[i], th, theta_1[i], theta_2[i], theta_12[i], LOS[i], occ)





    sys.exit()
