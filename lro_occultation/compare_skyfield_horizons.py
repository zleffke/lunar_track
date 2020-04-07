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

import utilities.utilities as utils

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
    cfg = utils.import_configs(args)
    #Import Horizons Vector and Pointing Data
    df_vec_lro = utils.HorizonsCSVToDataFrame(cfg['data']['path'], cfg['data']['lro']['vec_file'])
    df_eph_lro = utils.HorizonsCSVToDataFrame(cfg['data']['path'], cfg['data']['lro']['eph_file'])
    df_vec_moon = utils.HorizonsCSVToDataFrame(cfg['data']['path'], cfg['data']['moon']['vec_file'])
    df_eph_moon = utils.HorizonsCSVToDataFrame(cfg['data']['path'], cfg['data']['moon']['eph_file'])
    #df = df_vec.merge(df_eph, how='inner', on='datetime_jd [TDB]')
    dfs = [df_vec_lro, df_eph_lro, df_vec_moon, df_eph_moon]
    #df = pd.concat(dfs, join='outer')
    df = df_vec_lro.merge(df_eph_lro, how='inner', on='datetime_utc')
    df = df.merge(df_vec_moon, how='inner', on='datetime_utc')
    df = df.merge(df_eph_moon, how='inner', on='datetime_utc')
    print(df.info())
    #print(df)
    #sys.exit()

    #set up Skyfield data Loader path
    load = utils.setup_skyfield_loader(cfg)
    ts = load.timescale()
    if args.offset:
        a = 'datetime_utc'
        b = 'MOON 1-Way Prop Delay [sec]'
        #print(df[[a,b]])
        df['offset_utc'] = df.apply(lambda x: x[a]+pd.Timedelta(seconds=x[b]/args.offset), axis=1)
        t = ts.utc(np.array(df['offset_utc']))
    else:
        t = ts.utc(np.array(df['datetime_utc']))
    df['utc_iso'] = t.utc_iso(places=6)
    #print(df.keys())
    #load almanac
    e = load('de421.bsp')
    earth = e['earth']
    moon = e['moon']
    sun = e['sun']
    #setup ground station
    gs = sf.Topos(latitude_degrees=cfg['gs']['latitude'],
                  longitude_degrees=cfg['gs']['longitude'],
                  elevation_m=cfg['gs']['altitude'])

    #print(e)

    #---GS ----
    r_site = np.array(gs.at(t).position.km).transpose()
    mag_site = np.linalg.norm(r_site, axis=1)
    print("\n#---EARTH TO GS, r_site, Skyfield Computed----")
    print("              GS X [km]:", r_site.transpose()[0])
    print("              GS y [km]:", r_site.transpose()[1])
    print("              GS Z [km]:", r_site.transpose()[2])
    print("          GS Range [km]:", mag_site)

    r_e_bary = np.array(earth.at(t).position.km).transpose() #barycentric
    #r_s_bary = np.array(sun.at(t).position.km).transpose() #barycentric
    r_m_bary = np.array(moon.at(t).position.km).transpose() #barycentric
    r_moon = r_m_bary - r_e_bary
    mag_moon = np.linalg.norm(r_moon, axis=1)
    print("\n#---MOON Geocentric, r_moon, Skyfield Computed, from BARYCENTERIC----")
    print("            MOON X [km]:", r_moon.transpose()[0])
    print("            MOON y [km]:", r_moon.transpose()[1])
    print("            MOON Z [km]:", r_moon.transpose()[2])
    print("        MOON Range [km]:", mag_moon)

    r_moon2 = (earth).at(t).observe(moon).position.km.transpose()
    mag_moon = np.linalg.norm(r_moon2, axis=1)
    print("\n#---MOON Geocentric, r_moon2, Skyfield Computed, from Observe Method----")
    print("            MOON X [km]:", r_moon2.transpose()[0])
    print("            MOON y [km]:", r_moon2.transpose()[1])
    print("            MOON Z [km]:", r_moon2.transpose()[2])
    print("        MOON Range [km]:", mag_moon)

    #---GS TO MOON ----
    print("\n#---GS TO MOON, Skyfield Computed Compared to HORIZONS----")
    moon_rv = (earth+gs).at(t).observe(moon) #returns astrometric position
    [el_moon,az_moon,rho_moon] = moon_rv.apparent().altaz()
    df['Lunar Elevation [deg]'] = el_moon.degrees
    df['Lunar Azimuth [deg]'] = az_moon.degrees
    df['Lunar Range [km]'] = rho_moon.km

    print("  Lunar Azimuth [deg]*:", np.array(df['Lunar Azimuth [deg]']).transpose())
    print("Lunar Elevation [deg]*:", np.array(df['Lunar Elevation [deg]']).transpose())
    print("     Lunar Range [km]*:", np.array(df['Lunar Range [km]']).transpose())

    print("   Moon Azimuth [deg]`:", np.array(df['MOON Azimuth [deg]']).transpose())
    print(" Moon Elevation [deg]`:", np.array(df['MOON Elevation [deg]']).transpose())
    print("      Moon Range [km]`:", np.array(df['MOON Range [km]']).transpose())

    df['Azimuth Delta [deg]'] = df['Lunar Azimuth [deg]'] - df['MOON Azimuth [deg]']
    df['Elevation Delta [deg]'] = df['Lunar Elevation [deg]'] - df['MOON Elevation [deg]']
    df['Range Delta [km]'] = df['Lunar Range [km]'] - df['MOON Range [km]']

    print("* - Skyfield Computed")
    print("` - HORIZONS Computed\n")

    print("  Azimuth Delta [deg]:", np.array(df['Azimuth Delta [deg]']).transpose())
    print("Elevation Delta [deg]:", np.array(df['Elevation Delta [deg]']).transpose())
    print("     Range Delta [km]:", np.array(df['Range Delta [km]']).transpose())

    print('\n-------Residual Error Extremes-------')
    max_az_del = max(df['Azimuth Delta [deg]'])
    min_az_del = min(df['Azimuth Delta [deg]'])
    max_el_del = max(df['Elevation Delta [deg]'])
    min_el_del = min(df['Elevation Delta [deg]'])
    max_r_del = max(df['Range Delta [km]'])
    min_r_del = min(df['Range Delta [km]'])

    print('  Max Azimuth Delta [deg]:', max_az_del)
    print('  Min Azimuth Delta [deg]:', min_az_del)
    print('  Mag Azimuth Delta [deg]:', max_az_del - min_az_del)
    print()
    print('Max Elevation Delta [deg]:', max_el_del)
    print('Min Elevation Delta [deg]:', min_el_del)
    print('Mag Elevation Delta [deg]:', max_az_del - min_az_del)
    print()
    print('     Max Range Delta [km]:', max_r_del)
    print('     Min Range Delta [km]:', min_r_del)
    print('     Mag Range Delta [km]:', max_r_del - min_r_del)
    print()
    print("I think the remaining residuals may have to do with use of apparent")
    print("altaz in the skyfield calculation.....not sure though because HORIZONS")
    print("should also be delivering the apparent az/el/range")

    sys.exit()
