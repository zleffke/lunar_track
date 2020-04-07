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

    #---Earth to GS ----
    #set up Skyfield data Loader path
    load = utils.setup_skyfield_loader(cfg)
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
    print("\n#---EARTH TO MOON, r_moon, HORIZONS Computed----")
    print("            Moon X [km]:", np.array(df['MOON x [km]']).transpose())
    print("            Moon y [km]:", np.array(df['MOON y [km]']).transpose())
    print("            Moon Z [km]:", np.array(df['MOON z [km]']).transpose())
    print(" mag_moon |r_moon| [km]:", mag_moon)

    #---EARTH TO LRO ----
    pos_keys = ['LRO x [km]','LRO y [km]','LRO z [km]']
    r_lro = np.array(df[pos_keys]) #geocentric position
    mag_lro = np.linalg.norm(np.array(r_lro), axis=1)#earth center to lro
    print("\n#---EARTH TO LRO, r_lro, HORIZONS Computed----")
    print("             LRO X [km]:", np.array(df['LRO x [km]']).transpose())
    print("             LRO y [km]:", np.array(df['LRO y [km]']).transpose())
    print("             LRO Z [km]:", np.array(df['LRO z [km]']).transpose())
    print("   mag_lro |r_lro| [km]:", mag_lro)

    #---GS to Moon ----
    print("\n#---GS TO MOON, rho_moon, HORIZONS(`) & Skyfield(*) Computed----")
    rho_moon = r_moon - r_site
    rho_moon_mag = np.linalg.norm(np.array(rho_moon), axis=1)
    print("           Moon Azimuth [deg]`:", np.array(df['MOON Azimuth [deg]']).transpose())
    print("         Moon Elevation [deg]`:", np.array(df['MOON Elevation [deg]']).transpose())
    print("              Moon Range [km]`:", np.array(df['MOON Range [km]']).transpose())
    print("  rho_moon_mag |rho_moon| [km]:", rho_moon_mag)

    #---GS to LRO ----
    print("\n#---GS TO LRO, rho_lro, HORIZONS(`) & Skyfield(*) Computed----")
    rho_lro = r_lro - r_site
    rho_lro_mag = np.linalg.norm(np.array(rho_lro), axis=1)
    print("          LRO Azimuth [deg]`:", np.array(df['LRO Azimuth [deg]']).transpose())
    print("        LRO Elevation [deg]`:", np.array(df['LRO Elevation [deg]']).transpose())
    print("             LRO Range [km]`:", np.array(df['LRO Range [km]']).transpose())
    print(" rho_lro_mag |rho_lro| [km]*:", rho_lro_mag)

    #---MOON TO LRO ----
    print("\n#---MOON TO LRO, rho_ms, HORIZONS Computed----")
    rho_ms = r_lro - r_moon #selenocentric
    mag_ms = np.linalg.norm(np.array(rho_ms), axis=1)
    lro_alt = mag_ms - utils.R_moon #altitude of LRO
    print("      mag_ms |r_ms| [km]:", mag_ms)
    print("       LRO Altitude [km]:", lro_alt)
    print("   Maximum Altitude [km]:", max(lro_alt))
    print("   Minimum Altitude [km]:", min(lro_alt))

    rho_ms2 = rho_lro - rho_moon
    mag_ms2 = np.linalg.norm(np.array(rho_ms2), axis=1)
    lro_alt2 = mag_ms2 - utils.R_moon
    print("     mag_ms2 |r_ms| [km]:", mag_ms2)
    print("       LRO Altitude [km]:", lro_alt2)
    print("   Maximum Altitude [km]:", max(lro_alt2))
    print("   Minimum Altitude [km]:", min(lro_alt2))

    #---Do some math -----
    print("\n#---Math----")
    print('Moon Radius [km]:', utils.R_moon)
    moon_range = np.array(df['MOON Range [km]'])
    theta_1 = np.arccos(np.divide(utils.R_moon, moon_range)) * utils.rad2deg
    theta_2 = np.arccos(np.divide(utils.R_moon, mag_ms)) * utils.rad2deg
    theta_12 = theta_1 + theta_2
    print('   theta_1 [deg]:', theta_1)
    print('   theta_2 [deg]:', theta_2)
    print('  theta_12 [deg]:', theta_12)
    print()
    print('mag_moon[0]:', mag_moon[0])
    print('  mag_ms[0]:', mag_ms[0])

    theta = np.arccos(np.divide(np.diag(np.dot(-1*rho_moon, rho_ms.transpose())),moon_range*mag_ms))* utils.rad2deg
    print('theta [deg]:', theta)

    # for i, th in enumerate(theta):
    #     print (i, th, theta_12[i], th >= theta_12[i])

    sys.exit()
    #load timescale object, configure time DF
    ts = load.timescale()
    t = ts.tdb_jd(jd=np.array(df['datetime_jd [TDB]']))
    df['utc_iso'] = t.utc_iso(places=6)
    print(df.keys())
    #load almanac
    e = load('de421.bsp')
    earth = e['earth']
    moon = e['moon']
    #setup ground station
    gs = sf.Topos(latitude_degrees=cfg['gs']['latitude'],
                  longitude_degrees=cfg['gs']['longitude'],
                  elevation_m=cfg['gs']['altitude'])
    # print(type(gs_rv))
    od = skyfield.vectorlib.ObserverData()
    gs._snag_observer_data(od, t)

    #---GS TO MOON ----
    print("#---GS TO MOON, rho_moon, Skyfield Computed----")
    moon_rv = (earth+gs).at(t).observe(moon)
    moon_rv_rev = moon.at(t).observe(earth+gs)
    [el_moon,az_moon,rho_moon] = moon_rv.apparent().altaz()
    df['Lunar Elevation [deg]'] = el_moon.degrees
    df['Lunar Azimuth [deg]'] = az_moon.degrees
    df['Lunar Range [km]'] = rho_moon.km


    #print(type(moon_rv))
    #print(moon_rv.position.km)
    r_moon = np.array(moon_rv.position.km)
    v_moon = np.array(moon_rv.velocity.km_per_s)
    mag_r_moon = np.linalg.norm(r_moon.transpose(), axis=1)
    mag_v_moon = np.linalg.norm(v_moon.transpose(), axis=1)

    print("  Lunar Azimuth [deg]*:", np.array(df['Lunar Azimuth [deg]']).transpose())
    print("Lunar Elevation [deg]*:", np.array(df['Lunar Elevation [deg]']).transpose())
    print("     Lunar Range [km]*:", np.array(df['Lunar Range [km]']).transpose())
    #print('  r moon', r_moon)
    #print('  v moon', v_moon)
    print('      mag_r moon [km]*:', mag_r_moon)

    #print('  mag_v moon', mag_v_moon)
    #---GS TO LRO ----
    print("#---GS TO LRO, rho_sat, HORIZONS Computed ----")
    pos_keys = ['x [AU]','y [AU]','z [AU]']
    vel_keys = ['x_dot [AU/d]', 'y_dot [AU/d]', 'z_dot [AU/d]']
    #print(np.array(df[pos_keys]))
    lro_rv = skyfield.positionlib.Geocentric(np.array(df[pos_keys]).transpose(),
                                             np.array(df[vel_keys]).transpose(),
                                             t=t,
                                             center=399,
                                             target=-85)
    gs_rv = (gs).at(t) # geocentric position of gs
    #rho_sat = (lro_rv-gs_rv)
    #print(np.array(lro_rv.position.km)[:3])
    #print(np.array(gs_rv.position.km)[:3])
    r_sat = np.array(lro_rv.position.km) - np.array(gs_rv.position.km)
    #v_sat = np.array(rho_sat.velocity.km_per_s)
    mag_r = np.linalg.norm(np.array(r_sat).transpose(), axis=1)
    #mag_r = lro_rv.distance().km
    #mag_v = np.linalg.norm(v_sat.transpose(), axis=1)

    print("    LRO Azimuth [deg]`:", np.array(df['LRO Azimuth [deg]']).transpose())
    print("  LRO Elevation [deg]`:", np.array(df['LRO Elevation [deg]']).transpose())
    print("       LRO Range [km]`:", np.array(df['LRO Range [km]']).transpose())
    #print('  r_sat', r_sat)
    #print('  v_sat', v_sat)
    print('       mag_r LRO [km]*:', mag_r)
    print('   delta LRO `to* [km]:', (np.array(df['LRO Range [km]']).transpose()-mag_r))
    #print('  mag_v', mag_v)

    #---Moon to LRO------
    print("#---MOON TO LRO----")
    rho_ms = lro_rv - moon_rv
    r_ms = np.array(rho_ms.position.km)
    v_ms = np.array(rho_ms.velocity.km_per_s)
    # r_ms = r_sat - r_moon
    # v_ms = v_sat - v_moon
    mag_r_ms = np.linalg.norm(r_ms.transpose(), axis=1)
    mag_v_ms = np.linalg.norm(v_ms.transpose(), axis=1)
    print('  r_ms', r_ms)
    print('  v_ms', v_ms)
    print('  mag_r', mag_r_ms)
    print('  mag_v', mag_v_ms)


    #r_dot = np.dot(r,v) / np.linalg.norm(r, axis=1)
    #print('  r_dot', r_dot)
    print("* - Skyfield Computed")
    print("` - HORIZONS Computed")
    sys.exit()
    print(type(dif))
    print(dif.position.km, dif.velocity)
    #state_vector = dif.at(t)

    #[el,az,rho]=dif.apparent()

    sys.exit()















    for i, t in enumerate(vec[t_key]):
        #t = ts.tdb(jd=np.array(vec[t_key])[0])
        t_step = ts.tdb(jd=t)
        print (t_step, t_step.tdb, t_step.utc_iso())
        lro = skyfield.positionlib.Geocentric(tuple(vec[pos_keys][i]),
                                              tuple(vec[vel_keys][i]),
                                              t=t_step,
                                              center=cfg['body']['origin_center'],
                                              target=-85)
        print('lro', lro)
        #print(lro.position.km)
        #rv = (e['earth']+gs).at(t)


        geo_rv = gs.at(t_step)
        od = skyfield.vectorlib.ObserverData()
        gs._snag_observer_data(od, t_step)
        od.gcrs_position = geo_rv.position

        print(od.altaz_rotation, od.ephemeris, od.gcrs_position)
        print('gs', geo_rv)
        dif = lro - geo_rv
        dif = skyfield.positionlib.Geometric(dif.position,
                                             dif.velocity,
                                             observer_data=od)
        print(dif)
        print(dif.altaz())
        #print(dif.apparent())
        #print(dif.altaz())

    #print(lro.distance())

    #print(lro.at(t))
    sys.exit()
