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

deg2rad = math.pi / 180
rad2deg = 180 / math.pi
c       = float(299792458)    #[m/s], speed of light
R_moon = 1736.0 #Lunar Polar Radius, km


def setup_skyfield_loader(cfg):
    ''' setup skyfield loader '''
    cwd = os.getcwd()
    fp_load = '/'.join([cwd, cfg['data']['sf_path']])
    if not os.path.exists(fp_load):
        print('Skyfield data path doesn\'t exist: {:s}'.format(fp_load))
        print("creating...")
        os.makedirs(fp_load)
    load = sf.Loader(fp_load, verbose=True)
    return load

def import_configs(args):
    ''' setup configuration data '''
    fp_cfg = '/'.join([args.cfg_path,args.cfg_file])
    print (fp_cfg)
    if not os.path.isfile(fp_cfg) == True:
        print('ERROR: Invalid Configuration File: {:s}'.format(fp_cfg))
        sys.exit()
    print('Importing configuration File: {:s}'.format(fp_cfg))
    with open(fp_cfg, 'r') as json_data:
        cfg = json.load(json_data)
        json_data.close()
    print (cfg)
    return cfg

def import_horizons_data(cfg):
    ''' import JPL HORIZONS data from CSV '''
    fp_lro_vec = '/'.join([cwd,
                           cfg['data']['lro_path'],
                           cfg['data']['lro_vec_file']])
    if not os.path.isfile(fp_lro_vec) == True:
        print('ERROR: Invalid LRO Vector File: {:s}'.format(fp_lro_vec))
        sys.exit()
    print('Importing LRO Vector File: {:s}'.format(fp_lro_vec))
    df_vec = pd.read_csv(fp_lro_vec)
    df_name = cfg['data']['lro_vec_file'].split("_")
    df_vec.name = "_".join([df_name[0], df_name[1]])

    fp_lro_pnt = '/'.join([cwd,
                           cfg['data']['lro_path'],
                           cfg['data']['lro_pnt_file']])
    if not os.path.isfile(fp_lro_vec) == True:
        print('ERROR: Invalid LRO Pointing File: {:s}'.format(fp_lro_pnt))
        sys.exit()
    print('Importing LRO Vector File: {:s}'.format(fp_lro_pnt))
    df_pnt = pd.read_csv(fp_lro_pnt)
    df_name = cfg['data']['lro_pnt_file'].split("_")
    df_pnt.name = "_".join([df_name[0], df_name[1]])

    return [df_vec,df_pnt]

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
                       default="lro_occult_config.json",
                       help="Configuration File",
                       action="store")
    args = parser.parse_args()
    #--------END Command Line argument parser------------------------------------------------------
    subprocess.run(["reset"])
    #setup configuration file info
    cfg = import_configs(args)
    #set up Skyfield data Loader path
    #load = setup_skyfield_loader(cfg)
    #Import Horizons Vector and Pointing Data
    [df_vec, df_pnt] = import_horizons_data(cfg)
    df = df_vec.merge(df_pnt, how='inner', on='datetime_jd [TDB]')
    df.name = df_pnt.name.split('_')[0]
    columns = list(df_pnt.keys())
    columns.remove('datetime_jd [TDB]')
    column_map = {}
    for i, key in enumerate(columns): column_map[key] = "{:s} {:s}".format(df.name, key)
    print(column_map)
    df=df.rename(columns=column_map)
    del df_vec #destroy temporary dataframe, free up memory
    del df_pnt #destroy temporary dataframe, free up memory
    print(df)













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

    np.set_printoptions(precision=10, linewidth=150)
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
