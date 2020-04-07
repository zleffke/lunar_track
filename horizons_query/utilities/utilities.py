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
from astroquery.jplhorizons import Horizons

deg2rad = math.pi / 180
rad2deg = 180 / math.pi
c       = float(299792458)    #[m/s], speed of light
au2km   = 149597870.700 #km/au
day2sec = 86400.0 #days to seconds

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
    cfg['cwd'] = os.getcwd()
    print (cfg)
    return cfg

def query_ephemerides(cfg):
    #create observer location
    gs = {'lon': cfg['gs']['longitude'],
          'lat': cfg['gs']['latitude'],
          'elevation': cfg['gs']['altitude_m']/1000.0}
    #create query object
    if 'id' in cfg['horizons'].keys(): id=cfg['horizons']['id']
    else: id=cfg['horizons']['name']
    obj = Horizons(id       = id,
                   location = gs,
                   id_type  = cfg['horizons']['id_type'],
                   epochs   = cfg['horizons']['epochs'])
    eph = obj.ephemerides(quantities=cfg['horizons']['quantities'])
    keys = ['datetime_str', 'datetime_jd',
            'AZ', 'AZ_rate', 'EL', 'EL_rate',
            'delta', 'delta_rate', 'lighttime']
    df = pd.DataFrame(np.array(eph[keys]), columns=keys)
    df['AZ_rate'] = df['AZ_rate'] * 4.62963e-6 # arcsec per min to degrees per sec
    df['EL_rate'] = df['EL_rate'] * 4.62963e-6 # arcsec per min to degrees per sec
    df['delta']   = df['delta'] * au2km # AU to km
    df['lighttime'] = df['lighttime'] * 60 #minutes to seconds
    column_map = {'datetime_str':'datetime_str [UTC]',
                  'datetime_jd':'datetime_jd [UTC]',
                  'AZ':'Azimuth [deg]', 'AZ_rate':'Azimuth Rate [deg/sec]',
                  'EL':'Elevation [deg]', 'EL_rate':'Elevation Rate [deg/sec]',
                  'delta':'Range [km]', 'delta_rate':'Range Rate [km/sec]',
                  'lighttime':'1-Way Prop Delay [sec]'}
    df=df.rename(columns=column_map)
    df['datetime [ISO UTC]'] = pd.to_datetime(df['datetime_str [UTC]'], format="%Y-%b-%d %H:%M")
    return df

def generate_ephemerides_file_name(cfg):
    data_path = '/'.join([cfg['cwd'], cfg['data_path']])
    if not os.path.exists(data_path):
        print('Data path doesn\'t exist, creating...: {:s}'.format(data_path))
        os.makedirs(data_path)
    of = "_".join([cfg['horizons']['name'].upper(),
                   cfg['gs']['name'].upper(),
                   cfg['horizons']['type'].upper(),
                   cfg['horizons']['epochs']['start'].replace(" ", "T")+"Z",
                   cfg['horizons']['epochs']['stop'].replace(" ", "T")+"Z",
                   cfg['horizons']['epochs']['step']])
    of = ".".join([of,"csv"])
    o_fp = "/".join([data_path, of])
    print("            Data Output File: {:s}".format(o_fp))
    return o_fp

def query_vectors(cfg):
    print("Querying Vectors")
    #create query object
    if 'id' in cfg['horizons'].keys(): id=cfg['horizons']['id']
    else: id=cfg['horizons']['name']
    obj = Horizons(id       = id,
                   location = cfg['horizons']['coord_origin'],
                   id_type  = cfg['horizons']['id_type'],
                   epochs   = cfg['horizons']['epochs'])
    eph = obj.vectors()
    #print(eph)
    print(eph.columns)
    #return eph
    keys = ['datetime_str', 'datetime_jd', 'x', 'y', 'z', 'vx', 'vy', 'vz']
    df = pd.DataFrame(np.array(eph[keys]), columns=keys)
    df['x'] = df['x'] * au2km
    df['y'] = df['y'] * au2km
    df['z'] = df['z'] * au2km
    df['vx'] = df['vx'] * au2km / day2sec
    df['vy'] = df['vy'] * au2km / day2sec
    df['vz'] = df['vz'] * au2km / day2sec
    column_map = {'datetime_str':'datetime_str [UTC]',
                  'datetime_jd':'datetime_jd [UTC]',
                  'x':'x [km]', 'y':'y [km]', 'z':'z [km]',
                  'vx':'x_dot [km/sec]', 'vy':'y_dot [km/sec]', 'vz':'z_dot [km/sec]'}
    df=df.rename(columns=column_map)
    df['datetime [ISO UTC]'] = pd.to_datetime(df['datetime_str [UTC]'],
                                              format="A.D. %Y-%b-%d %H:%M:%S.%f")
    return df

def generate_vectors_file_name(cfg):
    data_path = '/'.join([cfg['cwd'], cfg['data_path']])
    if not os.path.exists(data_path):
        print('Data path doesn\'t exist, creating...: {:s}'.format(data_path))
        os.makedirs(data_path)
    of = "_".join([cfg['horizons']['name'].upper(),
                   cfg['horizons']['coord_origin'].upper(),
                   cfg['horizons']['type'].upper(),
                   cfg['horizons']['epochs']['start'].replace(" ", "T")+"Z",
                   cfg['horizons']['epochs']['stop'].replace(" ", "T")+"Z",
                   cfg['horizons']['epochs']['step']])
    of = ".".join([of,"csv"])
    o_fp = "/".join([data_path, of])
    print("            Data Output File: {:s}".format(o_fp))
    return o_fp
