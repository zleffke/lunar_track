#! /usr/bin/python3
#########################################
#    Title: JPL Horizons Data Downloader
#  Project: Lunar Track
#     Date: Mar 2020
#   Author: Zach Leffke, KJ4QLP
# Comments:
# - uses astroquery to download data from JPL HORIZONS Database and stores in CSV format.
# - downloads spacecraft observation data as opposed to state vectors or orbital elements.
#########################################
import sys
import os
import math
import string
import argparse
import json
import numpy as np
import pandas as pd
from astroquery.jplhorizons import Horizons

deg2rad = math.pi / 180
rad2deg = 180 / math.pi
c       = float(299792458)    #[m/s], speed of light
au2km   = 149597870.700 #km/au

if __name__ == "__main__":
    #--------START Command Line argument parser------------------------------------------------------
    parser = argparse.ArgumentParser(description="Download HORIZONS data and stores as CSV",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cwd = os.getcwd()
    cfg_fp_default = '/'.join([cwd, 'config'])
    cfg = parser.add_argument_group('Configuration File')
    cfg.add_argument('--cfg_path',
                       dest='cfg_path',
                       type=str,
                       default='/'.join([os.getcwd(), 'config']),
                       help="Configuration File Path",
                       action="store")
    cfg.add_argument('--cfg_file',
                       dest='cfg_file',
                       type=str,
                       default="horizons_query_config_lro.json",
                       help="Configuration File",
                       action="store")
    args = parser.parse_args()
    #--------END Command Line argument parser------------------------------------------------------

    # Config File Import
    fp_cfg = '/'.join([args.cfg_path,args.cfg_file])
    #print (fp_cfg)
    if not os.path.isfile(fp_cfg) == True:
        print('ERROR: Invalid Configuration File: {:s}'.format(fp_cfg))
        sys.exit()
    print('Importing configuration File: {:s}'.format(fp_cfg))
    with open(fp_cfg, 'r') as json_data:
        cfg = json.load(json_data)
        json_data.close()
    #print (cfg)
    #set up output data path
    data_path = '/'.join([cwd, cfg['data_path']])
    if not os.path.exists(data_path):
        print('Data path doesn\'t exist, creating...: {:s}'.format(fp_load))
        os.makedirs(data_path)
    of = "_".join([cfg['gs']['name'],
                   cfg['horizons']['name'],
                   cfg['horizons']['epochs']['start'].replace(" ", "T")+"Z",
                   cfg['horizons']['epochs']['stop'].replace(" ", "T")+"Z",
                   cfg['horizons']['epochs']['step']])
    of = ".".join([of,"csv"])
    o_fp = "/".join([data_path, of])
    print("            Data Output File: {:s}".format(o_fp))
    #load = sf.Loader(fp_load, verbose=True)

    #create observer location
    gs = {'lon': cfg['gs']['longitude'],
          'lat': cfg['gs']['latitude'],
          'elevation': cfg['gs']['altitude_m']/1000.0}

    #create query object
    obj = Horizons(id       = cfg['horizons']['name'],
                   location = gs,
                   id_type  = cfg['horizons']['id_type'],
                   epochs   = cfg['horizons']['epochs'])
    #execute query
    #print(obj)
    eph = obj.ephemerides(quantities=cfg['horizons']['quantities'])
    print(eph)
    print(eph.columns)
    #print(vec.columns)
    keys = ['datetime_str', 'AZ', 'AZ_rate', 'EL', 'EL_rate',
            'delta', 'delta_rate', 'lighttime']
    df = pd.DataFrame(np.array(eph[keys]), columns=keys)
    df['AZ_rate'] = df['AZ_rate'] * 4.62963e-6 # arcsec per min to degrees per sec
    df['EL_rate'] = df['EL_rate'] * 4.62963e-6 # arcsec per min to degrees per sec
    df['delta']   = df['delta'] * au2km # AU to km
    df['lighttime'] = df['lighttime'] * 60 #minutes to seconds
    column_map = {'AZ':'Azimuth [deg]', 'AZ_rate':'Azimuth Rate [deg/sec]',
                  'EL':'Elevation [deg]', 'EL_rate':'Elevation Rate [deg/sec]',
                  'delta':'Range [km]', 'delta_rate':'Range Rate [km/sec]',
                  'lighttime':'1-Way Prop Delay [sec]'}
    df=df.rename(columns=column_map)
    print(df)
    print("Exporting to: {:s}".format(o_fp))
    #Export Data
    df.to_csv(o_fp, index=False)
    sys.exit()
