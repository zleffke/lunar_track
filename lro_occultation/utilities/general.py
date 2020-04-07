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

deg2rad = math.pi / 180
rad2deg = 180 / math.pi
c       = float(299792458)    #[m/s], speed of light
#R_moon = 1736.0 #Lunar Polar Radius, km
R_moon = 1737.4 #Lunar Mean Radius, km
#R_moon = 1738.1 #Lunar Equatorial Radius, km

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
    #print (fp_cfg)
    if not os.path.isfile(fp_cfg) == True:
        print('ERROR: Invalid Configuration File: {:s}'.format(fp_cfg))
        sys.exit()
    print('Importing configuration File: {:s}'.format(fp_cfg))
    with open(fp_cfg, 'r') as json_data:
        cfg = json.load(json_data)
        json_data.close()
    cfg['cwd'] = os.getcwd()
    print(json.dumps(cfg, indent=4))
    return cfg

def HorizonsCSVToDataFrame(path, file):
    ''' import JPL HORIZONS data from CSV '''
    fp = '/'.join([os.getcwd(),path,file])
    if not os.path.isfile(fp) == True:
        print('ERROR: Invalid File: {:s}'.format(fp))
        sys.exit()
    print('Importing HORIZONS CSV File: {:s}'.format(fp))
    df = pd.read_csv(fp)
    df_name = file.split("_")
    name = "_".join([df_name[0], df_name[1]])
    drop_keys = ['datetime_str [UTC]', 'datetime_jd [UTC]', 'datetime [ISO UTC]']
    df['datetime_utc'] = pd.to_datetime(df['datetime [ISO UTC]'], utc=True)
    df = df.drop(labels=drop_keys, axis=1)
    columns = list(df.keys())
    columns.remove('datetime_utc')
    column_map = {}
    for i, key in enumerate(columns): column_map[key] = "{:s} {:s}".format(df_name[0], key)
    print(column_map)
    df=df.rename(columns=column_map)

    df.name = name
    print(df.info())
    return df
