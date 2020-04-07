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
import subprocess
import numpy as np
import pandas as pd

import utilities.utilities as utils

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
                       default="horizons_query_lro_observe.json",
                       help="Configuration File",
                       action="store")
    cfg.add_argument('--export',
                       dest='export',
                       type=int,
                       default=0,
                       help="If set, export to CSV",
                       action="store")
    args = parser.parse_args()
    #--------END Command Line argument parser------------------------------------------------------
    subprocess.run(["reset"])
    # Config File Import
    cfg = utils.import_configs(args)
    #set up output data path
    data_path = '/'.join([cfg['cwd'], cfg['data_path']])
    if not os.path.exists(data_path):
        print('Data path doesn\'t exist, creating...: {:s}'.format(data_path))
        os.makedirs(data_path)

    if cfg['horizons']['type'] == "vectors":
        df = utils.query_vectors(cfg)
        o_fp = utils.generate_vectors_file_name(cfg)
    if cfg['horizons']['type'] == "ephemerides":
        df = utils.query_ephemerides(cfg)
        o_fp = utils.generate_ephemerides_file_name(cfg)
    if cfg['horizons']['type'] == "elements":
        pass

    print(df)
    #Export Data
    print("Export Path: {:s}".format(o_fp))
    if args.export:
        print("Exporting...")
        df.to_csv(o_fp, index=False, date_format="%Y-%m-%dT%H:%M:%S.%fZ")
    else:
        print("Export flag not set")
    sys.exit()
    #print(vec.columns)

    print(df)
