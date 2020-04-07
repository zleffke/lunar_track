#! /usr/bin/python3
#########################################
#    Title: Observe
#  Project: Lunar Track
#     Date: Mar 2020
#   Author: Zach Leffke, KJ4QLP
# Comments:
# - reads in downloaded HORIZONS data from CSV file
# - does math.....
#########################################
import sys
import os
import math
import string
import argparse
import json
import datetime
import time
import numpy as np
import pandas as pd

import utilities.plotting
import utilities.satellite


deg2rad = math.pi / 180
rad2deg = 180 / math.pi
c       = float(299792458)    #[m/s], speed of light

if __name__ == "__main__":
    """ Main entry point to start the service. """
    #--------START Command Line argument parser------------------------------------------------------
    parser = argparse.ArgumentParser(description="Lunar Spacecraft Observations",
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
                       default="observe_config_lro.json",
                       help="Configuration File",
                       action="store")
    args = parser.parse_args()
    #--------END Command Line argument parser------------------------------------------------------

    #subprocess.run(["reset"])

    #print(sys.path)
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

    #set up input data path & file
    in_f = '/'.join([cwd, cfg['data_path'], cfg['data_file']])
    if not os.path.isfile(in_f) == True:
        print('ERROR: Invalid Input Data File: {:s}'.format(in_f))
        sys.exit()
    print('Importing Data File: {:s}'.format(in_f))
    #import to pandas dataframe
    df = pd.read_csv(in_f)
    df.name = cfg['data_file'].split("_")[1]

    df['datetime'] = pd.to_datetime(df['datetime_str'])
    doppler = utilities.satellite.Doppler_Shift(cfg['satellite']['freq'], df['Range Rate [km/sec]']*1000.0)
    df['Doppler Center [Hz]'] = doppler['center']
    df['Doppler Offset [Hz]'] = doppler['offset']

    print(df.keys().to_list())
    lam = utilities.satellite.Freq_2_Lambda(cfg['satellite']['freq'])
    df['Path Loss [dB]'] = utilities.satellite.Path_Loss(df['Range [km]']*1000.0, lam)

    print(df)
    #utilities.plotting.plot_offset_idx(0, df, None, 0)
    #utilities.plotting.plot_offset_ts(0, df, None, 0)
    utilities.plotting.plot_az_el_ts(0, df, None, 0)
    utilities.plotting.plot_az_el_polar_ts(0, df, None, 0)

    sys.exit()
