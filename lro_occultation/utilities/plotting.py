#!/usr/bin/env python3
#########################################
#   Title: pyephem helper utilities     #
# Project: TLE Match                    #
#    Date: Jan 2018                     #
#  Author: Zach Leffke, KJ4QLP          #
#########################################
import sys
import os
import math
import string
import struct
import numpy as np
import scipy
import datetime as dt
import pytz
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.dates import DateFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.dates as mdates
import matplotlib.colors as colors


deg2rad = math.pi / 180
rad2deg = 180 / math.pi
c       = float(299792458)    #[m/s], speed of light



def plot_moon_view2(idx, df, o_path=None, save = 0, R_moon_deg=0.26):
    a = df.name
    x = df['dAz [deg]'].values.tolist() #independent variable, input
    y1 = df['dEl [deg]'].values.tolist()

    df_los = df.loc[df['LOS'] == True]
    df_nlos = df.loc[df['LOS'] == False]

    x_los = df_los['dAz [deg]'].values.tolist()
    y_los = df_los['dEl [deg]'].values.tolist()
    x_nlos = df_nlos['dAz [deg]'].values.tolist()
    y_nlos = df_nlos['dEl [deg]'].values.tolist()

    limit = 1

    moon_circle = plt.Circle((0, 0), 0.26, color='green', fill = 0)
    #---- START Figure 1 ----
    xinch = 10
    yinch = 10
    fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
    ax1 = fig1.add_subplot(1,1,1)
    ax1.add_artist(moon_circle)
    #ax2 = ax1.twinx()
    ax1.set_xlim(xmin = -limit, xmax = limit)
    ax1.set_ylim(ymin = -limit, ymax = limit)
    #Configure Grids
    ax1.xaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')
    ax1.yaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')

    #Configure Labels and Title
    ax1.set_xlabel('dAz [deg]')
    ax1.set_ylabel('dEl [deg]')
    title = '{:s}'.format(a)
    ax1.set_title(title)

    #Plot Data
    #ax1.plot(x, y1, linestyle = '-', label="{:s}".format(a), markersize=1, markeredgewidth=0)
    #ax1.scatter(x, y1, c=df['datetime_utc'], s=2, cmap='plasma', alpha=0.75)

    ax1.plot(x_los, y_los, color='red', linewidth=1)
    ax1.plot(x_nlos, y_nlos, color='black', linewidth=1)
    #ax1.scatter(x, y1, c=df['datetime_utc'], s=30, cmap='YlOrRd', alpha=0.75)
    #ax1.plot(x, y1, c=df['datetime_utc'], cmap='plasma', alpha=0.75)

    #Save Figure
    if ((o_path != None) and save):

        fig_f = '{:s} Doppler Offset.png'.format(a)
        fig_fp = '/'.join([o_path,fig_f])
        print("Output Path:", fig_fp)
        print("Saving Figure {:2d}: {:s}".format(idx, fig_fp))
        fig1.savefig(fig_fp)

    plt.show()
    #plt.close(fig1)
    idx += 1
    return idx


def plot_moon_view(idx, df, o_path=None, save = 0):
    a = df.name
    x = df['dAz [deg]'].values.tolist() #independent variable, input
    y1 = df['dEl [deg]'].values.tolist()


    x_max = max(abs(df['dAz [deg]']))
    print(x_max)
    print(round(x_max,1))
    y_max = max(abs(df['dEl [deg]']))
    print(y_max)
    print(round(y_max,1))
    limit = max([round(y_max,1)+.1, round(x_max,1)+.1])
    print(limit)
    moon_circle = plt.Circle((0, 0), 0.26, color='black', fill = 0)
    #---- START Figure 1 ----
    xinch = 10
    yinch = 10
    fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
    ax1 = fig1.add_subplot(1,1,1)
    ax1.add_artist(moon_circle)
    #ax2 = ax1.twinx()
    ax1.set_xlim(xmin = -limit, xmax = limit)
    ax1.set_ylim(ymin = -limit, ymax = limit)
    #Configure Grids
    ax1.xaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')
    ax1.yaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')

    #Configure Labels and Title
    ax1.set_xlabel('dAz [deg]')
    ax1.set_ylabel('dEl [deg]')
    title = '{:s}'.format(a)
    ax1.set_title(title)

    #Plot Data
    #ax1.plot(x, y1, linestyle = '-', label="{:s}".format(a), markersize=1, markeredgewidth=0)
    #ax1.scatter(x, y1, c=df['datetime_utc'], s=2, cmap='plasma', alpha=0.75)

    ax1.plot(x, y1, color='black', linewidth=1)
    ax1.scatter(x, y1, c=df['datetime_utc'], s=30, cmap='YlOrRd', alpha=0.75)
    #ax1.plot(x, y1, c=df['datetime_utc'], cmap='plasma', alpha=0.75)

    #Save Figure
    if ((o_path != None) and save):

        fig_f = '{:s} Doppler Offset.png'.format(a)
        fig_fp = '/'.join([o_path,fig_f])
        print("Output Path:", fig_fp)
        print("Saving Figure {:2d}: {:s}".format(idx, fig_fp))
        fig1.savefig(fig_fp)

    plt.show()
    #plt.close(fig1)
    idx += 1
    return idx


def plot_deltas(idx, df, o_path=None, save = 0):
    a = df.name
    x = mdates.datestr2num(df['datetime_str'].values.tolist())
    dAz  = df['Azimuth Delta [deg]'].values.tolist()
    dEl  = df['Elevation Delta [deg]'].values.tolist()
    dRho = df['Range Delta [km]'].values.tolist()
    #---- START Figure 1 ----
    xinch = 14
    yinch = 7
    fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
    ax1 = fig1.add_subplot(1,1,1)
    ax2 = ax1.twinx()

    #Configure Grids
    ax1.xaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')
    ax1.yaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')

    #Configure Labels and Title
    ax1.set_xlabel('Time [UTC]')
    ax1.set_ylabel('Delta Azimuth & Delta Elevation [deg]')
    ax2.set_ylabel('Delta Range [km]')
    #ax2.set_ylim(-1500,1500)
    title = 'Lunar to {:s} Pointing and Range Deltas'.format(a)
    ax1.set_title(title)

    #Plot Data
    ln1 = ax1.plot(x, dAz, color='r', linestyle = '-', label="$\Delta$ Azimuth [deg]", markersize=1, markeredgewidth=0)
    ln2 = ax1.plot(x, dEl, color='b', linestyle = '-', label="$\Delta$ Elevation [deg]", markersize=1, markeredgewidth=0)
    ln3 = ax2.plot(x, dRho, color='g', linestyle = '-', label="$\Delta$ Range [km]", markersize=1, markeredgewidth=0)

    #Formate X Axis Timestamps
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d\n %H:%M:%S"))
    #ax1.set_xlim(x2[0] - dt.timedelta(minutes = 30), x2[-1] +  dt.timedelta(minutes=30))
    for label in ax1.xaxis.get_ticklabels():
        label.set_rotation(45)
    fig1.subplots_adjust(bottom=0.2)

    #Configure Legend
    lines = ln1+ln2+ln3
    labels = [l.get_label() for l in lines]
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width*0.90, box.height])
    #ax1.legend(lines, labels, loc='center left', numpoints = 1, bbox_to_anchor=(0.8, 0.95))
    ax1.legend(lines, labels, loc='center left', bbox_to_anchor=(1.05, 0.95))
    #ax1.legend(lines, labels, numpoints = 1)
    #ax1.legend(lns, labs, loc=0)

    #Save Figure
    if ((o_path != None) and save):

        fig_f = '{:s} Doppler Offset.png'.format(a)
        fig_fp = '/'.join([o_path,fig_f])
        print("Output Path:", fig_fp)
        print("Saving Figure {:2d}: {:s}".format(idx, fig_fp))
        fig1.savefig(fig_fp)

    plt.show()
    #plt.close(fig1)
    idx += 1
    return idx

def plot_az_el_polar_ts(idx, df, o_path=None, save = 0):
    a = df.name
    #x = [dt.datetime.utcfromtimestamp(element*1e-9) for element in df['datetime_str'].values.tolist()] #independent variable, input
    #x = mdates.num2date(df['datetime'].values.tolist())
    x = mdates.datestr2num(df['datetime_str'].values.tolist())
    az = (df['Azimuth [deg]']*deg2rad).values.tolist()
    el = df['Elevation [deg]'].values.tolist()
    #---- START Figure 1 ----
    xinch = 14
    yinch = 7
    fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
    ax1 = fig1.add_subplot(1,1,1, projection='polar')
    ax1.set_rlim(90, 0, 1)
    ax1.set_theta_offset(0.5*np.pi)
    ax1.set_theta_direction(-1)
    title = '{:s} Azimuth and Elevation'.format(a)
    ax1.set_title(title)
    #ax1.plot(az, el)
    ax1.scatter(az, el, c=az, s=1, cmap='plasma', alpha=0.75)

    plt.show()
    #plt.close(fig1)
    idx += 1
    return idx

def plot_az_el_ts(idx, df, o_path=None, save = 0):
    a = df.name
    #x = [dt.datetime.utcfromtimestamp(element*1e-9) for element in df['datetime_str'].values.tolist()] #independent variable, input
    #x = mdates.num2date(df['datetime'].values.tolist())
    x = mdates.datestr2num(df['datetime_str'].values.tolist())
    y1 = df['Azimuth [deg]'].values.tolist()
    y2 = df['Elevation [deg]'].values.tolist()
    #---- START Figure 1 ----
    xinch = 14
    yinch = 7
    fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
    ax1 = fig1.add_subplot(1,1,1)
    ax2 = ax1.twinx()

    #Configure Grids
    ax1.xaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')
    ax1.yaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')

    #Configure Labels and Title
    ax1.set_xlabel('Time [UTC]')
    ax1.set_ylabel('Azimuth [deg]')
    ax2.set_ylabel('Elevation [deg]')
    title = '{:s} Azimuth and Elevation'.format(a)
    ax1.set_title(title)

    #Plot Data
    ln1 = ax1.plot(x, y1, color='r', linestyle = '-', label="Azimuth", markersize=1, markeredgewidth=0)
    ln2 = ax2.plot(x, y2, color='b', linestyle = '-', label="Elevation", markersize=1, markeredgewidth=0)

    #Formate X Axis Timestamps
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d\n %H:%M:%S"))
    #ax1.set_xlim(x2[0] - dt.timedelta(minutes = 30), x2[-1] +  dt.timedelta(minutes=30))

    for label in ax1.xaxis.get_ticklabels():
        label.set_rotation(45)

    fig1.subplots_adjust(bottom=0.2)

    #Configure Legend
    box = ax1.get_position()
    #h1, l1 = ax1.get_legend_handles_labels()
    #h2, l2 = ax2.get_legend_handles_labels()
    #ax1.legend(h1, l1, loc='center left', numpoints = 1, bbox_to_anchor=(0.8, 0.85))
    #ax2.legend(h2, l2, loc='center left', numpoints = 1, bbox_to_anchor=(0.8, 0.85))

    # added these three lines
    lines = ln1+ln2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='center left', numpoints = 1, bbox_to_anchor=(0.1, 0.85))
    #ax1.legend(lns, labs, loc=0)

    #Save Figure
    if ((o_path != None) and save):

        fig_f = '{:s} Doppler Offset.png'.format(a)
        fig_fp = '/'.join([o_path,fig_f])
        print("Output Path:", fig_fp)
        print("Saving Figure {:2d}: {:s}".format(idx, fig_fp))
        fig1.savefig(fig_fp)

    plt.show()
    #plt.close(fig1)
    idx += 1
    return idx

def plot_offset_idx(idx, df, o_path=None, save = 0):
    a = df.name
    x = df.index.values.tolist() #independent variable, input
    y = df['Doppler Offset [Hz]'].values.tolist()
    #---- START Figure 1 ----
    xinch = 14
    yinch = 7
    fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
    ax1 = fig1.add_subplot(1,1,1)
    #ax2 = ax1.twinx()

    #Configure Grids
    ax1.xaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')
    ax1.yaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')

    #Configure Labels and Title
    ax1.set_xlabel('Index')
    ax1.set_ylabel('Doppler Offset [Hz]')
    title = '{:s} Doppler Offset [Hz]'.format(a)
    ax1.set_title(title)

    #Plot Data
    ax1.plot(x, df['Doppler Offset [Hz]'], linestyle = '-', label="{:s}".format(a), markersize=1, markeredgewidth=0)

    #Save Figure
    if ((o_path != None) and save):

        fig_f = '{:s} Doppler Offset.png'.format(a)
        fig_fp = '/'.join([o_path,fig_f])
        print("Output Path:", fig_fp)
        print("Saving Figure {:2d}: {:s}".format(idx, fig_fp))
        fig1.savefig(fig_fp)

    plt.show()
    #plt.close(fig1)
    idx += 1
    return idx

def plot_offset_ts(idx, df, o_path=None, save = 0):
    a = df.name
    #x = [dt.datetime.utcfromtimestamp(element*1e-9) for element in df['datetime_str'].values.tolist()] #independent variable, input
    #x = mdates.num2date(df['datetime'].values.tolist())
    x = mdates.datestr2num(df['datetime_str'].values.tolist())

    y = df['Doppler Offset [Hz]'].values.tolist()
    #---- START Figure 1 ----
    xinch = 14
    yinch = 7
    fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
    ax1 = fig1.add_subplot(1,1,1)
    #ax2 = ax1.twinx()

    #Configure Grids
    ax1.xaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')
    ax1.yaxis.grid(True,'major', linewidth=1)
    ax1.yaxis.grid(True,'minor')

    #Configure Labels and Title
    ax1.set_xlabel('Time [UTC]')
    ax1.set_ylabel('Doppler Offset [Hz]')
    title = '{:s} Doppler Offset [Hz]'.format(a)
    ax1.set_title(title)

    #Plot Data
    ax1.plot(x, df['Doppler Offset [Hz]'], linestyle = '-', label="{:s}".format(a), markersize=1, markeredgewidth=0)

    #Formate X Axis Timestamps
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d\n %H:%M:%S"))
    #ax1.set_xlim(x2[0] - dt.timedelta(minutes = 30), x2[-1] +  dt.timedelta(minutes=30))

    for label in ax1.xaxis.get_ticklabels():
        label.set_rotation(45)

    fig1.subplots_adjust(bottom=0.2)

    #Configure Legend
    box = ax1.get_position()
    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, loc='center left', numpoints = 1, bbox_to_anchor=(0.8, 0.85))

    #Save Figure
    if ((o_path != None) and save):

        fig_f = '{:s} Doppler Offset.png'.format(a)
        fig_fp = '/'.join([o_path,fig_f])
        print("Output Path:", fig_fp)
        print("Saving Figure {:2d}: {:s}".format(idx, fig_fp))
        fig1.savefig(fig_fp)

    plt.show()
    #plt.close(fig1)
    idx += 1
    return idx

# def plot_2poly_ts(idx, reg_x, p1,p2, o_path, save = 0):
#     a = p1['name'] + '_' + p2['name']
#     x = reg_x
#     y1 = p1['pf']['equation']
#     y2 = p2['pf']['equation']
#     #---- START Figure 1 ----
#     xinch = 14
#     yinch = 7
#     fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
#     ax1 = fig1.add_subplot(1,1,1)
#     #ax2 = ax1.twinx()
#
#     #Configure Grids
#     ax1.xaxis.grid(True,'major', linewidth=1)
#     ax1.yaxis.grid(True,'minor')
#     ax1.yaxis.grid(True,'major', linewidth=1)
#     ax1.yaxis.grid(True,'minor')
#
#     #Configure Labels and Title
#     ax1.set_xlabel('offset')
#     ax1.set_ylabel('Doppler Offset [Hz]')
#     title = '{:s} Doppler Offset [Hz]\n TCA Delta: {:3.3f}'.format(a, p2['tca_delta'])
#     ax1.set_title(title)
#
#     #Plot Data
#     ax1.plot(x, y1, color='r', linestyle = '-', label="{:s}".format(p1['name']), markersize=1, markeredgewidth=0)
#     ax1.plot(x, y2, color='b', linestyle = '-', label="{:s}".format(p2['name']), markersize=1, markeredgewidth=0)
#
#     #Formate X Axis Timestamps
#     #ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d\n %H:%M:%S"))
#     #ax1.set_xlim(x2[0] - dt.timedelta(minutes = 30), x2[-1] +  dt.timedelta(minutes=30))
#
#     #for label in ax1.xaxis.get_ticklabels():
#     #    label.set_rotation(45)
#
#     fig1.subplots_adjust(bottom=0.2)
#
#     #Configure Legend
#     box = ax1.get_position()
#     h1, l1 = ax1.get_legend_handles_labels()
#     ax1.legend(h1, l1, loc='center left', numpoints = 1, bbox_to_anchor=(0.8, 0.85))
#
#     #Save Figure
#     if save:
#
#         fig_f = '{:s} Doppler Offset.png'.format(a)
#         fig_fp = '/'.join([o_path,fig_f])
#         print "Output Path:", fig_fp
#         print "Saving Figure {:2d}: {:s}".format(idx, fig_fp)
#         fig1.savefig(fig_fp)
#
#     plt.show()
#     #plt.close(fig1)
#     idx += 1
#     return idx
#
#
# def plot_multi_doppler_ts(idx, dfs, o_path, save=0):
#     #---- START Figure 1 ----
#     xinch = 14
#     yinch = 7
#     fig1=plt.figure(idx, figsize=(xinch,yinch/.8))
#     ax1 = fig1.add_subplot(1,1,1)
#
#     #Configure Grids
#     ax1.xaxis.grid(True,'major', linewidth=1)
#     ax1.yaxis.grid(True,'minor')
#     ax1.yaxis.grid(True,'major', linewidth=1)
#     ax1.yaxis.grid(True,'minor')
#
#     #Configure Labels and Title
#     ax1.set_xlabel('Index')
#     ax1.set_ylabel('Doppler Offset [Hz]')
#     title = 'Multi Doppler Offset [Hz]'
#     ax1.set_title(title)
#
#     #Plot doppler curves
#     for df in dfs:
#         x = [dt.datetime.utcfromtimestamp(element*1e-9) for element in df['timestamp'].values.tolist()] #independent variable, input
#         if 'FOX' in df.name: col = 'r'
#         else: col = 'b'
#         ax1.plot(x, df['doppler_offset'], \
#                     color=col, linestyle = '-', \
#                     label="{:s}".format(df.name), markersize=1, markeredgewidth=0)
#
#     #Formate X Axis Timestamps
#     ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d\n %H:%M:%S"))
#     #ax1.set_xlim(x2[0] - dt.timedelta(minutes = 30), x2[-1] +  dt.timedelta(minutes=30))
#
#     for label in ax1.xaxis.get_ticklabels():
#         label.set_rotation(45)
#
#     fig1.subplots_adjust(bottom=0.2)
#
#     #Configure Legend
#     box = ax1.get_position()
#     h1, l1 = ax1.get_legend_handles_labels()
#     #ax1.legend(h1, l1, loc='center left', numpoints = 1, bbox_to_anchor=(0.8, 0.85))
#     ax1.legend(prop={'size': 9},ncol=2, bbox_to_anchor=(1,1))#numpoints = 1, bbox_to_anchor=(0.8, 0.85))
#     #Save Figure
#     if save:
#
#         fig_f = 'Multi_Doppler Offset.png'
#         fig_fp = '/'.join([o_path,fig_f])
#         print "Output Path:", fig_fp
#         print "Saving Figure {:2d}: {:s}".format(idx, fig_fp)
#         fig1.savefig(fig_fp)
#
#     plt.show()
#     #plt.close(fig1)
#     idx += 1
#     return idx
