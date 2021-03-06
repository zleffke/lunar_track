# lunar_track
Scripts for Lunar tracking

General idea is to use a combination of [Astroquery](https://pypi.org/project/astroquery/) to pull observer data from [JPL HORIZONS](https://ssd.jpl.nasa.gov/horizons.cgi) database as well as [Skyfield](https://pypi.org/project/skyfield/).  Scripts will also provide plotting utilities.  Also hope to provide real time tracking scripts in order to feed pointing programs for telescope and antenna pointing control.  General philosophy is to build a modular combination of scripts and execute them in particular order as opposed to building one monolithic script.

## Dependencies
All scripts are written for Python 3+.
```
sudo pip3 install astroquery skyfield matplotlib pandas
```

## Directory Structures
Each directory will house different scripts as they become more completed.

* sandbox - junk scripts for initial testing.
* horizons_query - astroquery based script.  Pulls data from JPL HORIZONS database and stores in CSV format.  Idea is to copy queried data into the 'data' directory for scripts that need the data.
* lro_occultation - scripts for simulating/computing Lunar occultation of the Lunar Reconnaissance Orbiter (LRO).

### Sub-Directories:
Each directory will contain some common sub directories.
* data - location for data download.
* config - location for configuration files.
* utilities - helper utilities for things like plotting.

## Script Descriptions:

#### jpl_horizons_to_csv.py
* Queries JPL HORIZONS database for observation information (time, az, el, range, az rate, el_rate, range rate, etc.)
* Config file is in config folder:  `horizons_query_config_<sc_id>.json`
* converts data from query response into pandas dataframe
* formats headers and converts values to more useful units in the dataframe
* exports data to data folder with filename syntax: `<GS>_<SC>_<START>_<STOP>_<STEP>.csv`

## Directories:

* data - location for data download.
* config - location for configuration files.
* utilities - helper utilities for things like plotting.
