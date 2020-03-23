# lunar_track
Scripts for Lunar tracking

General idea is to use a combination of [Astroquery](https://pypi.org/project/astroquery/) to pull observer data from [JPL HORIZONS](https://ssd.jpl.nasa.gov/horizons.cgi) database as well as [Skyfield](https://pypi.org/project/skyfield/).  Scripts will also provide plotting utilities.  Also hope to provide real time tracking scripts in order to feed pointing programs for telescope and antenna pointing control.  General philosophy is to build a modular combination of scripts and execute them in particular order as opposed to building one monolithic script.

### Dependencies
All scripts are written for Python 3+.

```
sudo pip3 install astroquery skyfield matplotlib pandas
```
### Directories:
* data - location for data download.
* config - location for configuration files.
* utilities - helper utilities for things like plotting.

### Script Description:

