#!/usr/bin/env python
# coding: utf-8

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from math import radians as rad
from math import degrees as deg

from datetime import timedelta
import xarray as xr


## Compute solar geometry..
# https://gml.noaa.gov/grad/solcalc/solareqns.PDF


## THIS PART IS FOR CALCULATING SOLAR NOON --- WE ULTIAMTELY JUST WANT THE #OF HOURS 
## THAT ARE +/- SOLAR NOON FOR A GIVEN DATETIME ... FOR NOW WE JUST ASSUME THAT 
## 12 IS THE SOLAR NOON
def eqtime(da):
    # Literally called "equation of time"   -- Minutes 
    #---------------------------------------
    ### equation of time - an astronomical term accounting for changes in the time of solar noon 
    ### for a given location over the course of a year. Earth's elliptical orbit and Kepler's law of 
    ### equal areas in equal times are the culprits behind this phenomenon. 
    ### Click here to see a plot of the equation of time vs. day of the year. 
    ### For more information on this phenomenon, see this offsite Analemma page.
    ### https://gml.noaa.gov/grad/solcalc/glossary.html#equationoftime
    #---------------------------------------
    eqt = 229.18*(0.000075 + 0.001868*np.cos(da) - 0.032077*np.sin(da) - 0.014615*np.cos(2*da) - 0.040849*np.sin(2*da))
    return eqt



def day_angle(doy):
    # doy is julian day
    return 2 * np.pi * (doy-1)/365  # RADIANS

def eccentricity(da):
    # da is the day_angle (radians)
    e = 1.000110 + 0.034221 * np.cos(da) + .001280*np.sin(da) + .000719*np.cos(2*da) + .000077*np.sin(2*da)
    return e

def declination(da):
    # da is the day angle
    d = .006918 - 0.399912*np.cos(da) + 0.070257*np.sin(da) - 0.006758*np.cos(2*da) # RADIANS
    return d

def zenith_angle(lat, dec, wt):
    # lat is latitude (radians)
    # dec is teh declination (radians)
    # w is the angular velocity of the earths rotation
    # t is the number of hours before/after solar noon
    
    # w is always the same...

    # w = .2618 # rad/hr == 15deg/hr

    # compute first part...
    cza = np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(wt)
    
    # take the inverse 
    return cza

def tr_ts(lat, dec):
    # return the sunrise and sunset times...
    w = .2618 # this is constant... earths rotation in radians
    Ts = np.arccos(-np.tan(dec)*np.tan(lat))/w  # time sunset
    Tr = -Ts # always symmetrical ...           # time sunrise
    return Ts, Tr


def equiv_hz_sfc(ks, h, lat):
    # this is the latitude of the "equivalent horizontal surface"
    # ks is the slope inclination
    # h is the slope azimuth
    # lat is the latitude
    return np.arcsin(np.sin(ks)*np.cos(h)*np.cos(lat) + np.cos(ks)*np.sin(lat))

def a_long_diff(ks,h,lat):
    # this is for the slope-normal terrestrial radiation 
    # difference in longitude between equivalent horizontal surface and slope
    return np.arctan(np.sin(h)*np.sin(ks)/(np.cos(ks)*np.cos(lat)-np.cos(h)*np.sin(ks)*np.sin(lat)))


def cza(lat,
        date):

    #lat:  lat in radians
    #long: long in radians
    #ks:   slope inclination in radians
    #h:    slope azimuth, radians clockwise from North (0 is North)
    #date: date (timesampt object)
    #t   : hour relative to solar noon (this will be updated later)
    
    # get the hour
    w = .2618           # rad/hr == 15deg/hr

    hour = date.hour # this is 0-24
    t = hour - 12
    
    # begin geometru calculations
    da = day_angle(date.dayofyear)
    ec = eccentricity(da)
    dec = declination(da)    
   
    # compute the zentith angle
    cza = zenith_angle(lat, dec, w*t)
    return cza



## Plot the irradiance
def cza_slp(lat,
            ks,
            h,
            date):

    #lat:  lat in radians
    #ks:   slope inclination in radians
    #h:    slope azimuth, radians clockwise from North (0 is North)
    #date: date (timesampt object)
    #t   : hour relative to solar noon (this will be updated later)
    
    # get the hour
    w = .2618           # rad/hr == 15deg/hr

    hour = date.hour # this is 0-24
    t = hour - 12
    
    # begin geometru calculations
    da = day_angle(date.dayofyear)
    ec = eccentricity(da)
    dec = declination(da)

    # do the adjustment for slopes
    a = equiv_hz_sfc(ks, h, lat)
    wta = w*t+a
    
    # now find the equiv latitude
    lat0 = np.arcsin(np.sin(ks)*np.cos(h)*np.cos(lat) + np.cos(ks)*np.sin(lat))
    
    # compute the zentith angle
    cza = zenith_angle(lat0, dec, wta)
    
    return cza


### WRF SPECIFIC RADIATION EQUATIONS ###
# SEE: https://github.com/wrf-model/WRF/blob/eed56d74b865af4ce9f98ea029226acfe52b0569/phys/module_ra_sw.F

# This computed the diffuse component of solar radiation, as implenented in WRF
# it is based on the ratio of the TOA atmospheric radiation and the
# radiation at the bottom atmospheric layer 
def diffuse_frac(swdown, date):
        # compute the diffuse fraction of radiation
        # this is based on what they do in WRF....
        # see github code.
        
        # FORTRAN CODE
        solcon =  1364.0 # this should probably be mult. by the eccentricity?
        da = day_angle(date.dayofyear)
        ec = eccentricity(da)
        solcon = solcon*ec #solar constant times eccentricity

        # in the WRF code, solcon is replaced by the SWDOWN at the highest model level
        # which is probably really close to the solar constant(times eccentricity...)

        diffuse_frac = np.minimum(1., 1./np.maximum(.1, 2.1 - 2.8*np.log(np.log(solcon/np.maximum(swdown, 1E-3)))))
        return diffuse_frac
    

# correction factor for SWDOWN .... 
def corr_fac(d, cza_slp0, cza0):
    # d is diffuse fract
    # cza_slp0 is the coszine zenith angle (rad) of sloping topo
    # cza0 is the cosine zenith angle of flat topography 
    corr_fac = d + (1-d)*cza_slp0/cza0
    
    ### apply some rules... fortran code:       ###
    # if ((slope.eq.0).or.(diffuse_frac.eq.1).or.(csza.lt.1.e-2)) 
    # then  ! no topographic effects when all radiation is diffuse or the sun is too close to the horizon
    ###
    
    corr_fac = corr_fac.where(d != 1, 1)
    corr_fac = corr_fac.where(cza0 > 1E-2, 1)
    return corr_fac


# make a function taht does it all 
def wrf_irrad_correct(lat, slope, aspect, date, swdown):
    # lat (rad)
    # lon (rad)
    # slope (rad)
    # azimith (rad)
    # date (pd.datetime)
    # swdown array
    
    # compute the diffuse fraction of swdown...
    d = diffuse_frac(swdown, date)
    
    # compute the flat and terrain-normal zenith angles 
    cza_slp0 = cza_slp(lat, slope, aspect, date)
    cza0 = cza(lat, date)
    
    # compute the corr fac..
    cf = corr_fac(d, cza_slp0, cza0)
    
    # now compute the radaition and return 
    return cf*swdown




### Equivalent slope stuff...
### This is so that we can compute the 

if __name__ == "__main__":
        ##example... 
        days_in_year = 365
        timezone = -7 # relative to GMT
        lat = rad(55)
        long = rad(-116.5)
        slope = rad(30)    # 30 Degree slope..
        azimuth = rad(145) # Clockwise from North(0)

        # check answer with book...
        #deg(equiv_hz_sfc(slope, azimuth, rad(lat)))

        # now do a...
        #deg(a_long_diff(slope, azimuth, rad(lat)))


        # it is done this way so that we can have more flexibility 
        # when doing the subgrid adjustment using the WRF fields 
        cza, scec  = irrad(lat, long, slope, azimuth, dr[10])
        ir = cza*scec


        aspect =  xr.open_dataset("500m_ModelRuns/aspect_regrid_500m.nc")
        slope = xr.open_dataset("500m_ModelRuns/slope_regrid_500m.nc")
        slope = slope.rename(dict(__xarray_dataarray_variable__="slope"))
        aspect = aspect.rename(dict(__xarray_dataarray_variable__="aspect"))


        dr = pd.date_range("2017-01-01 12:00", "2018-01-01 12:00", freq="h")

        # look at some different times...
        d1e = pd.to_datetime("2018-03-21 03:00")
        d2e = pd.to_datetime("2018-03-21 12:00")
        d3e = pd.to_datetime("2018-03-21 17:00")



        lat = slope.lat*np.pi/180
        lon = slope.lon*np.pi/180
        slope_rad = slope.slope*np.pi/180
        aspect_rad = aspect.aspect*np.pi/180


        #swnorm = [irrad(lat, lon, slope_rad, aspect_rad, d) for d in dr]
        may31_6am_rad = irrad(lat, lon, slope_rad, aspect_rad, d1e)
        may31_12am_rad = irrad(lat, lon, slope_rad, aspect_rad, d2e)
        may31_6pm_rad = irrad(lat, lon, slope_rad, aspect_rad, d3e)

