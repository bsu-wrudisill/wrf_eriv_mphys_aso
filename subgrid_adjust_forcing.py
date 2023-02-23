import xarray as xr
import numpy as np 
import pathlib as pl
import sys
from compute_srad_dingman import wrf_irrad_correct 
from datetime import timedelta
from metpy.calc import relative_humidity_from_mixing_ratio, mixing_ratio_from_relative_humidity
import metpy.calc as mpcalc
from metpy.units import units


def rad(xarr):
        return xarr*np.pi/180.


def adjust_precip(rainvar, dz, opg):
        # rainvar is the xarraay array 
        # that has the precip data in it
        # opg is the parameter to adjust the precipiitation...

        ## ADJUST THE PRECIP AND TEMPERATURE 
        updated_rainvar  = rainvar + rainvar*dz*opg

        # ADJUST THE TEMPERATURE 
        return updated_rainvar


def adjust_temperature(tempvar, dz, lr=-.0065):
    # tempvar is the temperature xarray 
    # dz is the array (same size of tempvar) of elevation differences 
    updated_tempvar = tempvar + dz*.0065
    return updated_tempvar 

#def adjust_pressure(pvar, dz):


def adjust_pressure(dz, pressure):
    pressure_adj = pressure*np.exp(-1*dz/8000) # from micromet paper; 8000m is atm. reference depth
    return pressure_adj 

def adjust_Q2_T2_PSFC(tempvar, pressure, q2, dz):
    # assume that the relative humidity from the original data 
    # is preserved everywhere... and that the 
    # Q2 gets adjusted to match the downscaled temperatures 
    rh = relative_humidity_from_mixing_ratio(pressure * units("Pa"), 
                                             tempvar * units("kelvin"), 
                                             q2* units("kg/kg") )
    
    # compute the new pressure 
    p_adj = adjust_pressure(dz, pressure)
    t_adj = adjust_temperature(tempvar, dz)
    
    # now compute the new q from the rh and adjusted temperature and pressure...
    q_adj = mixing_ratio_from_relative_humidity(p_adj * units("Pa"), t_adj * units("kelvin"),  rh)
    return q_adj, t_adj, p_adj

            
def adjust_swrad(swdown, lat, slope, aspect, date, tz=-7):
        # adjust the slope-normal radiation
        # this is mirroring the way that swnorm is calculated 
        # in WRF... ! NEEDS TO BE VERIFIED !
        # the SWDOWN is multiplied by the equivalent surface 
        # sun-zenith angle (takes into account the aspect 
        # and the slope)

        # INPUTS
        #swdown is xarray array
        #lat -- grids
        #lon  -- grids
        #slope is the same
        #azimuth is the same
        #date is a datetime object .... ASSUMES THAT 12:00 is TRUE SOLAR NOON
        
        # adjust theh date... the one that's listed is in UTC
        # need to adjust to local time.
        date = date + timedelta(hours=tz) 
        
        #cosz,scec = irrad(lat, lon, slope, aspect, date)
        #subgrid_swnorm = cosz*swdown
        subgrid_swnorm = wrf_irrad_correct(lat, slope, aspect, date, swdown)
        return subgrid_swnorm.where(subgrid_swnorm>0., 0.)
