#!/usr/bin/env python
# coding: utf-8
import xarray as xr
import xesmf as xe
import numpy as np 
import pathlib as pl
import pandas as pd
import matplotlib.pyplot as plt
import sys 


### regrid ASO data (that is in lat lon format) to the NoahMP/WRF dem 
### so that comparisions can be made 


# this is the target file -- must be provided
wrfinfile = pl.Path(sys.argv[1])
destfolder = pl.Path(sys.argv[2])

on=False 

# get the location of the processsed dems 
if on==True:
    aso_dat = pl.Path("./ASOdata/processed_dems")
    aso_swe_dat = pl.Path("./ASOdata/raw_swe_data")

    # open up the dem products 
    dem=xr.open_rasterio(aso_dat.joinpath("wgs84_latlon_elevation_50m.tif"))
    dem = dem.rename({'x': 'longitude', 'y': 'latitude'})

    aspect=xr.open_rasterio(aso_dat.joinpath("wgs84_latlon_aspect_50m.tif"))
    aspect = aspect.rename({'x': 'longitude', 'y': 'latitude'})

    slope=xr.open_rasterio(aso_dat.joinpath("wgs84_latlon_slope_50m.tif"))
    slope = slope.rename({'x': 'longitude', 'y': 'latitude'})


    # open up the wrfinput file -- target grid  
    target_wrfinput_file = xr.open_dataset(wrfinfile)
    target_wrfinput_file = target_wrfinput_file.rename({'XLONG': 'lon', 'XLAT': 'lat'})
    target_wrfinput_file['lat'] = target_wrfinput_file['lat'][0,:,:]
    target_wrfinput_file['lon'] = target_wrfinput_file['lon'][0,:,:]



    # do the regirdding 
    regridder = xe.Regridder(dem, target_wrfinput_file, 'bilinear')

    # regrid the dem/slope/aspect 
    dem_regrid = regridder(dem)
    dem_regrid = dem_regrid.isel(band=0)

    aspect_regrid = regridder(aspect)
    aspect_regrid = aspect_regrid.isel(band=0)

    slope_regrid = regridder(slope)
    slope_regrid = slope_regrid.isel(band=0)


    # write everything out 
    dem_regrid.to_netcdf(destfolder.joinpath("ASO_elevation_regridded.nc"))
    aspect_regrid.to_netcdf(destfolder.joinpath("ASO_aspect_regridded.nc"))
    slope_regrid.to_netcdf(destfolder.joinpath("ASO_slope_regridded.nc"))



target_wrfinput_file = xr.open_dataset(wrfinfile)
target_wrfinput_file = target_wrfinput_file.rename({'XLONG': 'lon', 'XLAT': 'lat'})
target_wrfinput_file['lat'] = target_wrfinput_file['lat'][0,:,:]
target_wrfinput_file['lon'] = target_wrfinput_file['lon'][0,:,:]


### NOW DO THE SAME FOR THE ASO SWE DATA:
for swedat in pl.Path("/home/wrudisill/scratch/EastLSM_Only/most_recent_aso_data/swedata/og_grid").glob("ASO_250M_SWE_bilin*.nc"):
        print(swedat)
        swe=xr.open_dataset(swedat)
        #swe = swe.rename({'Band1':'swe', 'lon': 'longitude', 'lat': 'latitude'})
        #swe = swe.drop(['time', 'band'])
        swe = swe.rename({"lat":"latitude", "lon":"longitude"})
        regridder = xe.Regridder(swe, target_wrfinput_file, 'bilinear')
        swe_regrid = regridder(swe.Band1)
        swe_regrid.to_netcdf(destfolder.joinpath("swe_regrid_250m_" + swedat.name))

        # done 


for swedat in pl.Path("/home/wrudisill/scratch/EastLSM_Only/most_recent_aso_data/snowdepths/og_grid").glob("ASO_SD_250m*.nc"):
        swe=xr.open_dataset(swedat)
        #swe = swe.rename({'__xarray_dataarray_variable__':'swe', 'lon': 'longitude', 'lat': 'latitude'})
        #swe = swe.drop(['time', 'band'])
        swe = swe.rename({'lon': 'longitude', 'lat': 'latitude'})
        regridder = xe.Regridder(swe, target_wrfinput_file, 'bilinear')
        swe_regrid = regridder(swe.Band1)
        swe_regrid.to_netcdf(destfolder.joinpath("sd_regrid_250m_" + swedat.name))

        # done 
