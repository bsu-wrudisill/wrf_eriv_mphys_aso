#!/usr/bin/env python
# coding: utf-8
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pathlib as pl
import sys

### Purpose -- patch the NoahMP grid with updated fields
### The model doesn't actually use these, but we 
### will use this file to adjust the forcings later 
### based on the updated DEM, aspect, and slope, etc

if len(sys.argv) < 6:
        print("order: geo_file, dem_layer, slope_layer, aspect_layer, destination")
        sys.exit()

geo_file = pl.Path(sys.argv[1])
dem_layer = pl.Path(sys.argv[2])
slope_layer = pl.Path(sys.argv[3])
aspect_layer = pl.Path(sys.argv[4])
dest_geo_file_path = pl.Path(sys.argv[5])

# open the dem
dem = xr.open_dataset(dem_layer)
dem = dem.drop("band")
dem = dem.rename(dict(__xarray_dataarray_variable__ = "HGT_M"))

# open up the aspect ...
aspect = xr.open_dataset(aspect_layer)
aspect = aspect.drop("band")

# open up the slope 
slope = xr.open_dataset(slope_layer)
slope = slope.drop("band")


# open the geofile
geo_file = xr.open_dataset(geo_file)

# get the patched dem... there mighte be nans in the high res layer
patched_DEM = np.where(dem.HGT_M <= 0, geo_file.HGT_M.isel(Time=0), dem.HGT_M)

# take the difference 
orig_topo = geo_file.HGT_M.isel(Time=0).copy()
hgt_diff = patched_DEM.copy() - orig_topo 


# update the values ...
geo_file.HGT_M.values[0,:,:] = patched_DEM
geo_file["DEMDIFF"] = geo_file.LU_INDEX    # this is just to make a var with the right dimensions...
geo_file["DEMDIFF"].values[0,:,:] = hgt_diff


# update the slope and aspect data...
geo_file["slope"] = geo_file.LU_INDEX
geo_file["slope"].values[0,:,:] = slope["__xarray_dataarray_variable__"]

geo_file["aspect"] = geo_file.LU_INDEX
geo_file["aspect"].values[0,:,:] = aspect["__xarray_dataarray_variable__"]

# write the file...
geo_file.to_netcdf(dest_geo_file_path.joinpath("geo_em_updatedDEM.d01.nc"))

