#!/usr/bin/env python
# coding: utf-8
import xarray as xr
import xesmf as xe
import numpy as np 
import pathlib as pl
import pandas as pd
import sys
from subgrid_adjust_forcing import *


#if (len(sys.argv) < 4):
#        print("Args: geog_file, wrfin_file, output_dir")
#        sys.exit()


subfile_path = pl.Path(sys.argv[1])
output_dir = pl.Path(sys.argv[2]) #folder location for the focing files

# open up the geog and wrfinput files

# open up the static data...
target_geog_file_loc = "/home/wrudisill/scratch/EastLSM_Only/most_recent_aso_data/geo_em_updatedDEM.d01.nc" 
target_wrfinput_file_loc = "/scratch/wrudisill/EastLSM_Only/most_recent_aso_data/wrfinput_d01.nc"

target_geog_file = xr.open_dataset(target_geog_file_loc)#"/home/wrudisill/scratch/EastLSM_Only/geog/geo_em.d01.nc")

veg_data = "/scratch/wrudisill/EastLSM_Only/most_recent_aso_data/daily_vegfrac.nc"
veg_data = xr.open_dataset(veg_data)


# open the target (wrfin) file and change the varnames so that we can regrid to it...
target_wrfinput_file = xr.open_dataset(target_wrfinput_file_loc)#"/home/wrudisill/scratch/EastLSM_Only/Run/DOMAIN/wrfinput_d01.nc")
target_wrfinput_file = target_wrfinput_file.rename({'XLONG': 'lon', 'XLAT': 'lat'})
target_wrfinput_file['lat'] = target_wrfinput_file['lat'][0,:,:]
target_wrfinput_file['lon'] = target_wrfinput_file['lon'][0,:,:]


# open up the subfiles....

#subfile_path = pl.Path("/scratch/wrudisill/EastLSM_Only/WRF_subset_files/Ishmael/WY2019").glob("*wsub*.nc")
#subfile_path = pl.Path("/scratch/wrudisill/EastLSM_Only/WRF_subset_files/Morrison/WY2019").glob("*wsub*.nc")
#subfile_path = pl.Path("/scratch/wrudisill/EastLSM_Only/WRF_subset_files/Thompson/WY2019").glob("*wsub*.nc")

subfile_path.mkdir(exist_ok=True)


# open up the geo_em file so that we can get the vegetation fraction...



# loop through them 
for subfile in subfile_path.glob("*wsub*.nc"):
#for subfile in [subfile_path.joinpath("Month12_SfcMet_wsub_WY2018.nc")]:
        #ds0 = xr.open_dataset("Month02_SfcMet_wsub_WY2018.nc")
        print(subfile)
        ds0 = xr.open_dataset(subfile)
        #ds = ds0.drop(["EAST_MASK","TAYLOR_MASK", "ACSWUPT", "LWUPT", "SWUPT", "HFX", "LH", "HR_PRCP", "SWDOWN"])
        ds = ds0.drop(["EAST_MASK","TAYLOR_MASK", "HFX", "LH", "HR_PRCP", "SWDOWN"])  # these are not in the "morrison" files
        ds = ds.rename({'XLONG': 'lon', 'XLAT': 'lat', 'SWNORM':'SWDOWN'})

        ## DO SOME REGRIDDING
        # do the regridding here 
        regridder = xe.Regridder(ds, target_wrfinput_file, 'bilinear')

        # just get the files that we want...
        varlist = [varname for varname in ds.keys() if varname != 'Times']  # huh... this is acceptable syntax
        newds = xr.Dataset(data_vars=None, attrs=ds.attrs)

        # go through the varlist...
        for var in varlist:
            var_regrid = regridder(ds[var])
            newds[var] = (['Time', 'south_north','west_east'], np.zeros_like(var_regrid))     
            newds[var] = var_regrid 

        # make temporary vegfrac var...
        newds["RAINC"] = (['XTIME', 'south_north','west_east'], np.zeros_like(var_regrid))
        newds["VEGFRA"] = (['XTIME', 'south_north','west_east'], np.ones_like(var_regrid))

        # regrid the mask... theres only one time for this one
        rgd_msk = regridder(ds0.EAST_MASK)
        
        # 
        #del newds["RAINNC"]
        
        # rename the rain
        newds = newds.rename(dict(ACCPRCP="RAINNC", lat="XLAT", lon="XLONG"))
    

        ## Now make one file per timestep...
        #output_dir = pl.Path("/scratch/wrudisill/EastLSM_Only/FORCING/hourly_files")

        demdiff = target_geog_file["DEMDIFF"].isel(Time=0)
        
        # loop through the times 
        lat = rad(target_geog_file["XLAT_M"].isel(Time=0))
        lon = rad(target_geog_file["XLONG_M"].isel(Time=0))
        slope = rad(target_geog_file["slope"].isel(Time=0))
        aspect = rad(target_geog_file["aspect"].isel(Time=0))

        # loop through a bunch of files...
        for i,XT in enumerate(newds["XTIME"]):
            xt = pd.to_datetime(XT.values)
            dtfmt = xt.strftime("%Y-%m-%d_%H")
    
            # make the output name...
            newname = "wrfout_d01_%s:00:00"%(dtfmt)
            
            # get the vegetation data... the year is hardcoded as 2010
            vgdate = pd.to_datetime("2010-%s"%(xt.strftime("%m-%d")))
            vgdata = veg_data.sel(time=vgdate)


            # write out the files..
            out_name = output_dir.joinpath(newname)
            out_dataset= newds.isel(XTIME=i)
            out_dataset["EAST_MASK"] = rgd_msk
            
            
            ## ADJUST THE PRECIP AND TEMPERATURE 
            #opg = .000 # 2x per km
            #out_dataset["RAINNC"] = adjust_precip(out_dataset["RAINNC"], demdiff, opg)
            
            # make some adjustments ...
            q_adj, t_adj, p_adj = adjust_Q2_T2_PSFC(out_dataset["T2"],
                                                    out_dataset["PSFC"],
                                                    out_dataset["Q2"],
                                                    demdiff)

            # add the veg in there..
            out_dataset["VEGFRA"][:,:] = vgdata.GREENFRAC.values 

            # ADJUST THE TEMPERATURE 
            out_dataset["Q2"][:,:]    = q_adj.values  #adjust_temperature(out_dataset["T2"], demdiff)
            out_dataset["T2"][:,:]   = t_adj.values + 3.5 # this is the bias at the snotel locations ... #adjust_temperature(out_dataset["T2"], demdiff)
            out_dataset["PSFC"][:,:] = p_adj.values  #adjust_temperature(out_dataset["T2"], demdiff)

            # ADJUST THE RADIATION...
            out_dataset["SWDOWN"] = adjust_swrad(out_dataset["SWDOWN"], lat, slope, aspect, xt)


            out_dataset.to_netcdf(out_name)

