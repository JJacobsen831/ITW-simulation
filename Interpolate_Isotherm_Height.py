#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 15:41:02 2023

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/island_wave/subroutines/')
import sys
sys.path.append('/home/jjacob2/python/island_wave/subroutines/')
import numpy as np
from netCDF4 import Dataset as nc4
import obs_depth_vs145 as dep


#file directory
root = '/home/jjacob2/runs/island_wave/lrgAmp/'
fname = 'output02/roms_his.nc'

#new file
dat_root = '/home/jjacob2/python/island_wave/data/'
newfile = 'TempHeight_lrgAmp2_5.nc'

#load file
ncfile = nc4(root+fname)
grid = nc4(root+'island_spongeNES_grid.nc')

#slicing
idx = {'lat' : slice(20,101),
       'lon' : slice(75,171)}

#grid
time = ncfile.variables['ocean_time'][:]/3600
xrho = grid.variables['x_rho'][idx['lat'],idx['lon']]
yrho = grid.variables['y_rho'][idx['lat'],idx['lon']]

#create file
f = nc4(dat_root+newfile, 'w', format = 'NETCDF4')

#dimensions
f.createDimension('x_rho', xrho.shape[1])
f.createDimension('y_rho', yrho.shape[0])
f.createDimension('time', time.shape[0])

#variables
xrho_ = f.createVariable('x_rho', 'f4', ('y_rho', 'x_rho'))
yrho_ = f.createVariable('y_rho', 'f4', ('y_rho', 'x_rho'))
time_ = f.createVariable('time', 'f4', 'time')
height_ = f.createVariable('t_height', 'f4', ('time', 'y_rho', 'x_rho'))

#depth
romsvars = {'Vstretching' : ncfile.variables['Vstretching'][0], \
            'Vtransform' : ncfile.variables['Vtransform'][0], \
            'theta_s' : ncfile.variables['theta_s'][0], \
            'theta_b' : ncfile.variables['theta_b'][0], \
            'N' : 200, \
            'h' : grid.variables['h'][:], \
            'hc': ncfile.variables['hc'][0]}
                
depth = dep._set_depth(ncfile, romsvars, 'rho', romsvars['h'])[:,idx['lat'],idx['lon']]

temp = ncfile.variables['temp'][:, :, idx['lat'],idx['lon']]

#height of 18 degree isotherm
t_height = np.empty((time.shape[0], depth.shape[1], depth.shape[2]))
for i in range(80) :
    for j in range(depth.shape[2]) :
        for t in range(time.shape[0]) :
            t_height[t,i,j] = np.interp(18, temp[t,:,i,j], depth[:,i,j])


#store data
xrho_[:] = xrho
yrho_[:] = yrho
time_[:] = time
height_[:] = t_height

f.close()

