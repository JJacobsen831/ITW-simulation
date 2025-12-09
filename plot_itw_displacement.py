#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 10:29:49 2023

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/island_wave/subroutines/')
import sys
sys.path.append('/home/jjacob2/python/island_wave/subroutines/')

import numpy as np
from netCDF4 import Dataset as nc4
import matplotlib.pyplot as plt

#interpolated file 
dataroot = '/home/jjacob2/python/island_wave/data/'
fname = 'BoundaryWave20degHeight.nc'

gridfile = '/home/jjacob2/runs/island_wave/rad_itw03/island_sponge_grid.nc'

#load files
ncfile = nc4(dataroot+fname)
grid = nc4(gridfile)
mask = np.array(grid.variables['mask_rho'][slice(20,101),slice(75,171)], 
                dtype = bool)


#grid
time = ncfile.variables['time'][200:]
xrho = ncfile.variables['x_rho'][:]/1000
yrho = ncfile.variables['y_rho'][:]/1000

#18 degree height
temp = ncfile.variables['t_height'][200:,:,:]

#range of height
hrange = np.ma.array(np.max(temp, axis = 0) - np.min(temp, axis = 0), mask = ~mask)


fig, (ax) = plt.subplots(1,1)
cf = ax.contourf(xrho, yrho, np.mean(temp, axis = 0))
ax.set_title('Average depth of 20 isotherm')
ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_xlim([50,70])
ax.set_ylim([20,40])
fig.gca().set_aspect('equal')
fig.colorbar(cf)

i = 30
j = 40
fig, (ax) = plt.subplots(1,1)
cf = ax.contourf(xrho, yrho, hrange)
ax.patch.set_facecolor('silver')
ax.scatter(xrho[i,j], yrho[i,j], c = 'white', marker = 'X')
ax.set_title('Range of 20 isotherm depth')
ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_xlim([50,70])
ax.set_ylim([20,40])
fig.gca().set_aspect('equal')
fig.colorbar(cf, label = '[m]')

fig, ax1 = plt.subplots(1,1)
ax1.plot(time, temp[:,i,j])
ax1.grid()
ax1.set_xlabel('Time [hr]')
ax1.set_ylabel('Depth of 20 Isotherm [m]')
