#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:29:59 2023

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/island_wave/subroutines/')
import sys
sys.path.append('/home/jjacob2/python/island_wave/subroutines/')

import numpy as np
from netCDF4 import Dataset as nc4
import matplotlib.pyplot as plt
import cmocean 

root = '/home/jjacob2/runs/island_wave/diurnal_npzd_itw01/'
flfile = nc4(root+'output02/roms_flt.nc')
grid = nc4(root+'island_spongeNES_grid.nc')

#grid
xrho = grid.variables['x_rho'][:]/1000 - 60
yrho = grid.variables['y_rho'][:]/1000 - 30
mask = np.array(grid.variables['mask_rho'][:], 
                dtype = bool)
h = np.ma.array(grid.variables['h'][:], mask = ~mask)

#lagrangian positions west spoke
xpos = flfile.variables['x'][11521:,:]/1000 - 60 
ypos = flfile.variables['y'][11521:,:]/1000- 30
zpos = flfile.variables['depth'][11521:,:]

#floats with initial depth at 20 m
i_west = slice(76,80)
i_north = slice(300,304)

fig, (ax1, ax2) = plt.subplots(1,2, sharey = True)
ax1.contourf(xrho, yrho,h)
ax1.plot(xpos[:,i_west], ypos[:,i_west])
ax1.set_xlim([-10,10])
ax1.set_ylim([-10,10])
ax1.set_aspect('equal', adjustable = 'box')
ax1.patch.set_facecolor('silver')
ax1.set_ylabel('y [km]')
ax1.set_xlabel('x [km]')

ax2.contourf(xrho, yrho,h)
ax2.plot(xpos[:,i_north], ypos[:,i_north])
ax2.set_xlim([-10,10])
ax2.set_ylim([-10,10])
ax2.set_aspect('equal', adjustable = 'box')
ax2.patch.set_facecolor('silver')
ax2.set_xlabel('x [km]')


xwest = [-7.5,-2]
fig, (ax,ax1) = plt.subplots(2,1, sharex = True, figsize=(4,6))
ax.plot(xpos[:,i_west], zpos[:,i_west])
ax.set_xlim(xwest)
ax.set_ylabel('depth [m]')

ax1.contourf(xrho, yrho,h)
ax1.plot(xpos[:,i_west], ypos[:,i_west])
ax1.set_xlim(xwest)
ax1.set_ylim([-6,6])
#ax1.set_aspect('equal', adjustable = 'box')
ax1.patch.set_facecolor('silver')
ax1.set_ylabel('y [km]')
ax1.set_xlabel('x [km]')


xnorth = [-1,11]
fig, (ax0,ax2) = plt.subplots(2,1, sharex = True, figsize=(4,6))
ax0.plot(xpos[:,i_north], zpos[:,i_north])
ax0.set_xlim(xnorth)

ax2.contourf(xrho, yrho,h)
ax2.plot(xpos[:,i_north], ypos[:,i_north])
ax2.set_xlim(xnorth)
ax2.set_ylim([-1,8])
#ax2.set_aspect('equal', adjustable = 'box')
ax2.patch.set_facecolor('silver')
ax2.set_xlabel('x [km]')

