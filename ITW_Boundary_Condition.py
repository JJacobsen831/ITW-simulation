#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 14:49:58 2023

ITW Boundary Condition - BAROTROPIC

@author: jjacob2
"""
import os
os.chdir('/home/jjacob2/python/island_wave/subroutines/')
import numpy as np
from netCDF4 import Dataset as nc4

#directory paths
root = '/home/jjacob2/runs/island_wave/rad_itw04/'
root_bry = root+'ncfiles/roms_island_LRG_01hrDT_WESN_bry.nc'
root_ini = root+'initial/roms_island_AnaTcline_LRG_ini.nc'
root_grd = root + 'island_spongeNES_LRGgrid.nc'

#time parameters
##NTime, DT
dur = (10)*86400 #n time steps
dt = 60*60 #seconds
ntime = int(np.ceil(dur/dt))
time = np.arange(0,ntime*dt,dt)
omega = (2*np.pi)/(24*3600)

#open bry file to be written
ncbry = nc4(root_bry, mode='w',format='NETCDF4_CLASSIC')

#open initial condition
ncini = nc4(root_ini, 'r')

#open grid file
grid = nc4(root_grd,'r')

#set numer of vertical levels
n_lev = 200
dims = (1, n_lev+1)

#grid dims
x_rho0 = grid.variables['x_rho'][:]
x_v0 = grid.variables['x_v'][:]
x_psi0 = grid.variables['x_psi'][:]
x_u0 = grid.variables['x_u'][:]

y_rho0 = grid.variables['y_rho'][:]
y_u0 =  grid.variables['y_u'][:]
y_psi0 = grid.variables['y_psi'][:]
y_v0 = grid.variables['y_v'][:]


#Define values 
# East West
salt = ncini.variables['salt'][0,:,:,0]
temp = ncini.variables['temp'][0,:,:,0]

#  repeat along western/eastern boundary
salt_we = np.repeat(salt[np.newaxis,:,:], repeats = ntime, axis = 0)
temp_we = np.repeat(temp[np.newaxis,:,:], repeats = ntime, axis = 0)
v_we = np.zeros((ntime, n_lev, x_v0.shape[0]))
vbar_we = np.array(v_we[:,0,:])
u_we = np.zeros((ntime, n_lev, x_u0.shape[0]))
ubar_we = np.array(u_we[:,0,:])

# North South
saltsn = ncini.variables['salt'][0,:,0,:]
tempsn = ncini.variables['temp'][0,:,0,:]

#  repeat along northern/southern boundaries
salt_sn = np.repeat(saltsn[np.newaxis,:,:], repeats = ntime, axis = 0)
temp_sn = np.repeat(tempsn[np.newaxis,:,:], repeats = ntime, axis = 0)
v_sn = np.zeros((ntime, n_lev, y_v0.shape[1]))
vbar_sn = np.array(v_sn[:,0,:])
u_sn = np.zeros((ntime, n_lev, y_u0.shape[1]))
ubar_sn = np.array(u_sn[:,0,:])

#Define time dependent with linear ramp
#Weighting to linearly ramp-up forcing over 17 hours
t_lim = 17*3600
ramp = np.empty(time.shape[0])
for i in range(time.shape[0]) :
    if (time[i] <= t_lim) :
        ramp[i] = 1/(t_lim)*time[i]
    else :
        ramp[i] = 1


zeta = 1.25*ramp*np.sin(omega*time)
zeta_we = np.repeat(zeta[:,np.newaxis], repeats = salt.shape[1], axis = 1)


#NetCDF writing
##Define dimensions in nc file
bry_time = ncbry.createDimension('bry_time',time.shape[0])

eta_rho = ncbry.createDimension('eta_rho',x_rho0.shape[0])
xi_rho = ncbry.createDimension('xi_rho', x_rho0.shape[1])

eta_psi = ncbry.createDimension('eta_psi',x_psi0.shape[0])
xi_psi = ncbry.createDimension('xi_psi', x_psi0.shape[1])

eta_v = ncbry.createDimension('eta_v', x_v0.shape[0])
xi_v = ncbry.createDimension('xi_v', x_v0.shape[1])

eta_u = ncbry.createDimension('eta_u',x_u0.shape[0])
xi_u = ncbry.createDimension('xi_u', x_u0.shape[1])

s_rho = ncbry.createDimension('s_rho',dims[1]-1)
s_w = ncbry.createDimension('s_w',dims[1])


##Define variables
btime = ncbry.createVariable('bry_time','f8',("bry_time"))

temp_west = ncbry.createVariable('temp_west','f8',
                                    ("bry_time", "s_rho", "eta_rho"))
temp_east = ncbry.createVariable('temp_east','f8',
                                    ("bry_time", "s_rho", "eta_rho"))
temp_south = ncbry.createVariable('temp_south','f8',
                                    ("bry_time", "s_rho", "xi_rho"))
temp_north = ncbry.createVariable('temp_north','f8',
                                    ("bry_time", "s_rho", "xi_rho"))

salt_west = ncbry.createVariable('salt_west','f8',
                                    ("bry_time", "s_rho", "eta_rho"))
salt_east = ncbry.createVariable('salt_east','f8',
                                    ("bry_time", "s_rho", "eta_rho"))
salt_south = ncbry.createVariable('salt_south','f8',
                                    ("bry_time", "s_rho", "xi_rho"))
salt_north = ncbry.createVariable('salt_north','f8',
                                    ("bry_time", "s_rho", "xi_rho"))


u_west = ncbry.createVariable('u_west','f8',("bry_time", "s_rho", "eta_u"))
u_east = ncbry.createVariable('u_east','f8',("bry_time", "s_rho", "eta_u"))
u_south = ncbry.createVariable('u_south','f8',("bry_time", "s_rho", "xi_u"))
u_north = ncbry.createVariable('u_north','f8',("bry_time", "s_rho", "xi_u"))

ubar_west = ncbry.createVariable('ubar_west','f8',("bry_time", "eta_u"))
ubar_east = ncbry.createVariable('ubar_east','f8',("bry_time", "eta_u"))
ubar_south = ncbry.createVariable('ubar_south','f8',("bry_time", "xi_u"))
ubar_north = ncbry.createVariable('ubar_north','f8',("bry_time", "xi_u"))

v_west = ncbry.createVariable('v_west','f8',("bry_time", "s_rho", "eta_v"))
v_east = ncbry.createVariable('v_east','f8',("bry_time", "s_rho", "eta_v"))
v_south = ncbry.createVariable('v_south','f8',("bry_time", "s_rho", "xi_v"))
v_north = ncbry.createVariable('v_north','f8',("bry_time", "s_rho", "xi_v"))

vbar_west = ncbry.createVariable('vbar_west','f8',("bry_time", "eta_v"))
vbar_east = ncbry.createVariable('vbar_east','f8',("bry_time", "eta_v"))
vbar_south = ncbry.createVariable('vbar_south','f8',("bry_time", "xi_v"))
vbar_north = ncbry.createVariable('vbar_north','f8',("bry_time", "xi_v"))

zeta_west = ncbry.createVariable('zeta_west','f8',("bry_time", "eta_rho"))
zeta_east = ncbry.createVariable('zeta_east','f8',("bry_time", "eta_rho"))
zeta_south = ncbry.createVariable('zeta_south','f8',("bry_time", "xi_rho"))
zeta_north = ncbry.createVariable('zeta_north','f8',("bry_time", "xi_rho"))

##Write variables
btime[:] = time
btime.field = 'ocean_time, scalar, series'
btime.long_name = 'boundary time'
btime.units = 'seconds since 0001-01-01 00:00:00'
btime.time_origin = '01-JAN-1900-01-01 00:00:00'

temp_west[:] = temp_we
temp_east[:] = temp_we
temp_south[:] = temp_sn
temp_north[:] = temp_sn

salt_west[:] = salt_we
salt_east[:] = salt_we
salt_south[:] = salt_sn
salt_north[:] = salt_sn

u_west[:] = u_we
u_east[:] = u_we
u_south[:] = u_sn
u_north[:] = u_sn

ubar_west[:] = ubar_we
ubar_east[:] = ubar_we
ubar_south[:] = ubar_sn
ubar_north[:] = ubar_sn

v_west[:] = v_we
v_east[:] = v_we
v_south[:] = v_sn
v_north[:] = v_sn

vbar_west[:] = vbar_we
vbar_east[:] = vbar_we
vbar_south[:] = vbar_sn
vbar_north[:] = vbar_sn

zeta_west[:] = zeta_we
zeta_east[:] = np.zeros(zeta_we.shape)
zeta_south[:] = np.zeros(temp_sn[:,0,:].shape)
zeta_north[:] = np.zeros(temp_sn[:,0,:].shape)

##close NetCDF bry file
ncbry.close()

