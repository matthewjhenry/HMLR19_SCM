import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from math import pi
import sigma_int

# plot meridional energy transports
# takes path to netcdf file as argument

def penergy_flux_vs_lat(path_to_file):
    
    # get variables
    f = netCDF4.Dataset(path_to_file)
    temp          = f.variables['temp'][:]
    sig           = f.variables['sigma'][:]
    lat           = f.variables['lat'][:]
    pot_temp      = f.variables['pot_temp'][:]
    v             = f.variables['v'][:]
    z_flux        = f.variables['mrdnl_z_flux'][:]
    temp_flux     = f.variables['mrdnl_tempflux'][:]
    sh_flux       = f.variables['vshum_avg'][:]
    shum_avg      = f.variables['shum_avg'][:]
    
    # set parameters
    del_sol         = 1.0
    start_xtropics  = 25.0
    deg             = pi/180
    l_cond          = 2500000
    radius          = 6371000
    p0              = 100000
    grav            = 9.8000
    cp              = 1.00464e+03
    
    lati, sigi = np.meshgrid(lat,sig)
    # temp and z fluxes include a cos(lat) factor
    temp_flux = temp_flux/np.cos(lati*deg)
    z_flux = z_flux/np.cos(lati*deg)
    
    # variables to be plotted
    energy_flux =  cp*temp_flux + grav*z_flux + l_cond*sh_flux
    dry_energy_flux =  cp*temp_flux+grav*z_flux
    hum_energy_flux =  l_cond*sh_flux
    
    # choose sigma levels
    # energy_flux     =  energy_flux[(sig>0.4)&(sig<0.6),:]
    # dry_energy_flux =  dry_energy_flux[(sig>0.4)&(sig<0.6),:]
    # hum_energy_flux =  hum_energy_flux[(sig>0.4)&(sig<0.6),:]
    # z_energy_flux   =  z_energy_flux[(sig>0.4)&(sig<0.6),:]
        
    # now integrate over sigma
    energy_flux_1d = sigma_int.my_integration(energy_flux, sig)
    dry_energy_flux_1d = sigma_int.my_integration(dry_energy_flux, sig)
    hum_energy_flux_1d = sigma_int.my_integration(hum_energy_flux, sig)
    # z_energy_flux_1d = sigma_int.my_integration(z_energy_flux, sig)
    
    # convert to power in PW
    rescale = 2.0*pi*radius*p0/grav*np.cos(lat*deg)/1e15
    energy_flux_1d = energy_flux_1d*rescale
    dry_energy_flux_1d = dry_energy_flux_1d*rescale
    hum_energy_flux_1d = hum_energy_flux_1d*rescale
    # z_energy_flux_1d = z_energy_flux_1d*rescale
    
    return energy_flux_1d #,dry_energy_flux_1d,hum_energy_flux_1d #,z_energy_flux_1d

def penergy_flux_2d(path_to_file):
    
    # get variables
    f = netCDF4.Dataset(path_to_file)
    temp          = f.variables['temp'][:]
    sig           = f.variables['sigma'][:]
    lat           = f.variables['lat'][:]
    pot_temp      = f.variables['pot_temp'][:]
    v             = f.variables['v'][:]
    z_flux        = f.variables['mrdnl_z_flux'][:]
    temp_flux     = f.variables['mrdnl_tempflux'][:]
    sh_flux       = f.variables['vshum_avg'][:]
    shum_avg      = f.variables['shum_avg'][:]
    
    # set parameters
    del_sol         = 1.0
    start_xtropics  = 25.0
    deg             = pi/180
    l_cond          = 2500000
    radius          = 6371000
    p0              = 100000
    grav            = 9.8000
    cp              = 1.00464e+03
    
    lati, sigi = np.meshgrid(lat,sig)
    # temp and z fluxes include a cos(lat) factor
    temp_flux = temp_flux/np.cos(lati*deg)
    z_flux = z_flux/np.cos(lati*deg)
    
    # variables to be plotted
    energy_flux =  cp*temp_flux + grav*z_flux + l_cond*sh_flux
    dry_energy_flux =  cp*temp_flux + grav*z_flux
    
    return energy_flux, dry_energy_flux

def MSE_vs_lat(path_to_file):
    
    # get variables
    f = netCDF4.Dataset(path_to_file)
    temp          = f.variables['temp'][:]
    sig           = f.variables['sigma'][:]
    lat           = f.variables['lat'][:]
    pot_temp      = f.variables['pot_temp'][:]
    v             = f.variables['v'][:]
    z_flux        = f.variables['mrdnl_z_flux'][:]
    temp_flux     = f.variables['mrdnl_tempflux'][:]
    sh_flux       = f.variables['vshum_avg'][:]
    shum_avg      = f.variables['shum_avg'][:]
    
    # set parameters
    del_sol         = 1.0
    start_xtropics  = 25.0
    deg             = pi/180
    l_cond          = 2500000
    radius          = 6371000
    p0              = 100000
    grav            = 9.8000
    cp              = 1.00464e+03
    
    lati, sigi = np.meshgrid(lat,sig)
    
    z_arr = 1-sig
    lati, z_arri = np.meshgrid(lat,z_arr)
    
    # variables to be plotted
    mse =  cp*temp + grav*z_arri + l_cond*shum_avg
    dry =  cp*temp + grav*z_arri
    moist = l_cond*shum_avg
        
    # now integrate over sigma
    #mse = sigma_int.my_integration(mse, sig)
    #dry = sigma_int.my_integration(dry, sig)
    #moist = sigma_int.my_integration(moist, sig)
    
    # convert to power in PW
    #rescale = 2.0*pi*radius*p0/grav*np.cos(lat*deg)/1e15
    #energy_flux_1d = energy_flux_1d*rescale
    #dry_energy_flux_1d = dry_energy_flux_1d*rescale
    
    return mse

def surf_MSE_vs_lat(path_to_file):
    
    # get variables
    f = netCDF4.Dataset(path_to_file)
    temp          = f.variables['temp'][29,:]
    sig           = f.variables['sigma'][:]
    lat           = f.variables['lat'][:]
    pot_temp      = f.variables['pot_temp'][:]
    v             = f.variables['v'][:]
    z_flux        = f.variables['mrdnl_z_flux'][:]
    temp_flux     = f.variables['mrdnl_tempflux'][:]
    sh_flux       = f.variables['vshum_avg'][:]
    shum_avg      = f.variables['shum_avg'][29,:]
    
    # set parameters
    del_sol         = 1.0
    start_xtropics  = 25.0
    deg             = pi/180
    l_cond          = 2500000
    radius          = 6371000
    p0              = 100000
    grav            = 9.8000
    cp              = 1.00464e+03
    
    lati, sigi = np.meshgrid(lat,sig)
    
    z_arr = 1-sig
    lati, z_arri = np.meshgrid(lat,z_arr)
    z_arri=z_arri[29,:]
    
    # variables to be plotted
    mse =  cp*temp + l_cond*shum_avg
    dry =  cp*temp + grav*z_arri
    moist = l_cond*shum_avg
    
    return mse,dry,moist
