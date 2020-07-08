from __future__ import print_function # make sure print behaves the same in 2.7 and 3.x
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from mean_lat import mean_lat
from plot_meridional_E_transport import penergy_flux_2d, penergy_flux_vs_lat, MSE_vs_lat, surf_MSE_vs_lat
import ClimateUtils as clim
import os.path
import sigma_int

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

class ncdf_obj:
    
    def __init__(self, dirname, ctl, test, af_1xCO2='', af_14xCO2=''):
        self.ctl = dirname + ctl + '/day_ave.nc'
        self.test = dirname + test + '/day_ave.nc'
        self.f_lin = netCDF4.Dataset(dirname + ctl + '/day_ave.nc')
        self.f_lin_14xCO2 = netCDF4.Dataset(dirname + test + '/day_ave.nc')
        self.f_lin.set_always_mask(False)
        self.f_lin_14xCO2.set_always_mask(False)
        
        self.sig_= self.f_lin.variables['sigma'][:]
        self.lat_= self.f_lin.variables['lat'][:]
        self.lat_mesh, self.sig_mesh = np.meshgrid(self.lat_,self.sig_)
        
        self.af_1xCO2 = af_1xCO2
        self.af_14xCO2 = af_14xCO2
        
        #self.trop_level_lin = self.f_lin.variables['trop_level1'][:]
        #self.trop_level_lin_14xCO2 = self.f_lin_14xCO2.variables['trop_level1'][:]
        
        # staggered latitude array
        self.lat_stag = np.zeros(self.lat_.shape[0]+1)
        self.lat_stag[0]=-90
        self.lat_stag[self.lat_.shape[0]]=90
        for i in range(self.lat_.shape[0]-1):
            self.lat_stag[i+1]=0.5*(self.lat_[i]+self.lat_[i+1])
            
        # staggered sig array
        self.sig_stag = np.zeros(self.sig_.shape[0]+1)
        self.sig_stag[0]=0
        self.sig_stag[self.sig_.shape[0]]=1
        for i in range(self.sig_.shape[0]-1):
            self.sig_stag[i+1]=0.5*(self.sig_[i]+self.sig_[i+1])
            
        self.phi = np.deg2rad( self.lat_ )
        self.phi_stag = np.deg2rad( self.lat_stag )
        
        self.temp_lin=self.f_lin.variables['temp'][:]
        self.temp_lin_14xCO2=self.f_lin_14xCO2.variables['temp'][:]
        self.temp_lin_diff = self.temp_lin_14xCO2 - self.temp_lin
        
        self.q_lin=self.f_lin.variables['shum_avg'][:]
        self.q_lin_14xCO2=self.f_lin_14xCO2.variables['shum_avg'][:]
        self.q_lin_diff = self.q_lin_14xCO2 - self.q_lin
        
        self.u_var_lin=self.f_lin.variables['u_var'][:]
        self.u_var_lin_14xCO2=self.f_lin_14xCO2.variables['u_var'][:]
        
        self.v_var_lin=self.f_lin.variables['v_var'][:]
        self.v_var_lin_14xCO2=self.f_lin_14xCO2.variables['v_var'][:]
        
        self.u=self.f_lin.variables['u'][:]
        self.u_14xCO2=self.f_lin_14xCO2.variables['u'][:]
        
        self.v=self.f_lin.variables['v'][:]
        self.v_14xCO2=self.f_lin_14xCO2.variables['v'][:]
        
        self.sat_lin=self.f_lin.variables['temp'][:][29,:]
        self.sat_lin_14xCO2=self.f_lin_14xCO2.variables['temp'][:][29,:]
        self.sat_lin_diff = self.sat_lin_14xCO2 - self.sat_lin
        
        self.s_lin = netCDF4.Dataset(dirname + ctl + '/surf_ave.nc')
        self.s_lin_14xCO2 = netCDF4.Dataset(dirname + test + '/surf_ave.nc')
        self.lon_= self.s_lin.variables['lon'][:]
        self.rhum_lin = self.f_lin.variables['rhum_avg'][:]
        self.rhum_lin_14xCO2 = self.f_lin_14xCO2.variables['rhum_avg'][:]
        
        # get extra variables
        self.pot_temp_lin = self.f_lin.variables['pot_temp'][:]
        self.d_dp_pot_temp_lin = self.f_lin.variables['d_dp_pot_temp'][:]
        self.w_lin = self.f_lin.variables['w_avg'][:]
        
        self.pot_temp_lin_14xCO2 = self.f_lin_14xCO2.variables['pot_temp'][:]
        self.d_dp_pot_temp_lin_14xCO2 = self.f_lin_14xCO2.variables['d_dp_pot_temp'][:]
        self.w_lin_14xCO2 = self.f_lin_14xCO2.variables['w_avg'][:]
        
        self.vert_tend_lin = self.temp_lin / self.pot_temp_lin * self.w_lin * self.d_dp_pot_temp_lin
        self.vert_tend_lin_14xCO2 = self.temp_lin_14xCO2 / self.pot_temp_lin_14xCO2 * self.w_lin_14xCO2 * self.d_dp_pot_temp_lin_14xCO2
        
        #self.fb_mse_div_diff_lin = -self.LW_flux(0, double = True)+self.LW_flux(0, double = False)
        '''
        self.lwflux_1xCO2 = self.s_lin.variables['flux_lw'][:]
        self.lwflux_1xCO2 = np.mean(self.lwflux_1xCO2,axis=0)
        self.lwflux_1xCO2 = np.mean(self.lwflux_1xCO2,axis=2)
        self.lwflux_14xCO2 = self.s_lin_14xCO2.variables['flux_lw'][:]
        self.lwflux_14xCO2 = np.mean(self.lwflux_14xCO2,axis=0)
        self.lwflux_14xCO2 = np.mean(self.lwflux_14xCO2,axis=2)
        '''
        
        #self.net_lw_surf = self.s_lin.variables['net_lw_surf'][:]
        #self.net_lw_surf = np.mean(self.net_lw_surf,axis=0)
        #self.net_lw_surf = np.mean(self.net_lw_surf,axis=1)
        
        #self.rf_net_lw_surf = self.s_lin.variables['rf_net_lw_surf'][:]
        #self.rf_net_lw_surf = np.mean(self.rf_net_lw_surf,axis=0)
        #self.rf_net_lw_surf = np.mean(self.rf_net_lw_surf,axis=1)
        
        #self.convec = np.mean(np.mean(self.s_lin.variables['dt_tg_convection'][:],axis=0),axis=2)
        #self.int_convec = sigma_int.my_integration(self.convec, self.sig_)
        
        self.sh_1xCO2    = self.f_lin.variables['shum_avg'][:]
        self.sh_14xCO2   = self.f_lin_14xCO2.variables['shum_avg'][:]
        self.sh_surf_1xCO2    = self.sh_1xCO2[29,:]
        self.sh_surf_14xCO2   = self.sh_14xCO2[29,:]
        self.z_1xCO2     = self.f_lin.variables['z'][:]
        self.z_14xCO2    = self.f_lin_14xCO2.variables['z'][:]
        self.mse_1xCO2   = clim.cp*self.temp_lin + clim.g*self.z_1xCO2 + clim.Lhvap*self.sh_1xCO2
        self.mse_14xCO2  = clim.cp*self.temp_lin_14xCO2 + clim.g*self.z_14xCO2 + clim.Lhvap*self.sh_14xCO2
        self.mse_surf_1xCO2 = self.mse_1xCO2[29,:]
        self.mse_surf_14xCO2 = self.mse_14xCO2[29,:]
        
        self.mse_calc_1xCO2   = clim.cp*self.temp_lin + clim.g*self.z_1xCO2 + clim.Lhvap * 0.8 * clim.qsat(self.temp_lin,self.sig_mesh*1000)
        self.mse_calc_14xCO2  = clim.cp*self.temp_lin_14xCO2 + clim.g*self.z_14xCO2 + clim.Lhvap * 0.8 * clim.qsat(self.temp_lin_14xCO2,self.sig_mesh*1000)
        self.mse_calc_surf_1xCO2 = self.mse_calc_1xCO2[29,:]
        self.mse_calc_surf_14xCO2 = self.mse_calc_14xCO2[29,:]
        
        self.mse_sat_1xCO2   = clim.cp*self.temp_lin + clim.g*self.z_1xCO2 + clim.Lhvap*clim.qsat(self.temp_lin,self.sig_mesh*1000)
        self.mse_sat_14xCO2  = clim.cp*self.temp_lin_14xCO2 + clim.g*self.z_14xCO2 + clim.Lhvap * clim.qsat(self.temp_lin_14xCO2,self.sig_mesh*1000)
        self.mse_sat_atm_1xCO2 = sigma_int.my_integration(self.mse_sat_1xCO2[self.sig_<0.5,:], self.sig_[self.sig_<0.5])
        self.mse_sat_atm_14xCO2 = sigma_int.my_integration(self.mse_sat_14xCO2[self.sig_<0.5,:], self.sig_[self.sig_<0.5])
        
        if os.path.isfile('/storage/mhenry/fixedsst_idealized_moist_forcing/' + self.af_1xCO2 + '/day_ave.nc'):
            tf_1xCO2 = netCDF4.Dataset('/storage/mhenry/fixedsst_idealized_moist_forcing/' + self.af_1xCO2 + '/day_ave.nc')
            tf_14xCO2 = netCDF4.Dataset('/storage/mhenry/fixedsst_idealized_moist_forcing/' + self.af_14xCO2 + '/day_ave.nc')
            self.dt_forcing = tf_14xCO2.variables['temp'][:] - tf_1xCO2.variables['temp'][:]
            
    def sig_change(self,field): # change field size from (sig_stag,lat) to (sig,lat)
        new_f = np.zeros((self.sig_.shape[0],self.lat_.shape[0]))
        for i in range(len(self.sig_)):
            for j in range(len(self.lat_)):
                new_f[i,j] = 0.5*(field[i,j]+field[i+1,j])
        return new_f

    def olr_contrib(self,lin=False, double=False):
        olr_contrib = np.zeros_like(self.temp_lin)
        for i in range(0,30):
            olr_contrib[i,:]= -self.LW_flux(i,double = double, lin = lin)+self.LW_flux(i+1,double = double, lin = lin)
        return olr_contrib

    def flux_cte(self, lin=True):
        if(lin):
            return clim.Lhvap*(0.07**2)*self.q_lin[29,:]*self.divergence_field(self.sat_lin)
        else:
            return clim.Lhvap*(0.07**2)*self.q_norm[29,:]*self.divergence_field(self.sat_norm)

    def eff_level(self,olr,T):
        # input : 1d OLR and 2d T
        # output : plot of Z_eff vs. latitude
        Teff = np.power(olr/clim.sigma,0.25)
        Zeff = np.ones_like(Teff)
        for i in range(T.shape[1]):
            Zeff[i]=self.sig_[find_nearest(Teff[i],T[:,i])]
        return Zeff
        
    def plot_pole(self, field_lin, lat):

        lr_lin=np.zeros((field_lin).shape[0],)

        for i in range(field_lin.shape[0]):
            lr_lin[i]=mean_lat(self.lat_[self.lat_>lat],field_lin[i,self.lat_>lat])
            
        return lr_lin
    
    def plot_eq(self, field_lin):

        lr_lin=np.zeros((field_lin).shape[0],)
        masks = [self.lat_>-10,self.lat_<10]
        mask = masks[0] & masks[1]

        for i in range(field_lin.shape[0]):
            lr_lin[i]=mean_lat(self.lat_[mask],field_lin[i,mask])
        
        return lr_lin
        
    def planck_feedback(self):
        surf = self.s_lin

        olr = surf.variables['olr'][:]
        olr = np.mean(olr,axis=0)

        olr_1deg = surf.variables['olr_pl'][:]
        olr_1deg = np.mean(olr_1deg,axis=0)

        forc = olr - olr_1deg 
        forc = np.mean(forc,axis=1)

        return forc
    
    def MSE_div_diff(self):
        div_MSE_1xCO2 = -self.LW_flux(0,double=False)
        div_MSE_14xCO2 = -self.LW_flux(0,double=True)
        return -div_MSE_14xCO2 + div_MSE_1xCO2
    
    def total_feedback(self, adjusted):
        self.fb_forc=self.forcing(adjusted)
        self.fb_sat_diff = self.temp_lin_14xCO2[29,:]-self.temp_lin[29,:]
        self.fb_mse_div_diff = -self.LW_flux(0, double = True)+self.LW_flux(0, double = False)     
        tot_fb = -(self.fb_forc + self.fb_mse_div_diff)/self.fb_sat_diff
        return tot_fb

        
    def SW_flux(self, level, double = False):
        # returns surface SW flux
        if(double):
            s = self.s_lin_14xCO2
        else:
            s = self.s_lin
            
        swgcm = s.variables['flux_sw'][:]
        swgcm = np.mean(swgcm,axis=0)
        swgcm = np.mean(swgcm,axis=2)
        swgcm = swgcm[level,:]
        swgcm = -swgcm
        return swgcm
    
    def LW_flux(self, level, double = False):
        # returns surface LW flux
        if(double):
            s = self.s_lin_14xCO2
        else:
            s = self.s_lin
            
        swgcm = s.variables['flux_lw'][:]
        swgcm = np.mean(swgcm,axis=0)
        swgcm = np.mean(swgcm,axis=2)
        swgcm = swgcm[level,:]
        swgcm = -swgcm
        return swgcm

    
    def forcing(self, adjusted):
        if(adjusted):
            f_1xCO2 = netCDF4.Dataset('/storage/mhenry/fixedsst_idealized_moist_forcing/' + self.af_1xCO2 + '/surf_ave.nc')
            f_14xCO2 = netCDF4.Dataset('/storage/mhenry/fixedsst_idealized_moist_forcing/' + self.af_14xCO2 + '/surf_ave.nc')
            
            lwflux_1xCO2 = f_1xCO2.variables['flux_lw'][:]
            lwflux_1xCO2 = np.mean(lwflux_1xCO2,axis=0)
            lwflux_1xCO2 = np.mean(lwflux_1xCO2,axis=2)
            lwflux_14xCO2 = f_14xCO2.variables['flux_lw'][:]
            lwflux_14xCO2 = np.mean(lwflux_14xCO2,axis=0)
            lwflux_14xCO2 = np.mean(lwflux_14xCO2,axis=2)
            
            return lwflux_1xCO2[0,:]-lwflux_14xCO2[0,:]
        else:
            surf=self.s_lin

            olr_1xCO2 = surf.variables['olr'][:]
            olr_1xCO2 = np.mean(olr_1xCO2,axis=0)

            olr_14xCO2 = surf.variables['inst_rf'][:]
            olr_14xCO2 = np.mean(olr_14xCO2,axis=0)

            forc = olr_1xCO2 - olr_14xCO2 
            forc = np.mean(forc,axis=1)

            return forc
    
    
    def plot2d_temp(self, field, title="Enter title", mini=0, maxi=0, save_dir=''):
        
        lat, sig = np.meshgrid(self.lat_,self.sig_)
        plt.pcolor(lat, sig, field, cmap='RdBu_r')
        plt.axis([lat.min(), lat.max(), sig.max(), sig.min()])
        plt.colorbar()
        if(mini != maxi):
            plt.clim(mini,maxi)
        plt.xticks(np.linspace(-90, 90, 7))
        plt.xlabel('Latitude (deg N)')
        plt.ylabel('Sigma level (p/ps)')
        plt.title(title)
        fig = plt.gcf()
        plt.show()
        
        if save_dir != "":
            fig.savefig(save_dir)
            
    def plot_pole_lr(self, title="Enter title",double=True,  save_dir=''):
        
        if(double == True):
            temp_norm_diff = self.temp_norm_2xCO2-self.temp_norm
            temp_lin_diff = self.temp_lin_2xCO2-self.temp_lin
        else :
            temp_norm_diff = self.temp_norm_14xCO2-self.temp_norm
            temp_lin_diff = self.temp_lin_14xCO2-self.temp_lin

        # Use 70-90 latitude band and take weighted average T change at each sigma level
        lr_norm=np.zeros((temp_norm_diff).shape[0],)
        lr_lin=np.zeros((temp_lin_diff).shape[0],)

        for i in range(temp_norm_diff.shape[0]):
            lr_norm[i]=mean_lat(self.lat_[self.lat_>70],temp_norm_diff[i,self.lat_>70])
            lr_lin[i]=mean_lat(self.lat_[self.lat_>70],temp_lin_diff[i,self.lat_>70])

        plt.plot(lr_lin,self.sig_,'b')
        plt.plot(lr_norm,self.sig_,'r')
        plt.gca().invert_yaxis()
        plt.grid()
        plt.xlabel('Temperature change (K)')
        plt.ylabel('Sigma level (p/p0)')
        plt.title(title)
        plt.legend(['linear','nonlinear'], loc="best")
        fig = plt.gcf()
        plt.show()
        
        if save_dir != "":
            fig.savefig(save_dir)
        
    def mse_flux(self, double=False):
        if double:
            [mse_flux,dry_flux,hum_flux]=penergy_flux_vs_lat(self.test)
            return mse_flux
        else:
            [mse_flux,dry_flux,hum_flux]=penergy_flux_vs_lat(self.ctl)
            return mse_flux
        
    def dry_flux(self, double=False):
        if double:
            [mse_flux,dry_flux,hum_flux]=penergy_flux_vs_lat(self.test)
            return dry_flux
        else:
            [mse_flux,dry_flux,hum_flux]=penergy_flux_vs_lat(self.ctl)
            return dry_flux
        
    def hum_flux(self, double=False):
        if double:
            [mse_flux,dry_flux,hum_flux]=penergy_flux_vs_lat(self.test)
            return hum_flux
        else:
            [mse_flux,dry_flux,hum_flux]=penergy_flux_vs_lat(self.ctl)
            return hum_flux
        
    def mse_flux_2d(self, double=False):
        if double:
            [mse_flux,dry_flux]=penergy_flux_2d(self.test)
            return mse_flux
        else:
            [mse_flux,dry_flux]=penergy_flux_2d(self.ctl)
            return mse_flux
        
    def dry_flux_2d(self, double=False):
        if double:
            [mse_flux,dry_flux]=penergy_flux_2d(self.test)
            return dry_flux
        else:
            [mse_flux,dry_flux]=penergy_flux_2d(self.ctl)
            return dry_flux
        
    def hum_flux_2d(self, double=False):
        if double:
            [mse_flux,dry_flux]=penergy_flux_2d(self.test)
            return mse_flux-dry_flux
        else:
            [mse_flux,dry_flux]=penergy_flux_2d(self.ctl)
            return mse_flux-dry_flux
    
    def mse_conv_1d(self, double=False):
        
        mse_flu = self.mse_flux(double)
        mse_flu_stag = np.zeros(mse_flu.shape[0]+1)
        mse_flu_stag[0]=0
        mse_flu_stag[self.lat_.shape[0]]=0
        for i in range(self.lat_.shape[0]-1):
            mse_flu_stag[i+1]=0.5*(mse_flu[i]+mse_flu[i+1])

        return ( -1./(2*np.math.pi*clim.a**2*np.cos(self.phi)) * np.diff( 1.E15*mse_flu_stag ) / np.diff(self.phi_stag) )
    
    def dry_conv_1d(self, double=False):
        
        mse_flu = self.dry_flux(double)
        mse_flu_stag = np.zeros(mse_flu.shape[0]+1)
        mse_flu_stag[0]=0
        mse_flu_stag[self.lat_.shape[0]]=0
        for i in range(self.lat_.shape[0]-1):
            mse_flu_stag[i+1]=0.5*(mse_flu[i]+mse_flu[i+1])

        return ( -1./(2*np.math.pi*clim.a**2*np.cos(self.phi)) * np.diff( 1.E15*mse_flu_stag ) / np.diff(self.phi_stag) )
    
    def hum_conv_1d(self, double=False):
        
        mse_flu = self.hum_flux(double)
        mse_flu_stag = np.zeros(mse_flu.shape[0]+1)
        mse_flu_stag[0]=0
        mse_flu_stag[self.lat_.shape[0]]=0
        for i in range(self.lat_.shape[0]-1):
            mse_flu_stag[i+1]=0.5*(mse_flu[i]+mse_flu[i+1])

        return ( -1./(2*np.math.pi*clim.a**2*np.cos(self.phi)) * np.diff( 1.E15*mse_flu_stag ) / np.diff(self.phi_stag) )
    
    def mse_conv_2d(self, double=False):
        mse_flu = self.mse_flux_2d(double)
        mse_flu_stag = np.zeros((self.sig_.shape[0],mse_flu.shape[1]+1))
        mse_conv = np.zeros((self.sig_.shape[0],mse_flu.shape[1]))
        mse_flu_stag[:,0]=0
        mse_flu_stag[:,self.lat_.shape[0]]=0
        for j in range(self.sig_.shape[0]):
            for i in range(self.lat_.shape[0]-1):
                mse_flu_stag[j,i+1]=0.5*(mse_flu[j,i]+mse_flu[j,i+1])
            mse_conv[j,:]=( -1./(2*np.math.pi*clim.a**2*np.cos(self.phi)) * np.diff( 1.E15*mse_flu_stag[j,:] ) / np.diff(self.phi_stag) )

        return mse_conv
    
    def dry_conv_2d(self, double=False):
        mse_flu = self.dry_flux_2d(double)
        mse_flu_stag = np.zeros((self.sig_.shape[0],mse_flu.shape[1]+1))
        mse_conv = np.zeros((self.sig_.shape[0],mse_flu.shape[1]))
        mse_flu_stag[:,0]=0
        mse_flu_stag[:,self.lat_.shape[0]]=0
        for j in range(self.sig_.shape[0]):
            for i in range(self.lat_.shape[0]-1):
                mse_flu_stag[j,i+1]=0.5*(mse_flu[j,i]+mse_flu[j,i+1])
            mse_conv[j,:]=( -1./(2*np.math.pi*clim.a**2*np.cos(self.phi)) * np.diff( 1.E15*mse_flu_stag[j,:] ) / np.diff(self.phi_stag) )

        return mse_conv
    
    def hum_conv_2d(self, double=False):
        mse_flu = self.hum_flux_2d(double)
        mse_flu_stag = np.zeros((self.sig_.shape[0],mse_flu.shape[1]+1))
        mse_conv = np.zeros((self.sig_.shape[0],mse_flu.shape[1]))
        mse_flu_stag[:,0]=0
        mse_flu_stag[:,self.lat_.shape[0]]=0
        for j in range(self.sig_.shape[0]):
            for i in range(self.lat_.shape[0]-1):
                mse_flu_stag[j,i+1]=0.5*(mse_flu[j,i]+mse_flu[j,i+1])
            mse_conv[j,:]=( -1./(2*np.math.pi*clim.a**2*np.cos(self.phi)) * np.diff( 1.E15*mse_flu_stag[j,:] ) / np.diff(self.phi_stag) )

        return mse_conv
    
    def mse_conv_diff(self, lin=True):
        return self.mse_convergence(lin, double=True)-self.mse_convergence(lin, double=False)
    
class ncdf_obj_tsurf:
    
    def __init__(self, dirname):
        self.dirn  = dirname
        self.ctl   = netCDF4.Dataset(dirname + '/ctl/day_ave.nc')
        self.trop  = netCDF4.Dataset(dirname + '/trop/day_ave.nc')
        self.poles = netCDF4.Dataset(dirname + '/poles/day_ave.nc')
        
        self.sig_ = self.ctl.variables['sigma'][:]
        self.lat_ = self.ctl.variables['lat'][:]
        
        self.temp_ctl  = self.ctl.variables['temp'][:]
        self.temp_trop  = self.trop.variables['temp'][:]
        self.temp_poles = self.poles.variables['temp'][:]
        self.sat_ctl   = self.ctl.variables['temp'][:][29,:]
        self.sat_trop   = self.trop.variables['temp'][:][29,:]
        self.sat_poles  = self.poles.variables['temp'][:][29,:]
        
    def plot2d_temp(self, field, title="Enter title", mini=0, maxi=0, save_dir=''):
        
        lat, sig = np.meshgrid(self.lat_,self.sig_)
        plt.pcolor(lat, sig, field, cmap='RdBu_r')
        plt.axis([lat.min(), lat.max(), sig.max(), sig.min()])
        plt.colorbar()
        if(mini != maxi):
            plt.clim(mini,maxi)
        plt.xticks(np.linspace(-90, 90, 7))
        plt.xlabel('Latitude (deg N)')
        plt.ylabel('Sigma level (p/ps)')
        plt.title(title)
        fig = plt.gcf()
        plt.show()
        
        if save_dir != "":
            fig.savefig(save_dir)
