import numpy as np
import numpy.matlib

def my_integration(f_,sigma_full_):

    # define sigma on half levels
    num_lev_    = sigma_full_.shape[0]
    sigma_half_ = np.zeros(num_lev_+1)

    #  lowest half level is ground
    sigma_half_[num_lev_] = 1
    
    sigma_half_[0]=0

    for j in range(num_lev_-1,-1,-1):
        sigma_half_[j] = 2*sigma_full_[j] - sigma_half_[j+1]
    
    sigma_half_[0]=0 
    d_sigma_  = np.diff(sigma_half_)
    d_temp_   = np.transpose(d_sigma_).reshape(f_.shape[0],1)
    d_sigma2_ = np.matlib.repmat(d_temp_,1,f_.shape[1])
    
    fin=np.sum(f_*d_sigma2_,0)
    
    return fin

def press_int_k(f, plev_full):
    plev_full_0 = np.concatenate([[0], plev_full],axis=0)
    d_plev = np.diff(plev_full_0)/100
    d_temp_ = np.transpose(d_plev).reshape(f.shape[0],1)
    d_plev2 = np.matlib.repmat(d_temp_, 1, f.shape[1])
    F      = np.sum(f*d_plev2,0)
    return F