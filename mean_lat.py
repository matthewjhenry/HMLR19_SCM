import math
import numpy as np

def mean_lat(lat,val):
    
    if(lat.shape[0] != val.shape[0]):
        print("Latitude and value vectors not of same size.")
        return 0
    else :
        # assume lat in degrees
        w = np.cos(lat*math.pi/180)
        return np.sum(val*w)/np.sum(w)