import numpy as np

def cal_area(lon,lat,dlon,dlat) :
# input : lon, lat (not array, but single element)
# output : area that matches the [lon,lat] coordinate 
#          units in cm^2

    # constant
    arad = np.float64(6.371E8) #cm
    dlon = np.float64(dlon*np.pi/180.)
    dlat = np.float64(dlat*np.pi/180.)
    lat = (90.-lat)*np.pi/180.
    area = arad*arad*np.sin(lat)*dlon*dlat

    return {'area':area}