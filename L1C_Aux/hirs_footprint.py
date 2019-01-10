#!/usr/bin/env python

import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def alpha2km(alpha, R, h):

    """ 
    Function to convert across-track scan angle (ie. scan mirror angle) to 
    across-track distance in km, given alpha (radians), Earth radius  and 
    orbital height.
    Inputs: Alpha: Along-track scan angle
            R: Earth Radius
            h: Orbital Height
    Returns: km: Across-track scan distance
    """

    #Gives the 1st quadrant angle (even though gamma is technically 2nd quadrant)
    gamma = math.asin(((R + h)/R)*math.sin(alpha))
    #Exploits the fact that returned gamma is in the 1st quadrant...
    beta  = gamma - alpha
    #Across-track position of centre of EFOV in km
    km    = R*beta

    return km

def getgacpix(hirs_spot, dtime, hirs_sensor):

    """
    Determine which GAC pixels fall into each HIRS footprint.  Pixels are
    included where the GAC lat/lon falls within the HIRS footrpint.
    Inputs: hirs_spot: HIRS pixel along scan line
            dtime: delta time of HIRS pixel along scan line
    Outputs:use_arr: output array of weights for pixels to use 
    """

    nx = 0
    ny = 0
    x0 = 0

    #----------------------------
    #Define a series of constants
    #----------------------------
    h = 819.      #Nominal METOP-A orbital height
    R = 6372.     #Earth Radius
    speed = 6.7   #Ground-track speed in km.s-1

    HIRS2 = ('noaa08','noaa09','noaa10','noaa12','noaa14')
    HIRS2i = ('noaa11','noaa14')
    HIRS3 = ('noaa15','noaa16','noaa17')
    HIRS4 = ('noaa18','noaa19','metopa','metopb')


    if hirs_sensor.lower() in HIRS2:   
        IFOV = math.radians(1.22)
    elif hirs_sensor.lower() in HIRS2i or hirs_sensor in HIRS3:
        IFOV = math.radians(1.4)
    elif hirs_sensor.lower() in HIRS4:
        IFOV = math.radians(0.69) 
                              # HIRS/2 IFOV ~= 0.0213 radians.
                              # HIRS/2i/3 is 1.4 degrees, HIRS/4 is 0.69 degrees.
    nu   = math.radians(1.8)  # 1.8 degree step angle ~=0.0314 radians
    nIFOVs = 56               #Number of HIRS FOVs
    half_scan = (nIFOVs-1)/2. #Number of HIRS FOVs between edge and nadir    
    Hz  = 6.                  #Raw AVHRR scan rate (lines/sec)
    anu = (Hz*2*math.pi)/39936. #AVHRR angle between pixels (39936 samples/sec @ 6Hz) 
    Hz  = 2.                  #Effective GAC scan rate (lines/sec - use 1 and skip 2 lines)
    anu = anu*5.              #Angle between GAC pixels (average 4, skip 1)    
    a = IFOV/(2.*anu)         #Semi-major (aka across-track) axis will always be the same in
                              #pixel units.
 
    #----------------
    #Main calculation
    #----------------
    alpha = (hirs_spot - 27.5)*nu   #Angle centre of HIRS spot from nadir (radians)
    x0    = (204.4*anu + alpha)/anu #Centre of HIRS FOV in GAC pixel units (starting
                                    #from zero).
    ix0   = round(x0)               #Closest integer HIRS spot
    dx0   = x0 - ix0                #Across-track offset of HIRS FOV centre from closest GAC
    dy0   = dtime * Hz              #Along-track offset of HIRS FOV from centre of closest GAC
    beta  = alpha2km(alpha, R, h) / R #Angle at centre of the Earth
    ay_km = speed*0.5*math.cos(beta)  #Spacing between GAC lines at given across-track position
    d     = R*math.sin(beta)/math.sin(alpha) #Distance from satellite to centre of HIRS spot (km)
    b     = ((d*IFOV)/2.0)/ay_km    #Semi-minor (aka along-track) axis of HIRS FOV in GAC 
                                    #pixel units.

    #---------------------------------------------------------------------
    #Run through the supplied array of GAC pixels and find those that fit.
    #---------------------------------------------------------------------
    nx = long(math.ceil(a))         #Number of pixels in sub-array (a is fixed so this is constant)
    ny = long(math.ceil(b))         #Number of rows in sub-array


    use_arr = np.zeros([(2*ny)+1,(2*nx)+1])  #Array of pixels centred on GAC pixel closest to
                                             #centre of HIRS FOV, inital weight = 0.


    for iy in np.arange(0,2*ny,1):
        
        dy = iy - ny + dy0    #Y-coordinate of GAC pixel from HIRS FOV centre (pixel units)
        
        for ix in np.arange(0,2*nx,1):

            dx = ix - nx -dx0 #X-coordinate of GAC pixel form HIRS FOV centre (pixel units)

            if ((dx**2/a**2)+(dy**2/b**2)) < 1.:
                use_arr[iy,ix] = 1.


    return use_arr,x0,a,b,nx,ny,ix0


if __name__ in '__main__':


    #Test code to check implementation of the functions above.
    plt.figure(figsize=[20,6])
    plt.ylim(-25,25)
    plt.xlabel('Across-track position / AVHRR GAC pixels')
    plt.ylabel('Along-track position / AVHRR GAC pixels')
    plt.title('HIRS Footprints')

    arr = np.arange(0,100,1)
    x_plot_arr = np.zeros([100])
    y_plot_arr = np.zeros([100])

    npnts = 100
    th = (math.pi*2.*arr) / npnts-1

    for i in np.arange(0,55,1):
        dtime = 0.2
        y0 = dtime*2.
    
        use_arr,x0,a,b,nx,ny,ix0 = getgacpix(i,dtime,hirs_sensor)

        for x in np.arange(0,th.shape[0],1):
            x_plot_arr[x] = x0+a*math.cos(th[x])
            y_plot_arr[x] = y0+b*math.sin(th[x])

        plt.plot(x_plot_arr,y_plot_arr,'r')

        x_arr = np.zeros([use_arr.shape[0],use_arr.shape[1]])
        y_arr = np.zeros([use_arr.shape[0],use_arr.shape[1]])
 
        for iy in np.arange(0,2*ny,1):

            y_arr[:,iy] = (ny - iy)
            
            for ix in np.arange(0,2*nx,1):

                x_arr[ix,iy] = ix0 - (nx - ix)



        ii = (use_arr == 1)

        plt.plot(x_arr[ii],y_arr[ii],'b.',markersize=1)

    plt.show()
