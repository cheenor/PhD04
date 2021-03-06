#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 11:00:28 2015

@author: jhchen
"""
import gdal,ogr    ### osheo can't be imported with pyhdf
import osgeo
import gdalinfo
import numpy as np
import matplotlib.pyplot as plt
### this is a fuction to obtain the lon and lat in gdal  
def pixel2coord(x, y,geoform):
    """Returns global coordinates from pixel x, y coords"""
    xoff=geoform[0]    
    a=geoform[1]
    b=geoform[2]
    yoff=geoform[3]
    d=geoform[4]
    e=geoform[5]
    xp = a * x + b * y + xoff
    yp = d * x + e * y + yoff
    return(xp, yp)
# xoff, a, b, yoff, d, e = ds.GetGeoTransform()
dirin="D:/MyPaper/PhD04/Cases/ETP/20100604_0704/MODIS/ETP20100604-32d/"
fnm="MYD08_D3.A2010154.006.2014264200556.hdf"
fpath=dirin+fnm
modis=gdal.Open(fpath)
modisinfo=gdalinfo(fpath)
listvar=modis.GetSubDatasets()
lc_data=gdal.Open(listvar[253][0])
lc=lc_data.ReadAsArray()
nx=lc_data.RasterXSize
ny=lc_data.RasterYSize
geoform=lc_data.GetGeoTransform()
lon=np.ndarray(shape=(nx,ny), dtype=float)
lat=np.ndarray(shape=(nx,ny), dtype=float)
for i in range(0,nx):
    for j in range(0,ny):
        lon[i,j],lat[i,j] =pixel2coord(i,j,geoform)
minlat=int(lat.min())
maxlat=int(lat.max())
minlon=int(lon.min())
maxlon=int(lon.max())
lonlabs=[]
latlabs=[]
for i in range(0,nx):
    itm="%d"%int(lon[i,0])
    lonlabs.append(itm)
for i in range(0,ny):
    itm="%d"%int(lat[0,i])
    latlabs.append(itm)    
fig,ax = plt.subplots(nrows=1,ncols=1)
#ax.set_ylim(35,80)
#ax.set_xlim(250,350)
plt.imshow(lc[:,:])
ax.set_xticks(range(0,nx,30))
xticklabels = [lonlabs[nn] for nn in range(0,nx,30)] 
ax.set_xticklabels(xticklabels, size=14)
ax.set_yticks(range(0,ny,15))
yticklabels = [latlabs[nn] for nn in range(0,ny,15)] 
ax.set_yticklabels(yticklabels, size=14)
plt.colorbar()
plt.show()
"""
fpath=dirin+"MYD08_D3.txt"
i=0
f=open(fpath,'w')
for nm in listvar:
    iband="%d "%i
    f.write(iband)
    itme="%s "%nm[1]
    f.write(itme)
    f.write('\n')
    i+=1
f.close()
"""    