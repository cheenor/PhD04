#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30  2015

@author: jhchen
"""
from netCDF4 import Dataset #用于读取netcdf文件
from mpl_toolkits.basemap import Basemap# 地图
import numpy as np  # 数组，数值计算
import matplotlib.pyplot as plt  # 画图 模块
import math  #数学公式模块
import string  #字符串处理
import os    # 系统模块，可以建立文件夹，删除文件，改文件名等等
os.system("cls") # 清除内存，类似与matlab的clc
etplat=[30,  37.5, 37.5, 30, 30]
etplon=[90, 90,  100,  100, 90]
wtplat=[30,  37.5, 37.5, 30, 30]
wtplon=[80, 80,  90,  90, 80]
neclat=[42.5, 50, 50, 42.5, 42.5]
neclon=[120, 120,  127.5, 127.5, 120]
npclat=[35,  42.5, 42.5, 35, 35]
npclon=[112.5, 112.5,  120,   120, 112.5]
mlyrlat=[27.5,  35, 35, 27.5, 27.5]
mlyrlon=[112.5, 112.5,  120,   120, 112.5]
prdlat=[20,  27.5, 27.5, 20, 20]
prdlon=[110, 110,  117.5,   117.5, 110]
fonttp = {'family' : 'serif',
        'color'  : 'w',
        'weight' : 'normal',
        'size'   : 18,
        }  
fontec = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 18,
        }  
dirout='D:/MyPaper/PhD04/Pics/'
dirin='D:/MyPaper/Phd01/data/'
filename='LAND_10MIN.nc'
fpath=dirin+filename
f=Dataset(fpath,'a')  # 打开netcdf文件
for a in f.variables:
    print a  # 这个循环将输出文件中所有变量的名字，如果你不知道经纬度等的名字的话，可以先执行这个循环
#
lon=f.variables['lon'][:]   # 得到经度 
lat=f.variables['lat'][:]    # 纬度
time=f.variables['time'][:]       # 时间 注意，时间是一个整数，这个以后在说
hight=f.variables['HT'][:]
nx=len(lon) # 
ny=len(lat)
fig = plt.figure(figsize=(10,8)) ### figure 8"X8" inches，图片大小
m=Basemap(projection='lcc', llcrnrlon=77.,llcrnrlat=15.,urcrnrlon=150.,urcrnrlat=53.,
          lat_1=15.,lat_2=45.,lat_0=35.,lon_0=105.,rsphere=6371200., resolution='l') 
m.drawcoastlines() # 显示地图海岸线
m.drawcountries() # 显示国界线
parallels = np.arange(-90,90,10.) #地图上的纬度线，间隔是10，范围是-90到90
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14) #地图上画出
# draw meridians
meridians = np.arange(0.,360.,10.) #地图上的经度线，间隔是10，范围是-90到90
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14) # 地图上画出
lon_grid=np.ndarray(shape=(ny,nx),dtype=float)
lat_grid=np.ndarray(shape=(ny,nx),dtype=float)
for i in range(0,ny):
    for j in range(0,nx):
        lon_grid[i,j]=lon[j]
        lat_grid[i,j]=lat[i]
x, y = m(lon_grid, lat_grid) # compute map proj coordinates. 
clevs= range(0,6000,500)  # 等值线间隔为10，从0开始
cs = m.contourf(x,y,hight[0,:,:],clevs,cmap='binary',extend='both')
cbar = m.colorbar(cs,location='bottom',pad="5%") # 色标
cbar.set_label('m',fontsize=16) # 色标下面标注单位
x, y = m(etplon, etplat) 
m.plot(x,y,c='w')
x, y = m(93, 33) 
plt.text(x,y,'ETP',fonttp)
x, y = m(wtplon, wtplat) 
m.plot(x,y,c='w')
x, y = m(83, 33) 
plt.text(x,y,'WTP',fonttp)
#
x, y = m(neclon, neclat) 
m.plot(x,y,c='k')
x, y = m(121, 45.5) 
plt.text(x,y,'NEC',fontec)
#
x, y = m(npclon, npclat) 
m.plot(x,y,c='k')
x, y = m(114, 38) 
plt.text(x,y,'NPC',fontec)
#
x, y = m(mlyrlon, mlyrlat) 
m.plot(x,y,c='k')
x, y = m(113.5, 31.5) 
plt.text(x,y,'MLYR',fontec)
#
x, y = m(prdlon, prdlat) 
m.plot(x,y,c='k')
x, y = m(112, 23.5) 
plt.text(x,y,'PRD',fontec)
#
plt.show() #显示图像
plt.savefig(dirout+'Land_10min.png',dpi=300)  #保存图片，这里保存为pdf，常用格式矢量格式 eps,ps等都是支持的
# 非矢量格式，png，jpg，tiff等
plt.close() # 关闭画图窗口

