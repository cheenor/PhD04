#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 09 09:50:12 2015

@author: jhchen
"""
#!/usr/bin/env python
"""
An animated image
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pylab import *

casenm="ETP2D1"
dirin='D:/MyPaper/PhD04/Cases/ETP/20100604_0704/Simulated/animate/'
dirout='D:/MyPaper/PhD04/Cases/ETP/20100604_0704/Simulated/'
varnm=["qc","qa","qb","qr","Total_Cloud_Water"]
filestr=varnm[4]
nt=2880
fig = plt.figure()

#def f(x, y):
#    return np.sin(x) + np.cos(y)
#
#x = np.linspace(0, 2 * np.pi, 120)
#y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
#im = plt.imshow(f(x, y), cmap=plt.get_cmap('jet'))
#it="before"
#print it

timstr="0000"
fpath=dirin+filestr+'_'+casenm+'_'+timstr+'.png'
image = plt.imread(fpath)
im=imshow(image)
plt.axis('off')
i=1
def updatefig(it): #(*args):
#    global x,y
#    print i
#    global it   
    timstr="%04d"%(it) 
    fpath=dirin+filestr+'_'+casenm+'_'+timstr+'.png'
    image = plt.imread(fpath)
    im.set_array(image)
#    global x,y
#    x += np.pi / 15.
#    y += np.pi / 20.
#    im.set_array(f(x,y))
    return im,
#
ani = animation.FuncAnimation(fig, updatefig,frames=nt, interval=5, blit=True)
#mywriter=animation.FFMpegWriter()
ani.save(dirout+casenm+filestr+'.mp4', fps=60)

plt.show()
