#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 28 20:19:28 2014

@author: Chenjh
"""
import matplotlib.pyplot as plt
def plot_contour_xdate(ax,xdat,ydat,zdat,tmlab,dx):
    ax=plt.contour(xdat,ydat,zdat,vmin=0,color='red',linestyle='-')
    ax=plt.contour(xdat,ydat,zdat, vmax=0,color='green',linestyle=':')
    xlen=len(xdat)
    ax.set_xticks(range(0,xlen,dx))
    xticklabels = [tmlab[nn] for nn in range(0,xlen,dx)]      
    ax.set_xticklabels(xticklabels, size=10)