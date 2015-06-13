#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 18:28:33 2015

@author: jhchen
"""
import urllib
import re
dirout='E:/Data/WebAQI/'
chinese=[]
pinyin=[]
def getHtml(url):
    page = urllib.urlopen(url)
    html = page.read()
    return html
def getlinenum(strlines,ns,ne,re):
    ll=-1
    for i in range(ns,ne):
        lnstr=strlines[i].strip()
        if lnstr==re:
            ll=i
            break
    return ll
def getname(lstr):
    l=len(lstr)
    cnm='' ; icnm=-1
    pnm='' ; ipnm=-1
    for i in range(1,l-1):
        sr1=lstr[i]
        sr0=lstr[i-1]
        sr2=lstr[i+1]
        if ipnm==1 :
            pnm=pnm+sr1 
            if sr2==r'"' :
                ipnm=-1
        if icnm==1 :
            cnm=cnm+sr1 
            if sr2==r'<' :
                icnm=-1
                break
        if sr1==r'/' and sr0==r'"' :
            ipnm=1
        if sr1==r'>' and sr0==r'"' :
            icnm=1
    return cnm,pnm
#
url='http://www.pm25.in/'
html = getHtml(url)
tmp=[]
tmp.append(html.split('\n'))
htmllines=[]
for sr in tmp[0]:
    htmllines.append(sr)   
nline=len(htmllines)
reg=r'<div><b>A.</b></div>'
k1=getlinenum(htmllines,0,nline,reg)
reg=r'</div> <!-- end of cities -->'
k2=getlinenum(htmllines,k1,nline,reg)
for i in range(k1,k2):
    linestr=htmllines[i].strip()
    reg='href='
    xstr=re.findall(reg,linestr)
    if xstr:
        cnm,pnm=getname(linestr)
        chinese.append(cnm)
        pinyin.append(pnm)
nl=len(chinese)
fpath=dirout+'AllCityName.txt'
fout=open(fpath,'w')
item='Index '+' Chinese  '+ 'Pinyin'+'\n'
fout.write(item)
for i in range(0,nl):
    item="%d "%i + chinese[i]+' '+ pinyin[i]+'\n'
    fout.write(item)
fout.close()


