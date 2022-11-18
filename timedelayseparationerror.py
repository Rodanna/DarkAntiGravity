#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 08:19:20 2022

@author: annarodriguez
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

arctomicrorad = 4.848137
pixelpermicrorad = 80/99.16
res = 240
rmax2 = 300 

x1 = np.zeros((5,5))  
y1 = np.zeros((5,5))
Timedelay = np.zeros((5,4,10))
extremum = np.zeros((5,5,2))

p = 0
for i in range (5,9):
    Timedelay[0][p] = np.loadtxt(f'Timedelay0.{i}(nr=700).txt')
    Timedelay[1][p] = np.loadtxt(f'Timedelay0.{i}(nr=1000).txt')
    Timedelay[2][p] = np.loadtxt(f'Timedelay0.{i}(nr=1500).txt')
    Timedelay[3][p] = np.loadtxt(f'Timedelay0.{i}(nr=2000).txt')
    Timedelay[4][p] = np.loadtxt(f'Timedelay0.{i}(nr=2238).txt')
    p +=1

extremum[0] = np.loadtxt('extremum700.txt')
extremum[1] = np.loadtxt('extremum1000.txt')  
extremum[2] = np.loadtxt('extremum1500.txt')  
extremum[3] = np.loadtxt('extremum2000.txt')  
extremum[4] = np.loadtxt('extremum2238.txt')    

for m in range(0,5):
    for kik in range (0,len(extremum)):
        y1[m][kik] = int(extremum[m][kik][0])
        x1[m][kik] = int(extremum[m][kik][1])

error700 = []
error1000 = []
error1500 = []
error2000 = []
for i in range (0,4): #differnet redshifts
    for j in range(10): #timdedelays
        for a in range (0,4): #different nr of Bessel roots
            if a == 0:
                error700.append(Timedelay[a][i][j]/Timedelay[4][i][j])
            if a == 1:
                error1000.append(Timedelay[a][i][j]/Timedelay[4][i][j])
            if a == 2:
                error1500.append(Timedelay[a][i][j]/Timedelay[4][i][j])
            if a == 3:
                error2000.append(Timedelay[a][i][j]/Timedelay[4][i][j])

TimeError = np.zeros((4,40))
TimeError[0] = sorted(np.asarray(error700))
TimeError[1] = sorted(np.asarray(error1000))
TimeError[2] = sorted(np.asarray(error1500))
TimeError[3] = sorted(np.asarray(error2000))
timemedian = np.zeros((4))
titles = [700,1000,1500,2000]

for i in range (0,4):
    timemedian[i] = TimeError[i][20]
    plt.figure()
    plt.grid()
    plt.hist(TimeError[i],bins = 100,rwidth = 0.6)
    plt.vlines(timemedian[i],0,7,'r',linestyles='dotted')
    plt.vlines(TimeError[i][int(len(TimeError[i])*0.16)],0,7,'r',linestyles='dotted')
    plt.vlines(TimeError[i][int(len(TimeError[i])*0.84)],0,7,'r',linestyles='dotted')
    plt.xlabel('Normalized error in time delay')
    plt.ylabel('Number count')
    plt.savefig(f'ErrorInTimedelay{i}.pdf')
    plt.show()

for i in range (0,4):
    print(i,timemedian[i],TimeError[i][int(len(TimeError[i])*0.16)],TimeError[i][int(len(TimeError[i])*0.84)])


dist0 = []
dist1 = []
dist2 = []
dist3 = []
dist4 = []
distances = np.zeros((5,10))
theta = np.zeros(6)
                                       
k = 0        

for n in  range(0,len(extremum)):
    for t in range(n+1,len(extremum)):
        dist0.append(np.sqrt((x1[0][n]-x1[0][t])**2+(y1[0][n]-y1[0][t])**2))
        dist1.append(np.sqrt((x1[1][n]-x1[1][t])**2+(y1[1][n]-y1[1][t])**2))
        dist2.append(np.sqrt((x1[2][n]-x1[2][t])**2+(y1[2][n]-y1[2][t])**2))
        dist3.append(np.sqrt((x1[3][n]-x1[3][t])**2+(y1[3][n]-y1[3][t])**2))
        dist4.append(np.sqrt((x1[4][n]-x1[4][t])**2+(y1[4][n]-y1[4][t])**2))

distances[0] = sorted(np.asarray(dist0)/pixelpermicrorad)
distances[1] = sorted(np.asarray(dist1)/pixelpermicrorad)
distances[2] = sorted(np.asarray(dist2)/pixelpermicrorad)
distances[3] = sorted(np.asarray(dist3)/pixelpermicrorad)
distances[4] = sorted(np.asarray(dist4)/pixelpermicrorad)

DistanceError = np.zeros((4,10))
distmedian = np.zeros(4)

for i in range (0,4):
    DistanceError[i] = sorted(distances[i]/distances[4])
    distmedian[i] = DistanceError[i][5]
    print(i,distmedian[i],DistanceError[i][int(len(DistanceError[i])*0.16)],DistanceError[i][int(len(DistanceError[i])*0.84)])

for i in range (0,4):
    plt.figure()
    plt.title(f'Error of Separations for nr = {titles[i]}')
    plt.grid()
    plt.hist(distances[i]/distances[4],bins = 100,rwidth = 0.6)
    plt.vlines(distmedian[i],0,7,'r',linestyles='dotted')
    plt.vlines(DistanceError[i][2],0,7,'r',linestyles='dotted')
    plt.vlines(DistanceError[i][8],0,7,'r',linestyles='dotted')
    plt.xlabel('Normalized error in Image Separations')
    plt.ylabel('Number Count')
    plt.savefig(f'ErrorInSeparations{i}.png')
    plt.show()

'''
df = pd.DataFrame(dict(nr700=distances[0],nr1000=distances[1],nr1500=distances[2],nr2000 = distances[3],nr2238 = distances[4]))

print(df.to_latex(index=False))  

df = pd.DataFrame(dict(nr700=DistanceError[0],nr1000=DistanceError[1],nr1500=DistanceError[2],nr2000 = DistanceError[3]))

print(df.to_latex(index=False))  '''
'''
for i in range(0,5):
    distances[i] = dist[i*10:(i+1)*10]'''
    
    
'''        
for m in range (0,6): #pixel to microradian
    theta[m] = distances[m]*arctomicrorad
    print(round(theta[m],2))'''
     
'''for m in range (0,6): #pixel to microradian
    theta[m] = rmax2/res*distances[m]
    print(round(theta[m],2))'''