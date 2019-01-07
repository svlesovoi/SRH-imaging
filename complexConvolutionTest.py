# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 04:51:08 2019

@author: Sergey
"""

import numpy as NP
import pylab as PL
import scipy.signal as SIG
import srhArray as SRH

N = 16
M = 40*N

qSun = NP.zeros((M,M))
aperEW = NP.zeros((M,M),dtype = 'complex')
aperS = NP.zeros((M,M),dtype = 'complex')
base = 5

telescope = SRH.SrhArray()
hourAngle = NP.deg2rad(30)
declination = NP.deg2rad(-23)

frequency = 7.5e9
uvPerPix = 8.0
radPerPix = 1./(uvPerPix * M)
qSunRadius = NP.deg2rad(2160/3600)/radPerPix/2
dishRadius = 1.8 / 2 / (3e8 / frequency) / uvPerPix

qSun[:,:] = 0.
for i in range(int(2*qSunRadius)):
    for j in range(int(2*qSunRadius)):
        x = i - qSunRadius
        y = j - qSunRadius
        if (NP.sqrt(x**2 + y**2) < qSunRadius):
            qSun[int(x) + M//2, int(y) + M//2] = 1.

uvQSun = NP.roll(NP.fft.fft2(qSun),(M//2,M//2),axis=(0,1))

corrPlot = NP.zeros(100)
hAdeg = NP.zeros(100)
for hA in range(100):
    hAdeg[hA] = 5 + hA*0.3
    print(hAdeg[hA])
    aperS[:,:] = 0. + 0j
    aperEW[:,:] = 0. + 0j
    for i in range(N):
        for j in range(N):
            x = i - N/2
            y = j - N/2
            if (NP.sqrt(x**2 + y**2) < dishRadius):
                for n in range(16):
                    uvw = telescope.baseline2uvw(NP.deg2rad(hAdeg[hA]), declination, 192 + n, 192.5)/(3e8 / frequency)/uvPerPix
                    aperS[int(NP.round(x + uvw[1] + M/2)), int(NP.round(y + uvw[0] + M/2))] = 1 + 0j
                for n in range(32):
                    uvw = telescope.baseline2uvw(NP.deg2rad(hAdeg[hA]), declination, 49 + n, 64.5)/(3e8 / frequency)/uvPerPix
                    aperEW[int(NP.round(x + uvw[1] + M/2)), int(NP.round(y + uvw[0] + M/2))] = 1 + 0j
    conv_SIG = SIG.fftconvolve(aperEW,aperS,mode='same')
    corrPlot[hA] = NP.abs(conv_SIG*uvQSun).mean()


#PL.imshow(NP.abs(uvQSun))
#PL.imshow(NP.abs(conv_SIG*uvQSun))
#PL.imshow(NP.abs(aperEW + aperS))
PL.plot(hAdeg, corrPlot)

