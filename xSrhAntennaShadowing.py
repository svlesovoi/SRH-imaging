#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:36:23 2018

@author: sergey
"""

import sys;
from PyQt5 import QtGui, QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import pylab as PL
import numpy as NP
import srhArray

class SrhPairCanvas(FigureCanvas):
    def arcsec_format(self, pix, pos):
        return '%3.1f' % (((pix - self.parent.M//2) * NP.rad2deg(self.parent.radPerPix)));

    def meter_format(self, pix, pos):
        meterPerPix = 0.05
        return '%3d' % (((pix - 125.) * meterPerPix));

    def uv_format(self, pix, pos):
        return '%3d' % ((pix - self.parent.M//2) * self.uvPerPix);

    def __init__(self, parent, width=5, height=4, dpi=100):
        self.parent = parent
        self.uvSize = 2*parent.N
        self.uvPerPix = parent.uvPerPix
        self.radSize = 2*parent.N
        self.radPerPix = parent.radPerPix
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.subplot = fig.add_subplot(111)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def setArcFormat(self):
        self.subplot.xaxis.set_major_formatter(PL.FuncFormatter(self.arcsec_format))
        self.subplot.xaxis.set_major_locator(PL.MultipleLocator(self.parent.M/10))
        self.subplot.yaxis.set_major_formatter(PL.FuncFormatter(self.arcsec_format))
        self.subplot.yaxis.set_major_locator(PL.MultipleLocator(self.parent.M/10))
        self.subplot.set_xlim(200., 400.)
        self.subplot.set_ylim(200., 400.)

    def setMeterFormat(self):
        self.subplot.xaxis.set_major_formatter(PL.FuncFormatter(self.meter_format))
        self.subplot.yaxis.set_major_formatter(PL.FuncFormatter(self.meter_format))

    def setUvFormat(self):
        self.subplot.xaxis.set_major_formatter(PL.FuncFormatter(self.uv_format))
        self.subplot.xaxis.set_major_locator(PL.MultipleLocator(self.parent.M/6))
        self.subplot.yaxis.set_major_formatter(PL.FuncFormatter(self.uv_format))
        self.subplot.yaxis.set_major_locator(PL.MultipleLocator(self.parent.M/6))
        self.subplot.set_xlim(100., 500.)
        self.subplot.set_ylim(100., 500.)

    def imshow(self, array):
        self.subplot.imshow(array, cmap='hot')
        self.draw()
    
class SrhAntennaShadowing(QtWidgets.QMainWindow):
    
    def onCorrect(self, value):
        self.correctBaseline = value
        
    def onCalc(self):
        frequency = float(self.frequency.toPlainText()) * 1e6
        diameter = float(self.diameter.toPlainText())
#        antennaA = int(self.antennaA.toPlainText())
#        antennaB = int(self.antennaB.toPlainText())
        self.fillDish(diameter, frequency)
        declination = NP.deg2rad(float(self.declination.toPlainText()))

        hourAngle = NP.deg2rad(float(self.hourAngle.toPlainText()))
        
        self.beamPattern[:,:] = 0.
        self.apertureCorr[:,:] = 0.
        self.uvPlaneB[:,:] = 0.
        u0 = self.M//2
        v0 = self.M//2
        self.uvPlaneB[u0 - self.N//2:u0 + self.N//2,v0 - self.N//2:v0 + self.N//2] = self.dish
        for antennaA in NP.linspace(49,80,32):
#        for antennaA in NP.linspace(antennaA,antennaA,1):
            print(antennaA)
            for antennaB in NP.linspace(192,177,16):
#            for antennaB in NP.linspace(antennaB,antennaB,1):
                uvw = self.SRH.baseline2uvw(hourAngle, declination, antennaA, antennaB) / (3e8 / frequency)
                baseline = NP.sqrt(uvw[0]**2 + uvw[1]**2)
                baselinePrime = (baseline + diameter)/2
                if (self.correctBaseline and baseline < diameter):
                    u1 = int(NP.round(self.M/2 + uvw[1]/self.uvPerPix*baselinePrime/baseline))
                    v1 = int(NP.round(self.M/2 - uvw[0]/self.uvPerPix*baselinePrime/baseline))
                else:
                    u1 = int(NP.round(self.M/2 + uvw[1]/self.uvPerPix))
                    v1 = int(NP.round(self.M/2 - uvw[0]/self.uvPerPix))
                
                self.uvPlaneA[:,:] = 0.
                self.uvPlaneA[u1 - self.N//2:u1 + self.N//2,v1 - self.N//2:v1 + self.N//2] = self.dish
                self.apertureCorr += self.fftConvolution(self.uvPlaneA, self.uvPlaneB).real
        self.baseline.setText('%2.1f, %2.1f, %2.1f'%(uvw[0], uvw[1], baseline))
        self.beamPattern = self.fftBeam(self.apertureCorr).real
        self.apertureCanvas.imshow(self.apertureCorr.real)
        self.response = self.halfFftConvolution(self.beamPattern, self.qSun).real
        self.patternCanvas.imshow(self.response)
        
    def fillDish(self, diameter, frequency):
        radius = diameter / 2 / (3e8 / frequency) / self.uvPerPix
        for i in range(self.N):
            for j in range(self.N):
                x=i - self.N/2
                y=j - self.N/2
                self.dish[i, j] = 0.
                if (NP.sqrt(x**2 + y**2) < radius):
                    self.dish[i, j] = 1.
                    if (x > 0):
                        self.shadowedDish[i, j] = 1.

    def fillQSun(self, diameter):
        for i in range(self.N):
            for j in range(self.N):
                x=i - self.N/2
                y=j - self.N/2
                self.qSun[i, j] = 0.
                if (NP.sqrt(x**2 + y**2) < diameter/2):
                    self.qSunSmall[i, j] = 1.
        self.qSun[self.M//2 - self.N//2:self.M//2 + self.N//2,self.M//2 - self.N//2:self.M//2 + self.N//2] = self.qSunSmall
        self.qSun[self.M//2, self.M//2] += 10.
        self.qSun[self.M//2+10:self.M//2+50, self.M//2-23:self.M//2-20] = 0.7

    def shift2D(self, arr):
        return NP.roll(arr, (arr.shape[0]//2, arr.shape[0]//2), axis=(0,1))
                       
    def fftConvolution(self, arr1, arr2):
        size = arr1.shape[0]
        corr = NP.roll((NP.fft.ifft2(NP.fft.fft2(arr1) * NP.conjugate(NP.fft.fft2(arr2)))),(size//2,size//2),axis=(0,1))
        corr1 = NP.flip(corr,0)
        corr2 = NP.flip(corr1,1)
        corr2 = NP.roll(corr2,(1,1),(0,1))
        return corr + corr2
        
    def halfFftConvolution(self, arr1, arr2):
        size = arr1.shape[0]
        return NP.roll((NP.fft.ifft2(NP.fft.fft2(arr1) * NP.conjugate(NP.fft.fft2(arr2)))),(size//2,size//2),axis=(0,1))
        
    def fftBeam(self, uvArr):
        return self.shift2D(NP.fft.fft2(self.shift2D(uvArr)))
        
    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self,parent)
        
        self.correctBaseline = False
        self.N = 100
        self.M = 6 * self.N
        self.uvPerPix = 12.0
        self.radPerPix = 1./(self.uvPerPix * self.M)
        self.dish = NP.zeros((self.N, self.N))
        self.shadowedDish = NP.zeros((self.N, self.N))
        self.SRH = srhArray.SrhArray()
        self.uvPlaneA = NP.zeros((self.M, self.M),dtype='complex')
        self.uvPlaneB = NP.zeros((self.M, self.M),dtype='complex')
        self.beamPattern = NP.zeros((self.M, self.M))
        self.apertureCorr = NP.zeros((self.M, self.M))
        self.qSun = NP.zeros((self.M, self.M))
        self.qSunSmall = NP.zeros((self.N, self.N))
        
        self.fillQSun(NP.deg2rad(2160/3600)/self.radPerPix)
        
        self.setGeometry(10,10,1350,800)
        self.hourAngleLabel = QtWidgets.QLabel('Hour Angle', self)
        self.declinationLabel = QtWidgets.QLabel('Declination', self)
        self.antennaALabel = QtWidgets.QLabel('AntennaA',self)
        self.antennaBLabel = QtWidgets.QLabel('AntennaB',self)
        self.frequencyLabel = QtWidgets.QLabel('Frequency [MHz]', self)
        self.diameterLabel = QtWidgets.QLabel('Dish diameter', self)
        self.calcButton = QtWidgets.QPushButton('Calc', self)
        self.baselineLabel = QtWidgets.QLabel('Baseline', self)
        self.correctButton = QtWidgets.QPushButton('Correct', self)
        self.correctButton.setCheckable(True)
        
        self.hourAngle = QtWidgets.QTextEdit(self)
        self.declination = QtWidgets.QTextEdit(self)
        self.frequency = QtWidgets.QTextEdit(self)
        self.diameter = QtWidgets.QTextEdit(self)
        self.antennaA = QtWidgets.QTextEdit(self)
        self.antennaB = QtWidgets.QTextEdit(self)
        self.baseline = QtWidgets.QLabel('0, 0, 0,', self)

        x0 = 10
        y0 = 10
        dX = 90
        dW = 85
        dY = 15
        dH = 25

        self.hourAngleLabel.setGeometry(x0,y0,dW,20)        
        self.declinationLabel.setGeometry(x0 + dX,y0,dW,20)
        self.frequencyLabel.setGeometry(x0 + 2*dX,y0,dW,20)
        self.diameterLabel.setGeometry(x0 + 3*dX,y0,dW,20)
        self.antennaALabel.setGeometry(x0 + 4*dX,y0,dW,20)
        self.antennaBLabel.setGeometry(x0 + 5*dX,y0,dW,20)
        self.calcButton.clicked.connect(self.onCalc)
        self.correctButton.clicked.connect(self.onCorrect)
        
        self.hourAngle.setGeometry(x0,y0 + dY,dW,dH)        
        self.declination.setGeometry(x0 + dX,y0 + dY,dW,dH)
        self.frequency.setGeometry(x0 + 2*dX,y0 + dY,dW,dH)
        self.diameter.setGeometry(x0 + 3*dX,y0 + dY,dW,dH)
        self.antennaA.setGeometry(x0 + 4*dX,y0 + dY,dW,dH)
        self.antennaB.setGeometry(x0 + 5*dX,y0 + dY,dW,dH)
        self.calcButton.setGeometry(x0 + 6*dX,y0,60,40)
        self.baselineLabel.setGeometry(x0 + 7*dX,y0,100,40)
        self.baseline.setGeometry(x0 + 7*dX,y0 + dY,100,40)
        self.correctButton.setGeometry(x0 + 8*dX,y0,40,40)
        
        self.hourAngle.setText('30.0')
        self.declination.setText('0.0')
        self.frequency.setText('6000.0')
        self.diameter.setText('1.8')
        self.antennaA.setText('64')
        self.antennaB.setText('192')

        self.apertureCanvas = SrhPairCanvas(self)
        self.apertureCanvas.setUvFormat()
        self.apertureCanvas.setGeometry(x0,100,700,700)
        self.patternCanvas = SrhPairCanvas(self)
        self.patternCanvas.setArcFormat()
        self.patternCanvas.setGeometry(x0 + 700,100,700,700)
        
        self.patternCanvas.imshow(self.qSun)
        
application = QtWidgets.QApplication.instance();
if not application:
    application = QtWidgets.QApplication(sys.argv);
    
if sys.platform == 'linux':
    font = QtGui.QFont()
    application.setFont(QtGui.QFont(font.defaultFamily(),8));

SRHantShad=SrhAntennaShadowing();
SRHantShad.setWindowTitle('SRH antenna shadowing')
SRHantShad.show();
sys.exit(application.exec_());
