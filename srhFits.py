# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 00:18:47 2018

@author: Sergey
"""
class srhFits():
    def __init__():
	self.isOpen = False

    def open(self, path)    
        try:
            self.hduList = fits.open(path)
            self.isOpen = True
            self.antennaNumbers = self.hduList[2].data['ANTENNA'][0]
            self.antennaA = self.hduList[2].data['ANTA'][0]
            self.antennaB = self.hduList[2].data['ANTB'][0]

            self.freqList = self.hduList[1].data['FREQUENCY'];
            self.freqListLength = self.freqList.size;
            self.dataLength = self.hduList[1].data['TIME'].size // self.freqListLength;
#            self.freqTime = NP.reshape(self.hduList[1].data['TIME'],(self.freqListLength,self.dataLength));
            self.freqTime = self.hduList[1].data['TIME']
            self.visListLength = self.hduList[1].data['VIS_LCP'].size // self.freqListLength // self.dataLength;
            self.visLcp = NP.reshape(self.hduList[1].data['VIS_LCP'],(self.freqListLength,self.dataLength,self.visListLength));
            self.visRcp = NP.reshape(self.hduList[1].data['VIS_RCP'],(self.freqListLength,self.dataLength,self.visListLength));
            self.ampLcp = NP.reshape(self.hduList[1].data['AMP_LCP'],(self.freqListLength,self.dataLength,self.antennaNumbers.size));
            self.ampRcp = NP.reshape(self.hduList[1].data['AMP_RCP'],(self.freqListLength,self.dataLength,self.antennaNumbers.size));
            self.visLcp.real = NP.sin(NP.pi/2*self.visLcp.real);
            self.visRcp.real = NP.sin(NP.pi/2*self.visRcp.real);
            self.visLcp.imag = NP.sin(NP.pi/2*self.visLcp.imag);
            self.visRcp.imag = NP.sin(NP.pi/2*self.visRcp.imag);

        except FileNotFoundError:
            print('File %s  not found'%name);
	