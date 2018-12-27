# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 00:18:47 2018

@author: Sergey
"""
import numpy as NP
from BadaryRAO import BadaryRAO

class srhArray():
    def __init__(self, fromFits = None):
        self.RAO = BadaryRAO()
        self.sizeOfUv = 256
        self.uvLcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex);
        self.uvRcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex);
        if fromFits is None:
            self.uvPlane = NP.zeros((self.uvSize, self.uvSize))
            