# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 00:18:47 2018

@author: Sergey
"""
import numpy as NP

class srhArray():
    def __init__(self, fromFits = None):
	self.uvSize = 256
        if fromFits is None:
            self.uvPlane = NP.zeros((self.uvSize, self.uvSize))
            