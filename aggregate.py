# -*- coding: utf-8 -*-
"""
Created on Thu May 18 13:22:24 2017

@author: Michelle

Editted on 6/21/2021
by Allison Wandling
"""

from location import Location 
import random

class Aggregate(Location):    
    def __init__(self, rad):
        Location.__init__(self)
        self.r = rad
        self.c = 1
        
###############################################################################
    def isOutOfBox(self, parameters):        
        [x, y, z] = Location.GetLocation(self)
        
        return((x < (0.0 + self.r * parameters.particleRadScale)) | (x > (parameters.boxSize - self.r * parameters.particleRadScale)) | \
               (y < (0.0 + self.r * parameters.particleRadScale)) | (y > (parameters.boxSize - self.r * parameters.particleRadScale)) | \
               (z < (0.0 + self.r * parameters.particleRadScale)) | (z > (parameters.boxSize - self.r * parameters.particleRadScale)))
        
###############################################################################
    def isOutTop(self, parameters):
        [x, y, z] = Location.GetLocation(self)
        return(z > (parameters.boxSize - self.r * parameters.particleRadScale))
        
###############################################################################
    def isOutBottom(self, parameters):
        [x, y, z] = Location.GetLocation(self)
        return(z < (0.0 + self.r * parameters.particleRadScale))
        
###############################################################################
    def SetRandomXLocation(self, parameters):
        [x, y, z] = Location.GetLocation(self)
        x = ((parameters.boxSize - (2.0 * self.r * parameters.particleRadScale)) * random.random()) + self.r * parameters.particleRadScale
        Location.SetLocation(self, parameters.latticeConstant, x, y, z)
        
###############################################################################
    def SetRandomYLocation(self, parameters):
        [x, y, z] = Location.GetLocation(self)
        y = ((parameters.boxSize - (2.0 * self.r * parameters.particleRadScale)) * random.random()) + self.r * parameters.particleRadScale
        Location.SetLocation(self, parameters.latticeConstant, x, y, z)
        
###############################################################################
    def SetRandomZLocation(self, parameters):
        [x, y, z] = Location.GetLocation(self)
        z = ((parameters.boxSize - (2.0 * self.r * parameters.particleRadScale)) * random.random()) + self.r * parameters.particleRadScale
        Location.SetLocation(self, parameters.latticeConstant, x, y, z)
        
###############################################################################
    def SetRandomLocation(self, parameters):
        self.SetRandomXLocation(parameters)
        self.SetRandomYLocation(parameters)
        self.SetRandomZLocation(parameters)
        
###############################################################################
    def SetLocation(self, parameters, x, y, z, checkForIsOutOfBox = False):        
        Location.SetLocation(self, parameters.latticeConstant, x, y, z)
        
        if (checkForIsOutOfBox == True):
            if (self.isOutOfBox(parameters)):
                raise ValueError("******\Error:  \nParticle is out of the box\n*****")
                
###############################################################################
    def GetRadius(self):
        return(self.r)
        
###############################################################################
    def GetParticleCount(self):
        return(self.c)
        
###############################################################################
