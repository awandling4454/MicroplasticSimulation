# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 13:19:59 2021

@author: allis
"""
from aggregate import Aggregate   
import numpy as np  
import math
from location import Location

from utilities import RecalculateBoxSize

from data import Data 

from runParameters import kappa
from runParameters import EDLConstants
from runParameters import allParticles
from runParameters import latticeConstant

from readFile import FileParameters



##############################################################################
def CreateParticles(Parameters, Data):
    
    for index in range(Parameters.particleCount):
        agg = Aggregate(Parameters.particleRad)
        
        if Parameters.graphForcesTest == False:
           loc = agg.SetRandomLocation(Parameters)
           Data.allParticles.append(loc) #puts particle locations into a list
           Data.allParticlesRadius.append(Parameters.particleRad)
        # else:
        #    agg.CheckLocation(Parameters, x, 0, 0)
            #x += Parameters.spacing   
        
    
    if Parameters.graphForcesTest == True:
        agg = Aggregate(FileParameters.particleRad)
        loc= agg.CheckLocation(latticeConstant, 0, 0, 0)
        Data.allParticles.append(agg)
        x = (Parameters.particleRad + Parameters.stericLayerThickness) * 2 + Parameters.spacing
   
    else:
        RecalculateBoxSize(Parameters, Data,allParticles, False)
    
    #[nudges, failed] = NudgeParticles(parameters, data, parameters.particleRad * .5, parameters.particleRad * .25)
    return([Data.allParticles,Data.allParticlesRadius])


#############################################################################
def CalcDistance(i, j):
    [x1, y1, z1] = i.GetLocation()
    [x2, y2, z2] = j.GetLocation()
    distance=math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    return(distance)
              
###############################################################################
def CalcSumRadiusMatrix(parameters, allParticlesRadius):
    rLst = np.full(len(allParticles),0) # dtype=parameters.doMathIn)
    for k in range(len(allParticles)):
        rLst[k] = allParticles[k].r
    #rLst = np.asarray(allParticlesRadius, dtype=np.longdouble) #turns allParticles lsit into an array
    #sumRadiusMatrix = np.add.outer(rLst, rLst) #creates a matrix of each radius added to another radius
    return(rLst)

##############################################################################
def CalcProductRadiusMatrix(parameters, allParticlesRadius):
    #multiplies the two radii together
    rLst = np.empty(len(allParticles), dtype = float, order = 'C') # dtype=parameters.doMathIn)
    for k in range(len(allParticles)):
      rLst.append(k,allParticles[k])
    #rLst = np.asarray(allParticlesRadius, dtype=np.float64)
    #productRadiusMatrix = np.multiply.outer(rLst, rLst)
    return(rLst)

##############################################################################
def CalcRadDivisionMatrix(parameters, allParticlesRadius):
        #(ri*rj)/(ri+rj)
        sumRadiusMatrix = CalcSumRadiusMatrix(parameters, allParticlesRadius)
        productRadiusMatrix = CalcProductRadiusMatrix(parameters, allParticlesRadius)
        
        radDivisionMatrix = np.divide(productRadiusMatrix,sumRadiusMatrix)
        return(radDivisionMatrix)            
        

###############################################################################
def CalcDistanceMatrix(parameters, allParticles):   
    # caclculates the distance of each particle center from center
    #xLst = np.full(len(allParticles), 0.0)# dtype=parameters.doMathIn)
    #yLst = np.full(len(allParticles), 0.0) # dtype=parameters.doMathIn)
    #zLst = np.full(len(allParticles), 0.0) # dtype=parameters.doMathIn)
    
    xLst = np.empty(0, dtype = float, order = 'C')
    yLst = np.empty(0, dtype=float, order = 'C')
    zLst = np.empty (0,dtype = float, order = 'C')
    
    for k in range(len(allParticles)):
        xyzParticles = allParticles[k]
        xLst= np.insert(xLst,k,xyzParticles[0])
        yLst = np.insert(yLst,k,xyzParticles[1])
        zLst = np.insert(zLst,k,xyzParticles[2])
        
    xSqDistMatrix = np.square(np.subtract.outer(xLst, xLst))
    ySqDistMatrix = np.square(np.subtract.outer(yLst, yLst))
    zSqDistMatrix = np.square(np.subtract.outer(zLst, zLst))
            
    distanceMatrix=   np.sqrt(np.add(xSqDistMatrix, np.add(ySqDistMatrix, zSqDistMatrix))) 
    return(distanceMatrix)

###############################################################################
def CalcSSDistanceMatrix(parameters, totalParticles,totalParticlesRadius):  
    #calculates the distance from the surface to surface of each particle
    #if number is negative that means that the particles are overlapping
    partRadSum = CalcSumRadiusMatrix(parameters, totalParticlesRadius)
    distanceMatrix=  CalcDistanceMatrix(parameters, totalParticles)
    SSDistanceMatrix = np.subtract(distanceMatrix, partRadSum)
    return(SSDistanceMatrix)    
###############################################################################
def CalcNaturalLogMatrix(parameters, allParticles, allParticlesRadius):
    #ln(1+exp(-kH))
    SSDistanceMatrix = CalcSSDistanceMatrix(parameters, allParticles, totalParticlesRadius)
    kappaSSDistanceMatrix = np.multiply(-kappa,SSDistanceMatrix)
    #creates numbers that are large to be exponentiated so change the units
    kappaSSDistanceMatrix_Mega = np.divide(kappaSSDistanceMatrix,1e6)
    #changes units to Megamenters
    exponentiated = np.exp(kappaSSDistanceMatrix_Mega, dtype = np.float64)
    plusOne = np.add(1e-6,exponentiated) #1e-6 because equation multiplied by 1e6
    NaturalLogMatrix = np.log(plusOne)
    
    return(NaturalLogMatrix)

##############################################################################
def EDLMatrix (parameters, allParticles, allParticlesRadius):
    #completing the full EDL equation 
    #only need to worry about kr>5 because the concentration would have to be
    #unlikely small for kr<5
    NaturalLogMatrix = CalcNaturalLogMatrix(parameters, allParticles, allParticlesRadius)
    radDivisionMatrix = CalcRadDivisionMatrix(parameters, allParticlesRadius)
    radDiviisionMatrix = np.multiply(radDivisionMatrix,1e-6) #because whole equation is multiplied by 1e-6
    NaturalLogAndRad =  np.multiply(NaturalLogMatrix,radDivisionMatrix) 
    EDLMatrix= np.multiply(EDLConstants,NaturalLogAndRad)
    EDLMatrix = np.multiply(EDLMatrix,1e6) #makes all the equation normal again
    #gives an answer in V #the more positive the number then the less likely to agglomerate
    
    return(EDLMatrix)

########################################################################
def CalcVdWMatrixWorker(parameters, partRadProd, partRadSum, partRadDiff, distanceSqMatrix, clearNaNErrors = False):
    A = parameters.A
    
    temp1 = np.multiply(2.0, partRadProd)   
    temp2 = np.square(partRadSum)
    temp3 = np.square(partRadDiff)
    temp4 = distanceSqMatrix

    x = np.where(np.subtract(temp4, temp2) == 0)
    y = np.where(np.subtract(temp4, temp3) == 0)
    
    if (len(x[0]) != 0): raise ValueError("******\Error:  \nabout to divide by 0.\n*****")
    if (len(y[0]) != 0): raise ValueError("******\Error:  \nabout to divide by 0.\n*****")
    
    temp5 = np.divide(temp1, np.subtract(temp4, temp2))
    temp6 = np.divide(temp1, np.subtract(temp4, temp3))
    
    
    temp7 = np.divide(np.subtract(temp4, temp2), np.subtract(temp4, temp3))
    
    # for this test some particles can overlap... We don't care about them
    # so just adjust the data so it works
    if clearNaNErrors:
        temp7[temp7 <= 0] = .1
        
    temp7 = np.log(temp7)
    
    VdW = np.multiply(A / 6.0, np.add(temp5, np.add(temp6, temp7)))
    VdW = np.negative(VdW)
    
    return(VdW)
    
###############################################################################
def CalcVdWMatrix(parameters, allParticles, distanceSqMatrix, clearNaNErrors = False):
    # Van der Waals equation
    # Where:
    #    R = center-to-center distance of two particles
    #    r_i = radius of particle 1
    #    r_j = radius of particle 2
    #    A = hamaker constant.  This is caculated in readFile.py and used here
    #
    # VdW= -(A/6) * \
    #           (((2*r_i*r_j) / (R**2-((r_i+r_j)**2))) + \
    #            ((2*r_i*r_j) / (R**2-((r_i-r_j)**2))) + \
    #          ln(R**2-((r_i+r_j)**2) / (R**2-((r_i-r_j)**2))))
    #
    # Now lets build code.  Build common parts first
    #
    # temp1 = 2*r_i*R_j
    # temp2 = (r_i+r_j)**2)
    # temp3 = (r_i-r_j)**2)
    # temp4 = R**2
    #
    # Substuting you get:
    # 
    # VdW= -(A/6) * \
    #           (((temp1) / (temp4-temp2)) + \
    #            ((temp1) / (temp4-temp3)) + \
    #          ln((temp4-temp2) / (temp4-temp3)))
    #
    # temp5 = ((temp1) / (temp4-temp2))
    # temp4 = ((temp1) / (temp4-temp3))
    # temp7 = ln((temp4-temp2) / (temp4-temp3)) 
    #
    # Final Subtuting:
    #
    # VdW= -(A/6) * (temp5 + temp6 + temp7)
    #
    # Since the diagnal is zero we will divide by zero, so ignore the exception  

    partRads = [part.r for part in allParticles]
    
    partRadDiff = np.subtract.outer(partRads, partRads)
    partRadSum = np.add.outer(partRads, partRads) 
    partRadProd = np.multiply.outer(partRads, partRads)
        
    VdW = CalcVdWMatrixWorker(parameters, partRadProd, partRadSum, partRadDiff, distanceSqMatrix, clearNaNErrors)

    return(VdW)
##############################################################################
[totalParticles, totalParticlesRadius] = CreateParticles(FileParameters, Data)
EDL = EDLMatrix(FileParameters, totalParticles, totalParticlesRadius)
distanceMatrix = CalcDistanceMatrix(FileParameters, privateParticleList)
np.fill_diagonal(distanceMatrix, 1.0)  # this will keep particles away from self
distanceSqMatrix = np.square(distanceMatrix) 
VDW = CalcVdWMatrix(FileParameters, totalParticles, distanceSqMatrix)
