# -*- coding: utf-8 -*-
"""
Created on Thu May 18 13:22:24 2017

@author: Michelle
"""
#from datetime import datetime
import os
import sys

     
from calcRoughness import CalculateRoughness
from readFile import FileParameters       
from data import Data       
from aggregate import Aggregate
###############################################################################
def ReadDump(parameters, data, fileName, displayDebug = False):
    dumpFile = open(fileName,'r')  # open file read only
    dumpFile.readline()                 # 'ITEM: TIMESTEP\n')    
    dumpFile.readline()                 # str(parameters.clock)+'\n')    
    dumpFile.readline()                 # 'ITEM: NUMBER OF ATOMS\n')    
    totalAtoms = int(dumpFile.readline())    #  
    #print("aggregates=", totalAtoms)
    dumpFile.readline()                 # 'ITEM: BOX BOUNDS pp pp ss\n') 
    
    line = dumpFile.readline()          # str(0)+' '+str(parameters.boxSize)+'\n')  
    [bitBucket, boxSizeStr] = line.split(' ')
    parameters.boxSize = float(boxSizeStr)
    dumpFile.readline()                 # str(0)+' '+str(parameters.boxSize)+'\n')    
    dumpFile.readline()                 # str(0)+' '+str(parameters.boxSize)+'\n')    
    dumpFile.readline()                 # 'ITEM: ATOMS x y z radius\n') 
              
    if (totalAtoms > 10000):
        raise ValueError("******\nMore then 10K particles. program Terminated\n*****")
        
    data.allParticles = []
    
    for j in range(totalAtoms):
        line = dumpFile.readline()
        
        [xStr, yStr, zStr, rStr] = line.split(' ') 
        x = float(xStr)
        y = float(yStr)
        z = float(zStr)
        r = float(rStr)
        #print("%20.14e %20.14e %20.14e %16.10e" % (x, y, z, r))
        aggregate = Aggregate(r)
        aggregate.SetLocation(parameters, x, y, z)
        data.allParticles.append(aggregate)
            
    dumpFile.close()

    
###############################################################################
def main():
    args = sys.argv[1:]
    print("Preforming roughness calculation")
    parameters = FileParameters()
    
    if len(args) == 0:
        raise ValueError("******\nNo dump files specified. program Terminated\n*****")


    for arg in args:  
        if (os.path.isfile(arg) == False):
            print("File %s does not exit. File Skipped" % (arg))
        else:
            #arg = "./Results/aggPoly_20vol/aggPoly_20vol.dump"
            parameters.baseDir = os.path.dirname(arg) + "/"
            #print("parameters.baseDir=", parameters.baseDir)
            
            [pathPlusBaseFileName, extention] = os.path.splitext(arg)
            
            if (extention != ".dump"):
                print("%s is not a .dump file.  File ignored" % (arg))
            else:
                data = Data()
                
                ReadDump(parameters, data, arg)
                fileName = pathPlusBaseFileName + "_roughness.csv"
                CalculateRoughness(parameters, data, fileName, range(95, 0, -5), False, True)
    
    print("Roughness calculation complete")
           
if __name__== "__main__":
    main()
