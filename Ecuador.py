# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:34:06 2020

@author: Ben
"""

"""
Import necessary packages
"""
import numpy as np
import pandas as pd
import functions as f
import scipy.io
import time
import seaborn as sns; sns.set()
from os import path, mkdir, remove
"""
Declare variables
"""

timeStepLength = 1            #Time step for the simulation in months
inputFileFolder = "C:/Users/Ben/OneDrive - UCL/GVP/"
inputFile = "GVP_Eruption_Results.xls.xlsx"          #File containing infos for the eruptions, from GVP website
outputFolder = "C:/Users/Ben/Desktop/results GVP/test/"
VEI4MAT = "atacazo_vei4.mat"
refVolcano = 'Atacazo'
probThreshold = 0.8     #Threshold of probability for the isopach
thresholdYear = 1842    #Year above which the record is complete
stopYear = 2020         #Year when the records end
startYear = -10000      #Start year of the holocène
cellSize = 0.009        #Size of the cells in degrees of Lat/Lon. DON'T CHANGE UNLESS IT FITS THE MATRIX FROM TEPHRAPROB
northLimit = None
southLimit = None
eastLimit = None
westLimit = -85
limits = [northLimit, southLimit, eastLimit, westLimit]

surfaceC = "yes"        # Compute the accumulation of C in the surface soil. "yes" if it needs to be computed. "no" if it doesn't, or "fuck off",
                        # or anything except "yes" actually.
    
mode = "sequential"     # Mode of inference of the eruptions. "stochastic" does a full stochastic approach
                        # for VEI4-5-6, "mixed" uses the historical data for VEI 5 and 6 and stochastic for
                        # VEI4. 
                        # "sequential" takes a more sequential approach : first computing
                        # if there's an eruption, then the VEI, then attribution to a volcano. "sequential" uses the stochastic
                        # approach.



if not path.exists(outputFolder):
            mkdir(outputFolder)
if path.exists(outputFolder + "parameters.txt"):
    remove(outputFolder + "parameters.txt")
paramLog = open(outputFolder + "parameters.txt", 'a')
paramLog.write("timeStepLength = " + str(timeStepLength) + "\n")
paramLog.write("inputFile = " + str(inputFile) + "\n")
paramLog.write("VEI4MAT = " + str(VEI4MAT) + "\n")
paramLog.write("refVolcano = " + str(refVolcano) + "\n")
paramLog.write("cellSize = " + str(cellSize) + "\n")
paramLog.write("probThreshold = " + str(probThreshold) + "\n")
paramLog.write("thresholdYear = " + str(thresholdYear) + "\n")
paramLog.write("stopYear = " + str(stopYear) + "\n")
paramLog.write("startYear = " + str(startYear) + "\n")
paramLog.write("limits = " + str(limits) + "\n")
paramLog.write("surfaceC = " + str(surfaceC) + "\n")
paramLog.write("mode = " + str(mode) + "\n")
paramLog.close()


"""
Import data
"""
fullPath = inputFileFolder + inputFile

data = pd.read_excel(fullPath, sheet_name="Eruption List")

mat = scipy.io.loadmat(inputFileFolder + VEI4MAT)

atacazoVEI4 = np.array(mat['atacazo_vei4'])
atacazoVEI5 = np.array(mat['atacazo_vei5'])
atacazoVEI6 = np.array(mat['atacazo_vei6'])

atacazoVEI4, atacazoVEI5, atacazoVEI6 = f.convert_VEI_to_latlon(atacazoVEI4, atacazoVEI5, atacazoVEI6, refVolcano, data)

refVEI = [atacazoVEI4,atacazoVEI5,atacazoVEI6]

"""
counterStart = time.perf_counter()

mode = "sequential"
outputFolder = "C:/Users/Ben/Desktop/results GVP/1842/OutputSequential/"
if not path.exists(outputFolder):
            mkdir(outputFolder)
if path.exists(outputFolder + "parameters.txt"):
    remove(outputFolder + "parameters.txt")
paramLog = open(outputFolder + "parameters.txt", 'a')
paramLog.write("timeStepLength = " + str(timeStepLength) + "\n")
paramLog.write("inputFile = " + str(inputFile) + "\n")
paramLog.write("VEI4MAT = " + str(VEI4MAT) + "\n")
paramLog.write("refVolcano = " + str(refVolcano) + "\n")
paramLog.write("probThreshold = " + str(probThreshold) + "\n")
paramLog.write("cellSize = " + str(cellSize) + "\n")
paramLog.write("thresholdYear = " + str(thresholdYear) + "\n")
paramLog.write("stopYear = " + str(stopYear) + "\n")
paramLog.write("startYear = " + str(startYear) + "\n")
paramLog.write("limits = " + str(limits) + "\n")
paramLog.write("surfaceC = " + str(surfaceC) + "\n")
paramLog.write("mode = " + str(mode) + "\n")
paramLog.close()
for i in range(100):
    counter0 = time.perf_counter()
    data = f.apply_coord_constraints(data, limits)
    refZone = f.get_ref_zone(data, refVolcano)
    probabilities = f.get_prob(data, startYear, stopYear, thresholdYear,timeStepLength, mode)
    eruptions = f.get_stoch_eruptions(data, probabilities, startYear, stopYear, thresholdYear, refZone, mode, timeStepLength)
    fileList = f.create_vei_files(inputFileFolder, refVolcano, eruptions, refVEI, refZone)
    grid, minLat, minLon = f.create_grid(inputFileFolder,fileList, cellSize, outputFolder, probThreshold)
    f.add_eruptions_to_grid(inputFileFolder,fileList, eruptions, grid, probThreshold, minLat, minLon)
    carbonGrid, logC = f.get_carbon_grid(grid, startYear, stopYear, surfaceC)
    count = f.save_results(outputFolder, carbonGrid, logC, eruptions)
    counter1 = time.perf_counter()
    print(str(count) + ": " + str(counter1-counter0) + " secondes")


"""


"""
Functions calls
"""

data = f.apply_coord_constraints(data, limits)

refCoords = f.get_ref_coords(data, refVolcano)

counter0 = time.perf_counter()

#Compute the probabilities for each volcano of erupting in a given time interval
probabilities = f.get_prob(data, startYear, stopYear, thresholdYear,timeStepLength, mode)
counter1 = time.perf_counter()
print('Probabilities :' + str(counter1-counter0))

#Get a list of eruptions for the volcanoes, based on the previously computed probabilities,
#on the full time interval defined between startYear and stopYear
eruptions = f.get_stoch_eruptions(data, probabilities, startYear, stopYear, thresholdYear, mode, timeStepLength)
counter2 = time.perf_counter()
print('Eruptions :' + str(counter2-counter1))

#Create VEI files which are lat/lon transpositions of the reference VEI file, and a reference fileList
#If the VEI files already exists, just generates the fileList for the sake of speed
fileList = f.create_vei_files(inputFileFolder, refVolcano, data, refVEI, refCoords)
counter3 = time.perf_counter()
print('Filelist :' + str(counter3-counter2))

#Create a grid containing the coordinates and in which the tephra deposits will be stored
grid, minLat, minLon = f.create_grid(inputFileFolder,fileList, cellSize, outputFolder, probThreshold)
counter4 = time.perf_counter()
print('grid :' + str(counter4-counter3))

#Append years of eruptions to the grid
f.add_eruptions_to_grid(inputFileFolder,fileList, eruptions, grid, probThreshold, minLat, minLon, cellSize)
counter5 = time.perf_counter()
print('gridErupt :' + str(counter5-counter4))

#Computes the accumulation of carbon for each cell from the startYear to the last eruption.
#Doesn't compute the surface carbon.
carbonGrid, logC = f.get_carbon_grid(grid, startYear, stopYear, surfaceC)
counter6 = time.perf_counter()
print('carbonGrid :' + str(counter6-counter5))

count = f.save_results(outputFolder, carbonGrid, logC, eruptions)


print('Total (' + str(count) + ') : ' + str(counter6-counter0))



"""
import matplotlib.pyplot as plt
import math
carbonDF = pd.DataFrame(carbonGrid)

plt.imshow(carbonDF, cmap=plt.cm.Reds, origin='lower', interpolation=None, extent=[math.ceil(minLon/cellSize)*cellSize + carbonDF.shape[1]*1000,math.ceil(minLon/cellSize)*cellSize,
                                               math.ceil(minLat/cellSize)*cellSize, math.ceil(minLat/cellSize)*cellSize + carbonDF.shape[0]*1000])
plt.colorbar()
plt.xlabel("Longitude (mètres)")
plt.ylabel("Latitude (mètres)")
plt.grid(False)
"""
