# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 18:42:15 2020

@author: Ben
"""
import UTMconv as utm
import numpy as np
import pandas as pd
import random as rd
from collections import Counter
import math
from os import path, mkdir,listdir
from pathlib import Path
from os.path import isfile, join
from scipy.integrate import quad


def get_index_positions(list_of_elems, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list

def convert_VEI_to_latlon(VEI4, VEI5, VEI6,refVolcano,data):
    index = get_index_positions(list(data.iloc[:,1]), refVolcano)[0]
    utmCoords = utm.project((data.iloc[index,23],data.iloc[index,22]))
    refZone = utmCoords[0]
    refZoneLetter = utmCoords[1]
    
    VEI4_conv = VEI4.copy()
    for i in range(len(VEI4_conv)):
        latlon = utm.unproject(refZone, refZoneLetter, VEI4_conv[i][0], VEI4_conv[i][1])
        VEI4_conv[i][0] = latlon[0]
        VEI4_conv[i][1] = latlon[1]
        
    VEI5_conv = VEI5.copy()
    for i in range(len(VEI5_conv)):
        latlon = utm.unproject(refZone, refZoneLetter, VEI5_conv[i][0], VEI5_conv[i][1])
        VEI5_conv[i][0] = latlon[0]
        VEI5_conv[i][1] = latlon[1]
        
    VEI6_conv = VEI6.copy()
    for i in range(len(VEI6_conv)):
        latlon = utm.unproject(refZone, refZoneLetter, VEI6_conv[i][0], VEI6_conv[i][1])
        VEI6_conv[i][0] = latlon[0]
        VEI6_conv[i][1] = latlon[1]
        
    return VEI4_conv, VEI5_conv, VEI6_conv
        
            
    
    
    
def apply_coord_constraints(data, coordConstraints):
    if not coordConstraints[0] == None:
        data = data.loc[(data['Latitude'] < coordConstraints[0])]
    if not coordConstraints[1] == None:
        data = data.loc[(data['Latitude'] > coordConstraints[1])]
    if not coordConstraints[2] == None:
        data = data.loc[(data['Longitude'] < coordConstraints[2])]
    if not coordConstraints[3] == None:
        data = data.loc[(data['Longitude'] > coordConstraints[3])]
        
    data = data[(data['VEI'] > 3)]
    data = data.reset_index(drop=True)

    return data

def get_ref_coords(data,refVolcano):
    indexLoc = data[data['Volcano Name'] == refVolcano].index.tolist()[0]
    refLat, refLon = (data.iloc[indexLoc,22], data.iloc[indexLoc,23])

    return (refLat, refLon)


def get_prob(data, startYear, stopYear,thresholdYear,timeStepLength, mode): #Compute the probabilities of eruption (VEI 3 & 4) for each volcano
  
    dataVEI4 = data.loc[(data['VEI'] > 3) & (data['VEI'] < 5) & (data['Start Year'] >= thresholdYear)]
    counterVolcanoesVEI4 = Counter(dataVEI4['Volcano Name'])
    probEruptions = []
    
    timeStepsProbs = (stopYear - thresholdYear) * (12 * timeStepLength)
    
    for i in counterVolcanoesVEI4.keys():
        newLine = [i, counterVolcanoesVEI4[i]/timeStepsProbs,
                                len(dataVEI4.loc[(dataVEI4['VEI'] == 4) 
                                                  & (dataVEI4['Volcano Name'] == i)]), 4]
        probEruptions.append(newLine)
    
    if not mode == "mixed":
        
        timeStepsProbsVEI56 = (stopYear - startYear) * (12 * timeStepLength)
        dataVEI5 = data.loc[(data['VEI'] == 5)]
        counterVolcanoesVEI5 = Counter(dataVEI5['Volcano Name'])
        
        
        for i in counterVolcanoesVEI5.keys():
            newLine = [i, counterVolcanoesVEI5[i]/timeStepsProbsVEI56, len(dataVEI5.loc[(dataVEI5['VEI'] == 5) & (dataVEI5['Volcano Name'] == i)]), 5]
            probEruptions.append(newLine)
            
        dataVEI6 = data.loc[(data['VEI'] == 6)]
        counterVolcanoesVEI6 = Counter(dataVEI6['Volcano Name'])
        
        
        for i in counterVolcanoesVEI6.keys():
            newLine = [i, counterVolcanoesVEI6[i]/timeStepsProbsVEI56, len(dataVEI6.loc[(dataVEI6['VEI'] == 6) & (dataVEI6['Volcano Name'] == i)]), 6]
            probEruptions.append(newLine)
            
        dataVEI4 = data.loc[(data['VEI'] == 4)]
        counterVolcanoesVEI4 = Counter(dataVEI4['Volcano Name'])
        for i in counterVolcanoesVEI4.keys():
            if len(dataVEI4.loc[(dataVEI4['Start Year'] >= thresholdYear) & 
                                (dataVEI4['Volcano Name'] == i)]) == 0:
                if len(dataVEI4.loc[(dataVEI4['Start Year'] < thresholdYear) & 
                                (dataVEI4['Volcano Name'] == i)]) > 0:
                    temp = dataVEI4.loc[(dataVEI4['Start Year'] < thresholdYear) & 
                                (dataVEI4['Volcano Name'] == i)]
                    newLine = [i, counterVolcanoesVEI4[i]/timeStepsProbsVEI56, len(temp), 4]
                    probEruptions.append(newLine)
        
    return probEruptions

def get_utm_zone(lat,long):
    x1 = long + 180
    x2 = x1 / 6
    numZone = round(x2)
    NS = None
    if lat >= 0:
        NS = "N"
    elif lat < 0:
        NS = "S"
    
    return numZone, NS

def get_stoch_eruptions(data, probabilities, startYear, stopYear,threshYear, mode, timeStepLength):
    
    timeSteps = int(round((stopYear - startYear) * (12 / timeStepLength)))
    eruptions = pd.DataFrame(columns=("Volcano","Year","VEI", "Lat","Lon"))    
    k=0
    
    if mode == "mixed":
        
        for j in range(len(probabilities)):
            x = rd.choices([1,0], [probabilities[j][1],1-probabilities[j][1]], k=(timeSteps))
            y = get_index_positions(x,1)
            z = [probabilities[j][3]]*len(y)
            x = [probabilities[j][0]]*len(z)
            locIndex = data.index[data['Volcano Name'] == x[0]].tolist()[0]

            latC = data.iloc[locIndex, 22]
            lonC = data.iloc[locIndex, 23]

            lat = [latC]*len(z)
            lon = [lonC]*len(z)
            
            for i in range(len(x)):
                eruptions.loc[k] = [x[i],math.floor(y[i]/12)+startYear,z[i], lat[i], lon[i]]
                k+=1
                
        dataVEI56 = data.loc[(data['VEI'] > 4) & (data['VEI'] < 7)]
        for i in range(len(dataVEI56)):
    
            latC = dataVEI56.iloc[i,22]
            lonC = dataVEI56.iloc[i,23]

            eruptions.loc[k] = [dataVEI56.iloc[i,1], dataVEI56.iloc[i,8], 
                                int(dataVEI56.iloc[i,5]), latC, lonC]
            k+=1
            
        dataVEI4 = data.loc[(data['VEI'] == 4)]
        counterVolcanoes = Counter(dataVEI4['Volcano Name'])
        for i in counterVolcanoes.keys():
            if len(dataVEI4.loc[(dataVEI4['Start Year'] >= threshYear) & 
                                (dataVEI4['Volcano Name'] == i)]) == 0:
                if len(dataVEI4.loc[(dataVEI4['Start Year'] < threshYear) & 
                                (dataVEI4['Volcano Name'] == i)]) > 0:
                    temp = dataVEI4.loc[(dataVEI4['Start Year'] < threshYear) & 
                                (dataVEI4['Volcano Name'] == i)]
                    for l in range(len(temp)):
                        
                        latC = temp.iloc[l,22]
                        lonC = temp.iloc[l,23]
                        
                        eruptions.loc[k] = [temp.iloc[l,1], temp.iloc[l,8], 
                                int(temp.iloc[l,5]),latC,lonC]
                        k+=1
                        
    elif mode == "stochastic":
        
        for j in range(len(probabilities)):
            x = rd.choices([1,0], [probabilities[j][1],1-probabilities[j][1]], k=(timeSteps))
            y = get_index_positions(x,1)
            z = [probabilities[j][3]]*len(y)
            x = [probabilities[j][0]]*len(z)
            locIndex = data.index[data['Volcano Name'] == probabilities[j][0]].tolist()[0]

            latC = data.iloc[locIndex,22]
            lonC = data.iloc[locIndex,23]
                
            lat = [latC]*len(z)
            lon = [lonC]*len(z)
            
            for i in range(len(x)):
                eruptions.loc[k] = [x[i],math.floor(y[i]/12)+startYear,z[i], lat[i], lon[i]]
                k+=1
                       
                
    elif mode == "sequential":
        
        probaEruptTot = 0
        probaEruptVEI4 = 0
        probaEruptVEI5 = 0
        probaEruptVEI6 = 0
        listVEI4 = []
        listVEI5 = []
        listVEI6 = []
        
        for i in range(len(probabilities)):
            probaEruptTot += probabilities[i][1]
            if probabilities[i][3] == 4:
                probaEruptVEI4 += probabilities[i][1]
                listVEI4.append(probabilities[i])
            if probabilities[i][3] == 5:
                probaEruptVEI5 += probabilities[i][1]
                listVEI5.append(probabilities[i])
            if probabilities[i][3] == 6:
                probaEruptVEI6 += probabilities[i][1]
                listVEI6.append(probabilities[i])

        x = rd.choices([1,0],[probaEruptTot, 1-probaEruptTot], k=timeSteps)
        y = get_index_positions(x,1)
        
        for i in range(len(y)):
            VEI = rd.choices([4,5,6], [probaEruptVEI4, probaEruptVEI5, probaEruptVEI6])

            vol = []
            prob = []
            if VEI == [4]:
                for j in range(len(listVEI4)):
                    vol.append(listVEI4[j][0])
                    prob.append(listVEI4[j][1])
            elif VEI == [5]:
                for j in range(len(listVEI5)):
                    vol.append(listVEI5[j][0])
                    prob.append(listVEI5[j][1])
            elif VEI == [6]:
                for j in range(len(listVEI6)):
                    vol.append(listVEI6[j][0])
                    prob.append(listVEI6[j][1])
                    
            volcano = rd.choices(vol, weights=prob)
            locIndex = data.index[data['Volcano Name'] == volcano[0]].tolist()[0]

            latC = data.iloc[locIndex,22]
            lonC = data.iloc[locIndex,23]
            
            eruptions.loc[k] = [volcano[0], math.floor(y[i]/12)+startYear,VEI[0], latC, lonC]
            k+=1
        
    eruptions = eruptions.sort_values(by=['Year'], ascending=False)
    return(eruptions)

def create_vei_files(inputFileFolder, refVolcano, data, refVEI, refCoords):
    refLat = refCoords[0]
    refLon = refCoords[1]
    
    savePath = inputFileFolder + "VEIs"
    for i in range(len(data)):
        
        if not path.exists(savePath):
            mkdir(savePath)
        if not path.exists(savePath / Path(str(data.loc[i,'Volcano Name']) + "VEI" + str(int(data.loc[i,'VEI'])) + ".csv")):
            latDecal = data.loc[i,'Latitude'] - refLat
            lonDecal = data.loc[i,'Longitude'] - refLon
            if(data.loc[i,'VEI']) == 4:
                matCopy = refVEI[0].copy()
            elif(data.loc[i,'VEI']) == 5:
                matCopy = refVEI[1].copy()
            elif(data.loc[i,'VEI']) == 6:
                matCopy = refVEI[2].copy()
            
            matCopy[:,1] += latDecal
            matCopy[:,0] += lonDecal
            
            np.savetxt((savePath / Path(data.iloc[i,1] + "VEI" + str(int(data.iloc[i,5])) + ".csv")), 
                       matCopy, delimiter=",")
    fileList = []
    uniqueID = (data["Volcano Name"] + "." + data["VEI"].astype(str)).unique()
    for i in range(len(uniqueID)):
        temp = uniqueID[i].split(".")
        fileList.append(temp[0] + "VEI" + temp[1] + ".csv")
    return fileList
    
def create_grid(inputFileFolder, fileList, cellSize, outputFolder, probThreshold):
    if not path.exists(outputFolder):
        mkdir(outputFolder)
    if not path.exists(outputFolder + "grid.npy"):
        minLat = None
        maxLat = None
        minLon = None
        maxLon = None
        savePath = inputFileFolder + "VEIs"
        
        for i in range(len(fileList)):
            temp1 = pd.read_csv(savePath / Path(fileList[i]))
            temp = pd.DataFrame()
            temp = temp1[temp1.iloc[:,2] >= probThreshold]
            if minLat == None:
                minLat = min(temp.iloc[:,1])
            elif minLat > min(temp.iloc[:,1]):
                minLat = min(temp.iloc[:,1])
                
            if maxLat == None:
                maxLat = max(temp.iloc[:,1])
            elif maxLat < max(temp.iloc[:,1]):
                maxLat = max(temp.iloc[:,1])
            
            if minLon == None:
                minLon = min(temp.iloc[:,0])
            elif minLon > min(temp.iloc[:,0]):
                minLon = min(temp.iloc[:,0])
                
            if maxLon == None:
                maxLon = max(temp.iloc[:,0])
            elif maxLon < max(temp.iloc[:,0]):
                maxLon = max(temp.iloc[:,0])
                
        maxLon = cellSize*(round(maxLon/cellSize))
        minLon = cellSize*(round(minLon/cellSize))
        maxLat = cellSize*(round(maxLat/cellSize))
        minLat = cellSize*(round(minLat/cellSize))
        
        mapMatrix = np.empty((int((maxLat - minLat)/cellSize),
                                  int((maxLon - minLon)/cellSize)), dtype=object)
        for i in range(mapMatrix.shape[0]):
            for j in range(mapMatrix.shape[1]):
                mapMatrix[i,j] = []
                mapMatrix[i,j].append(minLat+(i*cellSize))
                mapMatrix[i,j].append(minLon+(j*cellSize))
        np.save(outputFolder + "grid", mapMatrix)
        coords = [minLat,minLon]
        np.save(outputFolder + "coords", coords)
        
    else:
        mapMatrix = np.load(outputFolder + "grid.npy", allow_pickle=True)
        coords = np.load(outputFolder + "coords.npy")
        minLat = coords[0]
        minLon = coords[1]
    
    return mapMatrix, minLat, minLon

def add_eruptions_to_grid(inputFileFolder,fileList, eruptions, grid, probThreshold, minLat, minLon, cellSize):
    veis = dict()
    savePath = inputFileFolder + "VEIs"
    for i in range(len(fileList)):
        veis[fileList[i].split("/")[-1].rsplit(".",1)[0]] = pd.read_csv(savePath / Path(fileList[i]))
        
    for i in range(len(eruptions)):
        vei = eruptions.iloc[i,0] + "VEI" + str(eruptions.iloc[i,2])
        values = veis[vei]
        valuesThresh = values[(values.iloc[:,2] >= probThreshold)]
        
        for j in range(len(valuesThresh)):
            grid[int(((cellSize*round(valuesThresh.iloc[j,1]/cellSize)) - minLat)/cellSize),
                 int(((cellSize*round(valuesThresh.iloc[j,0]/cellSize)) - minLon)/cellSize)].append(eruptions.iloc[i,1])
        

def carbonAccumulation(x):
    res = 501.4*(x**(-0.55))
    return res

def get_carbon_grid(grid, startYear, stopYear, surfaceC):
    carbonGrid = np.empty((len(grid[:]),len(grid[0])))
    logC = [0] * (stopYear - startYear)
    for i in range(len(grid[:])):
        for j in range(len(grid[0])):
            sumC = 0
            if len(grid[i,j]) > 2:
                if not grid[i,j][-1] == startYear:
                    grid[i,j].append(startYear)
                k = len(grid[i,j]) -1
                while k > 2:
                    timeDif = grid[i,j][k-1] - grid[i,j][k]
                    if timeDif > 0:
                        amountC, error = quad(carbonAccumulation, 0, timeDif)
                        logC[grid[i,j][k-1]-startYear] += amountC/2
                        sumC += amountC/2
                    k -= 1
                if surfaceC == "yes":
                    amountC, error = quad(carbonAccumulation, 0, stopYear - grid[i,j][2])
                    sumC += amountC
            carbonGrid[i,j] = sumC
            
    return carbonGrid, logC

def save_results(outputFolder, carbonGrid, logC, eruptions):
    if not path.exists(outputFolder + "Runs/") :
        mkdir(outputFolder + "Runs/")
    np.savetxt(outputFolder + "Runs/run" + str(len(listdir(outputFolder + "Runs/"))+1) + ".csv",carbonGrid, delimiter=",")
    if not path.exists(outputFolder + "LogC/"):
        mkdir(outputFolder + "LogC/")
    np.savetxt(outputFolder + "LogC/logCrun" + str(len(listdir(outputFolder + "LogC/"))+1) + ".csv",logC, delimiter=",")
    if not path.exists(outputFolder + "LogErupt/"):
        mkdir(outputFolder + "LogErupt/")
    eruptions[['Volcano','Year','VEI']].to_csv(outputFolder + "LogErupt/eruptRun" + str(len(listdir(outputFolder + "LogErupt/"))+1) + ".csv", header=False, index=False)  
    count = len(listdir(outputFolder + "LogC/"))+1
    
    return count
        


