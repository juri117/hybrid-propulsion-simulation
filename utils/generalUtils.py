'''
Created on 15.12.2015

@author: Juri Bieler
'''

import numpy as np
import csv
import math
import os
import itertools
import operator
import matplotlib.pyplot as plt
from matplotlib import rc
import shutil
import const as c

def getMissionPath(mission):
    return 'data/missions/' + mission + "/"
    #return fixPath('data/missions/' + mission + "/")

#DEPRECATED -> moved to Mission.py
# def getMissionResPath(mission):
#     resPath = getMissionPath(mission) + "results/"
#     if(not os.path.isdir(resPath)):
#         os.makedirs(resPath)
#     return resPath


def calcNeededConfig(flight):
    Pneeded = max(flight.Pshaft)
    Fuel_needed = flight.calcUsedFuel()
    #Cbat_need = (-1. * (min(flight.SoC)) + 1.) * config["cbat"]
    #if(config["cbat"] <= 0.):
    #    Cbat_need = (-1. * (min(flight.SoC)) + 1.) * 0.1
    Cbat_needed = -1 * np.sum(flight.We)
    Cbat_needed = max(Cbat_needed, (-1. * (min(flight.SoC)) + 1.) * flight.config["cbat"])
    Pbat_needed = 0.
    #kleine Toleranz für rundungsfehler
    if(Cbat_needed < 1.):
        Cbat_needed = 0.
    else:
        Pbat_needed = max(abs(min(flight.chargePower)), abs(max(flight.chargePower)))
    return Pneeded, Fuel_needed, Cbat_needed, Pbat_needed

def copyIniToResults(mission):
    src = mission.getMissionPath() + "config.ini"
    dst = mission.getResultsPath() + "config.ini"
    shutil.copy(src, dst)

def getAllFolders(pref="mission", root="data/missions"):
    root = fixPath(root)
    dirs = next(os.walk(root))[1]
    missions = []
    for d in dirs:
        if(d.startswith(pref)):
            missions.append(d)
    return missions

def writeFlights2CSV(propulsionType, mission, flights):
    if(len(flights) > 0):
        filePath = mission.getResultsPath()
        with open(filePath + "runningConfigs_" + propulsionType + ".csv", 'w', newline='') as csvfile:
            heads = ["config", "pmax", "pice", "cbat", "fuel", "hybDegree", "hybDegreeMass", "mto"]
            
            #add changed (prefix) mass keys
            masses = addPrefix2Dict(flights[0].mass, "mass_")
            for k, _ in masses.items():
                heads.append(k)
            
            writer = csv.DictWriter(csvfile, fieldnames=heads, delimiter=c.CSV_DELIMITER)
            writer.writeheader()
            for fli in flights:
                dic = {}
                dic["pmax"] = fli.config["pmax"]
                dic["pice"] = fli.config["pice"]
                dic["cbat"] = fli.config["cbat"]
                dic["fuel"] = fli.config["fuel"]
                dic["hybDegree"] = 1 - (fli.config["pice"] / fli.config["pmax"])
                dic["hybDegreeMass"] = 1 - (fli.mass["ICE"] / (fli.mass["ICE"] + fli.mass["EMot"]))
                dic["mto"] = fli.m[0]
                masses = addPrefix2Dict(fli.mass, "mass_")
                dic.update(masses)
                dic["config"] = str(fli.config)
                #dic["mto"] = fli.mass
                
                writer.writerow(dic)
            

#data is a list of dicts
def write2CSV(filePath, data, heads=None):
    if(heads == None):
        heads = list(data[0].keys())
    with open(filePath, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=heads, delimiter=c.CSV_DELIMITER)
        writer.writeheader()
        for d in data:
            writer.writerow(d)
            
def addPrefix2Dict(dic, prefix):
    return dict([(prefix + k, v) for k, v in dic.items()])

#fix if path is needed in subcomponent (only Debug should need this)
def fixPath(path):
    if(not os.path.exists(path)):
        path = '../' + path
    return path

def readCsv(filename):
    filename = fixPath(filename)
    M = np.genfromtxt(filename, delimiter=';', skip_header=1)
    #M = np.genfromtxt(filename, delimiter=';', dtype=None, skip_header=1)
    return M

def readStrForNaN(filename, col, vec):
    filename = fixPath(filename)
    out = []
    with open(filename, 'rt') as f:
        reader = csv.reader(f, delimiter=';')
        next(reader)
        i=0
        for row in reader:
            if(math.isnan(vec[i])):
                out.append(row[col])
            else:
                out.append('')
            i+=1
    return out

def read1dAmeFile(filename):
    filename = fixPath(filename)
    row = 0
    rowCount = sum(1 for line in open(filename))
    fileRow = np.genfromtxt(filename, delimiter=' ', skip_header=row, skip_footer=rowCount-row-1)

    while(len(np.atleast_1d(fileRow)) <= 1):
        row += 1
        fileRow = np.genfromtxt(filename, delimiter=' ', skip_header=row, skip_footer=rowCount-row-1)

    M = np.genfromtxt(filename, delimiter=' ', skip_header=row, skip_footer=0)
   
    return M
            
def read2dAmeFile(filename):
    filename = fixPath(filename)
    row = 0
    rowCount = sum(1 for line in open(filename))
    fileRow = np.genfromtxt(filename, delimiter=' ', skip_header=row, skip_footer=rowCount-row-1)

    while(len(np.atleast_1d(fileRow)) <= 1):
        row += 1
        fileRow = np.genfromtxt(filename, delimiter=' ', skip_header=row, skip_footer=rowCount-row-1)
    
        
    x = np.genfromtxt(filename, delimiter=' ', skip_header=row, skip_footer=rowCount-row-1)
    row += 1
    y = np.genfromtxt(filename, delimiter=' ', skip_header=row, skip_footer=rowCount-row-1)
    row += 1
    M = np.genfromtxt(filename, delimiter=' ', skip_header=row, skip_footer=0)
   
    return x, y, M
    
def exportData(time, values, filename):
    output = '# Table format: XY'
    for i in range(0,len(time),1):
        out1 = formatToScience(time[i])
        out2 = formatToScience(values[i])
        output = output + '\n' + ' ' + out1 + ' ' + out2;
            
    f = open(filename, 'w')
    f.write(output) # python will convert \n to os.linesep
    f.close()
    
def formatToScience(num):
    strOut = '{:.14e}'.format(float(num))
    strOut = strOut.replace('e+', 'e+0')
    strOut = strOut.replace('e-', 'e-0')
    return strOut

def minRealToSecInt(min):
    return int(int(min) * 60 + round((min - int(min)) * 100))

def secIntToMinReal(sec):
    m, s = divmod(sec, 60)
    #h, m = divmod(m, 60)
    return m + (float(s) / 100.)

def secToTimeStr(sec):
    m, s = divmod(sec, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

def custom_div_cmap(numcolors=11, name='custom_div_cmap', mincol='blue', midcol='white', maxcol='red'):
    """ Create a custom diverging colormap with three colors
    Default is blue to white to red with 11 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    """
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(name=name, colors =[mincol, midcol, maxcol], N=numcolors)
    return cmap

def monotone_increasing(lst):
    pairs = zip(lst, lst[1:])
    return all(itertools.starmap(operator.le, pairs))

def monotone_decreasing(lst):
    pairs = zip(lst, lst[1:])
    return all(itertools.starmap(operator.ge, pairs))

def monotone(lst):
    return monotone_increasing(lst) or monotone_decreasing(lst)

def pltCosmetics(multiPlot=False):
    font = {'family' : 'sans-serif',
            #'weight' : 'bold',
            'size'   : 18}

    rc('font', **font)

    if(multiPlot):
        plt.subplots_adjust(left=0.09, bottom=0.072, right=0.7, top=0.95, wspace=None, hspace=0.25)
    else:
        plt.subplots_adjust(left=0.14, bottom=0.13, right=0.95, top=0.95, wspace=None, hspace=None)
    #rc('font', family='serif', serif='cm10')

    #rc('text', usetex=True)