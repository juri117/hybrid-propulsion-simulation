'''
Created on 01.12.2015
 
@author: Juri Bieler
'''
 
from Pilot import Pilot
from Logger import Logger
from Optimizer import Optimizer
from components.Mission import Mission
from Mesh import Mesh
from Ini import Ini
import utils.generalUtils as generalUtils
import sys
import const as c
from meshPlot import *
import time
from queue import Queue
from threading import Thread
import os

print('start...')
log = Logger()
t0 = time.time()

'''normale Missionen ausfürhren'''
if(True):
    missions = generalUtils.getAllFolders(pref="mission")
    resultsPostfix = ""
    
    #nur die nicht ausgefuehrten
    mis = []
    for m in missions:
        if(not os.path.isdir(generalUtils.getMissionPath(m) + "results")):
            mis.append(m)

    ''''wenn nicht auskommentiert werden NUR die Missionen, die noch kein Ergebnisordner haben ausgeführt'''
    #missions = mis
    
    ''''wen nicht auskommentiert werden NUR die in der Liste eingetragenen Missionen ausgeführt'''
    #missions = ["mission02_1", "mission03", "mission04"]
    
    print(missions)
    
    for m in missions:
        print("---------- " + m + " ----------")
        verbose = False
        mis = Mission(m, resultsPostfix)
        opt = Optimizer(mis, verbose=verbose)
        opt.saveFlights()
        ini = Ini(mis.NAME)
        ini.copyConfig(mis.getResultsPath() + "config.ini")
        del ini


'''Geschwindigkeit-Reichweite-Raster ausführen'''
if(False):
    mission = Mission('mission05')
    mesh = Mesh(mission)
    mesh.runMesh()
    plotMesh(mission, saveCSV=True)


t1 = time.time()
print('[main] Runtime [s]: %f' %(t1-t0))

log.terminate()
print('end!')