'''
Created on 31.12.2015

@author: Juri Bieler
'''

from Ini import Ini
from Optimizer import Optimizer
from Mission import Mission
from Plane import Plane
import const as c
import generalUtils
import time
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from meshPlot import *


class Mesh(object):

    ### CONSTRUCKTOR #####################################
    '''
    initialisiert die Mission namens "mission" und liest die Anweisungen wie das Geschwindigkeit-Reichweite-Raster zu erstellen ist aus der config.ini
    die erste Mission des Rasters (mit min. Geschw. und Reichweite) wird initialisiert
    '''
    def __init__(self, mission):
        self.mis = mission
        ini = Ini(mission.NAME)
        
        self.stratchRow = int(ini.get("mesh", "stratchRow"))
        
        self.Vmin = float(ini.get("mesh", "vmin"))
        self.Vmax = float(ini.get("mesh", "vmax"))
        self.Vadd = float(ini.get("mesh", "vadd"))
            
        self.Vact = self.Vmin
        
        self.Tadd = None
        if(ini.keyExists("mesh", "tmax") and ini.keyExists("mesh", "tadd")):
            self.Tmax = float(ini.get("mesh", "tmax"))
            self.Tadd = float(ini.get("mesh", "tadd"))
            self.Toffset = 0
            if(ini.keyExists("mesh", "toffset")):
                self.Toffset = float(ini.get("mesh", "toffset"))
                self.addTime(self.Toffset)
            self.mis.RESULTS_FOLDER = "results_t" + str(self.Tact).zfill(5) + "_v" + str(self.Vact).zfill(3)
        self.Torigin = self.mis.sec.copy()
        self.Tact = 0
        
            
        self.Radd = None
        self.Ract = 0
        if(ini.keyExists("mesh", "rmax") and ini.keyExists("mesh", "radd")):
            self.Rmax = float(ini.get("mesh", "rmax")) * 1000.
            self.Radd = float(ini.get("mesh", "radd")) * 1000.
            if(ini.keyExists("mesh", "roffset")):
                self.Roffset = float(ini.get("mesh", "roffset")) * 1000.
                self.Ract = self.Roffset
            self.mis.RESULTS_FOLDER = "results_r" + str(str(self.Ract/1000.)).zfill(6) + "_v" + str(self.Vact).zfill(3)

        self.Vact -= self.Vadd
        self.nextMission()


    '''
    Führt "Optimizer.saveFLights" für jeden Punkt im Raster aus, für jeden Punkt im Raster wird in 
    dem Missionsverzeichnis ein Ordner angelegt in dem die Ergebnisse in Textdateien gespeichert werden
    '''
    def runMesh(self):
        t0 = time.time()
        while(self.Ract <= self.Rmax):
            print("----------------------------------")
            print("actT: " + str(self.Tact) + " - " + str(self.mis.sec))
            #print("Torig: " + str(self.Torigin))
            print("actV: " + str(self.Vact) + " - " + str(self.mis.tas))
            print("actR: " + str(self.Ract / 1000.) + " - " + str(self.mis.calcDistance() / 1000.))
            
            self.mis.saveMission()
            
            verbose = False
            opt = Optimizer(self.mis, verbose=verbose, savePlots=False)
            opt.saveFlights(verbose=False)
            del opt
            
            ini = Ini(self.mis.NAME)
            ini.copyConfig(self.mis.getResultsPath() + "config.ini")
            del ini
            
            self.nextMission()
            #self.mis.plotMission()
        t1 = time.time()
        print('[Mesh] Runtime [s]: %f' %(t1-t0))


    '''
    initialisiert die nächste Mission im Raster
    erst wird die Geschwindigkeit hochgezählt, dann die Reichweite erhäht und die Geschwindigkeit wieder aufs minimum gesetzt.
    '''
    def nextMission(self):
        self.Vact += self.Vadd
        if(self.Vact > self.Vmax):
            self.Vact = self.Vmin
            if(self.Tadd != None):
                self.addTime(self.Tadd)
            if(self.Radd != None):
                self.Ract += self.Radd
        self.changeSpeed(self.Vact)
        
        if(self.Radd != None):
            T = int(self.Ract / self.Vact)
            self.changeTime(T)
            self.mis.RESULTS_FOLDER = "results_r" + str(str(self.Ract/1000.)).zfill(6) + "_v" + str(self.Vact).zfill(3)
        else:
            self.mis.RESULTS_FOLDER = "results_t" + str(self.Tact).zfill(5) + "_v" + str(self.Vact).zfill(3)
        self.mis.refraeshMission()
    
    
    '''
    addiert die Zeit "dT" zu der Reiseflugzeit der Mission
    '''
    def addTime(self, dT):
        self.Tact += dT
        for i in range(0, len(self.mis.sec), 1):
            if(i > self.stratchRow):
                self.mis.sec[i] += int(dT * 60.)
    
    
    '''
    setzt die Reiseflugzeit der Mission auf die Zeit "T"
    '''
    def changeTime(self, t):
        for i in range(0, len(self.mis.sec), 1):
            if(i > self.stratchRow):
                self.mis.sec[i] = int(self.Torigin[i] + t)
    
    
    '''
    setzt die Reiseflug-Geschwindigkeit der Mission auf die Geschwindigkeit "V"
    '''
    def changeSpeed(self, V):
        self.mis.cas[self.stratchRow] = V
        self.mis.cas[self.stratchRow + 1] = V


### MAIN #####################################

if __name__ == "__main__":
    mission = Mission('mission06_01')
    mesh = Mesh(mission)
    #mesh.runMesh()
    plotMesh(mission, saveCSV=False)
