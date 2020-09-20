'''
Created on 15.12.2015

@author: Juri Bieler
'''

from Ini import Ini as Ini

import numpy as np
import generalUtils as generalUtils
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import math
import scipy.interpolate as inter
import std_atm as SA
import os

class Mission:

    ### CONSTRUCKTOR #####################################
    '''
    liest Mission aus csv aus und berechnet Steigraten
    liest Startrollstrecke aus config.ini aus
    führt "refraeshMission" aus
    '''
    def __init__(self, missionFile, resultsPostfix = ""):
        self.NAME = missionFile
        self.PATH = generalUtils.getMissionPath(self.NAME)
        self.CSV_PATH = self.PATH + 'mission.csv'
        
        self.RESULTS_FOLDER = "results" + resultsPostfix
        
        #load mission file
        csv = generalUtils.readCsv(self.CSV_PATH)
        time = csv[:,0]
        self.sec = list(map(generalUtils.minRealToSecInt, time))
        self.height = csv[:,1]
        self.cas = csv[:,2]
        self.vSpecial = generalUtils.readStrForNaN(self.CSV_PATH, 2, self.cas)
        self.curveR = csv[:,3]
        
        ini = Ini(self.NAME)
        self.takeOffDistance = int(ini.get("general", "takeoffdistance"))   #[m]

        self.refraeshMission()        
        #print("Mission wurde initialisiert")
        
    def getResultsPath(self):
        resPath = self.PATH + self.RESULTS_FOLDER + "/"
        if(not os.path.isdir(resPath)):
            os.makedirs(resPath)
        return resPath
    
    def getMissionName(self):
        return self.NAME
    
    def getMissionPath(self):
        return self.PATH

    '''
    berechnet "tas", Steigrate "w_g" und erstellt funktionen mit der linearen interpolation aller Missionsparameter über der Zeit
    '''
    def refraeshMission(self):
        self.tas = np.array(self.CAS_to_TAS(self.cas, self.height))
        self.w_g = np.zeros(len(self.sec))
        for i in range(0, len(self.sec), 2):
            self.w_g[i] = self.calcW_g(i)
            self.w_g[i+1] = self.w_g[i]
        
        self.u_g = list(map(self.calcU_g, range(0, len(self.sec),1)))

        self.f_h = inter.interp1d(self.sec, self.height, kind='linear')
        self.f_R = inter.interp1d(self.sec, self.curveR, kind='linear')
        self.f_tas = inter.interp1d(self.sec, self.tas, kind='linear')
        self.f_u_g = inter.interp1d(self.sec, self.u_g, kind='linear')
        self.f_w_g = inter.interp1d(self.sec, self.w_g, kind='linear')

    ### FUNCTIONS #####################################

    def CAS_to_TAS(self, cas, h):
        rho0 = SA.alt2density(0, 'm', 'kg/m**3')
        last_val = 0.
        i = 0
        tas = []
        for v in cas:
            rho = SA.alt2density(h[i], 'm', 'kg/m**3')
            tas_val = v * math.sqrt(rho0/rho)
            if(math.isnan(tas_val)):
                tas_val = last_val
            else:
                last_val = tas_val
            tas.append(tas_val)
            i += 1
        return tas
    
    def calcDistance(self):
        s = 0
        for i in range(0, len(self.sec)-1, 1):
            s += (self.sec[i+1] - self.sec[i]) * ((self.u_g[i+1] + self.u_g[i]) / 2)
        return s
    
    
    '''
    alle folgenden Funktionen geben Missionsparameter zur Zeit t in sec. wieder (es wird linear interpoliert)
    '''
    def getHeight(self, t):
        return np.interp(t, self.sec, self.height)
    
    def getRadius(self, t):
        return np.interp(t, self.sec, self.curveR)
    
    def getTas(self, t):
        return np.interp(t, self.sec, self.tas)
    
    def getU_g(self, t):
        return np.interp(t, self.sec, self.u_g)
    
    def getW_g(self, t):
        return np.interp(t, self.sec, self.w_g)

    '''
    Rückgabe: horizontalfluggeschwindigkeit im Missionsabschnitt "step"
    '''
    def calcU_g(self, step):
        if(self.tas[step]**2 - self.calcW_g(step)**2 < 0.):
            return float("nan")
        return math.sqrt(self.tas[step]**2 - self.calcW_g(step)**2)
    
    '''
    Rückgabe: vertikalfluggeschwindigkeit im Missionsabschnitt "step"
    '''
    def calcW_g(self, step):
        if step + 1 >= len(self.sec):
            return 0
        deltH = self.height[step+1] - self.height[step]
        deltT = self.sec[step+1] - self.sec[step]
        if(deltT == 0):
            return 0
        return deltH / deltT
    
    
    ### PLOTS #####################################
    
    def plotMission(self):
        seconds = range(0, max(self.sec), 1)
        h = self.getHeight(seconds)
        v = self.getTas(seconds)

        plt.figure("Missionsprofil")
        
        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.75)
    
        par1 = host.twinx()
        
        host.set_xlabel("Zeit [s]")
        host.set_ylabel("H [m]")
        par1.set_ylabel("TAS [m/s]")

        p1, = host.plot(seconds, h, label="Flughoehe")
        p2, = par1.plot(seconds, v, label="Geschwindigkeit")
        #host.legend()
        
        host.axis["left"].label.set_color(p1.get_color())
        par1.axis["right"].label.set_color(p2.get_color())

        for i in range(0, len(self.sec), 2):
            if(self.w_g[i] != 0.):
                txt = "{:.1f}".format(self.w_g[i])
                x = self.sec[i]
                y = self.height[i]
                if(len(self.sec) > i + 1):
                    x = (self.sec[i] + self.sec[i+1]) / 2
                    y = (self.height[i] + self.height[i+1]) / 2
                host.annotate(txt, xy=(x, y), xycoords="data", va="center", ha="center", bbox=dict(boxstyle="round", fc="w"))
        
        host.annotate("w_g", xy=(0.95, 0.95), xycoords="axes fraction", va="center", ha="center", bbox=dict(boxstyle="round", fc="w"))
        
        host.axes.set_ylim([0., 1.1  * max(h)])
        #par1.axes.set_ylim([0., 1.2  * max(v)])
        
        plt.draw()
        plt.show()
        print('plot...')
    
    
    ### SAVE IT #####################################
     
    '''
    speichert Mission in csv
    '''
    def saveMission(self):
        csvData = []
        for i in range(0, len(self.sec), 1):
            row = {"Zeit [m.ss]":generalUtils.secIntToMinReal(self.sec[i]),
                   "Höhe":self.height[i],
                   "CAS":self.cas[i],
                   "R":self.curveR[i]}
            if(math.isnan(row["R"])):
                row["R"] = ""
            csvData.append(row)
        csvHead = ["Zeit [m.ss]", "Höhe", "CAS", "R"]
        generalUtils.write2CSV(self.getResultsPath() + "mission.csv", csvData, heads=csvHead)

if __name__ == "__main__":
    m = Mission('mission05')
    #m.plotMission()
    print(str(m.calcDistance()))
    
        