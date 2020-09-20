'''
Created on 30.12.2015

@author: Juri Bieler
'''

import numpy as np
import const as c

class Flight:

    ### CONSTRUCKTOR #####################################
    '''
    initialisiert Listen für Parameter die über der Flugzeit gespeichert werden
    '''
    def __init__(self, stepWidth, startMass, config):
        self.config = config
        self.stepWidth = stepWidth
        self.startMass = startMass
        self.distance = 0.
        
        self.t = []   #sec
        self.tas = [] #m/s
        self.H = []   #m
        self.r = []   #m
        self.F = []   #N
        self.m = []   #kg
        self.Pshaft = [] #Watt
        self.Nprop = []
        self.Pice = [] #Watt
        self.Nice = []
        self.Pemot = [] #Watt
        self.SoC = []   #%
        self.SoF = []   #%
        self.We = []
        self.chargePower = []   #W
        self.mass = {}
        
        self.abort = False
        self.fails = []
        self.runTime = 0
    
    ### FUNCTIONS #####################################
    
    '''
    speichert, dass ein Feheler mit der Beschreibung "failStr" aufgetreten ist
    '''
    def fail(self, failStr):
        self.abort = True
        self.fails.append(failStr)
        
    '''
    Rückgabe: Falls vorhanden Feherlmeldungen, Kommaseperiert
    '''
    def getErrorMessage(self):
        return ", ".join(self.fails)
    
    '''
    Fügt die Daten eines Simulationsschrittes hinzu
    '''
    def addStep(self, t, m, H, tas, r, F, Pshaft, Nprop, Pice, Nice, Pemot, SoC, SoF, We, chargePower):
        self.t.append(t)
        self.tas.append(tas)
        self.H.append(H)
        self.r.append(r)
        self.F.append(F)
        self.m.append(m)
        self.Pshaft.append(float(Pshaft))
        self.Nprop.append(float(Nprop))
        self.Pice.append(float(Pice))
        self.Nice.append(float(Nice))
        self.Pemot.append(float(Pemot))
        self.SoC.append(SoC)
        self.SoF.append(SoF)
        self.We.append(We)
        self.chargePower.append(float(chargePower))


#    def calcP(self, t):
#        return self.F[t] * self.tas[t]
#    
#    def calcP_kw(self, t):
#        return self.calcP(t) / 1000
    
    '''
    Rückgabe: Flugweg, flugzeugfest
    '''
    def calcDistance(self):
        m = np.trapz(self.tas, dx=self.stepWidth)
        return m
    
    '''
    Rückgabe: der tatsächlich benötigte Brennsotff in l
    '''
    def calcUsedFuel(self):
        return float((self.m[0] - self.m[-1]) / c.DENS_FUEL)
    
    '''
    Rückgabe: die tatsächlich benötigte elek. Energie in kW
    '''
    def calcUsedElec(self):
        return float((self.SoC[0] - self.SoC[-1]) * self.config["cbat"])
    
    '''
    Rückgabe: die Masse im letzten gespeicherten Simulationsschritt
    '''
    def getLastMass(self):
        if(len(self.m) > 0):
            return self.m[-1]
        return self.startMass
    
    '''
    prüft ob Masse, Batterie, Maximalleistung und Brennstoffverbrauch für die Mission ausgereicht haben
    wenn nicht werden Fehler gespeichert
    '''
    def checkFlight(self):
        if(max(self.m) > c.MTOM):
            self.fail("MTOM ueberschritten")
        if(min(self.SoC) < -0.1):
            self.fail("Batterie tiefentladen")
        if(max(self.SoC) > 1.1):
            self.fail("Batterie ueberladen")
        if(min(self.SoF) < -0.1):
            self.fail("Mehr Sprit gebraucht als mitgenommen")
        if(max(self.Pshaft) > 1.1 * self.config["pmax"]):
            self.fail("Mehr Sprit mitgenommen als in die Tanks passt")
        