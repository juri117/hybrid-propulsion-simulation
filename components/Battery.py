'''
Created on 07.12.2015

@author: Juri Bieler
'''

import const as c
from units import *
import generalUtils as generalUtils

import numpy as np
import scipy.interpolate as inter
import matplotlib.pyplot as plt

class Battery(object):
    ### CONSTRUCKTOR #####################################
    '''
    liest Batteriekennfeld (Leistungsdichte über Energiedichte) ein und interpoliert
    bestimmt idealen Punkt im Kennfeld für verhältnis von Kappazität "C" und max. Leistung "Pmax"
    '''
    def __init__(self, C, Pmax):
        self.verbose = True
        self.C = C
        self.Pmax = Pmax
        self.SoC = 1.
        
        self.etaCharge = 0.9
        self.etaDischarge = 0.9
        
        Pd_Ed = generalUtils.read1dAmeFile(c.BAT_FILE)
        Pd = Pd_Ed[:,0]
        Ed = Pd_Ed[:,1]
        Crate = np.divide(Pd, Ed)
        self.f_Pd_Ed = inter.interp1d(Ed, Pd, kind='linear')
        self.f_Crate = inter.interp1d(Crate, Ed, kind='quadratic')
        
        self.Ed_ideal = np.NaN
        self.Pd_ideal = np.NaN
        
        if(self.C > 0):
            Crate_bat = float(Pmax) / float(self.C)
            Crate_bat = max(Crate_bat, min(Crate))
            Crate_bat = min(Crate_bat, max(Crate))
            
            self.Ed_ideal = self.f_Crate(Crate_bat)
            self.Pd_ideal = self.f_Pd_Ed(self.Ed_ideal)
        
        #print("ED: " + str(self.Ed_ideal))
        #print("PD: " + str(self.Pd_ideal))


    '''
    Rückgabe: Batteriemasse in kg
    '''
    def calcMass(self):
        #Masse um Energie bereit zu stellen
        mass_E = self.C / self.Ed_ideal
        #Masse um Leistung bereit zu stellen
        mass_P = float(self.Pmax) / self.Pd_ideal
        return max(mass_P, mass_E)
    
    '''
    setzt SoC (state of charge) auf eins zurück, sollte es höher sein
    '''
    def checkSoC(self):
        if(self.SoC > 1.):
            self.SoC = 1.        
            
    '''
    wird ausgeführt egal ob die Arbeit "We" positiv oder negativ ist
    berechnet neuen Ladezustand der Batterie (SoC)
    wird mehr als der max. zulässige (in "const.py") Strom geladen, so wird "We" entsprechend reduziert
    Rückgabe: Ladearbeit "We", Ladeleistung "Pcharge"
    '''
    def charge(self, We, dt):
        if(self.C > 0.):
            #discharge
            if(We < 0.):
                We = We / self.etaDischarge
                self.SoC += We / self.C
            #charge
            else:
                if(We * 60**2 / dt > self.getMaxChargeRate()):
                    We = self.C * dt * c.CHARGE_RATE / 60**2
                We = We / self.etaCharge
                self.SoC += We / self.C
        if(self.full()):
            self.SoC = 1.
        if(self.verbose):
            ChargePower = (We * 60**2 / dt)
            print("[Battery] ChargePower: " + str(ChargePower))
        Pcharge = We * 60**2 / dt 
        return We, Pcharge
    
    '''
    Rückgabe: max. Laderate in W
    '''
    def getMaxChargeRate(self):
        return self.C * c.CHARGE_RATE
    
    '''
    Rückgabe: max. Ladearbeit in Wh
    '''
    def getMaxWeCharge(self, dt):
        return self.getMaxChargeRate() / (60**2 / dt)
    
    '''
    Rückgabe: Batterie voll? [True|False]
    '''
    def full(self):
        return (self.SoC >= 1.)

    '''
    Rückgabe: Batterie leer? [True|False]
    '''
    def empty(self):
        return (self.SoC < 0.)
    
    '''
    berechnet die Lade- und Entladeraten die während dem Flug "flight" auftraten
    '''
    def calcCharegeRates(self, flight):
        C = []
        P = []
        
        for i in range(1, len(flight.t), 1):
            dSoC = flight.SoC[i] - flight.SoC[i-1]
            dt = flight.t[i] - flight.t[i-1]
            #We = dSoC * self.C
            chargeRate = 0.
            Prate = 0.
            if(dt > 0):
                chargeRate = (dSoC * 60**2 / dt)
                Prate = (dSoC * 60**2 / dt) * self.C
            C.append(chargeRate)
            P.append(Prate)
        return C, P
    
    def plotChargeRate(self, flight):
        C = self.calcCharegeRates(flight)
        plt.plot(flight.t[:-1], C)
        plt.show()
    
    def plotPD_ED(self):
        x = range(10,270,10)
        #x = range(340,430,10)
        plt.plot(x, self.f_Pd_Ed(x))
        #plt.plot(Ed, Pd, "o")
        plt.show()
        
        x = range(1,110,1)
        #x = range(1,3,1)
        plt.plot(x, self.f_Crate(x))
        plt.show()

        print("done")

if __name__ == "__main__":
    print("[Battery] nothing")
    bat = Battery(8000., 20000.)
    bat.plotPD_ED()
    #print(str(bat.getMaxWeCharge(5)))
