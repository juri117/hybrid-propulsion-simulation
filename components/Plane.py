'''
Created on 01.12.2015

@author: Juri Bieler
'''

from Wing import Wing as Wing
import Propulsion as Propulsion
from Ini import Ini as Ini
import const as c
from Mission import Mission
import planePlot

import std_atm as SA
import numpy as np
import math
import scipy.interpolate as inter
import scipy.optimize


class Plane:

    ### CONSTRUCKTOR #####################################
    '''
    initialisiert "Wing", "Propulsion"
    Definiert Antriebskonfiguration entsprechend der Parameter in "config"
    falls "config" nicht definiert ist werden die Parameter aus der config.ini der Mission ausgelsen
    '''
    def __init__(self, mission, propulsionType=c.COMB, config={}):
        self.verbose = True

        self.w = Wing()

        self.propulsionType = propulsionType

        if(len(config) == 0):
            ini = Ini(mission.NAME)
            config["pmax"] = float(ini.get(self.propulsionType, "pmax"))
            config["fuel"] =  float(ini.get(self.propulsionType, "fuel"))
            if(self.propulsionType != c.COMB):
                config["cbat"] = float(ini.get(self.propulsionType, "cbat"))
                config["pice"] = float(ini.get(self.propulsionType, "pice"))
            if(ini.keyExists("general", "mpl")):
                config["mpl"] = float(ini.get("general", "mpl"))
            else:
                config["mpl"] = c.M_PL
            if(ini.keyExists(self.propulsionType, "pbat")):
                config["pbat"] = float(ini.get(self.propulsionType, "pbat"))
            
        #init für Verbrennungsmotor
        if(not "cbat" in config):
            config["cbat"] = 0.
        if(not "pbat" in config):
            config["pbat"] = 0.
        
        if(len(config) > 0):
            if(self.propulsionType == c.COMB):
                self.propulsion = Propulsion.Combustion(config)
            if(self.propulsionType == c.SERI):
                self.propulsion = Propulsion.SeriellHybrid(config)
            if(self.propulsionType == c.PARA):
                self.propulsion = Propulsion.ParallelHybrid(config)

        config["pice"] = self.propulsion.ice.Pmax
        self.config = config
        self.P_max = config["pmax"]
     
    ### FUNCTIONS #####################################
    
    ### Weight ###
    
    '''
    Rückgabe: aktuelles Fluggewicht
    '''
    def calcActWeight(self):
        weight = self.calcDryWeight() + (self.config["fuel"] * c.DENS_FUEL)
        return weight
    
    '''
    Rückgabe: Leergewicht (mit Antrieb)
    '''
    def calcDryWeight(self):
        weight = c.M_STRUCT + self.config["mpl"] + self.propulsion.calcPropulsionMasses()
        return weight
    
    
    ### Power ###
    
    '''
    Rückgabe: Wellenleistung und Drehzahl die in der Flughöhe "h" mit der Geschwindigkeit "tas" für den Schub "F" benötigt wird
    '''
    def calcPshaft(self, h, tas, F):
        return self.propulsion.calcPshaft(h, tas, F)

    '''
    Berechnet Verbrauch und Triebwerksparameter die in der FLughöhe "h" mit der Geschwindigkeit "tas" für den Schub "F" benötigt wird
    Rückgabe: Pshaft, Nprop, m_f, We, Pice, Nice, Pshaft, Wcharge
    '''
    def calcConsumption(self, h, tas, F, dt):
        return self.propulsion.calcConsumption(h, tas, F, dt)


    ### Thrust ###

    '''
    Führt Startlaufsimulation aus
    Rückgabe: P_to, t_abh, v_abh, F_0
    '''
    def calcTakeOffParameter(self, m, h, s):
        return self.w.calcTakeOffParameter(m, h, s, self.propulsionType, self.propulsion)
    
    '''
    Berechnet benötigten Schub mit Schubformel wie in Bachelorarbeit beschrieben
    Rückgabe: Schub
    '''
    def calcNeededThrust(self, m, h, u_g, w_g, curveR):
        if(u_g < c.V_MIN):
            return 0.
        
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        #Bahnneigungswinkel [deg]
        gamma = math.atan(float(w_g) / u_g)
        tas = math.sqrt(u_g**2. + w_g**2.)
        
        if(curveR == 0 or math.isnan(curveR)):
            phi = 0.
        else:
            phi = math.atan(tas**2. * math.cos(gamma) / (c.g * curveR))

        n = math.cos(gamma) / math.cos(phi)
        
        F1 = (math.sin(gamma) * m * c.g)
        F2 = ((self.w.C_W0 * c.S * rho * tas**2) / 2)
        F3 = (self.w.k * 2 * (n * m * c.g)**2 / (rho * c.S * tas**2))

        F = F1 + F2 + F3

        if(F < 0):
            return 0.
        return F
    
    
    ### Speed ###

    '''
    Berechnet Geschwindigkeit für steilstes Steigen
    Rückgabe: Fluggeschwindigkeit
    '''
    def calcV_Fmin(self, m, h, w_g, R):
        u_g = np.array(range(20, 75, 1))
        P = self.P_max * 0.85
        Diff = []
        Nprop = self.propulsion.calcNProp(P)
        
        for v in u_g:
            eta, Nprop = self.propulsion.prop.calcEfficiency(h, v, Nprop, P)
            F = eta * P / v
            Diff.append(F - self.calcNeededThrust(m, h, v, w_g, R))
        
        Diff = (np.array(Diff).T)[0]

        f_Diff = inter.interp1d(u_g, -1 * Diff, kind='cubic')
        u_gamma_max = scipy.optimize.fmin(f_Diff, 35, disp=False)
        return u_gamma_max

    '''
    Berechnet Fahrt minimalen Widerstands
    Rückgabe: Fluggeschwindigkeit
    '''
    def calcV_W_min(self, m, h, u_g, R):
        return self.w.calcV_W_min(m, h)

    '''
    Berechnet Fahrt mit geringstem Sinken
    Rückgabe: Fluggeschwindigkeit
    '''
    def calcV_dH_min(self, m, h, u_g, R):
        return self.w.calcV_dH_min(m, h)
    
    '''
    Plottet Grafik mit Ürinzip der Berechnung der Geschwindigkeit für steilstes Steigen
    '''
    def plotThrustOverflow(self):
        planePlot.plotThrustOverflow(self)

if __name__ == "__main__":
    p = Plane(Mission('mission04'), "combustion")
    print (p.calcNeededThrust(600, 1000, 25, 4, 0))
    print (p.calcNeededThrust(600, 1000, 25, 4, 100))
    print (p.calcNeededThrust(600, 1000, 25, 4, 0))
    p.plotThrustOverflow()