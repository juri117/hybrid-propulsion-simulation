'''
Created on 01.12.2015

@author: Juri Bieler
'''

import const as c
from units import * 
from Propeller import Propeller as Propeller
from ICE import ICE as ICE
from Emotor import Emotor as Emotor
from Battery import Battery as Battery
import numpy as np

class Propulsion(object):

    ### CONSTRUCKTOR #####################################
    '''
    genereller Konstruktor
    initialisiert "Propeller"
    '''
    def __init__(self):
        self.verbose = True
        self.mass = {}
        self.Pshaft = self.Pmax/2 #Watt
        self.prop = Propeller()
        self.SoF = 1. #StateOfFuel
        self.V_tank = self.calcV_tank()
        if(self.fuel <= 0.):
            self.fuel = 0.01

    '''
    Rückgabe: benötigte Wellenleistung
    '''
    def calcNeededP(self, tas, F, etaProp, gearEff):
        Pshaft = float(tas * F) / (etaProp * gearEff)
        return Pshaft

    '''
    berechnet die nötige Wellenleistung
    Rückgabe: Pshaft, Nprop
    '''
    def calcPshaft(self, h, tas, F):
        Pold = self.Pshaft
        Nprop = self.calcNProp(Pold)
        Pshaft = self.Pshaft  #first guess
        eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=self.ice.NiceToNprop(self.ice.getNmax()), fixN=True)
        Pshaft = self.calcNeededP(tas, F, eta, c.GEAR_EFF_COMB)
        
        itegrationLimit = 0.1
        i = 0
        while(abs(Pshaft - Pold) > itegrationLimit):
            Pold = Pshaft
            etaOld = eta
            Nprop = self.calcNProp(Pshaft)
            eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=self.ice.NiceToNprop(self.ice.getNmax()), fixN=True)
            Pshaft = self.calcNeededP(tas, F, eta, c.GEAR_EFF_COMB)
            i += 1
            if(i > 100 and Pshaft > Pold):
                if(self.verbose):
                    print("endlosschleife -> " + str(Pshaft) + " != [" + str(Pold) + "]")
                Pshaft = (Pshaft + Pold) / 2
                eta = (eta + etaOld) / 2
                Nprop = self.calcNProp(Pshaft)
                eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=self.ice.NiceToNprop(self.ice.getNmax()), fixN=True)
                Pshaft = self.calcNeededP(tas, F, eta, c.GEAR_EFF_PARA)
                break;
                
        #SHOW FINAL DEBUG INFO
        self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=self.ice.NiceToNprop(self.ice.getNmax()), verbose=True)
        
        Nice = self.ice.NpropToNice(Nprop)
        self.Pshaft = Pshaft
        if(self.verbose):
            print('[Propulsion] needed Thrust = ' + str(F))
            print('[Propulsion] TAS = ' + str(tas))
            print('[Propulsion] Nice = ' + str(Nice))
            print('[Propulsion] Pice = ' + str(Pshaft))
            print('[Propulsion] Pshaft = ' + str(Pshaft))
            print('[Propulsion] Nprop = ' + str(Nprop))
            print('[Propulsion] etaProp = ' + str(eta))
            print('[Propulsion] needed Iterations for calculation of Pshaft: ' + str(i))
        return Pshaft, Nprop

    '''
    Rückgabe: Propellerdrehzahl, passen zu Leistung "P"
    '''
    def calcNProp(self, P):
        return self.ice.calcNProp(P)
    
    '''
    berechnet erst benötigte Wellenleistung und dann motorverbrauch und Betriebsparameter 
    Rückgabe: Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge
    '''
    def calcConsumption(self, h, tas, F, dt):
        Pshaft, Nprop = self.calcPshaft(h, tas, F)
        return self.calcMotorConsumption(h, Pshaft, Nprop, dt)
    
    '''
    berechnet motorverbrauch und Betriebsparameter
    Rückgabe: Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge
    '''
    def calcMotorConsumption(self, h, Pshaft, Nprop, dt):
        Nice = self.ice.NpropToNice(Nprop)
        Mshaft = Pshaft / N2omega(Nice)
        m_f, FF = self.ice.calcConsumption(Nice, Mshaft, dt)
        self.SoF -= self.ice.calcConsumtionToL(m_f, dt) / self.fuel #self.V_tank
        self.bat.SoC = 0.
        if(self.verbose):
            print('[Propulsion] FF = ' + str(FF))
            print('[Propulsion] SoC = ' + str(self.bat.SoC))
            print('[Propulsion] SoF = ' + str(self.SoF))
        return Pshaft, self.ice.NiceToNprop(Nice), m_f, 0, Pshaft, Nice, 0., 0.
    
    '''
    berechnet Tankvolumen, entweder begrenzt durch MTOM oder durch max Tankvolumen in "const.py"
    Rückgabe: Tankvolumen in l
    '''
    def calcV_tank(self):
        m_propulsion = self.calcPropulsionMasses()
        maxFuel = (c.MTOM - (m_propulsion + c.M_STRUCT + self.m_pl)) / c.DENS_FUEL
        return min(c.V_TANK_MAX, maxFuel)
    
    '''
    defaultmethode
    Berechnet Antriebskomponentenmassen und speichert sie in Klassenvariable "mass" (dict)
    Rückgabe: Gesamtantriebsmasse
    '''
    def calcPropulsionMasses(self):
        return float('inf')


### - Propulsion: Combustion - #######################################################

class Combustion(Propulsion):
    
    '''
    Konstruktor für konventionelles Flugzeug
    initialisiert "ICE", "Battery"
    '''
    def __init__(self, config):
        self.Pmax = config["pmax"]
        self.m_pl = config["mpl"]
        self.ice = ICE(Pmax=self.Pmax)
        self.bat = Battery(0., 0.)
        self.bat.SoC = 0.
        self.fuel = config["fuel"]
        super(self.__class__, self).__init__()
        
        #if(self.fuel > 0):
            #self.SoF = 1.
            #self.SoF = min(float(self.fuel) / self.V_tank, 1.)

    '''
    Berechnet Antriebskomponentenmassen und speichert sie in Klassenvariable "mass" (dict)
    Rückgabe: Gesamtantriebsmasse
    '''
    def calcPropulsionMasses(self):
        self.mass["ICE"] = (self.Pmax / c.PD_ice) + c.PD_ice_add
        self.mass["Bat"] = 0.
        return sum(self.mass.values())


### - Propulsion: SeriellHybrid - ####################################################
class SeriellHybrid(Propulsion):

    '''
    Konstruktor für seriellen Hybrid
    initialisiert "ICE", "Battery"
    '''
    def __init__(self, config):
        self.Pemot_max = config["pmax"]
        self.m_pl = config["mpl"]
        self.Pice_max = config["pice"]
        self.ice = ICE(Pmax=self.Pice_max)
        self.Pice_opt = self.ice.Popt
        
        self.m_pl = config["mpl"]
        self.Nice = self.ice.N_opt
        self.Mice = self.Pice_opt / N2omega(self.Nice)
        self.emot = Emotor(self.Pemot_max)
        self.Pmax = self.Pemot_max
        self.bat = Battery(config["cbat"], config["pbat"])
        self.fuel = config["fuel"]
        super(self.__class__, self).__init__()
        
    '''
    Rückgabe: Propellerdrehzahl, passen zu Leistung "P"
    '''
    def calcNProp(self, P):
        #return self.ice.calcNProp(self.Pmax)
        return self.ice.getNforP(P, self.Pmax, verbose=False) / c.GEAR_RATIO
    
    '''
    berechnet die nötige Wellenleistung
    Rückgabe: Pshaft, Nprop
    '''
    def calcPshaft(self, h, tas, F):
        Pold = self.Pshaft
        Nprop = self.calcNProp(Pold)
        Pshaft = self.Pshaft  #first guess
        eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=c.PROP_MAX_N, fixN=True)
        Pshaft = self.calcNeededP(tas, F, eta, c.GEAR_EFF_SERI)
        
        itegrationLimit = 0.1
        i = 0
        while(abs(Pshaft - Pold) > itegrationLimit):
            Pold = Pshaft
            etaOld = eta
            Nprop = self.calcNProp(Pshaft)
            eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=c.PROP_MAX_N, fixN=True)
            Pshaft = self.calcNeededP(tas, F, eta, c.GEAR_EFF_SERI)
            i += 1
            if(i > 100 and Pshaft > Pold):
                if(self.verbose):
                    print("endlosschleife -> " + str(Pshaft) + " != [" + str(Pold) + "]")
                Pshaft = (Pshaft + Pold) / 2
                eta = (eta + etaOld) / 2
                Nprop = self.calcNProp(Pshaft)
                eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=c.PROP_MAX_N, fixN=True)
                Pshaft = self.calcNeededP(tas, F, eta, c.GEAR_EFF_PARA)
                break;
            
        #SHOW FINAL DEBUG INFO
        self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=c.PROP_MAX_N, verbose=True)
        
        self.Pshaft = Pshaft
        if(self.verbose):
            print('[Propulsion] needed Thrust = ' + str(F))
            print('[Propulsion] TAS = ' + str(tas))
            print('[Propulsion] Pshaft = ' + str(Pshaft))
            print('[Propulsion] Nprop = ' + str(Nprop))
            print('[Propulsion] etaProp = ' + str(eta))
            print('[Propulsion] needed Iterations for calculation of Pshaft: ' + str(i))
        return Pshaft, Nprop
    
    '''
    berechnet erst benötigte Wellenleistung und dann motorverbrauch und Betriebsparameter 
    Rückgabe: Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge
    '''
    def calcConsumption(self, h, tas, F, dt):
        if(F <= 0 and self.bat.full()):
            Pshaft = 0
            Nprop = 0
        else:
            Pshaft, Nprop = self.calcPshaft(h, tas, F)
        return self.calcMotorConsumption(h, Pshaft, Nprop, dt)
    
    '''
    berechnet motorverbrauch und Betriebsparameter
    Rückgabe: Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge
    '''
    def calcMotorConsumption(self, h, Pshaft, Nprop, dt):
        if(Pshaft <= 0 and self.bat.full()):
            Nprop = 0
            Mprop = 0
            Nice = self.ice.getNidle()
            m_f = 0
            FF = 0
            Pice = 0
            We_ice = 0
            We = 0
            We_bat = 0
        else:
            Mprop = Pshaft / N2omega(Nprop)
            Nice = self.Nice
            Mice = self.Mice
            Pice = self.Pice_opt
            
            #needed elec-Power
            We = self.emot.calcConsumption(Nprop, Mprop, dt)
            We_ice = We
            if(self.bat.empty()):
                We_ice = We + self.bat.getMaxChargeRate()
            Pice = self.ice.calcNeededICEPower(We_ice, dt)
            Pice = min(Pice, self.ice.Pmax)
            Nice = self.ice.calcNIce(Pice)
            Mice = Pice / N2omega(Nice)
            m_f, FF = self.ice.calcConsumption(Nice, Mice, dt)
            We_ice = self.ice.calcGeneratorWe(Pice, dt)
            
            We_bat = We_ice - We

        self.SoF -= self.ice.calcConsumtionToL(m_f, dt) / self.fuel #self.V_tank
        We, Wcharge = self.bat.charge(We_bat, dt)
        
        if(self.verbose):
            print('[Propulsion] Nice = ' + str(Nice))
            print('[Propulsion] Pice = ' + str(Pice))
            print('[Propulsion] FF = ' + str(FF))
            print('[Propulsion] SoC = ' + str(self.bat.SoC))
            print('[Propulsion] SoF = ' + str(self.SoF))
        return Pshaft, Nprop, m_f, We, Pice, Nice, Pshaft, Wcharge
    
    
    '''
    Berechnet Antriebskomponentenmassen und speichert sie in Klassenvariable "mass" (dict)
    Rückgabe: Gesamtantriebsmasse
    '''
    def calcPropulsionMasses(self):
        self.mass["ICE"] = 0.
        self.mass["Gen"] = 0.
        self.mass["Umr2"] = 0.
        if(self.Pice_max > 0):
            self.mass["ICE"] = (self.Pice_max / c.PD_ice) + c.PD_ice_add
            self.mass["Gen"] = (self.Pice_max / c.PD_emot) + c.PD_emot_add
            self.mass["Umr2"] = (self.Pice_max / c.PD_inv) + c.PD_inv_add
        self.mass["EMot"] = (self.Pmax / c.PD_emot) + c.PD_emot_add
        self.mass["Umr"] = (self.Pmax / c.PD_inv) + c.PD_inv_add
        
        self.mass["Bat"] = 0.
        if(self.bat.C > 0.):
            self.mass["Bat"] = self.bat.calcMass()
        return sum(self.mass.values())


### - Propulsion: ParallelHybrid - ###################################################

class ParallelHybrid(Propulsion):

    '''
    Konstruktor für parallelen Hybrid
    initialisiert "ICE", "Battery"
    '''
    def __init__(self, config):
        self.Pmax = config["pmax"]
        self.m_pl = config["mpl"]
        self.Pice_max = config["pice"]
        self.ice = ICE(Pmax=self.Pice_max)
        self.Pice_opt = self.ice.Popt
        self.Nice = self.ice.N_opt
        self.Mice = self.Pice_opt / N2omega(self.Nice)
        self.Pemot_max = max(self.Pmax - self.Pice_max, 0)
        self.emot = Emotor(self.Pemot_max)
        self.Pmax = self.ice.Pmax + self.emot.Pmax
        self.bat = Battery(config["cbat"], config["pbat"])
        self.fuel = config["fuel"]
        super(self.__class__, self).__init__()

    '''
    Rückgabe: Propellerdrehzahl, passen zu Leistung "P"
    '''
    def calcNProp(self, P):
        if(self.Pice_max == 0):
            return self.ice.getNforP(P, self.Pmax, verbose=False) / c.GEAR_RATIO
        return self.ice.calcNProp(P)

    '''
    berechnet die nötige Wellenleistung
    Rückgabe: Pshaft, Nprop
    '''
    def calcPshaft(self, h, tas, F):
        gearEff = c.GEAR_EFF_PARA
        #Getriebewirkungsgrad bei Wegfall des Getriebes
        if(self.Pice_max == 0):
            gearEff = c.GEAR_EFF_SERI
        Pold = self.Pshaft
        Nprop = self.calcNProp(Pold)
        Pshaft = self.Pshaft  #first guess
        eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=self.ice.NiceToNprop(self.ice.getNmax()), fixN=True)

        Pshaft = self.calcNeededP(tas, F, eta, gearEff)
        
        itegrationLimit = 0.1
        i = 0
        eta_s = []
        Nprop_s = []
        Pshaft_s = []
        while(abs(Pshaft - Pold) > itegrationLimit):
            Pold = Pshaft
            etaOld = eta
            Nprop = self.calcNProp(Pshaft)
            #print("------- " + str(i) + " -------")
            eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=self.ice.NiceToNprop(self.ice.getNmax()), fixN=True)
            Pshaft = self.calcNeededP(tas, F, eta, gearEff)
            i += 1
            if(i > 100 and Pshaft > Pold):
                if(self.verbose):
                    print("keine konvergenz -> " + str(Pshaft) + " != [" + str(Pold) + "]")
                Pshaft = (Pshaft + Pold) / 2
                eta = (eta + etaOld) / 2
                Nprop = self.calcNProp(Pshaft)
                eta, Nprop = self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=c.PROP_MAX_N, fixN=True)
                Pshaft = self.calcNeededP(tas, F, eta, gearEff)
                break;

        #SHOW FINAL DEBUG INFO
        self.prop.calcEfficiency(h, tas, Nprop, Pshaft, NpropMax=self.ice.NiceToNprop(self.ice.getNmax()), verbose=True)
        
        Nice = self.ice.NpropToNice(Nprop)
        self.Pshaft = Pshaft
        if(self.verbose):
            print('[Propulsion] needed Thrust = ' + str(F))
            print('[Propulsion] TAS = ' + str(tas))
            print('[Propulsion] Nice = ' + str(Nice))
            print('[Propulsion] Pshaft = ' + str(Pshaft))
            print('[Propulsion] Nprop = ' + str(Nprop))
            print('[Propulsion] etaProp = ' + str(eta))
            print('[Propulsion] needed Iterations for calculation of Pshaft: ' + str(i))
        return Pshaft, Nprop
    
    '''
    berechnet erst benötigte Wellenleistung und dann motorverbrauch und Betriebsparameter 
    Rückgabe: Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge
    '''
    def calcConsumption(self, h, tas, F, dt):
        if(F <= 0 and self.bat.full()):
            Nprop = 0
            Pshaft = 0
        else:
            Pshaft, Nprop = self.calcPshaft(h, tas, F)
        return self.calcMotorConsumption(h, Pshaft, Nprop, dt)
    
    '''
    berechnet motorverbrauch und Betriebsparameter
    Rückgabe: Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge
    '''
    def calcMotorConsumption(self, h, Pshaft, Nprop, dt):
        Nice = self.ice.NpropToNice(Nprop)
        if(Pshaft <= 0 and self.bat.full()):
            Pice = 0
            Pemot = 0
            Nice = self.ice.getNidle()
            Mice = 0
            Memot = 0
        else:
            Pice = Pshaft
            Pemot = 0
            if(Pshaft > self.Pice_max):
                Pice = self.Pice_max
                Pemot = Pshaft - Pice
            Mice = Pice / N2omega(Nice)
            Memot = Pemot / N2omega(Nice)

        m_f, FF = self.ice.calcConsumption(Nice, Mice, dt)
        We_ice = 0
        We = self.emot.calcConsumption(Nice, Memot, dt)

        self.SoF -= self.ice.calcConsumtionToL(m_f, dt) / self.fuel #self.V_tank
        We, Wcharge = self.bat.charge(We_ice - We, dt)
        
        if(self.verbose):
            print('[Propulsion] Pice = ' + str(Pice))
            print('[Propulsion] Pemot = ' + str(Pemot))
            print('[Propulsion] FF = ' + str(FF))
            print('[Propulsion] SoC = ' + str(self.bat.SoC))
            print('[Propulsion] SoF = ' + str(self.SoF))
        return Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge


    '''
    Berechnet Antriebskomponentenmassen und speichert sie in Klassenvariable "mass" (dict)
    Rückgabe: Gesamtantriebsmasse
    '''
    def calcPropulsionMasses(self):
        self.mass["ICE"] = 0.
        if(self.Pice_max > 0):
            self.mass["ICE"] = (self.Pice_max / c.PD_ice) + c.PD_ice_add
        
        self.mass["Umr"] = 0.
        self.mass["EMot"] = 0.
        if(self.Pemot_max > 0):
            self.mass["Umr"] = (self.Pemot_max / c.PD_inv) + c.PD_inv_add
            self.mass["EMot"] = (self.Pemot_max / c.PD_emot) + c.PD_emot_add
        
        self.mass["Bat"] = 0.
        if(self.bat.C > 0.):
            self.mass["Bat"] = self.bat.calcMass()
        return sum(self.mass.values())


### - MAIN - ####################################################
if __name__ == "__main__":
    ser = SeriellHybrid()
    par = ParallelHybrid()
    #Pshaft, Nshaft, m_f, We = ser.calcConsumption(2000, 50, 200, 60**2)
    kon = Combustion()
    #print('Nprop: ' + str(kon.calcNProp(40000)))
    print('done!')