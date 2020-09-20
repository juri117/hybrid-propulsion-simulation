'''
Created on 01.12.2015

@author: Juri Bieler
'''

from Flight import Flight
from Pilot import Pilot
from Plane import Plane
from Battery import Battery
from Ini import Ini as Ini
from Mission import Mission
import const as c
import generalUtils as generalUtils
from optiPlot import *

import math
import pylab as plt
import time
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy as np
from matplotlib import rc
from matplotlib.patches import Rectangle
from queue import Queue
from threading import Thread

class Optimizer(object):

    ### CONSTRUCKTOR #####################################
    '''
    initialisiert "Ini"
    Definiert Iterationsgrenzen
    '''
    def __init__(self, mission, verbose=True, savePlots=True):
        self.verbose = verbose
        self.savePlots = savePlots
        self.COMB = "combustion"
        self.SERI = "serialHybrid"
        self.PARA = "parallelHybrid"
        
        self.mission = mission
        self.ini = Ini(self.mission.NAME)
        self.m_pl = c.M_PL
        if(self.ini.keyExists("general", "mpl")):
            self.m_pl = float(self.ini.get("general", "mpl"))
            
        self.P_increment = c.P_INCREMENT

        #Liste von lauffähigen Flügen
        self.flights = []
        
        self.NOTRUNNABLE = 0
        self.RUNNABLE = 1
        self.OPTIMIZED = 2


    '''
    Optimiert iterativ einen Flug mit gegebener Konfiguration, bis genau die richtige Menga an Brennstoff mitgeführt wird
    Rückgabe: Instanz von "Flight"
    '''
    def flyWithNeededFuel(self, mission, propulsionType, config):
        limit = 0.1
        
        pic = Pilot(mission, propulsionType=propulsionType, config=config)
        pic.shutUp(savePlot=False)
        fl = pic.fly()
        while(abs(fl.calcUsedFuel() - config["fuel"]) > limit):
            config["fuel"] = fl.calcUsedFuel()
            pic = Pilot(mission, propulsionType=propulsionType, config=config)
            pic.shutUp(savePlot=False)
            fl = pic.fly()
            fl.checkFlight()
        
        if(self.verbose):
            print("Erfolg: " + str(fl.fails))
            print(propulsionType + " => " + str(config))
        return fl


    '''
    Führt die Optimierung für ein konventionelles Flugzeug durch
    Rückgabe: Erfolg [True|False]
    '''
    def runFlight(self):
        t0 = time.time()
        self.flights = []
        Pmax = 60000. #float(self.ini.get(c.COMB, "pmax"))
        fuel = 0. #float(self.ini.get(c.COMB, "fuel"))
        config = {"pmax":Pmax,"fuel":fuel,"mpl":self.m_pl}
        config, flight, suc = self.optimizeConfig(self.COMB, config)
        if(not suc):
            self.saveConfig(self.COMB, None)
        else:
            self.saveConfig(self.COMB, flight)
        t1 = time.time()
        print('[Optimizer] Runtime [s]: %f' %(t1-t0))
        return suc


    '''
    Fügt eine Instanz der Klasse "Flight" sortiert nach der Große des Verbrennungsmotors einer Liste von Flügen hinzu
    '''
    def addFlight(self, flight):
        print("config = " + str({k:round(v, 2) for k, v in flight.config.copy().items()}))
        i = 0
        for i in range(0, len(self.flights), 1):
            if(flight.config["pice"] < self.flights[i].config["pice"]):
                self.flights.insert(i, flight)
                return;
        self.flights.append(flight)


    '''
    Prüft ob eine Konfiguration optimiert bzw. flugfähig ist
    Vergleicht Eingangsparameter mit den tatsächlich benötigten Resourcen
    Rückgabe: [NOTRUNNABLE = 0|RUNNABLE = 1|OPTIMIZED = 2]
    '''
    def flightIsOptimized(self, config, flight):
        Pneeded, Fuel_needed, Cbat_needed, Pbat_needed = generalUtils.calcNeededConfig(flight)

        if(config["pmax"] * (1. + c.ITERATION_LIMIT) >= Pneeded
           and config["fuel"] * (1. + c.ITERATION_LIMIT) >= Fuel_needed
           and config["cbat"] * (1. + c.ITERATION_LIMIT) >= Cbat_needed
           and config["pbat"] * (1. + c.ITERATION_LIMIT) >= Pbat_needed):
            if(self.checkDeviation(config["pmax"], Pneeded) <= c.ITERATION_LIMIT
               and self.checkDeviation(config["fuel"], Fuel_needed) <= c.ITERATION_LIMIT
               and self.checkDeviation(config["cbat"], Cbat_needed) <= c.ITERATION_LIMIT
               and self.checkDeviation(config["pbat"], Pbat_needed) <= c.ITERATION_LIMIT):
                return self.OPTIMIZED
            return self.RUNNABLE
        return self.NOTRUNNABLE


    '''
    Führt eine Flugsimulation für die Antriebskonfiguration "propulsionType" mit den Parameter aus der "config" (Instanz von Dict) durch,
    iteriert so lange über max. Triebwerksleistung, Brennstoff, Batteriekappazität und Batterieleistung, bis "flightIsOptimized" "OPTIMIZED" zurück gibt
    Rückgabe: Dict mit Konfiguration, Instanz von "Flight", Erfolg [True|False]
    '''
    def optimizeConfig(self, propulsionType, config):
        if(propulsionType == self.COMB):
            config["cbat"] = 0.
            config["pbat"] = 0.
            config["fuel"] = 0.
            config["cbat"] = 1.
            config["pbat"] = 1.

        pic = Pilot(self.mission, propulsionType=propulsionType, config=config)
        pic.shutUp(savePlot=False)
        flight = pic.fly()
        
        if(flight.abort or config["cbat"] < 0. or config["cbat"] == float('inf')):
                print("[optimizeConfig] FEHLER: Flug abgebrochen: " + flight.getErrorMessage())
                print("[optimizeConfig] Batterie wird zurueck gesetzt!")
                #imaginized reset value
                config["cbat"] = 100.
                config["fuel"] = 10.
                pic = Pilot(self.mission, propulsionType=propulsionType, config=config)
                pic.shutUp(savePlot=False)
                flight = pic.fly()
                if(flight.abort):
                    print("[optimizeConfig] FEHLER: Flug abgebrochen: " + flight.getErrorMessage())
                    return config, flight, False

        configs = []
        iterations = 0
        alternating = False
        
        while(self.flightIsOptimized(config, flight) != self.OPTIMIZED):
            #print("->pmax: " + '{:<10}'.format(round(config["pmax"], 1)) + "\t" + "fuel: " + '{:<10}'.format(round(config["fuel"], 1)) + "\t" + "cbat: " + '{:<10}'.format(round(config["cbat"], 1)) + "\t" + "pbat: " + '{:<10}'.format(round(config["pbat"], 1)) + "\t" + str(config))

            if(self.flightIsOptimized(config, flight) == self.RUNNABLE):
                print("->pmax: " + str(round(config["pmax"], 1)).ljust(10) + "\t" + "fuel: " + str(round(config["fuel"], 2)).ljust(10) + "\t" + "cbat: " + str(round(config["cbat"], 1)).ljust(10) + "\t" + "pbat: " + str(round(config["pbat"], 1)).ljust(10) + "\t" + str(config))
                #auf Wiederholung der Ergebnisse prüfen
                if(self.configContains(configs, config)):
                    print("[optimizeConfig] FEHLER 111: Loesung alterniert")
                    alternating = True
                configs.append(config.copy())
            
            iterations += 1
            if(iterations > 150):
                print("[optimizeConfig] FEHLER 112: Flug abgebrochen: IterationsLimit: " + str(iterations))
                bestConf = self.findBestConfig(configs)
                if(bestConf != None):
                    pic = Pilot(self.mission, propulsionType=propulsionType, config=bestConf)
                    pic.shutUp(savePlot=False)
                    flight = pic.fly()
                    self.addFlight(flight)
                    return bestConf, flight, True
                return bestConf, flight, False

            #Batterie zu klein
            if(config["cbat"] < 0.):
                flight.fail("[optimizeConfig] FEHLER: Batterie zu klein -> " + str(config["cbat"]))
                print("[optimizeConfig] FEHLER 113: Flug abgebrochen: " + flight.getErrorMessage())
                return config, flight, False
            
            #zu schwer
            if(max(flight.m) > c.MTOM * (1 + c.MTOM_TOLERANCE)):
                flight.fail("[optimizeConfig] FEHLER: mtom ueberschritten")
                print("[optimizeConfig] FEHLER 114: Flug abgebrochen: " + flight.getErrorMessage())
                return config, flight, False
            
            #MTOM erreicht aber nicht genug Sprit
            if(max(flight.m) >= c.MTOM * (1 + c.MTOM_TOLERANCE)):
                if(pic.p.propulsion.V_tank * max(flight.SoF) < config["fuel"] and config["cbat"] <= 1.):
                    flight.fail("[optimizeConfig] FEHLER: mtom erreicht, aber Sprit fehlt")
                    print("[optimizeConfig] FEHLER 115: Flug abgebrochen: " + flight.getErrorMessage())
                    return config, flight, False

            Pneeded, Fuel_needed, Cbat_needed, Pbat_needed = generalUtils.calcNeededConfig(flight)
            
            if(alternating):
                #ALTERNIERT - die eingangsparameter werden nun nicht mehr mit dem gesamten Delta zur ausgabe addiert
                dela = 1./2.
                if(iterations > 50):
                    dela = 1./4.
                if(iterations > 100):
                    dela = 1./10.
                print("DELTA * " + str(dela) + " (" + str(iterations) + ")")
                
                Pmax_diff = self.checkDeviation(config["pmax"], Pneeded)
                Fuel_diff = self.checkDeviation(config["fuel"], Fuel_needed)
                Cbat_diff = self.checkDeviation(config["cbat"], Cbat_needed)
                Pbat_diff = self.checkDeviation(config["pbat"], Pbat_needed)
                
                if(Fuel_diff > c.ITERATION_LIMIT):
                    config["fuel"] = config["fuel"] + ((Fuel_needed - config["fuel"]) * dela)
                    print("fuel " + str(round(Fuel_diff, 4)))
                else:
                    if(Pmax_diff >= Fuel_diff and Pmax_diff >= Cbat_diff and Pmax_diff >= Pbat_diff):
                        config["pmax"] = config["pmax"] + ((Pneeded - config["pmax"]) * dela)
                        print("pmax " + str(round(Pmax_diff, 4)))
                        
                    if(Cbat_diff >= Pmax_diff and Cbat_diff >= Fuel_diff and Cbat_diff >= Pbat_diff):
                        config["cbat"] = config["cbat"] + ((Cbat_needed - config["cbat"]) * dela)
                        print("cbat " + str(round(Cbat_diff, 4)))
                    
                    if(Pbat_diff >= Pmax_diff and Pbat_diff >= Fuel_diff and Pbat_diff >= Cbat_diff):
                        config["pbat"] = config["pbat"] + ((Pbat_needed - config["pbat"]) * dela)
                        print("pbat " + str(round(Pbat_diff, 4)))
            else:
                config["pmax"] = Pneeded
                config["cbat"] = Cbat_needed
                config["fuel"] = Fuel_needed
                config["pbat"] = Pbat_needed
            
            pic = Pilot(self.mission, propulsionType=propulsionType, config=config)
            pic.shutUp(savePlot=False)
            flight = pic.fly()
            if(flight.abort):
                print("[optimizeConfig] FEHLER 116: Flug abgebrochen: " + flight.getErrorMessage())
                return config, flight, False

        if((pic.p.propulsion.V_tank * max(flight.SoF)) - flight.calcUsedFuel() < (-1) * c.FUELTANK_TOLERANCE):
            flight.fail("[optimizeConfig] FEHLER: mtom erreicht aber zu wenig Sprit")
            print("[optimizeConfig] FEHLER 117: Flug abgebrochen: " + flight.getErrorMessage())
            return config, flight, False
        
        if(max(flight.m) > c.MTOM):
            flight.fail("[optimizeConfig] FEHLER: mtom ueberschritten")
            print("[optimizeConfig] FEHLER 118: Flug abgebrochen: " + flight.getErrorMessage())
            return config, flight, False

        Pneeded, Fuel_needed, Cbat_needed, Pbat_needed = generalUtils.calcNeededConfig(flight)
        config["pmax"] = Pneeded
        config["fuel"] = Fuel_needed
        config["cbat"] = Cbat_needed
        self.addFlight(flight)
        return config, flight, True

    
    '''
    führt seriell Hybriden Flug mit variierendem Hybrdidisierungsgrad aus und findet die Minima des Energieverbrauches
    Rückgabe: minimum exisiterit [True|False]
    '''
    def runSeriellFlight(self):
        self.flights = []
        return self.optimizeHybridFlight(self.SERI)


    '''
    führt parallel Hybriden Flug mit variierendem Hybrdidisierungsgrad aus und findet die Minima des Energieverbrauches
    Rückgabe: minimum exisiterit [True|False]
    '''
    def runParallelFlight(self):
        self.flights = []
        return self.optimizeHybridFlight(self.PARA)


    '''
    führt hybride Flüge mit der Antriebskonfiguration "propulsionType" aus, wobei die Größe des Verbrennungsmotors
    entsprechend der Vorgabe in "Range" variiert wird
    die Ergebnisse der Simulationen werden bei erfolg in eine Liste in dieser Klasse gespeichert
    '''
    def runhybridFlights(self, propulsionType, Range=range(0, 80000, 1000), plotIt=True, multiTask=False):
        config = {"pmax":50000, "pice":0, "cbat":100, "fuel":50, "mpl":self.m_pl}

        for i in Range:
            print("--- Pice = " + str(i) + " ---")
            config = config.copy()
            config["pice"] = i
            configNew, _, suc = self.optimizeConfig(propulsionType, config)
            if(suc):
                config = configNew


    '''
    führt für einen Hybridantrieb der Antriebskonfiguration "propulsionType" erst Simulationen mit verschiedenen
    Verbrennungsmotorgrößen aus, wobei diese von 0 bis zur Maximalleistung des konventionellen Antriebs in
    "P_increment"-Schritten (in Konstruktor definiert) erhäht wird,
    anschließend werden die entstehenden lokalen Minima mit "optimizeHybridFlight_rec" genauer untersucht
    Rückgabe: Minimum gefunden [True|False]
    '''
    def optimizeHybridFlight(self, propulsionType):
        P_increment = self.P_increment

        Pice_max = c.P_ICE_MAX
        
        Pconv_max = float(self.ini.get(c.COMB, "pmax"))
        if(not math.isnan(Pconv_max)):
            Pice_max = Pconv_max
        
        Prange = list(range(0, int(Pice_max), int(P_increment)))
        Prange.append(int(Pice_max))
        
        self.runhybridFlights(propulsionType, Range=Prange, plotIt=False)
        
        P_increment = self.P_increment / 10.
        
        flightMins = self.findLocalMinima(self.flights)
        configBest = {"fuel":float("inf")}
        for f in flightMins:
            config = self.optimizeHybridFlight_rec(propulsionType, Pice_max, f.config, P_increment)
            if(config["fuel"] < configBest["fuel"]):
                configBest = config

        bestFlight = self.saveBestFlight(propulsionType)
        
        if(bestFlight != None):
            return True
        return False
      
            
    '''
    untersucht die Antriebskonfiguration mit minimalem Energiebedarf in "config" der Antriebart "propulsionType"
    da zuvor in relativ grobem Raster vorgegangen wurde, wird hier der Verlauf rechts und links des Minimums untersucht,
    um dessen exakte Postition zu bestimmen
    Rückgabe: Instanz von Dict "config" mit der Konfiguration mit dem minimalen Energiebedarf
    '''   
    def optimizeHybridFlight_rec(self, propulsionType, Pice_max, config, P_increment):
        if(self.verbose):
            print("[optimizeHybridFlight_rec] ----------------------")
            print("[optimizeHybridFlight_rec] P_increment = " + str(P_increment))
            print("[optimizeHybridFlight_rec] Pice = " + str(config["pice"]))
            print("[optimizeHybridFlight_rec] fuel = " + str(config["fuel"]))
        if(P_increment < c.P_INCREMENT_LIMIT or config["pice"] >= Pice_max or config["pice"] == 0):
            return config

        #check lower
        lowerPice = config["pice"] - P_increment
        if(lowerPice < c.P_MAX_ICE_MIN):
            lowerPice = 0
        if(not self.flightHasBeenRun(lowerPice)):
            configLow = {"pmax":config["pmax"], "pice":lowerPice, "cbat":config["cbat"], "fuel":config["fuel"], "mpl":self.m_pl}
            configLow, flight, suc = self.optimizeConfig(propulsionType, configLow)
            if(suc):
                if(self.calcEnergy(config) > self.calcEnergy(flight.config)):
                    return self.optimizeHybridFlight_rec(propulsionType, Pice_max, configLow, P_increment)
        
        #check higher
        higherPice = config["pice"] + P_increment
        if(higherPice > Pice_max):
            higherPice = Pice_max
        if(not self.flightHasBeenRun(higherPice)):
            configHigh = {"pmax":config["pmax"], "pice":higherPice, "cbat":config["cbat"], "fuel":config["fuel"], "mpl":self.m_pl}
            configHigh, flight, suc = self.optimizeConfig(propulsionType, configHigh)
            if(suc):
                if(self.calcEnergy(config) > self.calcEnergy(flight.config)):
                    return self.optimizeHybridFlight_rec(propulsionType, Pice_max, configHigh, P_increment)
        
        #reduce P_increment
        return self.optimizeHybridFlight_rec(propulsionType, Pice_max, config, P_increment/10.)


    '''
    Führt die Untersuchung und Optimierung für konventionellen Antrieb, den seriellen und den parallelen Hybridantrieb aus
    und speichert die Ergebnisse im results Ordner der Mission
    '''
    def saveFlights(self, verbose=True):
        t0 = time.time()
        self.saveFlight(c.COMB)
        
        self.saveFlight(c.SERI)
        generalUtils.writeFlights2CSV(c.SERI, self.mission, self.flights)
        flightsSeri = self.flights.copy()
        
        self.saveFlight(c.PARA)
        generalUtils.writeFlights2CSV(c.PARA, self.mission, self.flights)
        
        #plot all
        if(self.savePlots):
            plotHybrids(self.mission, flightsSeri, self.flights, self.P_increment, verbose=self.verbose)         

        t1 = time.time()
        print('[Optimizer] Runtime [s]: %f' %(t1-t0))


    '''
    Führt die Untersuchung und Optimierung für die Antriebsart "propulsionType" aus
    und speichert die Ergebnisse im results Ordner der Mission
    '''
    def saveFlight(self, propulsionType, verbose=True):
        self.flights = []
        #combustion
        if(propulsionType == self.COMB):
            suc = self.runFlight()
        #seriell
        if(propulsionType == self.SERI):
            suc = self.runSeriellFlight()
        #parallel
        if(propulsionType == self.PARA):
            suc = self.runParallelFlight()
        #saveFlightPlot
        if(suc):
            pic = Pilot(self.mission, propulsionType=propulsionType, config={})
            pic.shutUp(savePlot=self.savePlots)
            flight = pic.fly()



    ### UTILS #####################################
    
        
    '''
    Rückgabe: prozentuale Abweichung des value1 von value2 [0. - 1.], wenn einer der Werte 0 und der andere ungleich null ist -> 1.
    '''
    def checkDeviation(self, value1, value2):
        if(value1 == 0):
            if(value2 == 0):
                return 0.
            return 1.
        return abs((value1 - value2) / value1)
    
    
    '''
    sucht aus der in dieser Klasse gespeicherten Liste aus erfolgreich durchgeführten Flücken den mit dem minimalen Energiebedarf heraus
    Rückgabe: Instanz von "Flight"
    '''
    def saveBestFlight(self, propulsionType):
        bestFlight = None
        if(len(self.flights) > 0):
            bestFlight = self.flights[0]
            for f in self.flights:
                if(self.calcEnergy(f.config) < self.calcEnergy(bestFlight.config)):
                    bestFlight = f
        self.saveConfig(propulsionType, bestFlight)
        return bestFlight
        

    '''
    speichert die Konfiguration des Fluges "flight" in der config.ini des Missionsverzeichnisses
    '''
    def saveConfig(self, propulsionType, flight):
        if(flight == None):
            config = {"pmax":"nan", "pice":"nan", "cbat":"nan", "pbat":"nan", "fuel":"nan"}
        else:
            config = flight.config
        for key, value in iter(config.items()):
            if((propulsionType != self.COMB or (key != "cbat") or (key != "pbat")) and key != "mpl"):
                self.ini.set(propulsionType, key, value)
        
        #save general info
        if(propulsionType == c.COMB):
            dist = self.mission.calcDistance()
            if(math.isnan(dist)):
                dist = flight.distance
            self.ini.set("output", "distance", dist)
            self.ini.set("output", "flighttime", generalUtils.secToTimeStr(max(self.mission.sec)))
        else:
            if(not math.isnan(float(config["cbat"])) and config["cbat"] > 0.):
                self.ini.set("output", "c_" + propulsionType, str(min(flight.chargePower) / flight.config["cbat"]) + ", " + str(max(flight.chargePower) / flight.config["cbat"]))


    '''
    prüft, ob eine Antriebskonfiguration "conf" in der Liste von Konfigurationen "Configs" entahlten ist
    Rückgabe: enthalten [True|False]
    '''
    def configContains(self, configs, conf):
        for c in configs:
            match = True
            for key, value in iter(c.items()):
                if(not math.isclose(conf[key], value)):
                    match = False
            if(match):
                return True
        return False
    
    '''
    findet lokale minima des Energiebedarfs in der in dieser Klasse gespeicherten Liste von Flügen
    Rückgabe: Liste von Flügen mit lokal-minimalem Energiebedarf
    '''
    def findLocalMinima(self, flights):
        left = float('inf')
        FlightMins = []
        for i in range(0, len(flights), 1):
            right = float('inf')
            if(i + 1 < len(flights)):
                right = self.calcEnergy(flights[i+1].config)
            cost = self.calcEnergy(flights[i].config)
            if(cost < left and cost < right):
                FlightMins.append(flights[i])
            left = cost
        return FlightMins

    '''
    prüft ob ein Flug mit der Verbrennungsmotorleistung "Pice" bereits in der Liste der erfolgreich durchgeführten Flüge ist
    '''
    def flightHasBeenRun(self, Pice):
        for f in self.flights:
            if(math.isclose(f.config["pice"], Pice)):
                return True
        return False
    
        '''
    Rückgabe: Energiebedarf in kWh der Antriebskonfiguration "config"
    '''
    def calcEnergy(self, config):
        if("cbat" in config):
            return (config["fuel"] * c.ENERGY_FUEL_L) + (config["cbat"] / 1000)
        else:
            return config["fuel"] * c.ENERGY_FUEL_L
    
    
    '''
    Rückgabe: Instanz von dict "config" mit dem geringsten Energiebedarf aus der Liste "configs"
    '''
    def findBestConfig(self, configs):
        if(len(configs) == 0):
            return None
        best = configs[0]
        for con in configs:
            best = self.getBetterConfig(best, con)
        return best

    '''
    Rückgabe: die config, die weniger Energie benötigt
    '''
    def getBetterConfig(self, config1, config2):
        cost1 = self.calcEnergy(config1)
        cost2 = self.calcEnergy(config2)
        if(cost1 < cost2):
            return config1
        return config2
    
    
### MAIN #####################################

if __name__ == "__main__":
    mission = Mission('mission02_01')
    verbose = False
    opt = Optimizer(mission, verbose=verbose)
    opt.saveFlights()
    #opt.saveFlight(c.PARA)
    #opt.optimizeConfig(c.SERI, conf)
    print("ende!")