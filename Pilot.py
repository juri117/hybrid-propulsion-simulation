'''
Created on 01.12.2015

@author: Juri Bieler
'''

from Plane import Plane
from Mission import Mission
from Flight import Flight
import pilotPlot as pPlot
from Ini import Ini as Ini
import const as c

import numpy as np
import generalUtils as generalUtils
import math
import time

class Pilot:

    ### CONSTRUCKTOR #####################################
    '''
    initialisiert die Klassen Ini, Mission, Plane
    propulsionType: combustion, serialHybrid, parallelHybrid
    '''
    def __init__(self, mission, propulsionType="combustion", config={}):
        self.verbose = True
        self.savePlot = True
        self.propulsionType = propulsionType
        ini = Ini(mission.NAME)
        self.stepWidth = int(ini.get("general", "step"))
        self.p = Plane(mission, propulsionType=propulsionType, config=config)
        self.m = mission


    ### FUNCTIONS #####################################
    '''
    deaktiviert die Konsolen-Text-Ausgabe für alle Instanzen in dieser Klasse
    savePlot definiert, ob die Plots für einen einzelnen Flug gespeichert werden sollen
    '''
    def shutUp(self, savePlot=False):
        self.verbose = False
        self.savePlot = savePlot
        self.p.propulsion.verbose = False
        self.p.propulsion.ice.verbose = False
        self.p.propulsion.bat.verbose = False
        self.p.w.verbose = False


    '''
    startet die Simulation eines Fluges
    Rückgabe: instanz von Flight
    '''
    def fly(self):
        t0 = time.time()
        
        startWeight = self.p.calcActWeight()
        
        # in flight werden alle wichtigen Parameter über dem Flugverlauf gespeichert
        flight = Flight(self.stepWidth, startWeight, self.p.config)
        flight.mass = self.p.propulsion.mass
        flight.distance = self.m.calcDistance()
        
        # prüft ob das MTOM + Toleranz eingehalten wird
        if(startWeight > c.MTOM * (1 + c.MTOM_TOLERANCE)):
            flight.fail("[Pilot] MTOM ueberschritten")
            if(self.verbose):
                print(flight.fail)
            return flight
        
        # prüft ob das maximale Tankvolumen eingehlaten wird
        if(self.p.config["fuel"] > (1 + c.FUELTANK_TOLERANCE) * self.p.propulsion.V_tank):
            flight.fail("[Pilot] Tankvolumen reicht nicht aus.")
            if(self.verbose):
                print(flight.getErrorMessage())
            return flight
        
        missionState = 0;
        timeRange = range(0, max(self.m.sec) + 1, self.stepWidth)
        
        # Berechnung des Starts
        flight = self.takeOff(flight)

        # Berechnung des Fluges in einer Schleife, der Zietschritt ist in der config.ini der Mission definiert
        for t in timeRange:
            while(t >= self.m.sec[min(missionState + 1, len(self.m.sec) - 1)] and missionState < len(self.m.sec) - 1):
                missionState += 1
            
            if(self.verbose):
                print ('---TIME: ' + str(t) + '------------------')

            if(t == -1):
                print("break")
                
            if(t == -1):
                print("break")
                
            if(t == -1):
                print("break")

            mLast = flight.getLastMass()
            SoCLast = float(self.p.propulsion.bat.SoC)
            SoFLast = float(self.p.propulsion.SoF)
            H = float(self.m.getHeight(t))
            r = float(self.m.getRadius(t))
            w_g = float(self.m.getW_g(t))

            if(self.m.vSpecial[missionState] == ''):
                tas = float(self.m.getTas(t))
            else:
                tas = self.findSpeedSetting(mLast, H, w_g, r, self.m.vSpecial[missionState])

            u_g = math.sqrt(abs(tas**2 - w_g**2))

            F = self.p.calcNeededThrust(flight.getLastMass(), H, u_g, w_g, r)

            Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge = self.p.calcConsumption(H, tas, F, self.stepWidth)
            
            m = mLast - m_f
                
            SoC = float(self.p.propulsion.bat.SoC)
            SoF = float(self.p.propulsion.SoF)
                
            if(t == 0):
                m = mLast
                SoC = SoCLast
                SoF = SoFLast
                
            if(self.verbose):
                print('[Pilot] Mass [kg]: ' + str(m))
                print('[Pilot] Consumption [l/h]: ' + str(self.p.propulsion.ice.calcConsumtionToL_h(m_f, self.stepWidth)))
            
            flight.addStep(t, m, H, tas, r, F, Pshaft, Nprop, Pice, Nice, Pemot, SoC, SoF, We, Wcharge)

        t1 = time.time()
        flight.runTime = t1-t0

        # Konsolen-Ausgabe
        if(self.verbose):
            self.printReport(flight)
        
        # Speichern der Plots
        if(self.savePlot or self.verbose):
            self.plotFlight(flight)
            sfcMapPath = None
            propMapPath = None
            if(self.savePlot):
                fileName = self.m.NAME + "_" + self.propulsionType
                if(self.propulsionType == c.SERI or self.propulsionType == c.PARA):
                    fileName = fileName + "_Pice_" + str(int(flight.config["pice"]))
                sfcMapPath = self.m.getResultsPath() + "sfcMap_" + fileName
                propMapPath = self.m.getResultsPath() + "propMap_" + fileName
            self.p.propulsion.ice.plotEfficiencyMap(flight, fileName=sfcMapPath, showPlot=self.verbose, propType=self.propulsionType)
            self.p.propulsion.prop.plotEfficiencyMap(flight, fileName=propMapPath , showPlot=self.verbose)
  
        return flight


    '''
    gibt eine Zusammenfassung der relevanten Parameter aus
    zusätzlich wird angezeigt, ob die Mission mit der in config vorgegeben Konfiguration geflogen werden konnte
    '''
    def printReport(self, flight):
        print ('----------------------------')
        print ('---FINISH!------------------')
        print('[Pilot] config: ' + str(self.p.config))
        print('---POWER:')
        print('[Pilot] Pmax [W]: ' + str(self.p.propulsion.Pmax))
        print('[Pilot] Pmax_ice [W]: ' + str(self.p.propulsion.ice.Pmax))
        print('[Pilot] Pmax needed [W]' + str(max(flight.Pshaft)))
        print('[Pilot] displacement [l]: ' + str(self.p.propulsion.ice.V_h))
        print('[Pilot] ICE Efficiency [-]: ' + str(self.p.propulsion.ice.eff_total))
        
        print('---MASS:')
        print('[Pilot] PayloadMass [kg]: ' + str(flight.config["mpl"]))
        print('[Pilot] TakeOffMass [kg]: ' + str(flight.m[0]))
        print('[Pilot] LandingMass [kg]: ' + str(flight.m[-1]))
        print('[Pilot] WeightLoss [kg]: ' + str(flight.m[0] - flight.m[-1]))
        print('---FUEL:')
        print('[Pilot] max Fuel [L]: ' + str(self.p.config["fuel"]))
        print('[Pilot] used Fuel [L]: ' + str(flight.calcUsedFuel()))
        print('[Pilot] start SoF [%]: ' + str(flight.SoF[0]))
        print('[Pilot] end SoF [%]: ' + str(flight.SoF[-1]))
        print('---BAT:')
        print('[Pilot] elekt. Arbeit aus Batterie [W]: ' + str(np.sum(flight.We)))
        print('[Pilot] start SoC [%]: ' + str(flight.SoC[0]))
        print('[Pilot] end SoC [%]: ' + str(flight.SoC[-1]))
        print('[Pilot] min. SoC [%]: ' + str(min(flight.SoC)))
        chargeRate, Pcharge = self.p.propulsion.bat.calcCharegeRates(flight)
        print('[Pilot] installed. ChargeRate [kW]: ' + str(flight.config["pbat"] / 1000.))
        print('[Pilot] min. ChargeRate [C]: ' + str(min(chargeRate)) + " (" + str(min(flight.chargePower) / 1000.) + " kW)")
        if(flight.mass["Bat"] > 0.):
            print('[Pilot] PowerDensity [W/kg]: ' + str(max(abs(min(flight.chargePower)), abs(max(flight.chargePower))) / flight.mass["Bat"]))
        else:
            print('[Pilot] PowerDensity [W/kg]: 0')
        print('[Pilot] max. ChargeRate [C]: ' + str(max(chargeRate)) + " (" + str(max(flight.chargePower) / 1000.) + " kW)")
        
        print("[Pilot] Bat ED: " + str(self.p.propulsion.bat.Ed_ideal))
        print("[Pilot] Bat PD: " + str(self.p.propulsion.bat.Pd_ideal))
        
        print('---FLIGHT:')
        print('[Pilot] takeOffDistance [m]: ' + str(self.m.takeOffDistance))
        print('[Pilot] Distance: [km]: ' + str(flight.calcDistance() / 1000))
        print('[Pilot] Distance ueber Grund: [km]: ' + str(self.m.calcDistance() / 1000))
        print('[Pilot] total Flighttime: ' + generalUtils.secToTimeStr(max(flight.t)))
        print('[Pilot] timeStep [s]: ' + str(self.stepWidth))
        print('[Pilot] Runtime [s]: %f' %(flight.runTime))
        
        print('---ITERATION-DEBUG (installed, needed, diff, running):')
        
        pmax_need, fuel_need, cbat_need, pbat_need = generalUtils.calcNeededConfig(flight)
        print('[Pilot] valu [u]: ' + 'installed'.ljust(10) + '\t' + 'needed'.ljust(10) + '\t' + 'diff[%]'.ljust(10) + '\t-> ' + 'running')
        print('------------------' + '----------'.ljust(10) + '\t' + '----------'.ljust(10) + '\t' + '----------'.ljust(10) + '\t---' + '----------')
        pmax_in = round(flight.config["pmax"], 0)
        print('[Pilot] Pmax [W]: ' + str(pmax_in).ljust(10) + '\t' + str(pmax_need).ljust(10) + '\t' + str(round(100*(pmax_in - pmax_need) / pmax_in, 2)).ljust(10) + '\t-> ' + str(pmax_in >= pmax_need))
        fuel_in = round(flight.config["fuel"], 2) + 0.00001
        print('[Pilot] fuel [l]: ' + str(fuel_in).ljust(10) + '\t' + str(fuel_need).ljust(10) + '\t' + str(round(100*(fuel_in - fuel_need) / fuel_in, 2)).ljust(10) + '\t-> ' + str(fuel_in >= fuel_need))
        cbat_in = round(flight.config["cbat"], 2)
        print('[Pilot] cbat [W]: ' + str(cbat_in).ljust(10) + '\t' + str(cbat_need).ljust(10) + '\t' + str(round(100*(cbat_in - cbat_need) / (cbat_in + 0.001), 2)).ljust(10) + '\t-> ' + str(cbat_in >= cbat_need))
        pbat_in = round(flight.config["pbat"], 2)
        print('[Pilot] pbat [W]: ' + str(pbat_in).ljust(10) + '\t' + str(pbat_need).ljust(10) + '\t' + str(round(100*(pbat_in - pbat_need) / (pbat_in + 0.001), 2)).ljust(10) + '\t-> ' + str(pbat_in >= pbat_need))


    '''
    simuliert den Startlauf des FLugzeuges und speichert die Ergebnisse in flight
    Rückgabe: instanz von Flight
    '''
    def takeOff(self, flight):
        m = flight.getLastMass()
        h = self.m.getHeight(0)
        s = self.m.takeOffDistance
        Pto, t_to, v_to, F_0 = self.p.calcTakeOffParameter(m, h, s)
        
        SoC_0 = float(self.p.propulsion.bat.SoC)
        SoF_0 = float(self.p.propulsion.SoF)

        Nprop = self.p.propulsion.ice.calcNProp(Pto)
        Pshaft, Nprop, m_f, We, Pice, Nice, Pemot, Wcharge = self.p.propulsion.calcMotorConsumption(h, Pto, Nprop, t_to)

        flight.addStep((-1 * t_to), m, h, 0, np.nan, F_0, Pshaft, Nprop, Pice, Nice, Pemot, SoC_0, SoF_0, We, Wcharge)
        flight.addStep(0, m - m_f, h, v_to, np.nan, F_0, Pshaft, Nprop, Pice, Nice, Pemot, float(self.p.propulsion.bat.SoC), float(self.p.propulsion.SoF), 0., Wcharge)
        return flight


    '''
    wenn für einen Zeitschritt keine konstante FLuggeschwindigkeit vorgegeben ist, sondern ein kürzel genutzt wird, berechnet diese Funktion die entsprechende Geschwindigkeit
    Rückgabe: TAS
    '''
    def findSpeedSetting(self, m, H, w_g, R, txt):
        speedSettings = {'vEmax' : self.p.calcV_W_min,
                         'vdHmin' : self.p.calcV_dH_min,
                         'vFmin' : self.p.calcV_Fmin,
                         }
        
        return speedSettings[txt](m, H, w_g, R)

    '''
    zeigt Plots an und/oder speichert diese
    '''
    def plotFlight(self, flight):
        pPlot.plotFlight(flight, self.m, self.p, self.propulsionType, saveFile=self.savePlot, showPlot=self.verbose)


### MAIN #####################################

if __name__ == "__main__":
    mission = Mission('mission02_01')
    
    pic1 = Pilot(mission, propulsionType="combustion")
    flight1 = pic1.fly()

    pic2 = Pilot(mission, propulsionType="serialHybrid")
    flight2 = pic2.fly()
    
    conf = {'pice': 10000.0, 'mpl': 110.0, 'pmax': 90000.,
            'fuel': 15.5, 'cbat': 26400., 'pbat': 9500.}
    pic3 = Pilot(mission, propulsionType="parallelHybrid", config=conf)
    flight3 = pic3.fly()