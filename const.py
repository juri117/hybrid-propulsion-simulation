'''
Created on 13.12.2015

@author: Juri Bieler
'''

### Ausgabe, Dateien ############################################

''''Strings sind nur fürs Speichern der Dateien wichtig, bitte einfach so lassen'''
COMB = "combustion"
SERI = "serialHybrid"
PARA = "parallelHybrid"

'''Trennzeichen für CSV-Dateien'''
CSV_DELIMITER = ";"

'''Linienbreite für Plots in Pixel'''
LINE_WIDTH = 3


### ITERATION LIMITS ############################################

'''Iterationslimit für alle Parameter der Optimierung (pmax, fuel, cbat, pbat) in % [0 - 1]'''
ITERATION_LIMIT = 0.1 / 100.

'''Toleranz in der das maximale Abfluggewicht MTOM (siehe weiter unten) eingehalten werden muss in % [0 - 1]'''
MTOM_TOLERANCE = 0.05

'''Toleranz in der das maximale Tankvolumen V_TANK_MAX (siehe weiter unten) eingehalten werden muss in % [0 - 1]'''
FUELTANK_TOLERANCE = .1

'''definiert um wie genau verschiedene Verbrennungsmotorgräßen von Hybridantrieben untersucht werden sollen, Schrittweite in Watt
Anwendung in: Optimizer.optimizeHybridFlight'''
P_INCREMENT = 3000.

'''definiert das Limit für die Feinuntersuchtung der Minima des Energieverbrauches in Watt
Anwendung in: Optimizer.optimizeHybridFlight_rec'''
P_INCREMENT_LIMIT = 1.

'''kleinster Verbrennnungsmotor der verbraut werden kann in W'''
P_MAX_ICE_MIN = 1000.



### AERODYNAMICS ############################################

'''Flügelfläche in m^2'''
S = 17.42

'''max. Gleitzahl'''
E_MAX = 29.

'''Geschwindigkeit für E_MAX in m/s'''
V_E_MAX = 32.

'''min. Fluggeschwindigkeit (grobe schätzung reicht)'''
V_MIN = 18.

'''Oswald Faktor'''
OSWALD_FAKT = 0.8

'''Flügelstreckung'''
STRETCH = 18.62

#stallspeed
#V_MIN = 10.

'''Startparameter'''
C_A_TAKEOFF = 1.45
C_W_TAKEOFF = 0.06
'''Rollreibungsbeiwert'''
WHEEL_MY = 0.045 #Dubbel Q 87

### MASS ############################################

'''Erdbeschleunigung in m s^-2'''
g = 9.80665

'''max. Ablfugmasse in kg'''
MTOM = 900

'''Leermasse des Flugzeuges'''
M_EMPTY = 660

'''Masse des Originalmotors des FLugzeuges (mit Kühlung und LiMa)'''
M_ORIG_ENGINE = 78.

'''Strukturmasse des FLugzeuges ohne Motor'''
M_STRUCT = M_EMPTY - M_ORIG_ENGINE

'''standart Nutzlastmasse kann in jeder config.ini überschrieben werden'''
M_PL = 110.


### Leistungsdichten ############################################

'''Leistungsdichte ICE in kW/kg'''
PD_ice = 1070.
'''Anschlussmasse ICE in kg'''
PD_ice_add = 1.

'''Leistungsdichte E-Motor in kW/kg'''
PD_emot = 3760.
'''Anschlussmasse E-Motor in kg'''
PD_emot_add = 2.0

'''Leistungsdichte Umrichter in kW/kg'''
PD_inv = 15000.
'''Anschlussmasse Umrichter in kg'''
PD_inv_add = 0.4


### Brennstoff und Strom ############################################

'''Dichte Avgas (15C) kg/l'''
DENS_FUEL = 0.75

''''max Tankvolumen in l'''
V_TANK_MAX = 500

'''Strompreis in euro/kWh'''
PRICE_ELEC = 0.30 #konservativ -> GreenpeaceEnergy, eig. ~28cent
#https://www.greenpeace-energy.de/oekostrom/der-stromtarif.html

'''Avgaspreis in euro/l'''
PRICE_FUEL = 2.

''''Energiedichte von Avgas in kWh/kg'''
ED_FUEL = 12.1
''''Energiedichte von Avgas in kWh/l'''
ENERGY_FUEL_L = ED_FUEL * DENS_FUEL
#http://www.benzinpreis-aktuell.de/


### Propulsion ############################################

'''kleine Verbrennungsmotoren sind ineffizient (FEV-Streuband) (True), nicht beachten (False)'''
FIX_ICE_EFF = True

'''Verbrennungsmotorkennfeld, sfc über Drehzahl und Mitteldruck'''
SFC_MAP = 'data/ICE/sfc_Rotax914.data'

'''max p_e Verlauf über Drehzahl von Verbrennungsmotor'''
P_E_MAX_MAP = 'data/ICE/peMax_Rotax914.data'

'''FEV-Streuband, verlauf des min. sfc von Verbrennungmotoren über Hubvolumen'''
FEV_MAP = 'data/ICE/FEV-Efficiency.data'

'''optimale Arbeitslinie des Verbrennungsmotors, Motorleistung in Prozent [0-1.15] über Drehzahl'''
OPT_WORKINGLINE_ICE = 'data/ICE/Leistung-RPM-Map_S6_Opti.data'

'''x-Takter, Umdrehungen pro Zündung'''
#4-Takter
STROKES_ICE = 2

'''Propellerkennfeld'''
PROP_MAP = 'data/propeller/MTV-7-A_170-51.data'

'''Getriebeübersetzung zwischen ICE und Propeller'''
GEAR_RATIO = 2.6708

#P_MAX = 84500.

'''maximale Triebwerksleistung die unterscuht werden soll in Watt'''
P_ICE_MAX = 150000.

'''Getriebewirkungsgrade für die drei Antriebskonzepte in % [0 - 1]'''
GEAR_EFF_COMB = 0.98
GEAR_EFF_SERI = 1.
GEAR_EFF_PARA = 0.98
GEAR_EFF = {COMB:GEAR_EFF_COMB, SERI:GEAR_EFF_SERI, PARA:GEAR_EFF_PARA}

#N_IDLE = 1200.
#N_MAX = 5800.

'''GeneratorWirkungsgrad'''
GENERATOR_EFF = 0.95

'''Umrichter Wirkungsgrad'''
INVERTER_EFF = 0.98

'''ElektromotorWirkungsgrad'''
EMOTOR_EFF = 0.95

'''Prpellerdurchmesser in m'''
PROP_D = 1.7

'''max. Propellerdrehzahl in 1/min'''
PROP_MAX_N = 2800.


### Battery ############################################

'''mac. Laderate in C (vielfache der Kappazität)'''
CHARGE_RATE = 2.

'''Batterie-Ladewirkungsgrad'''
ETA_CHARGE = 0.9
'''Batterie-Entladewirkungsgrad'''
ETA_DISCHARGE = 0.9

'''Batterie-Kennfeld (Leistungsdichte über Energiedichte)'''
BAT_FILE = 'data/Battery/ED_PD_LiPo.data'
#BAT_FILE = 'data/Battery/ED_PD_LiS.data'

