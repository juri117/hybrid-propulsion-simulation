'''
Created on 02.12.2015

@author: Juri Bieler
'''

from units import * 
import generalUtils as generalUtils
import scipy.interpolate as inter
import scipy.optimize
import icePlot
#import pylab as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import colors, ticker, cm
import math
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import const as c
from Propeller import Propeller
import matplotlib.cbook as cbook
from scipy import interpolate

class ICE(object):
    '''
    liest Motorkennfeld aus und die max. Mitteldruckkurve
    berechnet benötigtes Hubvolumen für Motorleistung "Pmax"
    skaliert das Kennfeld entsprechend des FEV-Streubandes (kleine Motoren sind schlechter)
    interpoliert das Motorkennfeld
    '''
    def __init__(self, Popt=0, Pmax=0, verbose=True):
        self.verbose = verbose

        self.i = c.STROKES_ICE

        #------- sfcMap-max-curve -------
        self.maxPeMap = c.P_E_MAX_MAP
        self.maxPe = generalUtils.read1dAmeFile(self.maxPeMap)
        self.Nmax = self.maxPe[:,1]
        self.p_emax = self.maxPe[:,0]
        self.f_p_emax = inter.interp1d(self.Nmax, self.p_emax, kind='linear')
        self.N_maxOpPoint, self.p_e_maxOpPoint = self.calcMaxOpPoint()


        #------- efficiency -------
        FEV_FILE = c.FEV_MAP
        fev = generalUtils.read1dAmeFile(FEV_FILE)
        fev_V_h = fev[:,0]
        fev_sfc = fev[:,1]
        f_fev = inter.interp1d(fev_V_h, fev_sfc, kind='cubic')

        #------- sfcMap --------
        #PeFF [pa]
        #NFF [min**-1]
        #FF [g/kWh]
        fuelFlowMap = c.SFC_MAP
        PeFF, NFF, FF = generalUtils.read2dAmeFile(fuelFlowMap)
        NFF = np.array(NFF)
        PeFF = PeFF * 100000  #convert bar to pascal
        PeFF = np.array(PeFF)
        self.FF = np.array(FF)
        self.N_opt, self.p_e_opt, FFmin = self.findOptimum(PeFF, NFF, self.FF)
        
        #------- displacement -------
        self.V_h = 0.
        if(Popt > 0):
            self.V_h = self.calcV_h_4opt(Popt)   #m**3
        if(Pmax > 0):
            self.V_h = self.calcV_h_4max(Pmax)   #m**3            

        #------- efficiency correction -------
        eff = f_fev(self.V_h) / FFmin

        if(not c.FIX_ICE_EFF):
            eff = 1.
        self.eff_total = 1 / (float(FFmin * eff) * c.ED_FUEL / 1000)

        self.FF = np.multiply(self.FF, eff)
        self.f_FF = inter.interp2d(NFF, PeFF, self.FF, kind='linear')


        #------- optimal workingLine -------
        self.N_Pmap = c.OPT_WORKINGLINE_ICE
        self.N_P = generalUtils.read1dAmeFile(self.N_Pmap)    #P in %
        self.f_N_P = inter.interp1d(self.N_P[:,0], self.N_P[:,1], kind='linear')

        self.Popt = self.calcPopt()
        self.Pmax = self.calcPmax()

        self.gearRatio = c.GEAR_RATIO
        
        if(self.Pmax < c.P_MAX_ICE_MIN):
            self.Pmax = 0.
            self.Popt = 0.


    '''
    findet arbeitspunkt mit min sfc
    Rückgabe: Drehzahl, Mitteldruck, sfc (an dieser Stelle)
    '''
    def findOptimum(self, PeFF, NFF, FF):
        Pemin = 0
        Nmin = 0
        FFmin = 999
        
        for Pei in range(0,len(PeFF), 1):
            for Ni in range(0,len(NFF), 1):
                if(FFmin > FF[Pei][Ni]):
                    FFmin = FF[Pei][Ni]
                    Nmin = NFF[Ni]
                    Pemin = PeFF[Pei]
        return Nmin, Pemin, FFmin
        
    '''
    Rückgabe: Arbeit des Generators in der Zeit "dt" in sec.
    '''
    def calcGeneratorWe(self, Pice, dt):
        We = (Pice * (c.GENERATOR_EFF * c.INVERTER_EFF)) * (float(dt)/60**2)
        return We
    
    '''
    berechnet nätige Verbrennungsmotorleistung um vorgegebene elek. Leistung von Generator "We" zu erhalten
    Rückgabe: Verbrennungsmotorleistung
    '''
    def calcNeededICEPower(self, We, dt):
        WeLoss = We / (c.GENERATOR_EFF * c.INVERTER_EFF)
        Pice = WeLoss / (float(dt)/60**2)
        return Pice

    '''
    berechnet passende Verbrennungsmotordrehzahl für geforderte Leistung "P" entlang der Arbeitslinie
    bei max Verbrennungsmotorleistung "Pmax"
    Rückgabe: Dreahzahl
    '''
    def getNforP(self, P, Pmax, verbose=False):
        if(Pmax == 0.):
            #passiert nur, wenn der Verbrennungsmorot nicht vorhanden ist und spielt dann keine Rolle
            Pmax = 84500.
        P_pro = (P / Pmax) * max(self.N_P[:,0])
        if(P_pro < min(self.N_P[:,0])):
            P_pro = min(self.N_P[:,0])
            if(verbose):
                print('[ICE] minimal Power is 0 Watt, requestet: ' + str(P))
        if(P_pro > max(self.N_P[:,0])):
            P_pro = max(self.N_P[:,0])
            if(verbose):
                print('[ICE] maximal Power is ' + str(Pmax) + ' Watt, requestet: ' + str(P))
        return self.f_N_P(P_pro)
    
    '''
    berechnet passende Verbrennungsmotordrehzahl für geforderte Leistung "P" entlang der Arbeitslinie
    Rückgabe: Dreahzahl
    '''
    def calcNIce(self, Pice):
        Nprop = self.getNforP(Pice, self.Pmax)
        return Nprop

    '''
    berechnet passende Propellerdrehzahl für geforderte Leistung "P" entlang der Arbeitslinie
    Rückgabe: Dreahzahl
    '''
    def calcNProp(self, Pshaft):
        return self.NiceToNprop(self.calcNIce(Pshaft))
    
    '''
    rechnet Getriebeübersetzung um
    Rückgabe Verbrennungsmotordrehzahl
    '''
    def NpropToNice(self, Nprop):
        return Nprop * self.gearRatio
    
    '''
    rechnet Getriebeübersetzung um
    Rückgabe Propellerdrehzahl
    '''
    def NiceToNprop(self, Nice):
        return Nice / self.gearRatio

    '''
    Berechnet Brennstoffbedarf für Zeit "dt" in sec.
    Rückgabe: Brensoffmasse in kg, sfc in g/kWh
    '''
    def calcConsumption(self, N, M, dt):
        FF = self.calcFuelFlow(N, M)
        W = N2omega(N) * M * (float(dt)/60**2)
        m_f = (FF * (W/1000)) / 1000
        return m_f, FF
    
    '''
    rechnet Brennsoffmasse in Liter pro Stunde um
    '''
    def calcConsumtionToL_h(self, m_f, dt):
        return self.calcConsumtionToL(m_f, dt) * (60**2) / dt
    
    '''
    rechnet Brennsoffmasse in Liter um
    '''
    def calcConsumtionToL(self, m_f, dt):
        return m_f / c.DENS_FUEL

    '''
    berechnet das Hubvolumen, wenn die Leistung im Punkt mit min. sfc des Motors "P_opt" betragen soll
    Rückgabe: Hubvolumen in l
    '''
    def calcV_h_4opt(self, P_opt):
        M = P_opt / N2omega(self.N_opt)
        V_h = M * 2 * math.pi * self.i / self.p_e_opt
        return V_h
    
    '''
    berechnet das Hubvolumen, wenn die Maximalleistung des Motors "P_max" betragen soll
    Rückgabe: Hubvolumen in l
    '''
    def calcV_h_4max(self, P_max):
        N, p_e = self.getMaxOpPoint()
        M = P_max / N2omega(N)
        V_h = M * 2 * math.pi * self.i / p_e
        return V_h
    
    '''
    berechnet die Leistung im Punkt mit min. sfc des Motors
    '''
    def calcPopt(self):
        Popt = N2omega(self.N_opt) * self.p_e_to_M(self.p_e_opt)
        return Popt
        
    '''
    berechnet Maximalleistung des Motors
    '''
    def calcPmax(self):
        N, p_e = self.getMaxOpPoint()
        Pmax = N2omega(N) * self.p_e_to_M(p_e)
        return Pmax

    '''
    Rückgabe: Drehzahl und Mitteldruck des Betriebspunkt für Maximalleistung zurück
    '''
    def getMaxOpPoint(self):
        return self.N_maxOpPoint, self.p_e_maxOpPoint
    
    '''
    findet Drehzahl und Mitteldruck des Betriebspunkt mit max. Motorleistung
    Rückgabe: Drehzahl, Mitteldurck
    '''
    def calcMaxOpPoint(self):
        N = range(int(min(self.Nmax)),int(max(self.Nmax) + 1), 1)
        p_e = self.f_p_emax(N)
        Prod = list(np.multiply(N,p_e))
        imax = Prod.index([max(Prod)])
        p_e_max = p_e[imax]
        N_max = N[imax]
        return N_max, p_e_max
    
    '''
    Rückgabe: max. Mitteldruck der mit der Drehzahl "N" möglich ist (laut Kennfeld)
    '''
    def getMaxP_e(self, N):
        return self.f_p_emax(N)
    
    '''
    Rückgabe: sfc mit Drehzahl "N" und Moment "M"
    '''
    def calcFuelFlow(self, N, M):
        p_e = self.M_to_p_e(M)
        return self.getFuelFlow(N, p_e)
    
    '''
    rechnet Moment "M" in Mitteldruck "p_e" um
    Rückgabe: Mitteldruck
    '''
    def M_to_p_e(self, M):
        if(self.V_h > 0.):
            return M * 2 * math.pi * self.i / self.V_h
        return 0
    
    '''
    rechnet Mitteldruck "p_e" in Moment "M" um
    Rückgabe: Moment
    '''
    def p_e_to_M(self, p_e):
        M = p_e * self.V_h / (2 * math.pi * self.i)
        return M
    
    '''
    Rückgabe: interpoliertem fc aus Kennfeld
    '''
    def getFuelFlow(self, N, p_e):
        return self.f_FF(N, p_e)
    
    '''
    Rückgabe: Leerlaufdrehzahl des Motors
    '''
    def getNidle(self):
        return min(self.Nmax)
    
    '''
    Rückgabe: max. Drehzahll des Motors
    '''
    def getNmax(self):
        return max(self.Nmax)
    
    def plotEfficiencyMap(self, flight=None, fileName=None, showPlot=True, propType=c.COMB):
        icePlot.plotEfficiencyMap(self, flight=flight, fileName=fileName, showPlot=showPlot, propType=propType)
        
    def plotEfficiencyMap3D(self):
        icePlot.plotEfficiencyMap3D(self)
        
    def plotPmax(self):
        icePlot.plotPmax(self)
        
    def plotFuelConsumption(self):
        icePlot.plotFuelConsumption(self)


if __name__ == "__main__":
    ice = ICE(Pmax = 84500)

    ice.plotEfficiencyMap(showPlot=True, fileName="test")
