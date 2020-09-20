'''
Created on 02.12.2015

@author: Juri Bieler
'''

import utils.generalUtils as generalUtils
import scipy.interpolate as inter
import pylab as plt
import matplotlib.cm as cm
import numpy as np
import std_atm as SA
import const as c

class Propeller(object):
    '''
    Liest Propellerkennfeld ein und erstellt lineare Interpolationsfunktion
    '''
    def __init__(self):
        self.verbose = False
        self.efficiencyMap = c.PROP_MAP
        self.Cp, self.J, self.eta = generalUtils.read2dAmeFile(self.efficiencyMap)
        self.f_eta = inter.interp2d(self.Cp, self.J, self.eta, kind='linear')
        #self.Pmax = 84500   #W

    '''
    Rückgabe: Fortschrittsgrad "J"
    '''
    def calcJ(self, tas, N):
        if(N <= 0):
            return float('nan')
        J = tas * 60 / (N * c.PROP_D)
        return J
    
    '''
    Rückgabe: Leistungsbeiwert "C_P"
    '''
    def calcCp(self, h, N, Pshaft):
        if(N <= 0):
            return float('nan')
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        Cp = 100 * Pshaft / (rho * (N / 60)**3 * c.PROP_D**5)
        return Cp
        
    '''
    berechnet den Propellerwirkungsgrad abhängig von Flughöhe "h", Fluggeschwindigkeit "tas", Propellerdrehzahl "N" und
    Triebwerksleistung "Pshaft"
    Liegt J oder C_P auserhalb des Kennfeldes, wird die Motordrehzahl erhöhr, bis beide Werte innerhalb liegen
    Rückgabe: Wirkungsgrad "eta", gegebenenfalls neue Propellerdrehzahl "N"
    '''
    def calcEfficiency(self, h, tas, N, Pshaft, NpropMax=2000, verbose=False, fixN=False, counter=0):
        #fängt null Leistung ab
        if(Pshaft == 0):
            return 1., N

        J = self.calcJ(tas, N)
        Cp = self.calcCp(h, N, Pshaft)
        eta = self.getEfficiency(Cp, J)

        if(fixN and (N > 0. and N < NpropMax) and counter < 100):
            if(Cp < 2 or J < 0.2):
                N = N - 100
                #print('[Propeller] FIX -> N+100: ' + str(N) + ' -> J: ' + str(J) + " Cp: " + str(Cp) + " eta: " + str(eta))
                return self.calcEfficiency(h, tas, N, Pshaft, fixN=True, counter = counter + 1)
            if(Cp > 22 or J > 2.2):
                N = N + 100
                #print('[Propeller] FIX -> N+100: ' + str(N) + ' -> J: ' + str(J) + " Cp: " + str(Cp) + " eta: " + str(eta))
                return self.calcEfficiency(h, tas, N, Pshaft, fixN=True, counter = counter + 1)
        
        if(counter >= 100):
            print("[Propeller] WARNING: Timeout in calcEfficiency(...)")        

        if(verbose):    
            if(Cp < 2. / 1.1 or Cp > 22 * 1.1):
                print("[Propeller] Warning: Cp is out of Range [2-22]! Cp = " + str(Cp) + " -> f(" + str(h) + "," + str(tas) + "," + str(N) + "," + str(Pshaft) + ")")
            if(J < 0.2 / 1.1 or J > 2.2 * 1.1):
                print("[Propeller] Warning: J is out of Range [0.2-2.2]! J = " + str(J) + " -> f(" + str(h) + "," + str(tas) + "," + str(N) + "," + str(Pshaft) + ")")
        #eta = 0.7
        return eta, N

    '''
    liest interpoliertes Propellerkennfeld aus
    Rückgabe: Wirkungsgrad "eta"
    '''
    def getEfficiency(self, Cp, J):
        return self.f_eta(Cp, J)

    '''
    Rückgabe: Propellerkreisfläche
    '''
    def calcA_prop(self):
        return (np.pi / 4) * c.PROP_D**2
    
    def plotEfficiencyMap(self, flights=None, fileName=None, showPlot=True):
        J = [x / 100.0 for x in range(20, 221, 1)]
        Cp = [x / 10.0 for x in range(20, 221, 1)]
        
        eta = np.zeros([len(J), len(Cp)])

        for iCp in range(0, len(Cp), 1):
            for iJ in range(0, len(J), 1):
                eta[iJ, iCp] = self.getEfficiency(Cp[iCp], J[iJ])
        
        #fig = plt.figure()
        fig, ax = plt.subplots(figsize=(11,7))
        fig.canvas.set_window_title('propellerMap')
        #ax = fig.add_subplot(111)
        map = plt.pcolormesh(np.array(Cp), np.array(J), np.array(eta), vmin=0., vmax=1., cmap=generalUtils.custom_div_cmap(numcolors=256, mincol='#FF0000', maxcol='#00CC00', midcol='#FFFF00'))
        ax.set_xlabel("Leistungskoeffizient $C_{P}$ $10^{-2}$")
        #ax.axes.set_xlim([2, 22])
        ax.set_ylabel("Fortschrittsgrad $J$")
        #ax.axes.set_ylim([0.2, 2.2])
        cb = plt.colorbar()
        cb.set_label("Propellerwrikungsgrad $\eta_{Prop}$")
        plt.gca().invert_yaxis()

        if(flights != None):
            CpF_vec = []
            JF_vec = []
            #t = []
            colors = iter(cm.rainbow(np.linspace(0, 1, len(flights.t))))
            for i in range(0, len(flights.t), 1):
                if(flights.Pshaft[i] > 0 and flights.tas[i] > 0):
                    CpF = self.calcCp(flights.H[i], flights.Nprop[i], flights.Pshaft[i])
                    JF = self.calcJ(flights.tas[i], flights.Nprop[i])
                    #t.append(float(flights.t[i]) / flights.t[-1])
                    plt.scatter(CpF, JF, c=next(colors))
                    CpF_vec.append(CpF)
                    JF_vec.append(JF)
                    
            ax.set_xlim([min(min(Cp), min(CpF_vec)), max(max(Cp), max(CpF_vec))])
            ax.set_ylim([min(min(J), min(JF_vec)), max(max(J), max(JF_vec))])
        else:
            ax.set_xlim([min(Cp), max(Cp)])
            ax.set_ylim([min(J), max(J)])
        
        generalUtils.pltCosmetics()
        
        if(fileName != None):
            fig.savefig(fileName + ".png", format='png', dpi=400)
        if(showPlot):
            plt.show()
        plt.close('all')
        
if __name__ == "__main__":
    pr = Propeller()
    pr.plotEfficiencyMap(fileName="test.png")
    #eta = pr.getEfficiency(4, 1)
    
    #pr.calcEfficiency(1000, 30, 1500, 3000, 20000)
    #print(str(eta))
    
    eta = pr.calcEfficiency(1809.52380952,74,679.425675676,2042.45918826)
    print(str(eta))