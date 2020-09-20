'''
Created on 07.12.2015

@author: Juri Bieler
'''

import const as c
from units import *
from generalUtils import *

import numpy as np

class Emotor(object):
    ### CONSTRUCKTOR #####################################
    '''
    setzt "Pmax"
    '''
    def __init__(self, Pmax = 84000):
        self.Pmax = Pmax
    
    '''
    berechnet die Arbeit bei Drehzahl "N", bei Moment "M" für die Zeit "dt" in sec.
    Rückgabe: Arbeit
    '''
    def calcConsumption(self, N, M, dt):
        P = N2omega(N) * M
        We = (P / (c.EMOTOR_EFF * c.INVERTER_EFF)) * (float(dt)/60**2)
        return We

    def plotWeights(self):
        MOTOR_FILE = 'data/E-Motor/E-Motors.csv'
        csv = readCsv(MOTOR_FILE)
        name = csv[:,0]
        P = csv[:,1]
        m = csv[:,2]
        name = readStrForNaN(MOTOR_FILE, 0, csv[:,0])
        Plin = 1.1 * max(P)
        mlin = c.PD_emot_add + (1.1 * max(P) / (c.PD_emot / 1000))
        plt.plot([0, Plin], [c.PD_emot_add, mlin])
        
        Pav = np.sum(P[:-2]) / (len(P)-2)
        mav = (np.sum(m[:-2]) / (len(m)-2)) - c.PD_emot_add
        
        print("Ration [kW/gk]: " + str(1000 * Pav / mav))
        
        plt.plot(P, m, 'o')
        plt.xlim([0, 1.1 * max(P)])
        plt.ylim([0, 1.1 * max(m)])
        plt.xlabel("$P [kW]$")
        plt.ylabel("$m [kg]$")
        for i in range(0, len(name), 1):
            xAdd = 10
            yAdd = -6
            if(i % 2 == 1):
                xAdd = -10
                yAdd = 6
            plt.annotate(name[i], xy=(P[i], m[i]), xycoords="data",
                  xytext=(P[i] + xAdd, m[i] + yAdd), textcoords="data",
                  va="top", ha="center", fontsize=14,
                  bbox=dict(boxstyle="round", fc="w"),
                  arrowprops=dict(arrowstyle="->"))
            
        
        pltCosmetics()
        plt.draw()
        plt.show()

if __name__ == "__main__":
    emot = Emotor(Pmax = 84000)
    print('Cnsumption= ' + str(emot.calcConsumption(4000, 14, 60**2) / 1000) + ' kWh')
    emot.plotWeights()