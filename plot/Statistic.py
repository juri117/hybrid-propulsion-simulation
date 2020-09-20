'''
Created on 07.12.2015

@author: Juri Bieler
'''

import const as c
import generalUtils

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter

class Statistic(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        
    def plotFevBand(self):
        fig, ax = plt.subplots(figsize=(10,6))
        ax.grid(True)
        
        SFCfakt = 0.01
        im = plt.imread('../data/ICE/FEV-Streuband_merc_bg.png')
        implot = plt.imshow(im, extent=[0, 6, 200 * SFCfakt, 300 * SFCfakt])
    
        #------- efficiency -------
        FEV_FILE = 'data/ICE/FEV-Efficiency.data'
        fev = generalUtils.read1dAmeFile(FEV_FILE)
        fev_V_h = fev[:,0]
        fev_V_h = np.multiply(fev_V_h, 1000.)
        fev_sfc = fev[:,1]
        
        f_fev = inter.interp1d(fev_V_h, fev_sfc, kind='cubic')
        
        fev_sfc = np.divide(fev_sfc, 1./SFCfakt)
        #f_fev = inter.interp1d(fev_V_h, fev_sfc, kind='cubic')
        plt.plot(fev_V_h, fev_sfc, linewidth=3)
        
        #rotax914
        V_h_R914 = 1.2112
#         plt.plot([V_h_R914], [239 * SFCfakt], "ro", markersize=10)
#         ax.annotate(r'Rotax 914', xy=(V_h_R914 + 0, (239 + 0) * SFCfakt), xycoords="data",
#               xytext=(V_h_R914 + 1, (239 + 25) * SFCfakt), textcoords="data",
#               va="top", ha="center", fontsize=16,
#               bbox=dict(boxstyle="round", fc="w"),
#               arrowprops=dict(facecolor='black', shrink=0.05))
        
        #plt.gca().invert_yaxis()
        plt.xlabel("Hubraum $h_v$ in l")
        plt.ylabel("min. spez. Brennstoffverbrauch $sfc$ in g/kWh")
        
        ax.axes.set_ylim([0 * SFCfakt, 450 * SFCfakt])
        
        plt.draw()
        
#         labels = [item.get_text() for item in ax.get_yticklabels()]
#         newLabel = []
#         for i in range(0, len(labels), 1):
#             if(labels[i] != ''):
#                 lab = float(labels[i])
#                 lab = (lab / SFCfakt)/1000
#                 newLabel.append(str(lab))
#             else:
#                 newLabel.append(labels[i])
    
        #'0', '50', '100', '150', 
        newLabel = ['0', '50', '100', '150', '200', '250', '300', '350', '400', '450']
        ax.set_yticklabels(newLabel)
        
            
        hub = [0.29, 0.57, 0.86, 0.95, 1.211]
        #hub = np.multiply(hub, 10)
        eta = [0.250, 0.291, 0.323, 0.331, 0.346]
        sfc = np.multiply(eta, c.ED_FUEL)
        sfc = np.divide(1000., sfc)
        sfc = np.divide(sfc, 1./SFCfakt)
        plt.plot(hub[0], sfc[0], "o", label="Variante 1", markersize=15)
        plt.plot(hub[1], sfc[1], "o", label="Variante 2", markersize=15)
        plt.plot(hub[2], sfc[2], "o", label="Variante 3", markersize=15)
        plt.plot(hub[3], sfc[3], "o", label="Variante 4", markersize=15)
        plt.plot(hub[4], sfc[4], "o", label="Rotax 914", markersize=15)
        # 1 / (float(FFmin * eff) * c.ED_FUEL / 1000)
        
        h = [0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2.]
        
        #plt.plot(h, np.divide(f_fev(h), 1./SFCfakt), "o")
        #print(str(f_fev(0.29)))
        #plt.xlim([0, 2])
        #ax.axes.set_ylim([200 * SFCfakt, 400 * SFCfakt])

        generalUtils.pltCosmetics()
        plt.legend(loc=4)
        
        plt.show()
        
        
        
    def plotEmotor(self):
        MOTOR_FILE = 'data/E-Motor/E-Motors.csv'
        add = c.PD_emot_add
        dens = c.PD_emot
        self.plotWeights(MOTOR_FILE, add, dens, skipLast=2)
        
    def plotICE(self):
        MOTOR_FILE = 'data/ICE/ICE.csv'
        add = 1.
        dens = c.PD_ice
        self.plotWeights(MOTOR_FILE, add, dens)
    
    def plotWeights(self, MOTOR_FILE, add, dens, skipLast=0):
        csv = generalUtils.readCsv(MOTOR_FILE)
        P = csv[:,1]
        m = csv[:,2]
        name = generalUtils.readStrForNaN(MOTOR_FILE, 0, csv[:,0])
        Plin = 1.1 * max(P)
        mlin = add + (1.1 * max(P) / (dens / 1000))
        plt.plot([0, Plin], [add, mlin], linewidth=c.LINE_WIDTH)
        
        if(skipLast > 0):
            Pav = np.sum(P[:-skipLast]) / (len(P)-skipLast)
            mav = (np.sum(m[:-skipLast]) / (len(m)-skipLast)) - add
        else:
            Pav = np.sum(P) / len(P)
            mav = (np.sum(m) / len(m)) - add
        
        print("Ration [kW/gk]: " + str(1000 * Pav / mav))
        
        #PlinAv = 1.1 * max(P)
        #mlinAv = 0 + (1.1 * max(P) / (Pav/mav))
        #plt.plot([0, PlinAv], [add, mlinAv])
        
        plt.plot(P, m, 'o', markersize=10)
        plt.xlim([0, 1.1 * max(P)])
        plt.ylim([0, 1.1 * max(m)])
        plt.xlabel("Leistung $P$ in kW")
        plt.ylabel("Masse $m$ in kg")
        for i in range(0, len(name), 1):
            xAdd = 10
            yAdd = -6
            if(i % 2 == 1):
                xAdd = -10
                yAdd = 6
            if(len(name[i]) > 0):
                if(name[i][0] != '-'):
                    plt.annotate(name[i], xy=(P[i], m[i]), xycoords="data",
                          xytext=(P[i] + xAdd, m[i] + yAdd), textcoords="data",
                          va="top", ha="center", fontsize=14,
                          bbox=dict(boxstyle="round", fc="w"),
                          arrowprops=dict(arrowstyle="->", linewidth=2))

        generalUtils.pltCosmetics()
        plt.draw()
        plt.show()

if __name__ == "__main__":
    stat = Statistic()
    #stat.plotEmotor()
    #stat.plotICE()
    stat.plotFevBand()