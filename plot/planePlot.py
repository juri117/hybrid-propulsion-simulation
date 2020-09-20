'''
Created on 05.12.2015

@author: Juri Bieler
'''

import const as c
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import numpy as np
import scipy.interpolate as inter
import scipy.optimize
from matplotlib import rc

def plotThrustOverflow(plane):
    m = 900
    h = 2000
    u_g = np.array(range(15, 75, 1))
    w_g = 3.5
    R = 0
    ones = np.ones(len(u_g))
    Fneed = list(map(plane.calcNeededThrust, ones * m, ones * h, u_g, ones * w_g, ones * R))
    
    P = plane.P_max
    Favailable = []
    Diff = []
    Nprop = plane.propulsion.calcNProp(P)
    for v in u_g:
        eta, Nprop = plane.propulsion.prop.calcEfficiency(h, v, Nprop, P)
        F = eta * P / v
        Favailable.append(F)
        Diff.append(F - plane.calcNeededThrust(m, h, v, w_g, R))
        #print(str(eta) + " @ " + str(Nprop))
    
    Diff = (np.array(Diff).T)[0]
    
    G = m * c.g
    
    ax = host_subplot(111)
    ax2 = ax.twinx()
    
    ax.plot(u_g, Fneed, "b", label="benötigt", linewidth=c.LINE_WIDTH)
    ax.plot(u_g, Favailable, "g", label="vorhanden", linewidth=c.LINE_WIDTH)
    
    f_Diff = inter.interp1d(u_g, -1 * Diff, kind='cubic')
    u_gamma_max = scipy.optimize.fmin(f_Diff, 35)
    ax.plot(u_g, Diff, "r", label="Überschuss", linewidth=c.LINE_WIDTH)
    
    ax2.plot(u_g, np.divide(Fneed, G), "b--")
    ax2.plot(u_g, np.divide(Favailable, G), "g--")
    ax2.plot(u_g, np.divide(Diff, G), "r--")
    
    #u_gamma_max = slef.calcV_Fmin
    
    #plt.plot(u_g, -1 * f_Diff(u_g), 'x', label="")
    ax.plot([u_gamma_max, u_gamma_max], [min(Diff), max(max(Fneed), max(Favailable))], 'k--', linewidth=c.LINE_WIDTH)
    ax.plot([min(u_g), max(u_g)], [0, 0], 'k-')
    ax.set_xlim([min(u_g), max(u_g)])
    ymax = max(max(Fneed), max(Favailable), max(Diff))
    ymin = min(min(Fneed), min(Favailable), min(Diff))
    ax.set_ylim([ymin, ymax])
    ax2.set_ylim([ymin / G, ymax / G])
    ax.set_ylabel("Schub $F$ in N")
    ax.set_xlabel("Horizontalfluggeschwindigkeit $u_g$ in m/s")
    ax2.set_ylabel("$F/G$")
    
    plt.legend()
    font = {'family' : 'normal', 'size'   : 18}
    rc('font', **font)
    plt.subplots_adjust(left=0.17, bottom=0.13, right=0.85, top=0.95, wspace=None, hspace=None)
    plt.show()