'''
Created on 26.12.2015

@author: Juri Bieler
'''

from Plane import Plane
from Mission import Mission
from Flight import Flight
import components.Propulsion as Propulsion
import const as c

import numpy as np
import utils.generalUtils as generalUtils
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import math
import time


def plotFlight(flight, mission, plane, PropulsionType, saveFile=False, showPlot=True):
    plt.figure("Missionsprofil", figsize=(14,9))

    timeVec = np.divide(flight.t, 60.)
    timeMin = min(timeVec)
    timeMax = max(timeVec)
    #timeMax = 10
    
    legendDist = 1.1
    
    host = host_subplot(311, axes_class=AA.Axes)
    host.ticklabel_format(useOffset=False)
    plt.subplots_adjust(right=0.75)
    
    par1 = host.twinx()
    par1.ticklabel_format(useOffset=False)

    host.set_xlabel("Zeit $t$ in min")
    host.set_ylabel("Flughöhe $H$ in m")
    par1.set_ylabel("$TAS$ in m/s")
    
    host.locator_params(nbins=6)
    par1.locator_params(nbins=6)
    

    p1, = host.plot(timeVec, flight.H, label="$h$", linewidth=c.LINE_WIDTH)
    p2, = par1.plot(timeVec, flight.tas, label="$TAS$", linewidth=c.LINE_WIDTH)
    
    host.axes.set_ylim([0, 1.05 * max(flight.H)])
    par1.axes.set_ylim([0, 1.05 * max(flight.tas)])
    
    if((~np.isnan(flight.r)).any()):
        par2 = host.twinx()
        new_fixed_axis = par2.get_grid_helper().new_fixed_axis
        par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(60, 0))
        par2.axis["right"].toggle(all=True)
        par2.set_ylabel("Kurvenradius $r$ in m")
    
        p3, = par2.plot(timeVec, flight.r, label="$R$", linewidth=c.LINE_WIDTH)
        par2.axes.set_ylim([0, 1.1 * np.nanmax(flight.r)])
        par2.axis["right"].label.set_color(p3.get_color())
        legendDist = 1.18

    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    host.set_xlim([timeMin, timeMax])
    

    for i in range(0, len(mission.sec), 2):
        if(mission.w_g[i] != 0.):
            txt = "{:.1f}".format(mission.w_g[i])
            x = float(mission.sec[i]) / 60.
            y = mission.height[i]
            if(len(mission.sec) > i + 1):
                x = float(mission.sec[i] + mission.sec[i+1]) / (2. * 60.)
                y = (mission.height[i] + mission.height[i+1]) / 2
                x = max(x, 0.05 * float(mission.sec[-1]) / 60.)
                x = min(x, 0.95 * float(mission.sec[-1]) / 60.)
            host.annotate(txt, xy=(x, y), xycoords="data", va="center", ha="center", bbox=dict(boxstyle="round", fc="w", alpha=0.5),
                fontsize=18)

    plt.legend(bbox_to_anchor=(legendDist, 1.1), loc=2, borderaxespad=1)
    generalUtils.pltCosmetics(multiPlot=True)

    #SECOND PLOT################

    host = host_subplot(312, axes_class=AA.Axes)
    host.ticklabel_format(useOffset=False)
    plt.subplots_adjust(right=0.75)
    
    host.set_xlabel("Zeit $t$ in min")
    host.set_xlim([timeMin, timeMax])
    host.set_ylabel("Leistung $P$ in kW")
    host.locator_params(nbins=6)
    host.axes.set_ylim([0, 1.05 * (plane.P_max)/1000])

    p1, = host.plot(timeVec, np.divide(flight.Pshaft, 1000.), label="$P_{Welle}$", linewidth=c.LINE_WIDTH)
    p1, = host.plot(timeVec, np.divide(flight.Pice, 1000.), label="$P_{ice}$", linewidth=c.LINE_WIDTH)
    p1, = host.plot(timeVec, np.divide(flight.Pemot, 1000.), label="$P_{E-Motor}$", linewidth=c.LINE_WIDTH)

    plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)
    generalUtils.pltCosmetics(multiPlot=True)


    #THIRD PLOT################

    host = host_subplot(313, axes_class=AA.Axes)
    host.ticklabel_format(useOffset=False)
    plt.subplots_adjust(right=0.75)

    par1.ticklabel_format(useOffset=False)
     
    host.set_xlabel("Zeit $t$ in min")
    host.set_xlim([timeMin, timeMax])
    host.set_ylabel("Energiereserve in %")
    host.locator_params(nbins=6)
    host.axes.set_ylim([0., 1.])

    p1 = host.plot(timeVec, flight.SoC, label="$SoC$", linewidth=c.LINE_WIDTH)
    p1 = host.plot(timeVec, flight.SoF, label="$SoF$", linewidth=c.LINE_WIDTH)

    plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)
    generalUtils.pltCosmetics(multiPlot=True)
    plt.draw()

    if(saveFile):
        fileName = mission.NAME + "_" + PropulsionType
        if(PropulsionType == c.SERI or PropulsionType == c.PARA):
            fileName = fileName + "_Pice_" + str(int(flight.config["pice"]))
        plt.savefig(mission.getResultsPath() + fileName + ".png", format='png', dpi=400)

    if(showPlot):
        plt.show()
    plt.close('all')

    