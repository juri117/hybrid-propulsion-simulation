'''
Created on 01.12.2015

@author: Juri Bieler
'''

from Flight import Flight
from Plane import Plane
from Ini import Ini
import const as c
import generalUtils as generalUtils

import pylab as plt
import time
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy as np
from matplotlib import rc
from matplotlib.patches import Rectangle
from cProfile import label


def getResVectors(flights):
    Pice = []
    fuel = []
    mtom = []
    Cbat = []
    Pmax = []
    
    for flg in flights:
        Pice.append(flg.config["pice"])
        fuel.append(flg.calcUsedFuel())
        mtom.append(max(flg.m))
        Cbat.append(flg.config["cbat"])
        Pmax.append(flg.config["pmax"])
        
    return Pice, fuel, mtom, Cbat, Pmax

def calcHybridDegree(config):
    hybridDegree = 1 - (config["pice"] / config["pmax"])
    return hybridDegree

def findBestConfig(flights):
    bestConfig = {"fuel":float('inf')}
    for f in flights:
        if(f.config["fuel"] < bestConfig["fuel"]):
            bestConfig = f.config
    return bestConfig
    
def plotHybrids(mission, flightsSeri, flightsPara, P_increment, saveFile=True, verbose=True):
    maxs = {"fuel":0., "bat":0., "price":0., "energy":0.}
    for f in flightsSeri:
        maxs["fuel"] = max(maxs["fuel"], f.config["fuel"])
        maxs["bat"] = max(maxs["bat"], f.config["cbat"])
        PriceFuel = f.config["fuel"] * c.PRICE_FUEL
        PriceElec = ((f.SoC[0] - f.SoC[-1]) * f.config["cbat"] / 1000) * c.PRICE_ELEC
        PriceTotal = PriceFuel + PriceElec
        maxs["price"] = max(maxs["price"], PriceTotal)
        Energy = ((f.SoC[0] - f.SoC[-1]) * f.config["cbat"] / 1000) + (f.config["fuel"] * c.ENERGY_FUEL_L)
        maxs["energy"] = max(maxs["energy"], Energy)
        
    for f in flightsPara:
        maxs["fuel"] = max(maxs["fuel"], f.config["fuel"])
        maxs["bat"] = max(maxs["bat"], f.config["cbat"])
        PriceFuel = f.config["fuel"] * c.PRICE_FUEL
        PriceElec = ((f.SoC[0] - f.SoC[-1]) * f.config["cbat"] / 1000) * c.PRICE_ELEC
        PriceTotal = PriceFuel + PriceElec
        maxs["price"] = max(maxs["price"], PriceTotal)
        Energy = ((f.SoC[0] - f.SoC[-1]) * f.config["cbat"] / 1000) + (f.config["fuel"] * c.ENERGY_FUEL_L)
        maxs["energy"] = max(maxs["energy"], Energy)
        
    if(len(flightsSeri) > 0):
        plotHybrid(c.SERI, mission, flightsSeri, P_increment, saveFile=saveFile, verbose=verbose, maxs=maxs)
    if(len(flightsPara) > 0):
        plotHybrid(c.PARA, mission, flightsPara, P_increment, saveFile=saveFile, verbose=verbose, maxs=maxs)

def plotHybrid(propulsionType, mission, flights, P_increment, saveFile=True, verbose=True, maxs={}):

    Pice, fuel, mto, Cbat, Pmax = getResVectors(flights)
    hybDeg = []
    for f in flights:
        hybDeg.append(calcHybridDegree(f.config))

    #Pice = np.divide(Pice, 1000.)
    optConfig = findBestConfig(flights)
    if(len(optConfig) > 0):
        hybDegOpt = calcHybridDegree(optConfig)
    
    ini = Ini(mission.NAME)
    fuelComb =  float(ini.get(c.COMB, "fuel"))
    PoptComb =  float(ini.get(c.COMB, "popt"))
    PoptComb = PoptComb / 1000.
    PmaxComb =  float(ini.get(c.COMB, "pmax"))
    m_pl = c.M_PL
    if(ini.keyExists("general", "mpl")):
        m_pl = float(ini.get("general", "mpl"))
    

    #plt.figure(mission + " - " + propulsionType)
    host = host_subplot(111, axes_class=AA.Axes)
    #plt.subplots_adjust(right=0.75)
    
    par1 = host.twinx()
    par2 = host.twinx()
    
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(60, 0))
    par2.axis["right"].toggle(all=True)
    
    host.set_xlabel("Hybridisierungsgrad")
    host.set_xlim([-0.05, 1.05])
    #host.axes.set_xlim([min(flight.t), max(flight.t)])
    host.set_ylabel("Brennstoff in l")
    yMax = max(fuel)
    if("fuel" in maxs):
        host.set_ylim([0, 1.1 * maxs["fuel"]])
        yMax = maxs["fuel"]
    par1.set_ylabel("Abflugmasse in kg")
    #par1.axes.set_ylim([0, 1.1 * max(mto)])
    par1.set_ylim([(c.M_STRUCT + m_pl), c.MTOM])
    par2.set_ylabel("Batteriekappazität in kWh")
    if("bat" in maxs):
        par2.axes.set_ylim([0, 1.1 * (maxs["bat"] / 1000.)])
    host.invert_xaxis()
    

    p1, = host.plot(hybDeg, fuel, "o", label="$B$")
    p2, = par1.plot(hybDeg, mto, "o", label="$m_{TO}$")
    p3, = par2.plot(hybDeg, np.divide(Cbat, 1000.), "o", label="$C_{Bat}$")
    
    if(len(optConfig) > 0):
        p1_2 = host.plot([hybDegOpt, hybDegOpt], [0, 1.1 * max(fuel)], "b-")
        host.annotate("min", xy=(hybDegOpt, float(yMax) / 2.), xycoords="data", va="center", ha="center", bbox=dict(boxstyle="round", fc="w", alpha=0.5))
    
    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right"].label.set_color(p3.get_color())
    #plt.legend()
    host.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=3)
    
    
    generalUtils.pltCosmetics()
    font = {'size'   : 16}
    rc('font', **font)
    
    host.locator_params(nbins=7)
    plt.subplots_adjust(left=0.12, bottom=0.13, right=0.74, top=0.8, wspace=None, hspace=None)
    #font = {'size'   : 14}
    #rc('font', **font)
    plt.draw()
    if(saveFile):
        plt.savefig(mission.getResultsPath() + "optiMap_" + mission.NAME + "_" + propulsionType +".png", format='png', dpi=400)
    if(verbose):
        plt.show()
    plt.close('all')
    
    plotMass = True
    if(plotMass):
        komponents = ["ICE", "Gen", "Umr2", "B", "EMot", "Umr", "Bat"]
        colors = {"ICE":"FF0000", "Umr":"4FFF4F", "EMot":"00A100", "Bat":"00A6FF", "Gen":"FF7300", "Umr2":"FFBB00", "B":"D1D1D1"}
        #fig = plt.figure()
        ax = host_subplot(111)

        m_pl = flights[0].config["mpl"]
        valuesOld = list(np.ones(len(flights)) * (c.M_STRUCT + m_pl))

        add = 0
        plots = []
        legendKeys = []
        
        #for key, _ in iter(flights[0].mass.items()):
#         for key, _ in iter(sorted(flights[0].mass.items())):
#             values = []
#             for i in range(0, len(flights), 1):
#                 add = valuesOld[i]
#                 values.append(flights[i].mass[key] + add)
# 
#             ax.fill_between(hybDeg, valuesOld, values, facecolor="#" + colors[key], alpha=.7, label=key)
#             plots.append(Rectangle((0, 0), 1, 1, color="#" + colors[key]))
#             legendKeys.append(key)
#             valuesOld = values
            
        for i in range(0, len(komponents), 1):
            key = komponents[i]
            if(key in flights[0].mass):
                values = []
                for i in range(0, len(flights), 1):
                    add = valuesOld[i]
                    values.append(flights[i].mass[key] + add)
    
                ax.fill_between(hybDeg, valuesOld, values, facecolor="#" + colors[key], alpha=.7, label="$m_{" + key + "}$")
                plots.append(Rectangle((0, 0), 1, 1, color="#" + colors[key]))
                legendKeys.append("$m_{" + key + "}$")
                valuesOld = values
            
        #fuel...
        values = np.multiply(fuel, c.DENS_FUEL)
        values = np.add(values, valuesOld)
        ax.fill_between(hybDeg, valuesOld, values, facecolor="#" + colors["B"], alpha=.7, label="B")
        plots.append(Rectangle((0, 0), 1, 1, color="#" + colors["B"]))
        legendKeys.append("$m_{B}$")
            
        plt.xlim([min(hybDeg), max(hybDeg)])
        plt.ylim([(c.M_STRUCT + m_pl), c.MTOM])
        plt.legend(plots, legendKeys, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        generalUtils.pltCosmetics()
        ax.locator_params(nbins=7)
        plt.subplots_adjust(left=0.14, bottom=0.13, right=0.73, top=0.95, wspace=None, hspace=None)
        plt.xlabel("Hybridisierungsgrad")
        ax.set_xlim([-0.05, 1.05])
        plt.ylabel("Masse $m$ in kg")
        ax.invert_xaxis()
        
        if(len(optConfig) > 0):
            plt.plot([hybDegOpt, hybDegOpt], [(c.M_STRUCT + m_pl), c.MTOM], "b-")
            plt.annotate("$m_{B_{min}}$", xy=(hybDegOpt, (c.MTOM + (c.M_STRUCT + m_pl))), xycoords="data", color="b", va="center", ha="center", bbox=dict(boxstyle="round", fc="w", alpha=0.5))
    
        #draw konv. Mass
        p = Plane(mission, propulsionType=c.COMB, config={"pmax":PmaxComb, "fuel":fuelComb, "mpl":m_pl})
        massComb = p.calcActWeight()
        plt.plot([min(min(hybDeg), PoptComb), max(max(hybDeg), PoptComb)], [massComb, massComb], "k-")
        plt.plot([PoptComb], [massComb], "o", label="Masse konv.")

        plt.draw()
        if(saveFile):
            plt.savefig(mission.getResultsPath() + "mass_" + mission.NAME + "_" + propulsionType +".png", format='png', dpi=400)
        if(verbose):
            plt.show()
        plt.close('all')
        
    plotCost = True
    if(plotCost):
        host = host_subplot(111, axes_class=AA.Axes)
        PriceFuel = np.multiply(fuel, c.PRICE_FUEL)
        usedElec = []
        for f in flights:
            usedElec.append((f.SoC[0] - f.SoC[-1]) * f.config["cbat"] / 1000)
            
        PriceElec = np.multiply(usedElec, c.PRICE_ELEC)
        
        PriceTotal = list(np.add(PriceFuel, PriceElec))
        
        host.plot([min(hybDeg), max(hybDeg)], [fuelComb * c.PRICE_FUEL, fuelComb * c.PRICE_FUEL], "k-", label="konv.")
        host.plot(hybDeg, PriceFuel, ".-", label="Brennst.")
        host.plot(hybDeg, PriceElec, ".-", label="Strom")
        host.plot(hybDeg, PriceTotal, ".-", label="Summe")

        
        #host.plot([PoptComb], [fuelComb * c.PRICE_FUEL], "o", label="Sprit konv.")
        
        plt.legend(loc=4)

        minIndex = PriceTotal.index([min(PriceTotal)])
        host.plot([hybDeg[minIndex], hybDeg[minIndex]], [0, max(PriceTotal)], "k-")
        #max(PriceTotal) * 2. / 3.
        
        plt.xlabel("Hybridisierungsgrad")
        plt.xlim([-0.05, 1.05])
        plt.ylabel("Flugkosten in Euro")
        yMax = max(PriceTotal)
        if("price" in maxs):
            host.set_ylim([0, 1.1 * maxs["price"]])
            yMax = maxs["price"]
        host.invert_xaxis()
        
        host.annotate("min", xy=(hybDeg[minIndex], yMax), xycoords="data", va="center", ha="center", bbox=dict(boxstyle="round", fc="w", alpha=0.5))
        
        if(len(optConfig) > 0):
            plt.plot([hybDegOpt, hybDegOpt], [0, yMax], "b-")
            #max(PriceTotal) / 3.
            plt.annotate("$m_{F_{min}}$", xy=(hybDegOpt, float(yMax) / 10.), color="b", xycoords="data", va="center", ha="center", bbox=dict(boxstyle="round", fc="w", alpha=0.5))

        ax.locator_params(nbins=7)
        generalUtils.pltCosmetics()

        plt.draw()
        if(saveFile):
            plt.savefig(mission.getResultsPath() + "costs_" + mission.NAME + "_" + propulsionType +".png", format='png', dpi=400)
        if(verbose):
            plt.show()
        plt.close('all')

    plotBarChart = True
    if(plotBarChart):
        ax = host_subplot(111)
        ax2 = ax.twinx()
        
        
        
        hDLimit = max(hybDeg)
        #fls = []
        hDs = []
        eE = []
        eF = []
        #test = []
        #actP = 0.
        minF = flights[0]
        maxF = flights[0]
        
        #P_increment = P_increment * 2
        barWidth = 0.02
        
        for f in flights:
            P = int(round(f.config["pice"]))
            
            #hD = 1 - (f.config["pice"] / f.config["pmax"])
            hD = calcHybridDegree(f.config)
            #if(hD <= hDLimit):
            if(P % P_increment == 0):
                #fls.append(f)
                hDs.append(hD)
                e = (f.SoC[0] - f.SoC[-1]) * f.config["cbat"] / 1000
                eE.append(e)
                fuel = f.config["fuel"]
                eF.append(fuel)
                hDLimit -= 0.1
                #test.append(f.config["pice"])

            #lastF = f
            if(hD < calcHybridDegree(minF.config)):
                minF = f
            if(hD > calcHybridDegree(maxF.config)):
                maxF = f
        
        #letzten und ersten überschreiben mit min, bzw. max hD-Lösung
        if(len(hDs) > 0):
            hDs[0] = calcHybridDegree(maxF.config)
            e = (maxF.SoC[0] - maxF.SoC[-1]) * maxF.config["cbat"] / 1000
            eE[0] = e
            fuel = maxF.config["fuel"]
            eF[0] = fuel
            hDLimit -= 0.1
            
            hDs[-1] = calcHybridDegree(minF.config)
            e = (minF.SoC[0] - minF.SoC[-1]) * minF.config["cbat"] / 1000
            eE[-1] = e
            fuel = minF.config["fuel"]
            eF[-1] = fuel
            hDLimit -= 0.1
        else:
            f = flights[-1]
            hDs.append(calcHybridDegree(f.config))
            e = (f.SoC[0] - f.SoC[-1]) * f.config["cbat"] / 1000
            eE.append(e)
            fuel = f.config["fuel"]
            eF.append(fuel)
            
        
        
        
        p1 = ax.bar(np.add(hDs, barWidth / 2.), np.multiply(eF, c.ENERGY_FUEL_L), barWidth, color='#FFBB61', label="$E_{chem.}$", zorder=2)
        p2 = ax.bar(np.add(hDs, barWidth / 2.), eE, barWidth, color='#7BFF00', bottom=np.multiply(eF, c.ENERGY_FUEL_L), label="$E_{elekr.}$", zorder=2)
        
        p3 = ax2.bar(np.subtract(hDs, barWidth / 2.), np.multiply(eF, c.PRICE_FUEL), barWidth, color='#BD0000', label="$K_{chem.}$", zorder=2)
        p4 = ax2.bar(np.subtract(hDs, barWidth / 2.), np.multiply(eE, c.PRICE_ELEC), barWidth, color='#008C10', bottom=np.multiply(eF, c.PRICE_FUEL), label="$K_{elekr.}$", zorder=2)
        
        eComb = fuelComb * c.ENERGY_FUEL_L
        p5 = ax.plot([-0.05, 1.05], [eComb, eComb], "k-", label="konv.", zorder=1)
        kComb = fuelComb * c.PRICE_FUEL
        p5 = ax2.plot([-0.05, 1.05], [kComb, kComb], "k--",  zorder=1) #label="$K_{konv.}$",
        
        ax.set_xlim([-0.05, 1.05])
        
        
        quotE = maxs["energy"] / eComb
        quotK = maxs["price"] / kComb
        quot = max(quotE, quotK) * 1.05
        
        ax.set_ylim([0., quot * eComb])
        ax2.set_ylim([0., quot * kComb])
        #ax.set_ylim([0., maxs["energy"] * 1.05])
        #ax2.set_ylim([0., maxs["price"] * 1.05])
        
        
        plt.xlabel("Hybridisierungsgrad")
        ax.set_ylabel("Energie $E$ in kWh")
        ax2.set_ylabel("Betriebskosten $K$ in Euro")
        ax.invert_xaxis()
        ax.locator_params(nbins=7)
        plt.legend(bbox_to_anchor=(1.2, 1), loc=2, borderaxespad=0.)
        generalUtils.pltCosmetics()
        #right=0.67
        plt.subplots_adjust(left=0.12, bottom=0.13, right=0.64, top=0.95, wspace=None, hspace=None)
        #plt.show()
        
        plt.draw()
        if(saveFile):
            plt.savefig(mission.getResultsPath() + "energyCosts_" + mission.NAME + "_" + propulsionType +".png", format='png', dpi=400)
        if(verbose):
            plt.show()
        plt.close('all')

