'''
Created on 13.12.2015

@author: Juri Bieler
'''

from Ini import Ini
from Mission import Mission
from Plane import Plane
from Optimizer import Optimizer
import const as c
import generalUtils
import time
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib import ticker

def insertValueInMesh(mesh, x, y, val):
    while(y >= len(mesh)):
        mesh.append([])
    while(x >= len(mesh[y])):
        mesh[y].append([])
    mesh[y][x] = val
    return mesh

def insertValueInList(vec, val):
    if(val in vec):
        return vec, vec.index(val)
    if(len(vec) == 0):
        vec.append(val)
        return vec, 0
    i = 0
    for i in range(0, len(vec), 1):
        if(val < vec[i]):
            vec.insert(i, val)
            return vec, i;
    vec.append(val)
    return vec, i+1
    
def plotMesh(mis, saveCSV=True):
    results = generalUtils.getAllFolders(pref="results_r", root="data/missions/" + mis.NAME)
    print(str(results))
    Rvec = []
    Vvec = []
    meshComb = []
    meshSeri = []
    meshPara = []
    diffMeshSeriCost = []
    diffMeshParaCost = []

    diffMeshSeriFuel = []
    diffMeshParaFuel = []

    diffMeshSeriElec = []
    diffMeshParaElec = []
    hybDegSeriMesh = []
    hybDegParaMesh = []
    
    csvData = []
    for res in results:
        if(res == "results_r0260.0_v35.0"):
            print("stop")
        r = float(res.split("_")[1][1:])
        Rvec, Ri = insertValueInList(Rvec, r)
        v = float(res.split("_")[2][1:])
        Vvec, Vi = insertValueInList(Vvec, v)
        ini = Ini(mis.NAME, resPath=res+"/")
        
        #csvSeri = generalUtils.readCsv(mis.getResultsPath() + "/runningConfigs_serialHybrid.csv")
        #csvPara = generalUtils.readCsv(mis.getResultsPath() + "/runningConfigs_parallelHybrid.csv")

        #fuel / 100km
        fuel100Comb = 100. * float(ini.get(c.COMB, "fuel")) / (float(ini.get("output", "distance")) / 1000.)

        dist = float(ini.get("output", "distance")) / 1000.
        fuel100Comb = 100. * float(ini.get(c.COMB, "fuel")) / dist
        fuel100Seri = 100. * float(ini.get(c.SERI, "fuel")) / dist
        fuel100Para = 100. * float(ini.get(c.PARA, "fuel")) / dist

        elec100Seri = 0.1 * float(ini.get(c.SERI, "cbat")) / dist
        elec100Para = 0.1 * float(ini.get(c.PARA, "cbat")) / dist

        #costs
        valComb = fuel100Comb  * c.PRICE_FUEL
        valSeri = (fuel100Seri  * c.PRICE_FUEL) + (elec100Seri * c.PRICE_ELEC)
        valPara = (fuel100Para  * c.PRICE_FUEL) + (elec100Para * c.PRICE_ELEC)
        
        #energy
        valComb = fuel100Comb  * c.ENERGY_FUEL_L
        valSeri = (fuel100Seri  * c.ENERGY_FUEL_L) + (elec100Seri)
        valPara = (fuel100Para  * c.ENERGY_FUEL_L) + (elec100Para)
        
        diffMeshSeriCost = insertValueInMesh(diffMeshSeriCost, Ri, Vi, (valSeri / valComb) * 100.)
        diffMeshParaCost = insertValueInMesh(diffMeshParaCost, Ri, Vi, (valPara / valComb) * 100.)

        diffMeshSeriFuel = insertValueInMesh(diffMeshSeriFuel, Ri, Vi, fuel100Seri)
        diffMeshParaFuel = insertValueInMesh(diffMeshParaFuel, Ri, Vi, fuel100Para)

        diffMeshSeriElec = insertValueInMesh(diffMeshSeriElec, Ri, Vi, elec100Seri)
        diffMeshParaElec = insertValueInMesh(diffMeshParaElec, Ri, Vi, elec100Para)
            
        hybDegSeri = 1. - (float(ini.get(c.SERI, "pice")) / float(ini.get(c.SERI, "pmax")))
        hybDegSeriMesh = insertValueInMesh(hybDegSeriMesh, Ri, Vi, hybDegSeri)
        
        hybDegPara = 1. - (float(ini.get(c.PARA, "pice")) / float(ini.get(c.PARA, "pmax")))
        hybDegParaMesh = insertValueInMesh(hybDegParaMesh, Ri, Vi, hybDegPara)
        
        meshComb = insertValueInMesh(meshComb, Ri, Vi, valComb)
        meshSeri = insertValueInMesh(meshSeri, Ri, Vi, valSeri)
        meshPara = insertValueInMesh(meshPara, Ri, Vi, valPara)
        
        

        if(saveCSV):
            plConv = Plane(None, propulsionType=c.COMB, config=ini.getConfig(c.COMB))
            plSeri = Plane(None, propulsionType=c.SERI, config=ini.getConfig(c.SERI))
            plPara = Plane(None, propulsionType=c.PARA, config=ini.getConfig(c.PARA))
            
            
            plConv = Plane(None, propulsionType=c.COMB, config=ini.getConfig(c.COMB))
            plSeri = Plane(None, propulsionType=c.SERI, config=ini.getConfig(c.SERI))
            plPara = Plane(None, propulsionType=c.PARA, config=ini.getConfig(c.PARA))
            csvRow = {"TravelDistance(km)":r,
                      "TravelSpeed(m/s)":v,
                      "TotalDistance(km)":dist,
                      "FlightTime":ini.get("output", "flighttime"),
                      
                      "CONV:TakeOffMass(kg)":plConv.calcActWeight(),
                      "CONV:NeededFuel(l)":float(ini.get(c.COMB, "fuel")),
                      #"CONV:NeededElectricity(Wh)":0.,
                      "CONV:Costs(Euro)":valComb,
                      
                      "SERI:TakeOffMass(kg)":plSeri.calcActWeight(),
                      "SERI:NeededFuel(l)":float(ini.get(c.SERI, "fuel")),
                      "SERI:NeededElectricity(Wh)":float(ini.get(c.SERI, "cbat")),
                      "SERI:Costs(Euro)":valSeri,
                      "SERI:HybDeg":hybDegSeri,
                      
                      "PARA:TakeOffMass(kg)":plPara.calcActWeight(),
                      "PARA:NeededFuel(l)":float(ini.get(c.PARA, "fuel")),
                      "PARA:NeededElectricity(Wh)":float(ini.get(c.PARA, "cbat")),
                      "PARA:Costs(Euro)":valPara,
                      "PARA:HybDeg":hybDegPara}             
            csvData.append(csvRow)
            print(res + " - DONE!")

    #write CSV
    if(saveCSV):
        csvHead = ["TravelDistance(km)",
                   "TravelSpeed(m/s)",
                   "TotalDistance(km)",
                   "FlightTime",
                   "CONV:TakeOffMass(kg)",
                   "CONV:NeededFuel(l)",
                   #"CONV:NeededElectricity(Wh)",
                   "CONV:Costs(Euro)",
                   "SERI:TakeOffMass(kg)",
                   "SERI:NeededFuel(l)",
                   "SERI:NeededElectricity(Wh)",
                   "SERI:Costs(Euro)",
                   "SERI:HybDeg",
                   "PARA:TakeOffMass(kg)",
                   "PARA:NeededFuel(l)",
                   "PARA:NeededElectricity(Wh)",
                   "PARA:Costs(Euro)",
                   "PARA:HybDeg"]
        generalUtils.write2CSV(mis.PATH + "meshOpti.csv", csvData, heads=csvHead)

    Rvec = np.array(Rvec)
    Vvec = np.array(Vvec)
#     if(False):
#         meshComb = np.array(meshComb)
#         meshPara = np.array(meshPara)
#         meshSeri = np.array(meshSeri)
#         fig = plt.figure(figsize=(16,10))        
#         ax = fig.add_subplot(1, 1, 1, projection='3d')
#         fig.canvas.set_window_title('mesh')  
#         X, Y = np.meshgrid(Rvec, Vvec)
#         ax.plot_wireframe(X, Y, meshComb, rstride=1, cstride=1)
#         
#         ax.plot_wireframe(X, Y, meshPara, rstride=1, cstride=1, color="r")
#         ax.plot_wireframe(X, Y, meshSeri, rstride=1, cstride=1, color="g")
# 
#         ax.set_ylabel("Fluggeschwindigkeit in m/s")
#         ax.set_xlabel("Reichweite in km")
#         ax.set_zlabel("Brennstoffverbrauch auf 100 km")
#         #plt.ylim([0, max(np.amax(meshComb), max(np.amax(meshSeri), np.amax(meshPara)))])
#         plt.legend(["konventionell", "parallel", "seriell"])
#         plt.show()
    
    preparePlots(Rvec, Vvec, diffMeshSeriCost, diffMeshParaCost, diffMeshSeriFuel, diffMeshParaFuel, diffMeshSeriElec, diffMeshParaElec, hybDegSeriMesh, hybDegParaMesh, mis.PATH + mis.NAME)

def plotOneConfig(mis, saveCSV=True):
    results = generalUtils.getAllFolders(pref="results_", root="data/missions/" + mis.NAME)
    print(str(results))
    
    #ini = Ini(mis.NAME, resPath=results[-1]+"/")
    ini = Ini(mis.NAME, resPath="results_r0380.0_v45.0/")
    configComb = ini.getConfig(c.COMB)
    configSeri = ini.getConfig(c.SERI)
    configPara = ini.getConfig(c.PARA)
    
    #configSeri = {'mpl': 110, 'pmax': 88649.19999999995, 'pice': 51640.0, 'cbat': 15396.5684111, 'fuel': 39.832351270753406}
    
    Rvec = []
    #Svec = []
    Vvec = []
    meshComb = []
    meshSeri = []
    meshPara = []
    diffMeshSeriCost = []
    diffMeshParaCost = []

    diffMeshSeriFuel = []
    diffMeshParaFuel = []

    diffMeshSeriElec = []
    diffMeshParaElec = []
    csvData = []
    for res in results:
        print(res)
        #if(res == "results_r0020.0_v75.0"):
        #    print("stop")
        m = Mission(mis.NAME + "/" + res)
        r = float(res.split("_")[1][1:])
        Rvec, Ri = insertValueInList(Rvec, r)
        v = float(res.split("_")[2][1:])
        Vvec, Vi = insertValueInList(Vvec, v)
        
        opt = Optimizer(m, verbose=True, savePlots=False)
        fComb = opt.flyWithNeededFuel(m, c.COMB, configComb)
        fSeri = opt.flyWithNeededFuel(m, c.SERI, configSeri)
        fPara = opt.flyWithNeededFuel(m, c.PARA, configPara)
        
        
        dist = m.calcDistance() / 1000.
        fuel100Comb = 100. * fComb.calcUsedFuel() / dist
        fuel100Seri = 100. * fSeri.calcUsedFuel() / dist
        fuel100Para = 100. * fPara.calcUsedFuel() / dist

        elec100Seri = 0.1 * fSeri.calcUsedElec() / dist
        elec100Para = 0.1 * fPara.calcUsedElec() / dist

        valComb = fuel100Comb  * c.PRICE_FUEL
        valSeri = (fuel100Seri  * c.PRICE_FUEL) + (elec100Seri * c.PRICE_ELEC)
        valPara = (fuel100Para  * c.PRICE_FUEL) + (elec100Para * c.PRICE_ELEC)
        

        meshComb = insertValueInMesh(meshComb, Ri, Vi, valComb)
        meshSeri = insertValueInMesh(meshSeri, Ri, Vi, valSeri)
        meshPara = insertValueInMesh(meshPara, Ri, Vi, valPara)
        
        if(len(fSeri.fails) > 0):
            diffMeshSeriCost = insertValueInMesh(diffMeshSeriCost, Ri, Vi, float("nan"))
            diffMeshSeriFuel = insertValueInMesh(diffMeshSeriFuel, Ri, Vi, float("nan"))
            diffMeshSeriElec = insertValueInMesh(diffMeshSeriElec, Ri, Vi, float("nan"))
        else:
            diffMeshSeriCost = insertValueInMesh(diffMeshSeriCost, Ri, Vi, (valSeri / valComb) * 100.)
            diffMeshSeriFuel = insertValueInMesh(diffMeshSeriFuel, Ri, Vi, fuel100Seri)
            diffMeshSeriElec = insertValueInMesh(diffMeshSeriElec, Ri, Vi, elec100Seri)
            
        if(len(fPara.fails) > 0):
            diffMeshParaCost = insertValueInMesh(diffMeshParaCost, Ri, Vi, float("nan"))
            diffMeshParaFuel = insertValueInMesh(diffMeshParaFuel, Ri, Vi, float("nan"))
            diffMeshParaElec = insertValueInMesh(diffMeshParaElec, Ri, Vi, float("nan"))
        else:
            diffMeshParaCost = insertValueInMesh(diffMeshParaCost, Ri, Vi, (valPara / valComb) * 100.)
            diffMeshParaFuel = insertValueInMesh(diffMeshParaFuel, Ri, Vi, fuel100Para)
            diffMeshParaElec = insertValueInMesh(diffMeshParaElec, Ri, Vi, elec100Para)
        

        if(saveCSV):
            plConv = Plane(None, propulsionType=c.COMB, config=ini.getConfig(c.COMB))
            plSeri = Plane(None, propulsionType=c.SERI, config=ini.getConfig(c.SERI))
            plPara = Plane(None, propulsionType=c.PARA, config=ini.getConfig(c.PARA))
            
            
            plConv = Plane(None, propulsionType=c.COMB, config=ini.getConfig(c.COMB))
            plSeri = Plane(None, propulsionType=c.SERI, config=ini.getConfig(c.SERI))
            plPara = Plane(None, propulsionType=c.PARA, config=ini.getConfig(c.PARA))
            csvRow = {"TravelDistance(km)":r,
                      "TravelSpeed(m/s)":v,
                      "TotalDistance(km)":dist,
                      "FlightTime":ini.get("output", "flighttime"),
                      
                      "CONV:TakeOffMass(kg)":plConv.calcActWeight(),
                      "CONV:NeededFuel(l)":float(ini.get(c.COMB, "fuel")),
                      "CONV:NeededElectricity(Wh)":0.,
                      
                      "SERI:TakeOffMass(kg)":plSeri.calcActWeight(),
                      "SERI:NeededFuel(l)":float(ini.get(c.SERI, "fuel")),
                      "SERI:NeededElectricity(Wh)":float(ini.get(c.SERI, "cbat")),
                      
                      "PARA:TakeOffMass(kg)":plPara.calcActWeight(),
                      "PARA:NeededFuel(l)":float(ini.get(c.PARA, "fuel")),
                      "PARA:NeededElectricity(Wh)":float(ini.get(c.PARA, "cbat"))}             
            csvData.append(csvRow)
            print(res + " - DONE!")
            
        
    #write CSV
    if(saveCSV):
        csvHead = ["TravelDistance(km)",
                   "TravelSpeed(m/s)",
                   "TotalDistance(km)",
                   "FlightTime",
                   "CONV:TakeOffMass(kg)",
                   "CONV:NeededFuel(l)",
                   "CONV:NeededElectricity(Wh)",
                   "SERI:TakeOffMass(kg)",
                   "SERI:NeededFuel(l)",
                   "SERI:NeededElectricity(Wh)",
                   "PARA:TakeOffMass(kg)",
                   "PARA:NeededFuel(l)",
                   "PARA:NeededElectricity(Wh)"]
        generalUtils.write2CSV(mis.PATH + "meshOneConfig.csv", csvData, heads=csvHead)

    #fig, ax = plt.subplots(figsize=(16,10))
      
    #cmap = generalUtils.custom_div_cmap(numcolors=30, mincol='#008000', maxcol='#DD0000', midcol='#FFFF00')
    #norm=LogNorm(vmin=meshComb.min(), vmax=meshComb.max()), 
    
    Rvec = np.array(Rvec)
    Vvec = np.array(Vvec)
#     if(False):
#         meshComb = np.array(meshComb)
#         meshPara = np.array(meshPara)
#         meshSeri = np.array(meshSeri)
#         fig = plt.figure(figsize=(16,10))        
#         ax = fig.add_subplot(1, 1, 1, projection='3d')
#         fig.canvas.set_window_title('mesh')  
#         X, Y = np.meshgrid(Rvec, Vvec)
#         ax.plot_wireframe(X, Y, meshComb, rstride=1, cstride=1)
#         
#         ax.plot_wireframe(X, Y, meshPara, rstride=1, cstride=1, color="r")
#         ax.plot_wireframe(X, Y, meshSeri, rstride=1, cstride=1, color="g")
# 
#         ax.set_ylabel("Fluggeschwindigkeit in m/s")
#         ax.set_xlabel("Reichweite in km")
#         ax.set_zlabel("Brennstoffverbrauch auf 100 km")
#         #plt.ylim([0, max(np.amax(meshComb), max(np.amax(meshSeri), np.amax(meshPara)))])
#         plt.legend(["konventionell", "parallel", "seriell"])
#         plt.show()

    preparePlots(Rvec, Vvec, diffMeshSeriCost, diffMeshParaCost, diffMeshSeriFuel, diffMeshParaFuel, diffMeshSeriElec, diffMeshParaElec, hybDegSeriMesh, hybDegParaMesh, mis.PATH + mis.NAME)
        
def preparePlots(Rvec, Vvec, diffMeshSeriCost, diffMeshParaCost, diffMeshSeriFuel, diffMeshParaFuel, diffMeshSeriElec, diffMeshParaElec, hybDegSeriMesh, hybDegParaMesh, path):
    diffMeshSeriCost = np.array(diffMeshSeriCost)
    diffMeshParaCost = np.array(diffMeshParaCost)
    minColorCost = min(np.nanmin(diffMeshSeriCost), np.nanmin(diffMeshParaCost))
    maxColorCost = max(np.nanmax(diffMeshSeriCost), np.nanmax(diffMeshParaCost))
    limitCost = [minColorCost, maxColorCost]
    limitCost = [50, 150]
    
    
    diffMeshSeriFuel = np.array(diffMeshSeriFuel)
    diffMeshParaFuel = np.array(diffMeshParaFuel)
    minColorFuel = min(np.nanmin(diffMeshSeriFuel), np.nanmin(diffMeshParaFuel))
    maxColorFuel = max(np.nanmax(diffMeshSeriFuel), np.nanmax(diffMeshParaFuel))
    limitFuel = [minColorFuel, maxColorFuel]
    
    diffMeshSeriElec = np.array(diffMeshSeriElec)
    diffMeshParaElec = np.array(diffMeshParaElec)
    minColorElec = min(np.nanmin(diffMeshSeriElec), np.nanmin(diffMeshParaElec))
    maxColorElec = max(np.nanmax(diffMeshSeriElec), np.nanmax(diffMeshParaElec))
    limitElec = [minColorElec, maxColorElec]
    
    plotDiffMap(Rvec, Vvec, diffMeshSeriFuel, diffMeshSeriElec, diffMeshSeriCost, limitFuel, limitElec, limitCost, hybDegSeriMesh, path, "seri")
    plotDiffMap(Rvec, Vvec, diffMeshParaFuel, diffMeshParaElec, diffMeshParaCost, limitFuel, limitElec, limitCost, hybDegParaMesh, path, "para")

def plotDiffMap(Rvec, Vvec, diffMeshFuel, diffMeshElec, diffMeshCost, limitFuel, limitElec, limitCost, hybDegMesh, path, propType):
    
    Vvec = list(Vvec)
    Vvec.append(max(Vvec) + Vvec[1] - Vvec[0])
    Vvec = np.array(Vvec)
    Vvec = np.add(Vvec, (Vvec[0] - Vvec[1]) / 2)
    
    diffMeshFuel = np.ma.masked_invalid(diffMeshFuel)
    diffMeshElec = np.ma.masked_invalid(diffMeshElec)
    diffMeshCost = np.ma.masked_invalid(diffMeshCost)
    
    fig = plt.figure(figsize=(13,9))
    ax = host_subplot(311)
    ax.set_title("Brennstoffverbrauch")
    #ax.set_ylabel("Fluggeschwindigkeit in m/s")
    plt.ylim([min(Vvec), max(Vvec)])
    plt.xlim([min(Rvec), max(Rvec)])
    ax.set_xlabel("")
    ax.set_xticklabels([])
    
    cm = generalUtils.custom_div_cmap(numcolors=256, mincol='#FFFFFF', maxcol='#DD0000', midcol='#00CC00')
    cm.set_bad(color = '#000000', alpha = 1.)
    
    bunt = ax.pcolormesh(Rvec, Vvec, diffMeshFuel, vmin=limitFuel[0], vmax=limitFuel[1], cmap=cm)#, norm=LogNorm(vmin=FF.min(), vmax=FF.max()), cmap=cmap)
    cbar = fig.colorbar(bunt)
    cbar.set_label(r"Brennstoff $B$ in $\frac{\mathrm{l}}{100\mathrm{km}}$")
    
    ax = host_subplot(312)
    ax.set_title("Stromverbrauch")
    ax.set_xlabel("")
    ax.set_xticklabels([])
    
    
    cm = generalUtils.custom_div_cmap(numcolors=256, mincol='#FFFFFF', maxcol='#DD0000', midcol='#00CC00')
    cm.set_bad(color = '#000000', alpha = 1.)
    
    bunt = ax.pcolormesh(Rvec, Vvec, diffMeshElec, vmin=limitElec[0], vmax=limitElec[1], cmap=cm)#, norm=LogNorm(vmin=FF.min(), vmax=FF.max()), cmap=cmap)
    cbar = fig.colorbar(bunt)
    cbar.set_label(r"Strom $E_{Bat.}$ in  $\frac{\mathrm{kWh}}{100\mathrm{km}}$")
    
    ax.set_ylabel(r"Fluggeschwindigkeit $v$ in m/s")
    plt.ylim([min(Vvec), max(Vvec)])
    plt.xlim([min(Rvec), max(Rvec)])

    ax = host_subplot(313)
    ax.set_title("Energie")
    
#     cm = generalUtils.custom_div_cmap(numcolors=256, mincol='#00CC00', maxcol='#FF0000', midcol='#FFFFFF')
#     cm.set_bad(color = '#000000', alpha = 1.)
#     
#     bunt = ax.pcolormesh(Rvec, Vvec, diffMeshCost, vmin=limitCost[0], vmax=limitCost[1], cmap=cm)#, norm=LogNorm(vmin=FF.min(), vmax=FF.max()), cmap=cmap)
#     cbar = fig.colorbar(bunt)
#     cbar.set_label("Kosten in % von Konv.")
#     
#     #ax.set_ylabel("Fluggeschwindigkeit in m/s")
#     plt.ylim([min(Vvec), max(Vvec)])
#     plt.xlim([min(Rvec), max(Rvec)])
#     ax.set_xlabel("Reichweite in km")

    plotCompare(fig, ax, Rvec, Vvec, diffMeshCost, limitCost, hybDegMesh)
    
    generalUtils.pltCosmetics()
    plt.savefig(path + "_" + propType + "_overview_energy.png", format='png', dpi=400)
    #plt.show()
    plt.close('all')
    
    plotOnlyComp(Rvec, Vvec, diffMeshCost, limitCost, hybDegMesh, path, propType)
    
    
def plotCompare(fig, ax, Rvec, Vvec, diffMeshCost, limitCost, hybDegMesh):
    cm = generalUtils.custom_div_cmap(numcolors=256, mincol='#00CC00', maxcol='#FF0000', midcol='#FFFFFF')
    cm.set_bad(color = '#000000', alpha = 1.)
    
    bunt = ax.pcolormesh(Rvec, Vvec, diffMeshCost, vmin=limitCost[0], vmax=limitCost[1], cmap=cm)
    cbar = fig.colorbar(bunt)
    cbar.set_label("Energie $E$ in %")
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    
    #ax.set_ylabel("Fluggeschwindigkeit in m/s")
    plt.ylim([min(Vvec), max(Vvec)])
    plt.xlim([min(Rvec), max(Rvec)])
    ax.set_xlabel("Reichweite $R$ in km")
    
    for iR in range(0, len(Rvec), 1):
        for iV in range(0, len(Vvec)-1, 1):
            hybDeg = hybDegMesh[iV][iR]
            if(np.isnan(hybDeg)):
                hybDeg = 1.
                
            x = Rvec[iR] + 10
            y = Vvec[iV] + 2.5
            plt.plot(x, y, "ko", alpha=hybDeg)


def removeNaN(mesh, setTo=0):
    for row in range(0, len(mesh), 1):
        for col in range(0, len(mesh[row]), 1):
            if(math.isnan(mesh[row][col])):
                mesh[row][col] = setTo
    return mesh

def plotOnlyComp(Rvec, Vvec, diffMeshCost, limitCost, hybDegMesh, path, propType):
    fig = plt.figure(figsize=(16,4))
    ax = host_subplot(111)
    plotCompare(fig, ax, Rvec, Vvec, diffMeshCost, limitCost, hybDegMesh)
    generalUtils.pltCosmetics(multiPlot=False)
    plt.subplots_adjust(left=0.04, bottom=0.08, right=1.2, top=0.95, wspace=None, hspace=None)
    plt.savefig(path + "_" + propType + "_comp.png", format='png', dpi=400)
    #plt.show()
    plt.close('all')
    