'''
Created on 13.12.2015

@author: Juri Bieler
'''

from units import * 
import generalUtils as generalUtils
import scipy.interpolate as inter
import scipy.optimize
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
import matplotlib.ticker as mtick


'''
Plottet das FEV-Streuband mit original Grafik im Hintergrund
'''
def plotFEV():
    FEV_FILE = 'data/ICE/FEV-Efficiency.data'
    fev = generalUtils.read1dAmeFile(FEV_FILE)
    fev_V_h = fev[:,0]
    fev_sfc = fev[:,1]
    f_fev = inter.interp1d(fev_V_h, fev_sfc, kind='cubic')
    f_fev_l = inter.interp1d(fev_V_h, fev_sfc, kind='linear')
    h = range(0, int(100000 * max(fev_V_h)), 1)
    plt.plot(h, f_fev(np.divide(range(0, int(100000 * max(fev_V_h)), 1),100000)))
    plt.plot(h, f_fev_l(np.divide(range(0, int(100000 * max(fev_V_h)), 1),100000)))    
    plt.show()

def plotFuelConsumption(ice):
    N = [2500., 3000., 3500., 4000., 4500., 5000., 5500.]
    P = [25000., 34000., 42000., 51000., 59000., 66000., 73000.]
    M = np.divide(P, N2omega(np.array(N)))
    
    m_f_vec = list(map(ice.calcConsumption, N, M, 60**2 * np.ones(7)))
    m_f = []
    for m in m_f_vec:
        m_f.append(m[0])
    
    l_h = list(map(ice.calcConsumtionToL_h, m_f, 60**2 * np.ones(7)))
    
    plt.figure("Leistung")
    plt.plot(N, np.divide(P, 1000 * np.ones(7)))
    plt.figure("Verbrauch")
    plt.plot(N, l_h)
    #plt.plot(N, N*p_e)
    plt.show()
    
    print('jo')
    
def plotPmax(ice):
    #datafile = cbook.get_sample_data('data/ICE/P_N.png')
    #img = imread(datafile)
    #plt.scatter(x,y,zorder=1)
    #plt.imshow(img, zorder=0, extent=[0.5, 8.0, 1.0, 7.0])
    
    fig, ax = plt.subplots()
    
    Pfakt = 0.04
    
    plt.xlim([0,5800])
    plt.ylim([0,85000 * Pfakt])
    
    im = plt.imread('../data/ICE/P_N.png')
    implot = plt.imshow(im, extent=[2500, 5800, 0, 85000 * Pfakt])
    
    
    N = range(int(min(ice.Nmax)),int(max(ice.Nmax)) + 1, 1)
    p_e = ice.f_p_emax(N)
    M = list(map(ice.p_e_to_M, p_e))
    Prod = list(np.multiply(list(map(N2omega, N)),M))
    

    imax = Prod.index([max(Prod)])
    Prod = np.divide(Prod, 1/Pfakt)
    P_max = Prod[imax]
    N_max = N[imax]
    
    print("max: " + str(N_max) + " - " + str(P_max))
    
    plt.plot(N, Prod)
    plt.plot([N_max, N_max], [min(Prod), max(Prod)])
    
    prop = Propeller()
    NpropLine = list(map(ice.getNforP, np.divide(Prod, Pfakt), np.ones(len(Prod)) * ice.Pmax))
    NpropLine = list(map(ice.NpropToNice, NpropLine))
    plt.plot(NpropLine, Prod, "-")
    
    plt.draw()
    
    labels = [item.get_text() for item in ax.get_yticklabels()]
    newLabel = []
    for i in range(0, len(labels), 1):
        if(labels[i] != ''):
            lab = float(labels[i])
            lab = (lab / Pfakt)/1000
            newLabel.append(str(lab))
        else:
            newLabel.append(labels[i])

    ax.set_yticklabels(newLabel)

    plt.show()

def plotEfficiencyMap3D(ice):
    p_e = np.array(range(0, 2005000, 5000))
    N = np.array(range(1200,5810,10))
    
    FF = np.zeros([len(p_e), len(N)])
    
    for iPe in range(0, len(p_e), 1):
        for iN in range(0, len(N), 1):
            FF[iPe, iN] = ice.getFuelFlow(N[iN], p_e[iPe])
    
    fig = plt.figure(figsize=plt.figaspect(0.5))        
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    X = N
    Y = p_e
    X, Y = np.meshgrid(X, Y)
    R = FF
    Z = FF
    cmap = generalUtils.custom_div_cmap(numcolors=20, mincol='#007000', maxcol='#FF0000', midcol='#FFFF00')
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cmap, linewidth=0, antialiased=False)
    ax.invert_zaxis()
    
    fig.colorbar(surf, shrink=0.5, aspect=10)
    
    #---- Second subplot
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
    ax.invert_zaxis()

    plt.show()
    
def plotEfficiencyMap(ice, flight=None, fileName=None, showPlot=True, propType=c.COMB):
    fig, ax = plt.subplots(figsize=(11,7))
    fig.canvas.set_window_title('sfcMap')
    
    
    p_e = np.array(range(0, 1810000, 10000))
    N = np.array(range(1200,5810,10))
    
    FF = np.zeros([len(p_e), len(N)])
    
    for iPe in range(0, len(p_e), 1):
        for iN in range(0, len(N), 1):
            if(ice.getMaxP_e(N[iN]) >= p_e[iPe]):
                FF[iPe, iN] = ice.getFuelFlow(N[iN], p_e[iPe])
            else:
                FF[iPe, iN] = FF.max()

    p_e = np.divide(p_e, 100000)

    cmap = generalUtils.custom_div_cmap(numcolors=30, mincol='#008000', maxcol='#DD0000', midcol='#FFFF00')
    bunt = ax.pcolormesh(N, p_e, FF, norm=LogNorm(vmin=FF.min(), vmax=FF.max()), cmap=cmap)

    cbar = fig.colorbar(bunt)        
    uniq = np.unique(ice.FF)
    uniq = uniq[::4]
    cbar.set_ticks(uniq)
    cbar.set_ticklabels(uniq.astype(int))
    cbar.set_label("spez. Kraftstoffverbrauch $sfc$ in g/kWh")
    
    ax.set_ylabel("Mitteldruck $p_e$ in bar")
    ax.set_xlabel("Drehzahl $N$ in 1/min")
    #plt.xlim([min(N), max(N)])

    ax.plot(ice.Nmax, np.divide(ice.p_emax, 100000), 'r--', linewidth=3.0)
    
    Nmax, p_emax = ice.getMaxOpPoint()
    p_emax = p_emax / 100000
    
    
    ax.annotate(r'$\mathbf{P_{max}}$', xy=(Nmax, p_emax), xycoords="data",
              xytext=(Nmax - 400, p_emax + 2), textcoords="data",
              va="top", ha="center", fontsize=20,
              bbox=dict(boxstyle="round", fc="w"),
              arrowprops=dict(arrowstyle="->"))
    
    p_e_opt = ice.p_e_opt / 100000
    ax.annotate(r'$\mathbf{P_{opt}}$', xy=(ice.N_opt, p_e_opt), xycoords="data",
              xytext=(ice.N_opt - 800, p_e_opt + 1), textcoords="data",
              va="top", ha="center", fontsize=20,
              bbox=dict(boxstyle="round", fc="w"),
              arrowprops=dict(arrowstyle="->"))
    
    
    #ax.plot(ice.N_opt, ice.p_e_opt, 'o')
    #ax.plot(Nmax, p_emax, 'o')
    
    #stimmt NUR fuer Konventionell!!!
    #if(propType == c.COMB):
    if(False):
        tas = 40.
        h = 1000.
        if(flight != None):
            tas = sum(flight.tas) / float(len(flight.tas))
            h = sum(flight.H) / float(len(flight.H))
        prop = Propeller()
        
        N_prop = []
        p_e_prop = []

        p_e_prop_min = []
        p_e_prop_max = []
        
        for iN in range(0, len(N), 1):
            for iPe in range(0, len(p_e), 1):
                Pshaft = ice.p_e_to_M(100000 * p_e[iPe]) * N2omega(N[iN])
                Nprop = ice.NiceToNprop(N[iN])
                J = prop.calcJ(tas, Nprop)
                Cp = prop.calcCp(h, Nprop, Pshaft)
                if (not ((Cp < 2 or Cp > 22) or (J < 0.2 or J > 2.2)) and (Nprop <= c.PROP_MAX_N)):
                    p_e_prop.append(p_e[iPe])
                    
            if(len(p_e_prop) > 0):
                p_e_prop_min.append(min(p_e_prop))
                p_e_prop_max.append(max(p_e_prop))
            else:
                p_e_prop_min.append(min(p_e))
                p_e_prop_max.append(min(p_e))
            N_prop.append(N[iN])
            #N_prop = []
            p_e_prop = []
                    
        
        ax.fill_between(N_prop, min(p_e) *  np.ones(len(p_e_prop_min)), p_e_prop_min, facecolor="#000000", alpha=.25, label="prop")
        ax.fill_between(N_prop, p_e_prop_max, max(p_e) * np.ones(len(p_e_prop_max)), facecolor="#000000", alpha=.25, label="prop")
    
    if(flight != None):
        P = flight.Pice
        Nice = flight.Nice
        M = np.divide(P, list(map(N2omega, Nice)))
        p_e_ice = list(map(ice.M_to_p_e, M))
        p_e_ice = np.divide(p_e_ice, 100000)
        
        p_e_ice = [round(x, 1) for x in p_e_ice]
        Nice = [round(x/50, 0)*50 for x in Nice]

        vecLen = len(Nice)
        Nx_vec = []
        p_ex_vec = []
        while(len(Nice) > 0):
            Nx = Nice[0]
            p_ex = p_e_ice[0]
            dens, Nice, p_e_ice = countApperance(Nice, p_e_ice, Nx, p_ex)
            ax.plot(Nx, p_ex, "o", ms=min(max(60 * dens / vecLen, 5), 30), color="black")
            Nx_vec.append(Nx)
            p_ex_vec.append(p_ex)
        
        ax.set_xlim([min(min(N), min(Nx_vec)), max(max(N), max(Nx_vec))])
        ax.set_ylim([min(min(p_e), min(p_ex_vec)), max(max(p_e), max(p_ex_vec))])
    else:
        ax.set_xlim([min(N), max(N)])
    
    stepWidth = max(int(int(ice.Pmax) / 10), 1)
    P = range(0, int(ice.Pmax) + stepWidth, stepWidth)
    N = list(map(ice.calcNIce, P))
    omega = list(map(N2omega, N))
    M = np.divide(P, omega)
    
    p_e = list(map(ice.M_to_p_e, M))
    p_e = np.divide(p_e, 100000)
    ax.plot(N, p_e, 'k-', alpha=0.6)#, linewidth=3.0)
    
    generalUtils.pltCosmetics()

    #plt.draw()
    if(fileName != None):
        fig.savefig(fileName + ".png", format='png', dpi=400)
    if(showPlot):
        plt.show()
    plt.close('all')


def countApperance(aL, bL, a, b):
    count = 0
    aNew = []
    bNew = []
    for i in range(0, len(aL), 1):
        if(aL[i] == a and bL[i] == b):
            count += 1
        else:
            aNew.append(aL[i])
            bNew.append(bL[i])
    return count, aNew, bNew