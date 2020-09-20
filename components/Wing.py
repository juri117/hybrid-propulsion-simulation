'''
Created on 15.12.2015

@author: Juri Bieler
'''

from Propeller import Propeller as Propeller
import generalUtils as generalUtils
import numpy as np
import math
import matplotlib.pyplot as plt
import std_atm as SA
import const as c

class Wing:

    ### CONSTRUCKTOR #####################################
    '''
    berechnet k und C_W0
    '''
    def __init__(self):
        self.verbose = True
        self.k = self.calcK()
        self.C_W0 = self.calcC_W0()
     
    ### FUNCTIONS #####################################

    ### Thrust ###
    
    '''
    Passt die Startleistung (1. Schätzung 60000 W) so lange an, bis sie so dimensioniert ist, dass nach der Strecke "s"
    abgehoben wird
    Rückgabe: P_to, t_abh, v_abh, F_0
    '''
    def calcTakeOffParameter(self, m, h, s, propulsionType, propulsion):
        P_increment = 10000.
        P_to = 60000.
        
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        v_abh = math.sqrt(2. * m * c.g / (rho * c.S * c.C_A_TAKEOFF))
        
        a_req = v_abh**2 / (2 * s)
        
        a_mean, F_0 = self.calcMeanAccelor(m, h, P_to, propulsionType, propulsion.ice.getNmax())
        
        
        while(abs(a_req - a_mean) > 0.000001):
            if(a_mean > a_req and P_increment > 0):
                P_increment = P_increment * (-0.1)
            if(a_mean < a_req and P_increment < 0):
                P_increment = P_increment * (-0.1)
            #print("P -> " + str(P_to))
            P_to += P_increment
            a_mean, F_0 = self.calcMeanAccelor(m, h, P_to, propulsionType, propulsion.ice.getNmax())
        
        t_abh = v_abh / a_mean
        s_abh = (a_mean / 2) * t_abh**2
        
        if(self.verbose):
            print("[Wing] a_mean = " + str(a_mean))
            print("[Wing] t_abh = " + str(t_abh))
            print("[Wing] s_abh = " + str(s_abh))
            print("[Wing] P_to = " + str(P_to))
        
        return P_to, t_abh, v_abh, F_0
    
    '''
    Simuliert einen Startlauf, bis die abhebegeschwindigkeit erreicht ist
    Rückgabe: mittlere Beschleunigung während dem Startlauf "a_mean", Standschub "F0"
    '''
    def calcMeanAccelor(self, m, h, Pmax, propulsionType, NiceMax):
        pr = Propeller()
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        v_abh = math.sqrt(2. * m * c.g / (rho * c.S * c.C_A_TAKEOFF))
        G = m * c.g
        
        #Getribe...
        P_to = Pmax * c.GEAR_EFF[propulsionType]
        A_prop = pr.calcA_prop()
        eta_prop_0 = 0.4 #geschaetzt
        #aus Standschub.pdf
        #F_0 = 0.4 * (2 * rho * A_prop * P_to**2)**(1./3.)
        F_0 = eta_prop_0 * (P_to**(2./3.) * (2 * rho * A_prop)**(1./3.))
        N_pr_to = NiceMax / c.GEAR_RATIO    #1/min

        a_vec = []
        
        dt = 1.
        t = 0.
        v = 0.
        a = 0.
        s = 0.
        
        while(v < v_abh):
            v = v + a
            #Thrust
            eta = 0.5
            if(v > 11):
                eta, N = pr.calcEfficiency(h, v, N_pr_to, P_to)
            F = min(F_0 , eta * (P_to / (v + 0.01)))
            #accelaration
            A = c.C_A_TAKEOFF * c.S * (rho / 2) * v**2
            W = c.C_W_TAKEOFF * c.S * (rho / 2) * v**2
            a = (F - W - (c.WHEEL_MY * (G - A))) / m
            a_vec.append(a)
            t += dt
            if(t > 1200):
                a_mean = np.inf
                a_vec.append(np.inf)
                break
            
        a_mean = np.mean(a_vec)
        return a_mean, F_0
    

    ### Speed ###
    
    '''
    Berechnet Fahrt minimalen Widerstands
    Rückgabe: Fluggeschwindigkeit
    '''
    def calcV_W_min(self, m, h):
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        V_W_min = math.sqrt(2. * m * c.g / (rho * c.S * math.sqrt(self.C_W0 / self.k)))
        return V_W_min
    
    '''
    Berechnet Fahrt mit geringstem Sinken
    Rückgabe: Horizontalanteil der Fluggeschwindigkeit
    '''
    def calcU_g_W_min(self, m, h):
        C_A = self.calcC_aFromC_w(2 * self.C_W0)
        return self.calcU_g(m, h, C_A)
    
    '''
    Berechnet Fahrt mit geringstem Sinken
    Rückgabe: Fluggeschwindigkeit
    '''
    def calcV_dH_min(self, m, h):
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        V_dH_min = math.sqrt(2. * m * c.g / (rho * c.S * math.sqrt(3. * self.C_W0 / self.k)))
        return V_dH_min

    '''
    Berechnet Fahrt mit geringstem Sinken
    Rückgabe: Horizontalanteil der Fluggeschwindigkeit
    '''
    def calcU_g_dH_min(self, m, h):
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        U_g_dH_min = math.sqrt(2. * m * c.g / (rho * c.S * math.sqrt(3. * self.C_W0 / self.k)))
        return U_g_dH_min

    '''
    Berechnet zu "C_A" passende Geschwindigkeit aus Gleitflugpolare
    Rückgabe: Fluggeschwindigkeit
    '''
    def calcV(self, m, h, C_A):
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        gamma = math.atan(self.calcC_WFromC_A(C_A) / C_A)
        V = math.sqrt(2 * m * c.g * math.cos(gamma) / (rho * C_A * c.S))
        return V

    '''
    Berechnet zu "C_A" passende Geschwindigkeit aus Gleitflugpolare
    Rückgabe: Horizontalanteil der Fluggeschwindigkeit
    '''
    def calcU_g(self, m, h, C_A):
        gamma = math.atan(self.calcC_WFromC_A(C_A) / C_A)
        V = self.calcV(m, h, C_A)
        return V * math.cos(gamma)
    
    '''
    Berechnet zu "C_A" passende Geschwindigkeit aus Gleitflugpolare
    Rückgabe: Vertikalanteil der Fluggeschwindigkeit
    '''
    def calcW_g(self, m, h, C_A):
        gamma = math.atan(self.calcC_WFromC_A(C_A) / C_A)
        V = self.calcV(m, h, C_A)
        return V * math.sin(gamma)
    
    
    
    ### Cw, Ca, Eps ###
    
    #source: FM Skript Luckner p. 127
    '''
    bestimmt passendes C_A für horizontalfluggeschwindigkeit "u_g"
    Rückgabe: "C_A"
    '''
    def calcCa(self, m, h, u_g):
        if(u_g <= 0.):
            return 0.
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        C_A = (math.sqrt(2. * m * c.g / (c.S * rho)) * (1. / u_g)) ** 2.
        return C_A
     
    #source: FM Skript Luckner p. 127
    '''
    bestimmt passendes C_W für horizontal und vertikalfluggeschwindigkeit "u_g" und "w_g"
    Rückgabe: "C_W"
    '''
    def calcCw(self, m, h, u_g, w_g):
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        c_w = w_g * self.calcCa(m, h, u_g)**(3./2.) / math.sqrt(2. * m * c.g / (c.S * rho))
        return c_w
    
    '''
    LilienthalPolare("C_A")
    Rückgabe: "C_W"
    '''
    def calcC_WFromC_A(self, C_A):
        C_W = self.C_W0 + (self.k * C_A**2.)
        return C_W
    
    '''
    LilienthalPolare("C_W")
    Rückgabe: "C_A"
    '''
    def calcC_aFromC_w(self, C_W):
        C_A = math.sqrt((C_W - self.C_W0) / self.k)
        return C_A
    
    '''
    berechnet k-Faktor der Lilienthalpolare aus Flügelgeometrie in "const.py"
    Rückgabe: "k"
    '''
    def calcK(self):
        k = 1. / (c.STRETCH * math.pi * c.OSWALD_FAKT)
        return k
     
    '''
    berechnet CW_0 anhand von bestem Gleiten in "const.py"
    Rückgabe: "C_W0"
    '''
    def calcC_W0(self):
        rho = SA.alt2density(0, 'm', 'kg/m**3')
        C_W0 = (2 * c.MTOM * c.g / (rho * c.S * c.V_E_MAX**2))**2 * self.k
        return C_W0
    
    '''
    Rückgabe: Gleitverhältnis
    '''
    def calcEps(self, m, h, u_g):
        C_A = self.calcCa(m, h, u_g)
        c_w = self.calcC_WFromC_A(C_A)
        eps = c_w / C_A
        return eps
    
    '''
    Rückgabe: Gleitzahl
    '''
    def calcE(self, m, h, u_g):
        return 1. / self.calcEps(m, h, u_g)
    

    ### PLOTS #####################################
    
    def plotAll(self):
        self.plotEps()
        self.plotSpeedPolar()
        self.plotLilienthalPolar()

    def plotEps(self):
        h = 0.
        m = 900.
        u_g = range(20, 90, 2)
        eps = list(map(self.calcE, m * np.ones(len(u_g)), h * np.ones(len(u_g)), u_g))
        plt.figure(r"Gleitzahl")
        pG = plt.plot(u_g, eps, '-', linewidth=c.LINE_WIDTH)
        plt.ylabel(r"Gleitzahl $\epsilon$")
        plt.xlabel(r"Fluggeschwindigkeit $TAS$ in m/s")
        #plt.legend(pG, ['Gleitzahl'], loc=4)
        generalUtils.pltCosmetics()
        plt.show()
        plt.close('all')
        
    
    def plotSpeedPolar(self):
        plt.figure("Geschwindigkeitspolare")
        h = 0.
        m = 900.
        c_a = [x / 100.0 for x in range(14, 260, 1)]
        
        u_g = list(map(self.calcU_g, m * np.ones(len(c_a)), h * np.ones(len(c_a)), c_a))
        w_g = list(map(self.calcW_g, m * np.ones(len(c_a)), h * np.ones(len(c_a)), c_a))
        
        #plt.plot(self.polar[:,0], self.polar[:,1], '-')
        plt.plot([0, 3 * c.E_MAX], [0, 3], '--', linewidth=c.LINE_WIDTH)
        plt.plot(u_g, w_g, '-', linewidth=c.LINE_WIDTH)
        C_W0 = self.calcC_W0()
        V_W_min = plt.plot(self.calcU_g_W_min(m, h), self.calcW_g(m, h, self.calcC_aFromC_w(2. * C_W0)), 'o', label="$V_{MW}$", markersize=12)
        plt.plot(self.calcU_g_dH_min(m, h), self.calcU_g_dH_min(m, h) * math.sin(math.atan(self.calcEps(m, h, self.calcU_g_dH_min(m, h)))), 'o', label="$V_{w_{g, min}}$", markersize=12)
        #plt.plot(c.V_E_MAX, c.V_E_MAX / c.E_MAX, 'o', markersize=12)

        plt.grid(True)
        plt.gca().invert_yaxis()
        plt.ylabel("Vertikalfluggeschwindigkeit $w_g$ in m/s")
        plt.xlabel("Horizontalfluggeschwindigkeit $u_g$ in m/s")
        plt.legend(loc=3)
        generalUtils.pltCosmetics()
        plt.show()
        plt.close('all')
        
    def plotLilienthalPolar(self):
        h = 0.
        c_a = [x / 10.0 for x in range(-1, 18, 1)]
        c_w = list(map(self.calcC_WFromC_A, c_a))
        
        plt.figure("Lilienthalpolare")
        #p1 = plt.plot(c_w_guess, c_a_guess, '-')
        p2 = plt.plot(c_w, c_a, '-', linewidth=c.LINE_WIDTH)
        p3 = plt.plot(2. * self.C_W0 , self.calcC_aFromC_w(2. * self.C_W0), 'o', label="$2 C_{W0}$", markersize=12)
        p4 = plt.plot(self.C_W0, 0, 'o', label="$C_{W0}$", markersize=12)
        
        plt.ylabel("Auftriebsbeiwert $C_A$")
        plt.xlabel("Widerstandsbeiwert $C_W$")
        #plt.legend(loc=4)
        axes = plt.gca()
        axes.set_xlim([0, 1.1 * max(c_w)])
        axes.set_ylim([-0.6, 2.])
        plt.legend(loc=4)
        generalUtils.pltCosmetics()
        plt.show()
        plt.close('all')
        
        
    '''
    Plottet parameter über dem Startlauf
    '''
    def plotTakeOffRun(self, m, h):
        pr = Propeller()
        rho = SA.alt2density(h, 'm', 'kg/m**3')
        v_abh = math.sqrt(2 * m * c.g / (rho * c.S * c.C_A_TAKEOFF))
        G = m * c.g
        
        #Getribe...
        P_to = c.GEAR_EFF_COMB * 84500    #kW
        A_prop = pr.calcA_prop()
        #aus Standschub.pdf
        eta_prop_0 = 0.4 #geschaetzt
        F_0 = eta_prop_0 * (P_to**(2./3.) * (2 * rho * A_prop)**(1./3.))
        #F_0 = 0.4 * (2 * rho * A_prop * P_to**2)**(1./3.)
        #F_0 = 1800 #MUELL!!!!!!!!!!!!!!!
        N_pr_to = 5800. / c.GEAR_RATIO #1/min

        #v_abh = 100
        t_vec = []
        v_vec = []
        F_vec = []
        A_vec = []
        W_vec = []
        a_vec = []
        
        t = 0
        v = 0
        a = 0
        
        while(v < v_abh):
            t_vec.append(t)
            #speed
            v = v + a
            v_vec.append(v)
            
            #Thrust
            N = N_pr_to
            eta = 0.5
            if(v > 11):
                eta, N = pr.calcEfficiency(h, v, N_pr_to, P_to)
            F = min(F_0 , eta * (P_to / (v + 0.01)))
            print(str(v) + " -> " + str(eta) + " -> " + str(F) + " -> " + str(N))
            F_vec.append(F)
            
            #accelaration
            A = c.C_A_TAKEOFF * c.S * (rho / 2) * v**2
            A_vec.append(A)
            W = c.C_W_TAKEOFF * c.S * (rho / 2) * v**2
            W_vec.append(W)
            a = (F - W - (c.WHEEL_MY * (G - A))) / m
            a_vec.append(a)
            
            t += 1
            
        a_mean = np.mean(a_vec)
        t_abh = v_abh / a_mean
        s_abh = (a_mean / 2) * t_abh**2
        
        print("t_abh = " + str(t_abh))
        print("s_abh = " + str(s_abh))
        #plt.plot(t_vec, F_vec, label="A")
        plt.plot(t_vec, A_vec, label="A")
        #plt.plot(t_vec, W_vec, label="W")
        #plt.plot(t_vec, a_vec, "o-", label="a")
        #plt.plot(t_vec, v_vec, "x-", label="v")
        plt.xlabel("$Zeit [s]$")
        plt.legend()
        plt.show()
        return s_abh

if __name__ == "__main__":
    w = Wing()
    #w.calcA(900, 0, 0.01)
    #w.plotTakeOffRun(900, 0)
    w.plotAll()
    #w.plotEps()
    #w.plotSpeedPolar()
    #plt.show()
    #w.calcTakeOffParameter(900, 0, 245, c.COMB)
    print("done")