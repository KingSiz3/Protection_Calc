# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 14:35:55 2020

@author: ASGRM
"""
import numpy as np
from scipy.optimize import fsolve
from materials import Material
from dataclasses import dataclass
from tabulate import tabulate


g = {
    'Alpha': u'\u0391',
    'Beta': u'\u0392',
    'Gamma': u'\u0393',
    'Delta': u'\u0394',
    'Epsilon': u'\u0395',
    'Zeta': u'\u0396',
    'Eta': u'\u0397',
    'Theta': u'\u0398',
    'Iota': u'\u0399',
    'Kappa': u'\u039A',
    'Lamda': u'\u039B',
    'Mu': u'\u039C',
    'Nu': u'\u039D',
    'Xi': u'\u039E',
    'Omicron': u'\u039F',
    'Pi': u'\u03A0',
    'Rho': u'\u03A1',
    'Sigma': u'\u03A3',
    'Tau': u'\u03A4',
    'Upsilon': u'\u03A5',
    'Phi': u'\u03A6',
    'Chi': u'\u03A7',
    'Psi': u'\u03A8',
    'Omega': u'\u03A9',
    'alpha': u'\u03B1',
    'beta': u'\u03B2',
    'gamma': u'\u03B3',
    'delta': u'\u03B4',
    'epsilon': u'\u03B5',
    'zeta': u'\u03B6',
    'eta' : u'\u03B7',
    'theta': u'\u03B8',
    'iota': u'\u03B9',
    'kappa': u'\u03BA',
    'lamda': u'\u03BB',
    'mu': u'\u03BC',
    'nu': u'\u03BD',
    'xi': u'\u03BE',
    'omicron': u'\u03BF',
    'pi': u'\u03C0',
    'rho': u'\u03C1',
    'sigma': u'\u03C3',
    'tau': u'\u03C4',
    'upsilon': u'\u03C5',
    'phi': u'\u03C6',
    'chi': u'\u03C7',
    'psi': u'\u03C8',
    'omega': u'\u03C9',
    'degree': '\u00B0C',
}

def calc_con_non_adi_factor(con_mat, ins_mat, S, t, F):
    '''
    Calculation of non-adiabatic factor for conductors and spaced wire screens
    ref: IEC60949/1988 page 11 section 5.1
    Parameters
    ----------
    con_mat : class Metal(Material) (package materials)
        conductor material.
    ins_mat : class Insulator(Material) (package materials)
        insulation material surrounding the conductor.
    S : float
        [mm] cross-section of conductor.
    t : float
        [s] time.
    F : float
        factor to account for imperfect thermal contact.

    Returns
    -------
    epsilon : float
        non-adiabatic factor for conductors.

    '''
    C1 = 2464 # [mm/m]
    C2 = 1.22 # [K*m*mm**2/J]
    sigma_c = con_mat.sigma # [J/K*m**3]
    sigma_i = ins_mat.sigma # [J/K*m**3]
    rho_i = ins_mat.thermal_res
    A = C1 / sigma_c * (sigma_i / rho_i)**(1 / 2) # [(mm**2/s)**(1/2)]
    B = C2 / sigma_c * (sigma_i / rho_i) # [mm**2/s]
    epsilon = (1 + F * A * (t / S)**(1 / 2) + F**2 * B * (t / S))**(1 / 2)
    return epsilon

def calc_screen_non_adi_factor(con_mat, ins_mat_1, ins_mat_2, delta, F, t):
    '''
    Calculation of non-adiabatic factor for sheaths, screens and wires
    ref: IEC60949/1988 page 13 section 6.1
    
    Parameters
    ----------
    con_mat : class Metal(Material) (package materials)
        wire screen material.
    ins_mat_1 : class Insulator(Material) (package materials)
        insulation material surrounding the wire screen.
    ins_mat_2 : class Insulator(Material) (package materials)
        insulation material surrounding he wire screen.
    delta : float
        [mm] thickness of screen, sheath or armour
    F : float
        factor to account for imperfect thermal contact.
    t : float
        [s] time.

    Returns
    -------
    epsilon : float
        non-adiabatic factor for sheaths, screens and wires.

    '''
    sigma_1 = con_mat.sigma
    sigma_2 = ins_mat_1.sigma
    sigma_3 = ins_mat_2.sigma
    rho_2 = ins_mat_1.thermal_res
    rho_3 = ins_mat_2.thermal_res
    
    M = (((sigma_2 / rho_2)**(1 / 2)+(sigma_3 / rho_3)**(1 / 2)) \
         / (2 * sigma_1 * delta * 10**(-3))) * F # s**(-1/2)
    
    epsilon = 1 + 0.61 * M * t**(1 / 2) \
              - 0.069 * (M * t**(1 / 2))**2 \
              + 0.0043 * (M * t**(1 / 2))**3
    return epsilon, M

def calc_thermal_scc_adibatic(con_mat, S, theta_f, theta_i, t):
    '''
    Calculation of adiabatic short-circuit current
    General expression
    ref: IEC60949/1988 page 9 section 3
    
    Parameters
    ----------
    con_mat : class Metal(Material) (package materials)
        wire screen material.
    S : float
        [mm] cross-section of conductor.
    theta_f : float
        [C] Final temperature.
    theta_i : float
        [C] Initial temperature.
    t : float
        time (Short-circuit clearance time).

    Returns
    -------
    I_AD : float
        Permissible adiabatic short-circuit current

    '''
    K = con_mat.K
    beta = con_mat.beta
    I_AD = K * S * (np.log( (theta_f + beta) / (theta_i + beta)))**(1/2) * t**(-1/2)   
    return I_AD

#done #not tested
def calc_con_thermal_scc_non_adiabatic(con_mat, ins_mat, S, theta_f, theta_i, t, F):
    
    epsilon = calc_con_non_adi_factor(con_mat, ins_mat, S, t, F)
    
    I_AD = calc_thermal_scc_adibatic(con_mat, S, theta_f, theta_i, t)
    
    return [I_AD, I_AD*epsilon, epsilon]

#done #not tested
def calc_screen_thermal_scc_non_adiabatic(con_mat, ins_mat_1, ins_mat_2, S, theta_f, theta_i, t, delta, F):
        
    epsilon, M = calc_screen_non_adi_factor(con_mat, ins_mat_1, ins_mat_2, delta, F, t)
    
    I_AD = calc_thermal_scc_adibatic(con_mat, S, theta_f, theta_i, t)
    
    return [I_AD, I_AD*epsilon, epsilon, M]
 
def calc_G(d, P):
    '''
    Geometric factor for xxx
    
    Parameters
    ----------
    d : float
        [mm] mean diameter of the wire screen
    P : float
        [mm] helix pitch of the wires

    Returns
    -------
    float
        Geometric factor
    '''
    return (1 + (np.pi * d / P)**2)**(1 / 2)
    
#mat1, S1, G1, mat2, S2, G2
def calc_limit(mat, S, G, theta_i, theta_f):
    K = mat.K
    rho_20 = mat.rho_20
    beta = mat.beta
    limit = ((K * G * rho_20 *10**6)/(20 + beta))**2 \
        *((theta_f + beta)**2 - (theta_i + beta)**2)
    return limit

def calc_R(mat, S, G, theta):
    rho_20 = mat.rho_20
    beta = mat.beta
    return ((rho_20 * ((theta + beta)/(20 + beta)))/(S * 10**(-6))) * G

def calc_R_f(comp_ref, comp_limit):
    return ((comp_limit.S**2*comp_limit.K**2)/(comp_ref.S**2*comp_ref.K**2)\
                    *(comp_limit.R_f**2-comp_limit.R_i**2)+comp_ref.R_i**2)**(1/2)

def calc_theta_f(R_f, S, G, rho_20, beta):
    theta_f = (R_f*S*10**(-6))/(G*rho_20)*(20+beta)-beta
    return theta_f

class sys(): #name solver?
    def __init__(self, components):
        
        # List of components
        self.components = components
                
        # Calculating permissible energy
        self.E = calc_E(components) 
                
    def calc_Iad(self, t):
        return (self.E/t)**(1/2)
    
    def calc_t(self, Iad):
        return self.E/Iad**2
    
    def calc_S(self, comp_ref, Eref):
        self.Sn = calc_opt_S(self.components, comp_ref, Eref)
        
    def print_results(self):
        print('System')
        print(tabulate([ ['Name' , 'Value',  'Unit'],
                         ['Permissible energy' , f'{self.E:.0f}', '[J]']]))
        
        for comp in self.components:
            print(f'Component: {comp.name}')
            print(tabulate([ ['Name' , 'Value', 'Unit'],
                           ['R\u1DA0' , f'{comp.R_f*1000:.4f}', f'{g["Omega"]}/km'],
                           [f'{g["theta"]}\u1DA0' , f'{comp.theta_f:.2f}', f'{g["degree"]}']]))


class _group():
    _allcomp = []
    
    def add_component(self):
        self._allcomp.append(self)
        self._allcomp.sort()

@dataclass()
class component(_group):
    name      : str
    mat       : Material
    theta_i   : float
    theta_max : float
    S         : float
    d         : float = 0.0
    P         : float = 1.0
    #mylist: list = field(default_factory=list)
    
    @property
    def limit(self):
        return calc_limit(self.mat, self.S, self.G, self.theta_i, self.theta_max)
    
    @property # Initial resistance
    def R_i(self):
        return calc_R(mat=self.mat, S=self.S, G=self.G, theta = self.theta_i)
    
    @property
    def G(self):
        return calc_G(self.d, self.P)
    
    @property
    def R_ui(self):
        return self.R_i*self.S
    
    @property
    def K(self):
        return self.mat.K
    
    @property
    def R_f(self):
        if self._allcomp[0] == self:
            return calc_R(mat=self.mat, S=self.S, G=self.G, theta = self.theta_max)
        else:
            return calc_R_f(self, self._allcomp[0])
    
    @property
    def R_uf(self):
        return self.R_f*self.S
    
    @property
    def theta_f(self):
        return calc_theta_f(R_f = self.R_f,
                            S  = self.S,
                            G  = self.G,
                            rho_20 = self.mat.rho_20,
                            beta = self.mat.beta)
    
    def __post_init__(self):
        self.add_component()
        #self.mylist.append(self)
        
    def __lt__(self, other):
        return self.limit < other.limit

def calc_opt_S(comp, comp_ref, Eref):
    # extracting variables
    R_ui, R_uf, S, K = list(zip(*[[c.R_ui, c.R_uf, c.S, c.K] for c in comp]))
    
    # finding index to component for optimization
    n = comp.index(comp_ref) 
    
    phi = 0
    xrange = list(range(len(comp)))
    del xrange[n]
    
    for j in xrange:
        phi += S[j]*K[j]*np.log((R_uf[n]*K[n]+R_uf[j]*K[j])/(R_ui[n]*K[n]+R_ui[j]*K[j]))
    phi = phi**2 + Eref*np.log(R_uf[n]/R_ui[n])
    temp1 = np.log(R_uf[n]/R_ui[n])
    temp2 = 0
    for j in xrange:
        for k in xrange:  
            temp2 +=  S[j]*K[j]*S[k]*K[k]*np.log((R_uf[j]*K[j]+R_uf[k]*K[k])/(R_ui[j]*K[j]+R_ui[k]*K[k]))
    phi -= temp1*temp2
    
    temp3 = 0
    for j in xrange:
        temp3 += S[j]*K[j]*np.log((R_uf[n]*K[n]+R_uf[j]*K[j])/(R_ui[n]*K[n]+R_ui[j]*K[j]))
    Sn = (phi**(1/2)-temp3)/(K[n]*np.log(R_uf[n]/R_ui[n]))
    return Sn
    
# Note: udtrykket vil altid maksimere energien.
# flere componenter = mere energi.
# men slut resistans/temperaturer forbliver den samme
def eq(S, K, R_i, R_f):
    return np.sum([
        S[j] * S[k] *K[j] *K[k]\
           * np.log(
             (R_f[j] * S[j] * K[j] + R_f[k] * S[k] * K[k])
           / (R_i[j] * S[j] * K[j] + R_i[k] * S[k] * K[k])
                   )      
                   for j in range(len(S))
                       for k in range(len(S))])


def eq2(S, K, R_i, R_f, Q, j, k):
     return S[j] * S[k] *K[j] *K[k]\
            * np.log(
             ((Q+R_i[j]**2 * S[j]**2 * K[j]**2)**(1/2) + (Q+ R_f[k]**2 * S[k]**2 * K[k]**2)**(1/2))
           / (R_i[j] * S[j] * K[j] + R_i[k] * S[k] * K[k])
                   )

def dsigma(func, frm, to, *args):
    #double sigma function
    result = 0
    for j in range(frm, to+1):
        for k in range(frm, to+1): 
            result += func(*args, j, k)
    return result

def calc_E(list_comp):
    # components shall be sorted and initialized for this function to work
    
    # extract variables
    R_i, R_f, S, K = list(zip(*[[c.R_i, c.R_f, c.S, c.K] for c in list_comp]))
    return eq(S, K, R_i, R_f)

def eval_Q(components, Q):
    R_i, R_f, S, K = list(zip(*[[c.R_i, c.R_f, c.S, c.K] for c in components]))

    #myeq = partial(eq, S, K, R_i, R_f)
    return dsigma(eq2, 0, 1, S, K, R_i, R_f, Q)

def cal_Q(components):
    E = calc_E(components)
    R_i, R_f, S, K = list(zip(*[[c.R_i, c.R_f, c.S, c.K] for c in components]))
    
    func = lambda Q : dsigma(eq2, 0, 1, S, K, R_i, R_f, Q ) - E
    Q_initial_guess = 0
    Q = fsolve(func, Q_initial_guess)
    
    return Q


def iterate(components, Eref, dt):
    
    # Initial temperatures
    theta = [comp.theta_i for comp in components]
    
    # Maximal permissible temperatures
    theta_max = [comp.theta_max for comp in components]
  
    for t in np.arange(0, 5+dt, dt):
        
        # Calculate resistance with current temperature
        R = [calc_R(comp.mat, comp.S, comp.G, temp) for comp, temp in zip(components, theta)]
        
        # Calculate resistance ratio
        RT = sum([1/Rj for Rj in R])
        
        # Calculate current in each conductor (Energy normalized at 1s)
        I = [(Eref**(1/2))*((1/Rj)/RT) for Rj in R]
        
        # Break if any temperature reached their maximal limit
        if any([th2 <= th1 for th1 , th2 in zip(theta, theta_max)]):
            #print('Maximal permissible energy!')
            break
        
        # Break if the energy excess the reference level Eref
        if (sum(I)**2*t >= Eref):
            #print('Temperature as a function of energy level!')
            break
        
        # delta T for each conductor
        dT = [np.exp((I_j**2*dt)/(c.mat.K**2*c.S**2))*(theta_j+c.mat.beta)-c.mat.beta-theta_j for I_j, theta_j, c in zip(I, theta, components)]
        
        # incremental of temperatures
        theta = [a + b for a, b in zip(theta, dT)]
        
    return theta, sum(I)**2*(t) #minus dt?

def bisection_search(comps, comp_opt, Eref):   
    E = calc_E(comps)
    a = 0
    b = 100
      
    S_init = comp_opt.S 
    
    if round(Eref) == round(E):
        #print('Cable cannot withstand the energy!')
        return comp_opt.S

    while (abs(E-Eref) or (b - a)/2) > 0.00001: 
        c = (a+b)/2
             
        comp_opt.S = c
        #comps = find_R_f(comps) #updating R
        
        E = calc_E(comps)
        
        if np.sign(E-Eref) == 1:
            b = c
        else:
            a = c
        
    comp_opt.S = S_init
    return c

    
#Tests

#input
'''
t = 0.5
screen = component(name = 'screen',
                  mat = cu,
                    S = 50.1197,
                    theta_i = 70,
                    theta_max = 250,
                    d = 65.25,
                    P = 430)

sheath = component( name = 'sheath',
                    mat = al,
                    S = 45,
                    theta_i = 70,
                    theta_max = 250)


'''
'''
t = 0.52
screen = component( name = 'screen',
                    mat = cu,
                    S = 17.223720290493482,
                    theta_i = 80,
                    theta_max = 250,
                    d = 48.4,
                    P = 387)

sheath = component( name = 'sheath',
                    mat = al,
                    S = 33,
                    theta_i = 80,
                    theta_max = 250)

#comp2.R_i = 7.5835*10**(-4)

s = sys([screen, sheath])
print(s.E)
print((s.E/t)**(1/2))
s.calc_S(screen, 7974**2*t)
print(s.Sn)
print(bisection_search(s.components, s.components[0], 7974**2*0.52))
'''