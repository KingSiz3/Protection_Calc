# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 08:48:40 2020

@author: ASGRM
"""

from math import pi, sin

####################### 60287-1-1 2006 #######################################
# Table 2 - Skin and proximity effects - Experimental values for the
# coefficients k_s and k_p.
##############################################################################
def b(d_c, d_i):
    # Not tested
    # d_c [mm] is the inside diameter of the conductor (central duct)
    # d_i [mm] is the outside diameter of equivalent solid conductor having the same central duct.
    return ((d_c - d_i) / (d_c + d_i)) * ((d_c + 2*d_i) / (d_c + d_i))**2

def c(n, b, c):
    # Not tested
    alpha = 1 / (1 + sin(pi /n))**2
    psi = (pi / n + 2 / 3) / (2 * (1 + pi / n))
    k_s =(12 * c * ((alpha * c - 0.5)**2 + (alpha * c - 0.5) * (psi - alpha)\
          * c + 0.33 * (psi - alpha)**2 * c**2 ) + b * (3 - 6 * b + 4 * b**2))**(0.5)
    return k_s 

Table_2 = {
    'cu' : {
        'Round, solid'             : {None  : {'k_s' : 1,     'k_p' : 1}},
        'Round, stranded'          : {'Yes' : {'k_s' : 1,     'k_p' : 0.8},
                                      'No'  : {'k_s' : 1,     'k_p' : 1  }},
        'Round, segmental'         : {None  : {'k_s' : 0.435, 'k_p' : 0.37}},
        'Hollow, helical stranded' : {'Yes' : {'k_s' : b,     'k_p' : 0.8}},
        'Sector-shaped'            : {'Yes' : {'k_s' : 1,     'k_p' : 0.8},
                                      'No'  : {'k_s' : 1,     'k_p' : 1  }},
        },
    'al' : {
        'Round, solid'             : {None     : {'k_s' : 1,     'k_p' : 1}},
        'Round, stranded'          : {'Either' : {'k_s' : 1,     'k_p' : 0.8}},
        'Round, segmental 4'       : {'Either' : {'k_s' : 0.28,  'k_p' : 0.37}},
        'Round, segmental 5'       : {'Either' : {'k_s' : 0.19,  'k_p' : 0.37}},
        'Round, segmental 6'       : {'Either' : {'k_s' : 0.12,  'k_p' : 0.37}},
        'Segmenta with \
         peripheral strands'       : {'Either' : {'k_s' : c,     'k_p' : 0.37}},
        }
    }
##############################################################################