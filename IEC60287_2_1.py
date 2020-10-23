# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 08:50:36 2020

@author: ASGRM
"""

from math import pi, log

####################### 60287-2-1 2015 #######################################
# Table 1 - Thermal resistivities of materials
##############################################################################

Table_1 = {
    'Insulating materials' : {
        'Paper insulation in solid type cables' : 6.0,
        'Paper insulation in oil-filled cables' : 5.0,
        'Paper insulation in cables with external gas pressure' : 5.5,
        'Paper insulation in cables with internal gas pressure, pre-impreganted' : 5.5,
        'Paper insulation in cables with internal gas pressure, mass-impreganted' : 6.0,
        'PE' : 3.5,
        'XLPE' : 3.5,
        'PPL' : 5.5,
        'Polyvinyl chloride, <= 3 kV cables' : 5.0,
        'Polyvinyl chloride, > 3 kV cables'  : 6.0,
        'EPR <= 3 kV cables' : 3.0,
        'EPR, > 3 kV cables'  : 5.0,
        'Butyl rubber' : 5.0,
        'Rubber' : 5.0
        },
    'Protective coverings' : {
        'Compounded jute and fibrous materials' : 6.0,
        'Rubber sandwich protection' : 6.0,
        'Polychloroprene' : 5.5,
        'PVC, <= 35 kV cables' : 5.0,
        'PVC, > 35 kV cables'  : 6.0,
        'PVC/bitumen on corrugated aluminium sheaths' : 6.0,
        'PE' : 3.5
        },
    'Materials for duct installations' : {
        'Concrete' : 1.0,
        'Fibre' : 4.8,
        'Asbestos' : 2.0,
        'Earthenware' : 1.2,
        'PVC' : 6.0,
        'PE' : 3.5        
        }
    }
##############################################################################


def cal_thermal_res(rho_T : float,
                    inner_diameter : float,
                    outer_diameter : float) -> float:
    """
    Calculates the thermal resistance between two components for a single core
    cable.
    Reference: IEC60287_2_1 page 10 and 14 (T_1 / T_2 / T_3)
    """
    T_1 = rho_T / (2 * pi) * log(outer_diameter / inner_diameter)
    return T_1













'''
# can be deleted
def cal_thermal_res_filling(rho_T, inner_diameter, outer_diameter):
    return rho_T / (2 * pi) * log(outer_diameter / inner_diameter)
    #return rho_T/(6*pi)
'''
'''
# can be deleted
def cal_thermal_res2(rho_T, inner_diameter, outer_diameter):
    # Ref page 14 60287
    t_2 = (outer_diameter-inner_diameter)/2
    return rho_T / (2 * pi) * log(1 + 2*t_2 / inner_diameter)
'''