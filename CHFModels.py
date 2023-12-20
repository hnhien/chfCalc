from iapws import IAPWS97
import numpy as np

""" Tong correlation
    P     : Pressure (MPa)
    T     : Temperature (K)
    dTsub : Subcooling (K)
    G     : Mass flow rate (kg/m2s)
    xe    : equivalent quality (-)
"""
def Tong(P=None, G=None, T=None, dTsub=None, xe=None, D=None, modified=None):
    # Fluid properties
    sat_liq = IAPWS97(P=P, x=0)
    sat_stm = IAPWS97(P=P, x=1)
    if not T:
        T = sat_liq.T - dTsub
    water = IAPWS97(P=P, T=T)
    # Latent heat (kJ/kg)
    hfg = sat_stm.h - sat_liq.h
    # Equivalent quality (-)
    if not xe:
        qual = (water.h - sat_liq.h) / hfg
    else:
        qual = xe
    # Dynamic viscousity (Pa.s)
    muf = water.mu
    # Reynolds number
    Re = G * D / muf
    # Considering modification
    if modified==0 or modified==None:
        C = 1.76 - 7.433*qual + 12.222 * qual**2.
        n = 0.6
    elif modified=="celata":
        C = 0.27 + 5.93e-2 * P
        if qual > -0.1:
            C = C * (0.825 + 0.986 * qual)
        n = 0.5
    # Critical heat flux (MW/m2) 
    res = C / Re**n * G * hfg / 1.0e3
    return res



""" Westinghouse
    P (pressure)       : MPa
    T (temperature)    : K
    dTsub (Subcooling) : K
    G (Mass flow rate) : kg/m2s
    h (enthalpy)       : kJ/kg
    q'' (heat flux)    : MW/m2
"""
def W3(P=None, T=None, dTsub=None, G=None, D=None, Pi=None, Ti=None, xe=None, LDh = None):
    # Fluid properties
    sat_liq = IAPWS97(P=P, x=0)
    sat_stm = IAPWS97(P=P, x=1)
    if not T:
        T = sat_liq.T - dTsub
    if not dTsub:
        dTsub = sat_liq.T - T
    hfg = sat_stm.h - sat_liq.h
    loc_liq = IAPWS97(P=P, T=T)
    if not xe:
        qual = (loc_liq.h - sat_liq.h) / hfg
    else:
        qual = xe
    inlet_liq = IAPWS97(P=Pi, T=Ti)
    # CHF
    a = (2.022 - 0.06238 * P) + (0.1722 - 0.01427 * P) * np.exp((18.177 - 0.5987*P)*qual)
    b = (0.1484 - 1.596*qual + 0.1729*qual*abs(qual)) * 2.326*G + 3271
    c = 1.157 - 0.869 * qual
    d = 0.2664 + 0.8357 * np.exp(-124.1*D)
    e = 0.8258 + 3.413e-4 * (sat_liq.h - inlet_liq.h)
    res = a*b*c*d*e / 1e3
    return res