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


""" Weisman and Pei model
    Dh : Hydrodynamic diameter
    Ph : Heated perimeter
    Ac : Cross-section area
"""
def WeiPei(P=None, T=None, dTsub=None, G=None, L=None, Dh=None, Ph=None, Ac=None):
    # Fluid properties
    sliq = IAPWS97(P=P, x=0)
    sstm = IAPWS97(P=P, x=1)
    if not T:
        T = sliq.T - dTsub
    if not dTsub:
        dTsub = sliq.T - T
    hfg = sstm.h - sliq.h  #kJ/kg
    liq = IAPWS97(P=P, T=T)
    # Calculate single phase heat transfer
    Prf = liq.Prandt
    Ref = G * Dh / liq.mu
    H1p = 0.023 * Ref**0.8 * Prf**0.4 * liq.k / Dh
    # Calculate yb+
    ybp = 0.1 * (liq.sigma * Dh * liq.rho)**0.5 / liq.mu
    # Calcualte friction factor
    rn = 1.0e-4
    ffric = 0.0055 * (1 + (2.0e4 * rn + 1.0e6/Ref)**(1./3.))
    # Calculate multiplier of hf - h1d
    mult = liq.c0 * 1.0e3 / H1p
    tmp = G * (ffric / 8.)**0.5
    if ybp > 0 and ybp <= 5.0:
        mult = mult - Prf * ybp / tmp
    elif ybp <= 30.0:
        tmp1 = Prf + np.log(1.0 + Prf(ybp/5.0 - 1.0))
        mult = mult - 5.0 * tmp1 / tmp
    else:
        tmp1 = Prf + np.log(1.0 + 5.0*Prf) + 0.5 *np.log(ybp / 30.0)
        mult = mult - 5.0 * tmp1 / tmp
    # Iteration
    qchf, err = 1.0e5, 9.9e9
    H0 = 0.075 #[(s.C)^-1]
    qc = 0.0
    while err > 1.0:
        hf_h1d = mult * qchf
        # Iteration to calculate xavg and hl
        xavg, err1 = 0, 9.9
        while err1 > 0.01:
            h1 = (liq.h - xavg * sstm.h) / (1.0 - xavg) * 1.0e3
            liq1 = IAPWS97(P=P, h=h1/1.0e3)
            # Calcualte epsilon = pumping heat flux / evaporative heat flux
            eps = liq1.rho * (sliq.h - h1) / sstm.rho / hfg
            #
            rho_avg = sstm.rho / (sstm.rho / liq1.rho * (1.0 - xavg) + xavg)
            alpha = xavg * rho_avg / sstm.rho
            #
            qb = qchf * (h1 -sliq.h) / hf_h1d * 1.0e3
            if h1 > sliq.h - hf_h1d:
                qc = H0 * hfg / (1.0 / sstm.rho - 1.0 / sliq.rho) * Ac / Ph *alpha * (sstm.T - liq1.T)
            else:
                qc = 0.0
            #
            xavg_new = (qb / (1.0 + eps) - qc) * Ph * L / G / Ac / hfg / 1.0e3
            err = abs(xavg_new - xavg)
        # Calculate ib and psi
        mu_avg = 0




            


