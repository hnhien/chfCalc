from iapws import IAPWS97
import numpy as np
import math

""" Common function
"""
def Reynolds(G=None, Dh=None, mu=None):
    return G * Dh / mu


def fricFactor(G=None, Dh=None, mu=None, roughness=None):
    rn = roughness or 1.0e-4
    Re = Reynolds(G=G,Dh=Dh,mu=mu)
    return 5.5e-3 * (1 + (2.0e4 * (rn/Dh) + 1.0e6/Re)**(1.0/3.0))

def wallShearStress(G=None, Dh=None, rho=None, mu=None, roughness=None):
    ff = fricFactor(G=G, Dh=Dh, mu=mu, roughness=roughness)
    return ff / 8.0 / rho * G * G

def dbLahey(G=None, Dh=None, rhol=None, rhog=None, mul=None, sigma=None, roughness=None):
    g = 9.81
    tau = wallShearStress(G=G, Dh=Dh, rho=rhol,mu=mul)
    return 0.015 * (sigma * Dh / tau)**0.5 / (1 + 0.1 * g * (rhol - rhog)/tau * Dh)**0.5

def States(P=None, T=None, dTsub=None, xe=None):
    sliq = IAPWS97(P=P, x=0)
    sstm = IAPWS97(P=P, x=1)
    if not T:
        T = sliq.T - dTsub
    if not dTsub:
        dTsub = sliq.T - T    
    liq = IAPWS97(P=P, T=T)
    hfg = sstm.h - sliq.h
    if not xe:
        xe = (liq.h - sliq.h) / hfg
    return (sliq, sstm, liq, T, dTsub, hfg, xe)


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


"""
    Weisman and Ileslamlou model
"""            
def WeisIles(P=None, T=None, dTsub=None, G=None, xe=None, L=None, Dh=None, Ph=None):
    # Liquid in bubbly layer is assumed saturated
    sliq, sstm, liq, T, dTsub, hfg, x2 = States(P=P, T=T, dTsub=dTsub)
    # Eqn 18 multiplier
    f18 = (sliq.h * (1.0 - x2) + sstm.h * x2 - liq.h) * 1.0e3
    # Db lahey
    Db = dbLahey(G=G, Dh=Dh, rhol=liq.rho, rhog=sstm.rho, mul=liq.mu, sigma=sstm.sigma)
    # calculate ib
    k = 2.4
    Gs2h = G * 3600.0
    if Gs2h > 9.7e6:
        a = 0.135 / (G/9.7E6)**0.3
    else:
        a = 0.135   
    Re = Reynolds(G=G, Dh=Dh, mu=liq.mu)
    ib = 0.462 * k**0.6 / Re**0.1 * (Db/Dh)**0.6 * (1 + a * (liq.rho - sstm.rho)/sstm.rho)
    # Iteratiion
    qw, err = 1.0e5, 1.0e9
    while err > 1.0:
        v11 = qw * x2 / f18 / sstm.rho
        rho_ = liq.rho # correct need?
        sv_ = ib * G / rho_
        x = v11 / sv_ / 2**0.5
        psi = np.exp(-x*x) / np.sqrt(2 * np.pi) - x / 2**0.5 * math.erfc(x)
        qw_cal = G * psi * ib * f18
        err = abs(qw_cal - qw)
        qw = qw_cal
    return qw
        

def newton_raphson(f, df, x0, tol):
    if abs(f(x0)) < tol:
        return x0
    return newton_raphson(f, df, x0 - f(x0)/df(x0), tol)
    
def Katto90(P=None, T=None, dTsub=None, G=None, xe=None, L=None, Dh=None, Ph=None):
    # Satutated state
    sliq, sstm, liq, T, dTsub, hfg, xe = States(P=P, T=T, dTsub=dTsub)
    xe_loc = (liq.h - sliq.h) / hfg
    hfg *= 1.0e3
    cpl = liq.cp * 1.0e3
    Re = Reynolds(G=G, Dh=Dh, mu=liq.mu)
    Pe = G * cpl * Dh / liq.k
    r_rho = sstm.rho / liq.rho
    Pe = G * cpl * Dh / liq.k
    # Forced convection
    hfc = 0.023 * Re**0.8 * liq.Pr**0.4 * liq.k / Dh
    # Iteration
    qw = 1.0e5
    err = 9.999e9
    while err > 100.0:
        # liquid sublayer thickness
        phi0 = 230. * (qw / G / hfg)**0.5
        dTwl = ((phi0 - 1.0) * dTsub + qw/hfc) / phi0
        qfc = hfc * dTwl
        qb = qw - qfc
        delta = (np.pi * 0.0584**2.0 / 2.0) * r_rho**0.4 * (1.0 + r_rho) \
            * (liq.sigma / sstm.rho) * (sstm.rho * hfg / qb)**2.0
        # true quality
        xen = qw / sliq.rho / hfg
        if Pe < 7.0e4:
            xen *= -0.0022 * Dh / sliq.k * cpl * sliq.rho
        else:
            xen *= -154 / G * sliq.rho
        if xe_loc <= xen:
            x_true = 0.0
        else:
            x_true = 1.0 + (xe_loc - 1.0) / (1.0 - xen * np.exp(xe_loc/xen - 1.0))
        # fluid density
        rho = 1.0 / (x_true/sstm.rho + (1.0 - x_true)/liq.rho)
        alpha = x_true / (x_true + (1.0 - x_true) * r_rho)
        mu = alpha * sstm.mu + (1.0 - alpha) * liq.mu * (1.0 + 2.5 * alpha)
        # K factor
        if alpha == 0.0:
            kc = 6.0e4 / (1.0 + 254.8 * r_rho**2.83) / Re**0.8
        elif alpha < 0.25:
            kc = 1.5e4 / (1.0 + 87.2 * r_rho**1.19) / Re**0.8
        elif alpha < 0.7:
            kc = 6.4e3 / (1.0 + 87.3 * r_rho**1.28) / Re**0.8
        else:
            kc = 5.0e3/Re**0.8
        # Velocity Ub
        f = lambda x: 2.0*x - 10.0**(-x) + 2.0 * np.log10(Re) - 0.8
        fp = lambda x: 2.0 + np.log(10.0) / 10.0**x
        f_ = newton_raphson(f, fp, 1.0e-3, 1.0e-4)
        fric = 10.0**(2.0*f_)
        tauw = fric * rho  * (G/rho)**2.0 / 8.0
        ybp = delta * (tauw/rho)**0.5 / (mu/rho)
        if ybp < 30:
            Ubp = 5.5 + 2.5 * np.log(ybp)
        elif ybp > 5:
            Ubp = 5.0 + 5.0 * np.log(ybp / 5.0)
        elif ybp > 0:
            Ubp = ybp
        else:
            raise TypeError("Error")
        Ub = kc * Ubp * tauw/rho
        # Length of vapor blanket
        Lb = 2.8 * np.pi * liq.sigma * (sstm.rho + liq.rho) / (sstm.rho * liq.rho * Ub * Ub)
        # wall heat flux
        qchf = delta * liq.rho * hfg * Ub / Lb
        # Error
        err = abs(qchf - qw)
        qw = (qw + qchf) / 2.0

    # Wall sheare
    

