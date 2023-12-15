from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from scipy import special
import os

from iapws import IAPWS97

"""Plot function
"""
def plot_scatter(title=None, data=None, range=None, cline=None):
    # Create output directory
    dir = "plots"
    if not os.path.exists(dir):
        os.mkdir(dir)
    path = os.path.join(dir, title)
    # Generate xy line
    x0, y0, x1, y1 = range
    xy = np.linspace(x0,x1)
    # Plot
    fig, ax = plt.subplots(figsize=(5,5))
    for key in data.keys():
        tmp = data[key]
        ax.scatter(x=tmp["x"], y=tmp["y"], marker=tmp["mark"], label=key)
    if cline:        
        ax.plot(xy,xy)
    else:
        tmp = np.ones(xy.size)
        ax.plot(xy,tmp)
    ax.set_xlim([x0, x1])
    ax.set_ylim([y0, y1])
    ax.set_xlabel("$q^{''}_{chf, exp} (MW/m^2)$")
    ax.set_ylabel("$q^{''}_{chf, Tong} (MW/m^2)$")

    #plt.savefig(path, dpi=300)
    plt.show()



""" Tong correlation
    P     : Pressure (MPa)
    T     : Temperature (K)
    dTsub : Subcooling (K)
    G     : Mass flow rate (kg/m2s)
"""
def Tong(P=None, G=None, T=None, dTsub=None, D=None, modified=None):
    # Fluid properties
    sat_liq = IAPWS97(P=P, x=0)
    sat_stm = IAPWS97(P=P, x=1)
    if not T:
        T = sat_liq.T - dTsub
    water = IAPWS97(P=P, T=T)
    # Latent heat (kJ/kg)
    hfg = sat_stm.h - sat_liq.h
    # Equivalent quality (-)
    xe = (water.h - sat_liq.h) / hfg
    # Dynamic viscousity (Pa.s)
    muf = water.mu
    # Reynolds number
    Re = G * D / muf
    # Considering modification
    if modified==0 or modified==None:
        C = 1.76 - 7.433*xe + 12.222 * xe**2.
        n = 0.6
    elif modified=="celata":
        C = 0.27 + 5.93e-2 * P
        if xe > -0.1:
            C = C * (0.825 + 0.986 * xe)
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
def W3(P=None, T=None, dTsub=None, G=None, D=None, Pi=None, Ti=None, LDh = None):
    # Fluid properties
    sat_liq = IAPWS97(P=P, x=0)
    sat_stm = IAPWS97(P=P, x=1)
    if not T:
        T = sat_liq.T - dTsub
    if not dTsub:
        dTsub = sat_liq.T - T
    hfg = sat_stm.h - sat_liq.h
    loc_liq = IAPWS97(P=P, T=T)
    xe = (loc_liq.h - sat_liq.h) / hfg
    inlet_liq = IAPWS97(P=Pi, T=Ti)
    # CHF
    a = (2.022 - 0.06238 * P) + (0.1722 - 0.01427 * P) * np.exp((18.177 - 0.5987*P)*xe)
    b = (0.1484 - 1.596*xe + 0.1729*xe*abs(xe)) * 2.326*G + 3271
    c = 1.157 - 0.869 * xe
    d = 0.2664 + 0.8357 * np.exp(-124.1*D)
    e = 0.8258 + 3.413e-4 * (sat_liq.h - inlet_liq.h)
    res = a*b*c*d*e / 1e3
    return res

""" Weisman & Pei
    P (local pressure) : MPa
    T (local temperature) : K
    G (total axial mass velocity) : kg/m2h
"""
def WeismanPei(P=None, T=None, G=None, Dh=None):
    # Fluid properties
    sat_liq = IAPWS97(P=P, x=0)
    sat_stm = IAPWS97(P=P, x=1)
    loc_liq = IAPWS97(P=P, T=T)
    rhof, muf = sat_liq.rho, sat_liq.mu
    rhog, sigma = sat_stm.rho, sat_stm.sigma    
    Hl, Cpl, kl, mul = loc_liq.h * 1.0e3, loc_liq.cp0 * 1.0e3, loc_liq.k, loc_liq.mu
    g, edh = 9.81, 1.0e-4
    # Reynold number
    Re = G * Dh / muf
    # Prandtl number
    Pr = Cpl * mul/ kl
    # Iteration
    qchf, err = 1.0, 999999.0
    while err > 1.0:
        # Friction factor
        ffric = 0.0055 * (1 + 20000 * edh + 1e6 / Re)**(1./3.)
        # Calculation of hld
        ybp = 0.1 * (sigma * g * Dh * rhof)**0.5 / muf
        tmp = Cpl * qchf / Hl
        tmp1 = qchf / G / (ffric/8.0)**0.5
        if ybp <= 5.0:
            hfhld = tmp - tmp1 * Pr * ybp
        elif ybp <= 30.0:
            tmp2 = tmp1 * (Pr + np.log(1.0 + Pr * (ybp/5.0 - 1.0)))
            hfhld = tmp - 5.0 * tmp2
        else:
            tmp3 = tmp1 * (Pr + np.log(1.0 + 5.0 * Pr) + 0.5 * np.log(ybp/30.0))
            hfhld = tmp - 5.0 * tmp3
    
    # Bubble diameter
    
    
    tauw = ffric / 8. / loc_liq.rho * G**2.
    Dp = 0.015 * (loc_liq.sigma * Dh / tauw)**0.5
    # Turbulent intensity
    alpha = 0.135 / (G/9.7e6)**0.3
    drho = (loc_liq.rho - sat_stm.rho) / sat_stm.rho
    ib = (Dp/Dh)**0.6 / Re**0.1 * (1. + alpha * drho)
    # Psi depends on qchf so that iteration is needed
    qchf, err = 1.0, 1.0
    while err < 1.0:
        v11 = qchf / sat_stm.rho / loc_liq.h / 1.0e3
        sigvp = G / loc_liq.rho * ib
        tmp = v11 / sigvp
        psi = 0.5 / np.pi * np.exp(-0.5 * tmp**2) - 0.5 *tmp*special.erfc(tmp / 2**0.5)
        xe = (loc_liq.h - sat_stm.h) / (sat_liq.h - sat_stm.h)


# Read data
celata1992 = pd.read_csv("celata1992.csv", header=0)

# Plot u=30, p = {0.8, 1.4, 2.5}
df = celata1992
df1 = df.loc[abs(df["u_out_m/s"] - 30) < 1.5][["P_exit_Mpa","Tsub_in_K", "qchf_MW/m2"]]
print(df1)

df08 = df1.loc[abs(df["P_exit_Mpa"] - 0.8) < 0.2][["Tsub_in_K", "qchf_MW/m2"]]
df14 = df1.loc[abs(df["P_exit_Mpa"] - 1.4) < 0.2][["Tsub_in_K", "qchf_MW/m2"]]
df25 = df1.loc[abs(df["P_exit_Mpa"] - 2.5) < 0.2][["Tsub_in_K", "qchf_MW/m2"]]

tmp = {
    "p = 0.8 MPa" : {"x" : df08["Tsub_in_K"], "y" : df08["qchf_MW/m2"], "mark" : "s"},
    "p = 1.4 MPa" : {"x" : df14["Tsub_in_K"], "y" : df14["qchf_MW/m2"], "mark" : "o"},
    "p = 2.5 MPa" : {"x" : df25["Tsub_in_K"], "y" : df25["qchf_MW/m2"], "mark" : "^"}
}
#plot_scatter(title="fig1", data=tmp)

# CHF correlations
df["qchf_Tong"] = None
df["qchf_WH"] = None

for i in df.index:
    P = df["P_exit_Mpa"][i]
    T = df["To_c"][i] + 273.15
    G = df["G_kg/m2s"][i]
    D = df["D_mm"][i] / 1000.
    dTsub = df["Tsub_out_K"][i]
    Pi = df["P_in_Mpa"][i]
    Ti = df["Ti_c"][i] + 273.15
    df["qchf_Tong"][i] = Tong(P=P, G=G, T=T, D=D, modified="celata")
    df["qchf_WH"][i] = W3(P=P, T=T, dTsub=dTsub,G=G, D=D, Pi=Pi, Ti=Ti, LDh=0.1/D)
print(df)

# Tong correlation
tmp = {
    "Tong" : {"x" : df["qchf_MW/m2"], "y" : df["qchf_Tong"], "mark" : "s"}
}
plot_scatter(title="Tong", data=tmp, range=[10,10,70,70], cline=True)

# Westinghouse correlation
tmp = {
    "WH" : {"x" : df["qchf_MW/m2"], "y" : df["qchf_WH"], "mark" : "s"}
}
plot_scatter(title="Westinghouse", data=tmp, range=[10,10,70,70], cline=True)

tmp = {
    "WH" : {"x" : df["P_exit_Mpa"], "y" : df["qchf_WH"]/df["qchf_MW/m2"], "mark" : "s"}
}
plot_scatter(title="Westinghouse", data=tmp, range=[0,0,3.0,10.0], cline=False)