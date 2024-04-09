import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from iapws import IAPWS97

# Set nan value
nan = -9999



"""Plot function
"""
def plot_scatter(title=None, data=None, fig_size=None, cline=None):
    # Create output directory
    dir = "plots"
    if not os.path.exists(dir):
        os.mkdir(dir)
    path = os.path.join(dir, title)

    # Plot
    fsiz = fig_size or (5,5)
    fig, ax = plt.subplots(figsize=fsiz)
    labels, range = None, None
    for key in data.keys():
        if key == "labels":
            labels = data[key]
        elif key == "range":
            range = data[key]
        else:
            tmp = data[key]
            ax.plot(tmp["x"], tmp["y"], tmp["ls"], markerfacecolor='None')
    if range != None:
        if cline:     
            ax.plot(range["x"], range["y"])
        ax.set_xlim(range["x"])
        ax.set_ylim(range["y"])
    if labels != None:
        ax.set_xlabel(labels["x"])
        ax.set_ylabel(labels["y"])
    plt.savefig(path, dpi=300)
    plt.show()


def Read_CSV(file=None):
    """ Read CSV file
    """
    df = pd.read_csv(file, header=[0], skiprows=[1])
    df = df.fillna(nan)
    #print(df.keys())
    return df


def Tong(P=None, G=None, Tl=None, dTl=None, dhl=None, xe=None, 
             Dh=None, zchf=None, fq=None, modified=None):
    """ Tong correlation
        P     : Pressure (MPa)
        G     : Mass flow rate (kg/m2s)
        Tl    : Temperature (K)
        dTl   : Subcooling (K)
        xe    : equilibrium quality (-)
        Dh    : Hydrodynamic diameter (m)
        Dr    : Rod diameter (m)
        Dc    : Shroud diameter (m)
        Wc    : Shroud width (m)
        s     : Spacing (m)
    """
    # Fluid properties
    P /= 1.0e6 
    liq_sat = IAPWS97(P=P, x=0)
    vap_sat = IAPWS97(P=P, x=1)
    # Latent heat
    hfg = vap_sat.h - liq_sat.h
    # Determine fluid at location of CHF
    if Tl == nan:
        if dTl == nan:
            if xe == nan:
                raise TypeError("Missing input parameters!")
            else:
                hl = liq_sat.h + xe * hfg
                liq = IAPWS97(P=P, h=hl)
        else:
            Tl = liq_sat.T + dTl
            liq = IAPWS97(P=P, T=Tl)
    else:
        liq = IAPWS97(P=P, T=Tl)
    # Equivalent quality (-)
    if xe == nan:
        qual = (liq.h - liq_sat.h) / hfg
    else:
        qual = xe
    # Dynamic viscousity (Pa.s)
    muf = liq.Liquid.mu
    # Reynolds number
    Re = G * Dh / muf
    # Considering modification
    if modified==0 or modified==None:
        C = 1.76 - 7.433*qual + 12.222 * qual**2.
        n = 0.6
    elif modified=="celata":
        C = 0.27 + 5.93e-2 * P
        if qual > -0.1:
            C = C * (0.825 + 0.986 * qual)
        n = 0.5

    # Non uniform
    cf = 185.6 * (1 - xe)**4.31 / G**0.478
    z = np.linspace(0,zchf,200)
    fz = (0.82 + 0.68 * np.cos(2.44 * z - 2.44))
    gz = fz * np.exp(cf*(z - zchf))
    Iz = np.trapz(gz, z)
    Fz = cf * Iz / (fz[-1] * (1 - np.exp(-cf * zchf)))

    # Critical heat flux (W/m2) 
    res = C / Re**n * G * hfg * 1.0e3 #/ Fz
    return res