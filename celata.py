from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os

from iapws import IAPWS97
from CHFModels import *

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


# Read data
celata1992 = pd.read_csv("celata1992.csv", header=0)
df = celata1992

#=======================================================================================#
# CHF correlations
df["qchf_Tong"] = None
df["qchf_WH"] = None
df["xe"] = None

for i in df.index:
    P = df["P_exit_Mpa"][i]
    T = df["To_c"][i] + 273.15
    G = df["G_kg/m2s"][i]
    D = df["D_mm"][i] / 1000.
    dTsub = df["Tsub_out_K"][i]
    Pi = df["P_in_Mpa"][i]
    Ti = df["Ti_c"][i] + 273.15
    # calcualte xe
    sat_liq = IAPWS97(P=P, x=0)
    sat_stm = IAPWS97(P=P, x=1)
    water = IAPWS97(P=P, T=T)
    df["xe"][i] = (water.h - sat_liq.h) / (sat_stm.h - sat_liq.h)
    # calcualte chf
    df["qchf_Tong"][i] = Tong(P=P, G=G, T=T, D=D, modified="celata")
    df["qchf_WH"][i] = W3(P=P, T=T, dTsub=dTsub,G=G, D=D, Pi=Pi, Ti=Ti, LDh=0.1/D)
    #df["qchf_Weis"][i] = WeisIles(P=P, T=T, dTsub=dTsub, G=G, Dh=D)
    #df["qchf_Katto"][i] = Katto90(P=P, T=T, dTsub=dTsub, G=G, Dh=D)
print(df)

#=======================================================================================#
# Tong correlation
tmp = {
    "Tong" : {"x" : df["qchf_MW/m2"], "y" : df["qchf_Tong"], "mark" : "s"}
}
plot_scatter(title="Tong", data=tmp, range=[10,10,70,70], cline=True)

#=======================================================================================#
# Westinghouse correlation
tmp = {
    "WH" : {"x" : df["qchf_MW/m2"], "y" : df["qchf_WH"], "mark" : "s"}
}
plot_scatter(title="Westinghouse", data=tmp, range=[10,10,70,70], cline=True)

tmp = {
    "WH" : {"x" : df["P_exit_Mpa"], "y" : df["qchf_WH"]/df["qchf_MW/m2"], "mark" : "s"}
}
plot_scatter(title="Westinghouse", data=tmp, range=[0,0,3.0,10.0], cline=False)

#=======================================================================================#
# Plot u=30, p = {0.8, 1.4, 2.5}
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

#=======================================================================================#
x = df["xe"]
y = df["G_kg/m2s"]
z = df['qchf_MW/m2']

ax = plt.axes(projection='3d')
ax.scatter(x,y,z, c=z, cmap='inferno', s=5, alpha=0.5)
plt.show()