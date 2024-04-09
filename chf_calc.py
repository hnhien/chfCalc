import chf_predictor as chf
import numpy as np

# Read data
file = 'chf_data.csv'
df = chf.Read_CSV(file=file)

# Greenwood et al. data
# Pressure measured at inlet & assumed constant along test section
dfs  = df[df["YEAR"]==2017]
print(dfs.describe())

# Calculate chf
shapes = {"C" : "circle", "S" : "square"}
dfs["CHF_TONG"] = None
dfs["K"], dfs["Re"] = None, None
for i in dfs.index:
    s = dfs["GEOM"][i]
    shape = shapes[s[0]]
    nrods = int(s[1:])
    if shape == "circle":
        if nrods == 0:
            dh = dfs["SHROUD_DIM"][i]
            fa = 1.0
        elif nrods == 1:
            dh = dfs["SHROUD_DIM"][i] - dfs["ROD_DIM"][i]
            fa = 1.0
    elif shape == "square":
        dh = dfs["ROD_DIM"][i]
        fa = (dfs["ROD_PITCH"][i]**2.0 - np.pi * dfs["ROD_DIM"][i]**2.0 / 4) / (dfs["SHROUD_DIM"][i]**2.0)
    G = df["G"][i] * fa
    P_e = dfs["P_IN"][i]
    T_e = df["T_EXIT"][i]
    dT_e = dfs["TSUB_EXIT"][i]
    xe = dfs["XE_EXIT"][i]
    zchf = dfs["Z_CHF"][i]
    # calcualte chf
    dfs["CHF_TONG"][i] = chf.Tong(P=P_e, G=G, Tl=T_e, dTl=dT_e, xe=xe, Dh=dh, 
                                  zchf=zchf, modified="celata")
print(dfs.describe())


tmp = {
    "labels" : {"x" : "$q^{''}_{chf, exp} (MW/m^2s)$",
               "y" : "$q^{''}_{chf, Tong} (MW/m^2s)$"},
    "range" : {"x" : [0, 2], "y" : [0, 2]},
    "center" : {"x" : [0, 2], "y" : [0, 2], "ls" : "--"},
    "Tong" : {"x" : dfs["CHF"] / 1e6, "y" : dfs["CHF_TONG"] / 1e6, "ls" : "rs"}
}
chf.plot_scatter(title="Tong_Greenwood_uniform", data=tmp, fig_size=(6,6))


tmp = {
    "labels" : {"x" : "$X_e (-)$",
               "y" : "$q^{''}_{chf, Tong} / q^{''}_{chf, exp} (-)$"},
    "range" : {"x" : [-0.5, 0.5], "y" : [0, 2]},
    "center" : {"x" : [-0.5, 0.5], "y" : [1, 1], "ls" : "--"},
    "Tong" : {"x" : dfs["XE_EXIT"], "y" : dfs["CHF_TONG"] / dfs["CHF"], "ls" : "rs"}
}
chf.plot_scatter(title="Tong_Greenwood_uniform_xe", data=tmp, fig_size=(9,6))