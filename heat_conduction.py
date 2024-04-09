import numpy as np

L = 0.005
qw = 1.0e6
Ti = 100
Tinf = 0

rho = 8000
k = 16.2
cp = 500
alpha =k / (rho * cp)

tmp = qw * L / k / (Ti - Tinf)

t = np.arange(0,1,0.001)
tau = t * alpha/L/L

for ti in t:
    wn, err = 0.0, 1.0
    while err > 0.0001:
        wn += 0.00001
        An = 4 * np.sin(wn) / (2 * wn + np.sin(2 * wn))
        err = An / np.exp(wn * wn) - tmp / wn / np.tan(wn)
    


