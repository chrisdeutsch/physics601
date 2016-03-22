import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf

### Fit Funktion
def lifetime(t, t0, sigma, A0, At, tau0, taut, bg):
    a0 = (sigma**2 + tau0 * t0) / (np.sqrt(2.0) * sigma * tau0)
    b0 = (tau0 * (t - t0) - sigma**2) / (np.sqrt(2.0) * sigma * tau0)
    at = (sigma**2 + taut * t0) / (np.sqrt(2.0) * sigma * taut)
    bt = (taut * (t - t0) - sigma**2) / (np.sqrt(2.0) * sigma * taut)
    
    print(a0, b0, at, bt)
    
    return A0 / (2.0 * tau0) * np.exp((sigma**2 - 2.0 * tau0 * (t- t0)) / (2.0 * tau0**2)) * (erf(a0) - erf(b0)) + At / (2.0 * taut) * np.exp((sigma**2 - 2.0 * taut * (t- t0)) / (2.0 * taut**2)) * (erf(at) - erf(bt)) + bg


data = pd.read_csv("data/temperatures/50_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)

t0 = 1000.0
sigma = 50.0
A0 = 10000
At = 1000
tau0 = 100e-9
taut = 100e-9
bg = 0.0


data = data[(data.index > 740) & (data.index < 2000)]
plt.plot(data.index, data.cnt, "-")
plt.plot(data.index, lifetime(data.index, t0, sigma, A0, At, tau0, taut, bg))
print(lifetime(1000.0, t0, sigma, A0, At, tau0, taut, bg))

plt.show()