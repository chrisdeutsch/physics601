import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import erf
from uncertainties import unumpy, correlated_values, ufloat
from lmfit import Model

plt.style.use("publication")

### Fit Funktion
def lifetime(t, t0, sigma, A0, At, tau0, taut, bg):
    a0 = (sigma**2 + tau0 * t0) / (np.sqrt(2.0) * sigma * tau0)
    b0 = (tau0 * (t - t0) - sigma**2) / (np.sqrt(2.0) * sigma * tau0)
    at = (sigma**2 + taut * t0) / (np.sqrt(2.0) * sigma * taut)
    bt = (taut * (t - t0) - sigma**2) / (np.sqrt(2.0) * sigma * taut)
    
    return A0 / (2.0 * tau0) * np.exp((sigma**2 - 2.0 * tau0 * (t- t0)) / (2.0 * tau0**2)) * (erf(a0) + erf(b0)) + At / (2.0 * taut) * np.exp((sigma**2 - 2.0 * taut * (t- t0)) / (2.0 * taut**2)) * (erf(at) + erf(bt)) + bg

lifetime_model = Model(lifetime)

### Fit Startwerte
t0 = 975
sigma = 37.8
A0 = 55680.
At = 40000.
tau0 = 24.7
taut = 108.0
bg = 0.0

# Zeitkalibration
m = 0.006721094760572852 # ns / Kanal
dm = 1.7194029982964498e-5


data = pd.read_csv("data/acrylglass/longtermmeasurement.txt", sep="\t",
                   index_col=0, names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)

data = data[800:1800]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt=".", markersize=0.2,
             label="Lebenszeitspektrum")

result = lifetime_model.fit(data.cnt.get_values(),
                            weights=data.dcnt.get_values()**-1,
                            t=data.index.get_values(),
                            t0 = 975, sigma = 37.8, A0 = 55680.,At = 40000.,
                            tau0 = 24.7, taut = 108.0, bg = 0.0)

plt.plot(data.index.get_values(), result.best_fit, "-", label="Anpassung")


plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")
plt.legend(loc=0)

plt.savefig("figures/acryl_glas.pdf")
plt.close()

# Umrechnung in ps
vals = result.params

tau0 = vals.get("tau0").value
dtau0 = vals.get("tau0").stderr
tau1 = vals.get("taut").value
dtau1 = vals.get("taut").stderr

tau0_ps = tau0 * m
dtau0_ps = np.sqrt(tau0**2 * dm**2 + m**2 * dtau0**2)
tau1_ps = tau1 * m
dtau1_ps = np.sqrt(tau1**2 * dm**2 + m**2 * dtau1**2)

print(result.fit_report())
print("tau0: {} +- {}".format(tau0_ps, dtau0_ps))
print("tau1: {} +- {}".format(tau1_ps, dtau1_ps))
