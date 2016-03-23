import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf

plt.style.use("publication")

### Fit Funktion
def lifetime(t, t0, sigma, A0, At, tau0, taut, bg):
    a0 = (sigma**2 + tau0 * t0) / (np.sqrt(2.0) * sigma * tau0)
    b0 = (tau0 * (t - t0) - sigma**2) / (np.sqrt(2.0) * sigma * tau0)
    at = (sigma**2 + taut * t0) / (np.sqrt(2.0) * sigma * taut)
    bt = (taut * (t - t0) - sigma**2) / (np.sqrt(2.0) * sigma * taut)
    
    return A0 / (2.0 * tau0) * np.exp((sigma**2 - 2.0 * tau0 * (t- t0)) / (2.0 * tau0**2)) * (erf(a0) + erf(b0)) + At / (2.0 * taut) * np.exp((sigma**2 - 2.0 * taut * (t- t0)) / (2.0 * taut**2)) * (erf(at) + erf(bt)) + bg


### Fit Startwerte
t0 = 975
sigma = 37.8
A0 = 55680.
At = 40000.
tau0 = 24.7
taut = 108.0
bg = 0.0

# Plotrange
x = np.linspace(700.0, 2000.0, 1000)

# Fit results
fit_res = []


### Lifetime Raumtemperatur
data = pd.read_csv("data/temperatures/raumtemperatur.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="+", label="Spektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([22.6, 22.8, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)

plt.savefig("figures/lifetimes/raumtemperatur.pdf")
plt.close()


### Lifetime 50°C
data = pd.read_csv("data/temperatures/50_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="+", label="Spektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([50.0, 51.3, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)

plt.savefig("figures/lifetimes/50.pdf")
plt.close()


### Lifetime 62°C
data = pd.read_csv("data/temperatures/62_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="+", label="Spektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([63.0, 62.6, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)

plt.savefig("figures/lifetimes/62.pdf")
plt.close()


### Lifetime 77°C
data = pd.read_csv("data/temperatures/77_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="+", label="Spektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([76.0, 76.8, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)

plt.savefig("figures/lifetimes/77.pdf")
plt.close()


### Lifetime 88°C
data = pd.read_csv("data/temperatures/88_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="+", label="Spektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([88.7, 88.9, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)

plt.savefig("figures/lifetimes/88.pdf")
plt.close()


### Lifetime 99°C
data = pd.read_csv("data/temperatures/99_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="+", label="Spektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([100.3, 99.8, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)

plt.savefig("figures/lifetimes/99.pdf")
plt.close()


### Lifetime 110°C
data = pd.read_csv("data/temperatures/110_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="+", label="Spektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([110.1, 110.0, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)

plt.savefig("figures/lifetimes/110.pdf")
plt.close()


### Lifetime 120°C
data = pd.read_csv("data/temperatures/120_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="+", label="Spektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([121.3, 121.1, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)

plt.savefig("figures/lifetimes/120.pdf")
plt.close()


t0, sigma, A0, At, tau0, taut, bg

### Lifetime table
lifetimes = pd.DataFrame(fit_res, columns=["T0", "T1", "t0", "sigma", "A0", "At", "tau0", "taut", "bg", "dt0", "dsigma", "dA0", "dAt", "dtau0", "dtaut", "dbg"])

# Intensitäten
lifetimes["I0"] = lifetimes.A0 / (lifetimes.A0 + lifetimes.At)
lifetimes["It"] = lifetimes.At / (lifetimes.A0 + lifetimes.At)

# Mittlere Lebenszeit
lifetimes["taubar"] = lifetimes.I0 * lifetimes.tau0 + lifetimes.It * lifetimes.taut
lifetimes["dtaubar"] = lifetimes.dtau0

### Mean lifetime plot
plt.errorbar(lifetimes.T0, lifetimes.taubar, yerr=lifetimes.dtaubar.get_values(), fmt="o")

def s_curve(T, tauf, taut, sigma, S, H):
    kB = 1.0
    a = sigma * np.exp(S / kB - H / (kB * T))
    return tauf * (1.0 + a * taut) / (1.0 + a * tauf)

popt, pcov = curve_fit(s_curve, lifetimes.T0, lifetimes.taubar)

x = np.linspace(20.0, 130.0, 1000)
plt.plot(x, s_curve(x, *popt), "-")

plt.savefig("figures/lifetime_s_curve.pdf")

plt.close()
