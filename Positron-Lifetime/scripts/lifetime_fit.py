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

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([22.6, 22.8, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)
plt.xlim(800, 1800)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/lifetimes/raumtemperatur.pdf")
plt.close()


### Lifetime 50°C
data = pd.read_csv("data/temperatures/50_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([50.0, 51.3, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)
plt.xlim(800, 1800)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/lifetimes/50.pdf")
plt.close()


### Lifetime 62°C
data = pd.read_csv("data/temperatures/62_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([63.0, 62.6, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)
plt.xlim(800, 1800)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/lifetimes/62.pdf")
plt.close()


### Lifetime 77°C
data = pd.read_csv("data/temperatures/77_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([76.0, 76.8, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)
plt.xlim(800, 1800)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/lifetimes/77.pdf")
plt.close()


### Lifetime 88°C
data = pd.read_csv("data/temperatures/88_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([88.7, 88.9, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)
plt.xlim(800, 1800)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/lifetimes/88.pdf")
plt.close()


### Lifetime 99°C
data = pd.read_csv("data/temperatures/99_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([100.3, 99.8, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)
plt.xlim(800, 1800)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/lifetimes/99.pdf")
plt.close()


### Lifetime 110°C
data = pd.read_csv("data/temperatures/110_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([110.1, 110.0, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)
plt.xlim(800, 1800)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/lifetimes/110.pdf")
plt.close()


### Lifetime 120°C
data = pd.read_csv("data/temperatures/120_grad.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)
data = data[700:2000]

plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
perr = np.sqrt(np.diag(pcov))
plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")

fit_res.append([121.3, 121.1, *popt.tolist(), *perr.tolist()])

plt.legend(loc=0)
plt.xlim(800, 1800)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/lifetimes/120.pdf")
plt.close()


t0, sigma, A0, At, tau0, taut, bg

### Lifetime table
lifetimes = pd.DataFrame(fit_res, columns=["T0", "T1", "t0", "sigma", "A0", "At", "tau0", "taut", "bg", "dt0", "dsigma", "dA0", "dAt", "dtau0", "dtaut", "dbg"])

# Zeitkalibration
m = 0.006721094760572852 # ns / Kanal
dm = 1.7194029982964498e-5

# -> ps / Kanal
m = 1000.0 * m
dm = 1000.0 * dm

# Mittlere Temperaturen
lifetimes["Tbar"] = (lifetimes.T0 + lifetimes.T1) / 2.0
lifetimes["dTbar"] = 0.0 * lifetimes.Tbar + 1.0

# Intensitäten
lifetimes["I0"] = lifetimes.A0 / (lifetimes.A0 + lifetimes.At)
lifetimes["It"] = lifetimes.At / (lifetimes.A0 + lifetimes.At)

# Mittlere Lebenszeit
lifetimes["taubar"] = lifetimes.I0 * lifetimes.tau0 + lifetimes.It * lifetimes.taut
lifetimes["dtaubar"] = lifetimes.dtau0

### Umrechnung von Kanal -> ps
picosec = pd.DataFrame()
picosec["Tbar"] = lifetimes.Tbar
picosec["dTbar"] = lifetimes.dTbar
picosec["tau0"] = lifetimes.tau0 * m
picosec["dtau0"] = np.sqrt(lifetimes.tau0**2 * dm**2 + m**2 * lifetimes.dtau0**2)
picosec["taut"] = lifetimes.taut * m
picosec["dtaut"] = np.sqrt(lifetimes.taut**2 * dm**2 + m**2 * lifetimes.dtaut**2)
picosec["taubar"] = lifetimes.taubar * m
picosec["dtaubar"] = np.sqrt(lifetimes.taubar**2 * dm**2 + m**2 * lifetimes.dtaubar**2)


### Mean lifetime plot
plt.errorbar(picosec.Tbar, picosec.taubar, xerr=picosec.dTbar.get_values(), yerr=picosec.dtaubar.get_values(), fmt="o", label="mittlere Lebenszeit (Anpassung)")

def s_curve(T, tauf, taut, sigma, S, H):
    kB = 1.0
    a = sigma * np.exp(S / kB - H / (kB * T))
    return tauf * (1.0 + a * taut) / (1.0 + a * tauf)

popt, pcov = curve_fit(s_curve, picosec.Tbar, picosec.taubar, p0=[4.5e+2, 3.8e+2, 2.9e-1, -7.6, -190.0])

x = np.linspace(20.0, 140.0, 1000)
plt.plot(x, s_curve(x, *popt), "-", label="Anpassung S-Kurve")

plt.xlabel(r"Temperatur~$T$ / \si{\degreeCelsius}")
plt.ylabel(r"mittlere Lebenszeit~$\bar{\tau}$ / \si{ps}")
plt.legend(loc="lower right")


plt.savefig("figures/lifetime_s_curve.pdf")
plt.close()


### Latex Tabellen
from scripts.tools import round

out = lifetimes[["Tbar", "t0", "dt0", "sigma", "dsigma", "A0", "dA0", "At", "dAt", "tau0", "dtau0", "taut", "dtaut"]]

out.columns = [r"{$\bar{T}$ / \si{\degreeCelsius}}",
               r"{$t_0$ / Kanal}", 
               r"{$\Delta t_0$ / Kanal}",
               r"{$\sigma$ / Kanal}",
               r"{$\Delta \sigma$ / Kanal}",
               r"{$A_0$}",
               r"{$\Delta A_0$}",
               r"{$A_t$}",
               r"{$\Delta A_t$}",
               r"{$\tau_0$ / Kanal}",
               r"{$\Delta \tau_0$ / Kanal}",
               r"{$\tau_t$ / Kanal}",
               r"{$\Delta \tau_t$ / Kanal}"]

out.to_latex("tables/lifetime_fit.tex", index = False,
             formatters=[round(1), round(2), round(2), round(2), round(2),
                         round(0), round(0), round(0), round(0), round(1),
                         round(1), round(1), round(1), ],
             column_format="S[table-format=2.1]S[table-format=3.2]S[table-format=1.2]S[table-format=2.2]S[table-format=1.2]S[table-format=5.0]S[table-format=4.0]S[table-format=5.0]S[table-format=4.0]S[table-format=2.1]S[table-format=1.1]S[table-format=2.1]S[table-format=1.1]", escape=False)

out = picosec[["Tbar", "tau0", "dtau0", "taut", "dtaut", "taubar", "dtaubar"]]

out.columns = [r"{$\bar{T}$ / \si{\degreeCelsius}}",
               r"{$\tau_0$ / \si{ps}}",
               r"{$\Delta \tau_0$ / \si{ps}}",
               r"{$\tau_t$ / \si{ps}}",
               r"{$\Delta \tau_t$ / \si{ps}}",
               r"{$\bar{\tau}$ / \si{ps}}",
               r"{$\Delta \bar{\tau_t}$ / \si{ps}}"]

out.to_latex("tables/mean_lifetimes.tex", index = False,
             formatters=[round(1), round(1), round(1), round(1), round(1), round(1), round(1)],
             column_format="S[table-format=2.1]S[table-format=3.1]S[table-format=2.1]S[table-format=3.1]S[table-format=2.1]S[table-format=3.1]S[table-format=2.1]", escape=False)
