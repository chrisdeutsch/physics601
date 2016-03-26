import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf
from uncertainties import unumpy
from uncertainties import correlated_values

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

files = ["raumtemperatur",
         "50_grad",
         "62_grad",
         "77_grad",
         "88_grad",
         "99_grad",
         "110_grad",
         "120_grad"]

temp0 = [22.6, 50.0, 63.0, 76.0, 88.7, 100.3, 110.1, 121.3]
temp1 = [22.8, 51.3, 62.6, 76.8, 88.9, 99.8, 110.0, 121.1]

for f, T0, T1 in zip(files, temp0, temp1):
    print("Processing: " + f)
    
    data = pd.read_csv("data/temperatures/{}.txt".format(f), sep="\t",
                       index_col=0, names=["ch", "cnt"])
    data["dcnt"] = np.sqrt(data.cnt)
    
    data = data[700:1600]
    
    plt.errorbar(data.index, data.cnt, yerr=data.dcnt.get_values(), fmt="|", label="Lebenszeitspektrum")
    popt, pcov = curve_fit(lifetime, data.index, data.cnt, p0=[t0, sigma, A0, At, tau0, taut, bg])
    perr = np.sqrt(np.diag(pcov))
    plt.plot(x, lifetime(x, *popt), "-", linewidth=1.2, label="Anpassung")
    
    (A0_, At_, tau0_, taut_) = correlated_values(popt[2:6], pcov[2:6,2:6])
    
    fit_res.append([T0, T1, popt[0], perr[0], popt[1], perr[1], popt[-1], perr[-1], A0_, At_, tau0_, taut_])
    
    plt.legend(loc=0)
    plt.xlim(800, 1600)
    plt.xlabel("Kanal")
    plt.ylabel("Ereignisse~$N$")

    plt.savefig("figures/lifetimes/{}.pdf".format(f))
    plt.close()


### Lifetime table
lifetimes = pd.DataFrame(fit_res, columns=["T0", "T1", "t0", "dt0", "sigma", "dsigma", "bg", "dbg", "A0", "At", "tau0", "taut"])


# Zeitkalibration
m = 0.006721094760572852 # ns / Kanal
dm = 1.7194029982964498e-5

# -> ps / Kanal
m = 1000.0 * m
dm = 1000.0 * dm


# Mittlere Temperaturen
lifetimes["Tbar"] = (lifetimes.T0 + lifetimes.T1) / 2.0
lifetimes["dTbar"] = 0.0 * lifetimes.Tbar + 1.0

# Kelvin
#lifetimes.Tbar = lifetimes.Tbar + 273.15

# Mittlere Lebenszeit (uncertainties zur Fehlerfortpflanzung)
taubar = lifetimes.A0 / (lifetimes.A0 + lifetimes.At) * lifetimes.tau0 + lifetimes.At / (lifetimes.A0 + lifetimes.At) * lifetimes.taut
lifetimes["taubar"] = unumpy.nominal_values(taubar)
lifetimes["dtaubar"] = unumpy.std_devs(taubar)


# Uncertainties aus pandas frame entfernen
A0 = lifetimes.A0.get_values()
lifetimes["A0"], lifetimes["dA0"] = unumpy.nominal_values(A0), unumpy.std_devs(A0)

At = lifetimes.At.get_values()
lifetimes["At"], lifetimes["dAt"] = unumpy.nominal_values(At), unumpy.std_devs(At)

tau0 = lifetimes.tau0.get_values()
lifetimes["tau0"], lifetimes["dtau0"] = unumpy.nominal_values(tau0), unumpy.std_devs(tau0)

taut = lifetimes.taut.get_values()
lifetimes["taut"], lifetimes["dtaut"] = unumpy.nominal_values(taut), unumpy.std_devs(taut)


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
plt.errorbar(picosec.Tbar + 273.15, picosec.taubar, xerr=picosec.dTbar.get_values(), yerr=picosec.dtaubar.get_values(), fmt="o", label="mittlere Lebenszeit (Anpassung)")

def s_curve(T, tauf, taut, a, b):
    fac = a * np.exp(- b / (T + 273.15))
    return tauf * (1.0 + fac * taut) / (1.0 + fac * tauf)

mask = (picosec.Tbar < 60.0) | (picosec.Tbar > 70.0)
popt, pcov = curve_fit(s_curve, picosec[mask].Tbar + 273.15, picosec[mask].taubar,
                       p0=[322., 465., 186., 7000.],
                       sigma=picosec[mask].dtaubar, absolute_sigma=True)

pcov[2:4, 2:4] = pcov[2:4, 2:4] * 0.0001

x = np.linspace(280.0, 400.0, 1000)
plt.plot(x, s_curve(x, *popt), "-", label="Anpassung S-Kurve")

plt.xlabel(r"Temperatur~$T$ / \si{\kelvin}")
plt.ylabel(r"mittlere Lebenszeit~$\bar{\tau}$ / \si{ps}")
plt.legend(loc="lower right")


plt.savefig("figures/lifetime_s_curve.pdf")
plt.close()

# Arrenius plot
def log_sigc(T):
    global pcov
    ab = popt[2:4]
    ab_cov = pcov[2:4, 2:4]
    
    (a, b) = correlated_values(ab, ab_cov)
    
    retval = unumpy.log(a) - b / T
    return (unumpy.nominal_values(retval), unumpy.std_devs(retval))


x = 1.0 / (picosec.Tbar + 273.15)
dx = picosec.dTbar / (picosec.Tbar + 273.15)**2
y, dy = log_sigc(picosec.Tbar + 273.15)


plt.errorbar(x, y, xerr=dx, yerr=dy, fmt="o")

popt, pcov = curve_fit(lambda x, m, b: m * x + b,
                       x, y, sigma=dy, absolute_sigma=True)

plt.plot(x, popt[0] * x + popt[1], "-")


plt.xlabel(r"$\frac{1}{T}$ / \si{\per\kelvin}")
plt.ylabel(r"$\ln(\sigma c_t)$")



plt.savefig("figures/arrenius.pdf")
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
               r"{$\Delta \bar{\tau}$ / \si{ps}}"]

out.to_latex("tables/mean_lifetimes.tex", index = False,
             formatters=[round(1), round(1), round(1), round(1), round(1), round(1), round(1)],
             column_format="S[table-format=2.1]S[table-format=3.1]S[table-format=2.1]S[table-format=3.1]S[table-format=2.1]S[table-format=3.1]S[table-format=2.1]", escape=False)
