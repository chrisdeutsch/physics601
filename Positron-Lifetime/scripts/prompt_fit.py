import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.style.use("publication")

data = pd.read_csv("data/prompt/prompt_curves_4ns.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)

data = data[(data.index > 500) & (data.index < 5500)]

plt.xlim(500, 5500)
#plt.plot(data.index, data.cnt)

t = []
mu = []
dmu = []

### Fit
def gauss(x, mu, sigma, norm, bg):
    return norm * np.exp(-0.5 * np.power((x - mu) / sigma, 2)) / np.sqrt(2.0 * np.pi) / sigma + bg

# Fit 1
fitr = data[(data.index > 850) & (data.index < 1050)]
popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[950.0, 50.0, 100000.0, 0.0], sigma=fitr.dcnt, absolute_sigma=True)

t.append(20.0)
mu.append(popt[0])
dmu.append(np.sqrt(np.diag(pcov))[0])

plt.plot(fitr.index, gauss(fitr.index, *popt), "-")


# Fit 2
fitr = data[(data.index > 1425) & (data.index < 1625)]
popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[1525.0, 50.0, 100000.0, 0.0], sigma=fitr.dcnt, absolute_sigma=True)

t.append(24.0)
mu.append(popt[0])
dmu.append(np.sqrt(np.diag(pcov))[0])

plt.plot(fitr.index, gauss(fitr.index, *popt), "-")


# Fit 3
fitr = data[(data.index > 2030) & (data.index < 2230)]
popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[2130.0, 50.0, 100000.0, 0.0], sigma=fitr.dcnt, absolute_sigma=True)

t.append(28.0)
mu.append(popt[0])
dmu.append(np.sqrt(np.diag(pcov))[0])

plt.plot(fitr.index, gauss(fitr.index, *popt), "-")


# Fit 4
fitr = data[(data.index > 2610) & (data.index < 2810)]
popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[2710.0, 50.0, 100000.0, 0.0], sigma=fitr.dcnt, absolute_sigma=True)

t.append(32.0)
mu.append(popt[0])
dmu.append(np.sqrt(np.diag(pcov))[0])

plt.plot(fitr.index, gauss(fitr.index, *popt), "-")


# Fit 5
fitr = data[(data.index > 3220) & (data.index < 3420)]
popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[3320.0, 50.0, 100000.0, 0.0], sigma=fitr.dcnt, absolute_sigma=True)

t.append(36.0)
mu.append(popt[0])
dmu.append(np.sqrt(np.diag(pcov))[0])

plt.plot(fitr.index, gauss(fitr.index, *popt), "-")


# Fit 6
fitr = data[(data.index > 3810) & (data.index < 4010)]
popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[3910.0, 50.0, 100000.0, 0.0], sigma=fitr.dcnt, absolute_sigma=True)

t.append(40.0)
mu.append(popt[0])
dmu.append(np.sqrt(np.diag(pcov))[0])

plt.plot(fitr.index, gauss(fitr.index, *popt), "-")


# Fit 7
fitr = data[(data.index > 4420) & (data.index < 4620)]
popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[4520.0, 50.0, 100000.0, 0.0], sigma=fitr.dcnt, absolute_sigma=True)

t.append(44.0)
mu.append(popt[0])
dmu.append(np.sqrt(np.diag(pcov))[0])

plt.plot(fitr.index, gauss(fitr.index, *popt), "-")


# Fit 8
fitr = data[(data.index > 5000) & (data.index < 5200)]
popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[5100.0, 50.0, 100000.0, 0.0], sigma=fitr.dcnt, absolute_sigma=True)

t.append(48.0)
mu.append(popt[0])
dmu.append(np.sqrt(np.diag(pcov))[0])

plt.plot(fitr.index, gauss(fitr.index, *popt), "-")

plt.close()


### Prompt-Kurve
prompt = pd.DataFrame(np.array([t, mu, dmu]).transpose(), columns=["t", "mu", "dmu"])
plt.errorbar(prompt.mu, prompt.t, xerr=prompt.dmu, yerr=np.zeros(8)+0.1, fmt="o")

f = lambda x, m, b: m * x + b
popt, pcov = curve_fit(f, prompt.mu, prompt.t)

fx = np.linspace(0.0, 6000.0, 1000)
fy = f(fx, *popt)

plt.plot(fx, fy, "-")

plt.xlim(0.0, 6000.0)
plt.xlabel("Kanal")
plt.ylabel("Verzögerung $\Delta t$ / ns")

plt.show()
