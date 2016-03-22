import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.style.use("publication")


data = pd.read_csv("data/prompt/prompt_curve_20min.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)

### Prompt plot
plt.errorbar(data.index, data.cnt, yerr=data.dcnt, fmt="+", label="Messpunkte")

### Fit
def gauss(x, mu, sigma, norm, bg):
    return norm * np.exp(-0.5 * np.power((x - mu) / sigma, 2)) / np.sqrt(2.0 * np.pi) / sigma + bg

mu = 963.0
fitr = data[(data.index > mu - 140.0) & (data.index < mu + 140.0)]

popt, pcov = curve_fit(gauss, fitr.index, fitr.cnt, p0=[960.0, 50.0, 1000000.0, 390], sigma=fitr.dcnt, absolute_sigma=True)

plt.plot(fitr.index, gauss(fitr.index, *popt), "-", label="Anpassung")

plt.xlim(750, 1200)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")
plt.legend(loc="best")

### Chi Quadrat
chisq = np.power(gauss(fitr.index, *popt) - fitr.cnt, 2) / fitr.dcnt**2
ndf = len(chisq) - 5
chisq = chisq.get_values().sum() / ndf

### FWHM
sigma = popt[1]
dsigma = np.sqrt(np.diag(pcov))[1]
a = 2.0 * np.sqrt(2.0 * np.log(2.0))
fwhm = a * sigma
dfwhm = a * dsigma


print("Chi^2/ndf: {}".format(chisq))
print("sigma: {} +- {}".format(sigma, dsigma))
print("FWHM: {} +- {}".format(fwhm, dfwhm))


plt.savefig("figures/fwhm_fit.pdf")
plt.close()
