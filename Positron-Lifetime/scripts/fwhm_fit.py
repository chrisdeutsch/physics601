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
plt.xlim(750, 1200)

### Fit
def gauss(x, mu, sigma, norm, bg):
    return norm * np.exp(-0.5 * np.power((x - mu) / sigma, 2)) / np.sqrt(2.0 * np.pi) / sigma + bg


mu = float(data.idxmax()[0])
fit_data = data[(data.index > mu - 150.0) & (data.index < mu + 150.0)]

popt, pcov = curve_fit(gauss, fit_data.index, fit_data.cnt, p0=[960.0, 50.0, 1000000.0, 390], sigma=fit_data.dcnt, absolute_sigma=True)

plt.plot(fit_data.index, gauss(fit_data.index, *popt), "-", label="Anpassung")

plt.legend(loc="best")

### Chi Quadrat
chisq = np.power(gauss(fit_data.index, *popt) - fit_data.cnt, 2) / fit_data.dcnt**2
ndf = len(chisq) - 5
chisq = chisq.get_values().sum()
chisqndf = chisq / ndf

print(chisqndf)

plt.savefig("figures/fwhm_fit.pdf")