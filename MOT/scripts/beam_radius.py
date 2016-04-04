import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Model
from scipy.special import erf

plt.style.use("publication")

data = pd.read_csv("data/beam_radius.csv", comment="#")
data.columns = ["d", "P", "dP"]

# Remove offset
data.d -= data.d.min()

# d error
data["dd"] = 0.02


# Fit beam
def erf_fit(d, P0, w, mu):
    return 0.5 * P0 * (1.0 - erf( np.sqrt(2.0) * (d - mu) / w))

beam_mod = Model(erf_fit)
beam_fit = beam_mod.fit(data.P, d=data.d, weights=data.dP**-1, P0=2.8, w=0.3, mu=0.7)

with open("data/beam_radius_fit.txt", "w") as f:
    f.write(beam_fit.fit_report())


# PLOT
plt.xlim(-0.1, 1.2)
plt.ylim(-0.1, 3.0)
plt.xlabel(r"Klingenposition $d$ / \si{cm}")
plt.ylabel(r"Leistung $P$ / \si{mW}")

plt.errorbar(data.d, data.P, data.dP, xerr=data.dd, fmt="o", label="Messpunkte", zorder=2)

x = np.linspace(-0.1, 1.2, 100)
plt.plot(x, beam_fit.eval(d=x), "-", label="Anpassung", zorder=1)

plt.legend(loc=0)

plt.tight_layout()
plt.savefig("figures/beam_radius_fit.pdf")
plt.close()