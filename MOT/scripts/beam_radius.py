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

# To millimeter
data.d *= 10


# d error
data["dd"] = 0.2


# Fit beam
def erf_fit(d, P0, w, mu):
    return 0.5 * P0 * (1.0 - erf( np.sqrt(2.0) * (d - mu) / w))

beam_mod = Model(erf_fit)
beam_fit = beam_mod.fit(data.P, d=data.d, weights=data.dP**-1, P0=2.8, w=3.0, mu=7.0)

with open("data/beam_radius_fit.txt", "w") as f:
    f.write(beam_fit.fit_report())


# PLOT
from scripts.tools import set_mpl_comma
set_mpl_comma()
plt.xlim(-1, 12)
plt.ylim(-0.1, 3.0)
plt.xlabel(r"Klingenposition $x$ / \si{mm}")
plt.ylabel(r"Leistung $P$ / \si{mW}")

plt.errorbar(data.d, data.P, data.dP, xerr=data.dd, fmt="o", label="Messpunkte", zorder=2)

x = np.linspace(-1, 12, 100)
plt.plot(x, beam_fit.eval(d=x), "-", label="Anpassung", zorder=1)

plt.legend(loc=0)

plt.tight_layout()
plt.savefig("figures/beam_radius_fit.pdf")
plt.close()


# Latex
from scripts.tools import round
out = data[["d", "dd", "P", "dP"]]
out.columns = ["{$x$ / \si{mm}}", "{$\Delta x$ / \si{mm}}",
               "{$P$ / \si{mW}}", "{$\Delta P$ / \si{mW}}"]

out.to_latex("tables/beam_radius.tex", index=False,
             formatters=[round(1), round(1), round(2), round(2)], 
             column_format="S[table-format=1.1]S[table-format=1.1]S[table-format=1.2]S[table-format=1.2]", 
             escape=False)


