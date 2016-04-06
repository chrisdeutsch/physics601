import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Model

plt.style.use("publication")

data = pd.read_csv("data/lambda_4.csv", comment="#")
data.columns = ["phi", "P_in", "P_out"]

# Subtract background
data.P_in -= 248.0
data.P_out -= 232.0

data["dP"] = 5.0


# Fit to incoming beam ???


# Fit to outgoing beam
def out_f(phi, c):
    return np.zeros(len(phi)) + c

out_model = Model(out_f)
out_fit = out_model.fit(data.P_out, phi=data.phi, weights=data.dP**-1, c=110.0)

with open("data/lambda_4_out_fit.txt", "w") as f:
    f.write(out_fit.fit_report())


### PLOT
# incoming beam
plt.xlim(0.0, 180.0)
# plt.ylim(-0.1, 3.0)
plt.xlabel(r"Winkel der Verzögerungsplatte $\varphi$ / \si{\degree}")
plt.ylabel(r"Fluoreszenzleistung $P$ / \si{nW}")

plt.errorbar(data.phi, data.P_in, data.dP, fmt="o", label="Messpunkte", zorder=2)

# x = np.linspace(-0.1, 1.2, 100)
# plt.plot(x, beam_fit.eval(d=x), "-", label="Anpassung", zorder=1)

# plt.legend(loc=0)

plt.tight_layout()
plt.savefig("figures/lambda_4_in.pdf")
plt.close()


# outgoing beam
plt.xlim(0.0, 180.0)
# plt.ylim(-0.1, 3.0)
plt.xlabel(r"Winkel der Verzögerungsplatte $\varphi$ / \si{\degree}")
plt.ylabel(r"Fluoreszenzleistung $P$ / \si{nW}")

plt.errorbar(data.phi, data.P_out, data.dP, fmt="o", label="Messpunkte", zorder=2)

x = np.linspace(0.0, 180.0, 100)
plt.plot(x, out_fit.eval(phi=x), "-", label="Anpassung", zorder=1)

plt.legend(loc=0)

plt.tight_layout()
plt.savefig("figures/lambda_4_out.pdf")
plt.close()


### LATEX Tabelle
from scripts.tools import round
out = data[["phi","P_in","P_out"]]
out.columns = [r"{$\varphi$ / \si{\degree}}",
               r"{$P_\mathrm{ein.}$ / \si{nW}}",
               r"{$P_\mathrm{ref.}$ / \si{nW}}"]

out.to_latex("tables/lambda_4.tex", index=False,
             formatters=[round(0), round(0), round(0)], 
             column_format="S[table-format=3.0]S[table-format=3.0]S[table-format=3.0]", 
             escape=False)






