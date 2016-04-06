import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("publication")

data = pd.read_csv("data/magnetic_field.csv", comment="#")
data.columns = ["I", "P"]

data["Pbg"] = data.P - 188
data["dPbg"] = 5.0


# Plot
plt.xlim(-0.2, 5.4)
plt.ylim(-10, 170)
plt.xlabel(r"Spulenstrom $I$ / \si{\ampere}")
plt.ylabel(r"Fluoreszenzleistung $P$ / \si{nW}")

plt.errorbar(data.I, data.Pbg, data.dPbg, fmt="o")

plt.tight_layout()
plt.savefig("figures/magnetic_field.pdf")
plt.close()

# Latex
from scripts.tools import round
out = data[["I", "Pbg"]]
out.columns = ["{Strom $I$ / \si{A}}", "{Leistung $P$ / \si{nW}}"]

out.to_latex("tables/magnetic_field.tex", index=False,
             formatters=[round(2), round(0)],
             column_format="S[table-format=1.2]S[table-format=3.0]",
             escape=False)
