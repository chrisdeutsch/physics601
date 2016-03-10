import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

plt.style.use("publication")

# raw data
raw_data = pd.read_csv("data/atlantis/muon_energy_loss.csv", comment="#")

# data for analysis
data = raw_data.drop(["p_in", "p_out"], axis=1)


### MOMENTUM
# p = pT cosh(eta)
data["p_in"] = data.p_T_in * np.cosh(data.eta_in)
data["p_out"] = data.p_T_out * np.cosh(data.eta_out)
# with error: dp = |cosh(eta) dpT|    (error in eta is neglected)
data["dp_in"] = np.abs(np.cosh(data.eta_in) * data.dp_T_in)
data["dp_out"] = np.abs(np.cosh(data.eta_out) * data.dp_T_out)

# comparison with atlantis (max. 1% deviation)
p_in_ratio = data.p_in / raw_data.p_in
p_out_ratio = data.p_out / raw_data.p_out

if ((p_in_ratio < 0.99) | (p_in_ratio > 1.01)).any():
    print("WARNING: deviation in p_in calculation")

if ((p_out_ratio < 0.99) | (p_out_ratio > 1.01)).any():
    print("WARNING: deviation in p_out calculation")


### MOMENTUM LOSS
data["p_loss"] = np.abs(data.p_in) - np.abs(data.p_out)

# Error: dp_loss = sqrt(dp_in^2 + dp_out^2)
data["dp_loss"] = np.sqrt(data.dp_in**2 + data.dp_out**2)


### ENERGY
muon_mass = 0.105658  # GeV
data["E_in"] = np.sqrt(muon_mass**2 + data.p_in**2)
data["E_out"] = np.sqrt(muon_mass**2 + data.p_out**2)

# Error: dE = sqrt(p^2 dp^2 / (m^2 + p^2))
data["dE_in"] = np.sqrt(data.p_in**2 * data.dp_in**2 / (muon_mass**2 + data.p_in**2))
data["dE_out"] = np.sqrt(data.p_out**2 * data.dp_out**2 / (muon_mass**2 + data.p_out**2))


### ENERGY LOSS
data["E_loss"] = data.E_in - data.E_out
data["dE_loss"] = np.sqrt(data.dE_in**2 + data.dE_out**2)


### MEAN ENERGY LOSS
df = data[["E_loss", "dE_loss"]].dropna()
f = lambda x, b: 0.0 * x + b
popt, pcov = curve_fit(f, df.index, df.E_loss, sigma=df.dE_loss)

chisquare = ((np.zeros(df.E_loss.size) + popt[0]) - df.E_loss)**2 / df.dE_loss**2
chisquare = chisquare.sum()
# ndf = N - n - 1 (N: data points; n: fit parameters)
chisquare_red = chisquare / (df.E_loss.size - 2)

print("Average energy loss: {:.2f} +- {:.2f}".format(popt[0], 
                                                     np.sqrt(np.diag(pcov))[0]))
print("chi/ndf: {:.3f}".format(chisquare_red))


### LATEX
from scripts.tools import round

# Momentum calculation
out = data[["event", "p_T_in", "dp_T_in", "eta_in", "p_in", "dp_in", 
                     "p_T_out", "dp_T_out", "eta_out", "p_out", "dp_out"]]
out.columns = ["{Ereignis}",
               "{$p_T$ / \si{\GeV}}",
               "{$\sigma_{p_T}$ / \si{\GeV}}",
               "{$\eta$}",
               "{$p$ / \si{\GeV}}",
               "{$\sigma_p$ / \si{\GeV}}",
               "{$p_T^\prime$ / \si{\GeV}}",
               "{$\sigma_{p_T^\prime}$ / \si{\GeV}}",
               "{$\eta^\prime$}",
               "{$p^\prime$ / \si{\GeV}}",
               "{$\sigma_{p^\prime}$ / \si{\GeV}}"]

out.to_latex("tables/muon_momenta.tex", index=False,
             formatters=[str, round(1), round(1), round(2), round(1), round(1),
                              round(1), round(1), round(2), round(1), round(1)],
             column_format="rSSSSSSSSSS", escape=False)


# Momentum and energy loss
out = data[["event", "E_loss", "dE_loss"]]
out = out.dropna()
out.columns = ["{Ereignis}",
               "{$\Delta E$ / \si{\GeV}}",
               "{$\sigma_{\Delta E}$ / \si{\GeV}}"]

out[out["{Ereignis}"] <= 11].to_latex("tables/muon_energy_loss_pt1.tex",
    index=False, formatters=[str, round(1), round(1)],
    column_format="rSS", escape=False)

out[out["{Ereignis}"] > 11].to_latex("tables/muon_energy_loss_pt2.tex",
    index=False, formatters=[str, round(1), round(1)],
    column_format="rSS", escape=False)

### PLOTS
def make_plots():
    # Pseudorapidity dependence
    plt.errorbar(np.abs(data.eta_in), data.E_loss, yerr=data.dE_loss,
                 fmt="o", zorder=2, label="Messpunkte")
    plt.plot(np.linspace(0.0, 3.0, 100), np.zeros(100) + popt[0], "-", zorder=1,
             label=r"$\langle \Delta E \rangle$")

    plt.xlabel("Pseudorapidität~$|\eta|$")
    plt.ylabel("Energieverlust~$\Delta E$ / GeV")
    plt.grid()
    plt.legend(loc="lower left")
    
    plt.savefig("figures/muon_energy_loss/eta.pdf")
    plt.close()

    # Momentum dependence
    plt.errorbar(np.abs(data.p_in), data.E_loss, xerr=data.dp_in,
                 yerr=data.dE_loss, fmt="o", zorder=2, label="Messpunkte")
    plt.plot(np.linspace(0.0, 250.0, 100), np.zeros(100) + popt[0], "-",
             zorder=1, label=r"$\langle \Delta E \rangle$")
    
    plt.xlabel("Impuls~$p$")
    plt.ylabel("Energieverlust~$\Delta E$ / GeV")
    plt.grid()
    plt.legend(loc="lower left")

    plt.savefig("figures/muon_energy_loss/momentum.pdf")
    plt.close()
