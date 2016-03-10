import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# raw data
raw_data = pd.read_csv("data/atlantis/invariant_mass_zee.csv", comment="#")

# data for analysis
data = raw_data.drop(["p_1", "p_2"], axis=1)
data = data.sort_values(by="event")


### MOMENTUM
# p = pT cosh(eta)
data["p_1"] = data.p_T_1 * np.cosh(data.eta_1)
data["p_2"] = data.p_T_2 * np.cosh(data.eta_2)
# with error: dp = |cosh(eta) dpT|    (error in eta is neglected)
data["dp_1"] = np.abs(np.cosh(data.eta_1) * data.dp_T_1)
data["dp_2"] = np.abs(np.cosh(data.eta_2) * data.dp_T_2)

# comparison with atlantis (max. 1% deviation)
p_1_ratio = data.p_1 / raw_data.p_1
p_2_ratio = data.p_2 / raw_data.p_2

if ((p_1_ratio < 0.99) | (p_2_ratio > 1.01)).any():
    print("WARNING: deviation in p_in calculation")

if ((p_1_ratio < 0.99) | (p_2_ratio > 1.01)).any():
    print("WARNING: deviation in p_out calculation")


### THETA CMS
data["theta_1"] = 2 * np.arctan(np.exp(-data.eta_1))
data["theta_2"] = 2 * np.arctan(np.exp(-data.eta_2))


### UNIT VECTORS
data["unitv_x_1"] = np.sin(data.theta_1) * np.cos(data.phi0_1 * np.pi / 180.0)
data["unitv_y_1"] = np.sin(data.theta_1) * np.sin(data.phi0_1 * np.pi / 180.0)
data["unitv_z_1"] = np.cos(data.theta_1)

data["unitv_x_2"] = np.sin(data.theta_2) * np.cos(data.phi0_2 * np.pi / 180.0)
data["unitv_y_2"] = np.sin(data.theta_2) * np.sin(data.phi0_2 * np.pi / 180.0)
data["unitv_z_2"] = np.cos(data.theta_2)

### COS
data["cos"] = data.unitv_x_1 * data.unitv_x_2 + \
              data.unitv_y_1 * data.unitv_y_2 + \
              data.unitv_z_1 * data.unitv_z_2


### ENERGIES
mass_electron = 0.000511  # GeV
# E = sqrt(m^2 + p^2)
data["E_1_from_p"] = np.sqrt(mass_electron**2 + data.p_1**2)
data["E_2_from_p"] = np.sqrt(mass_electron**2 + data.p_2**2)
# dE = sqrt( p^2 dp^2 / (m^2 + p^2) )
data["dE_1_from_p"] = np.sqrt( data.p_1**2 * data.dp_1**2 /
                               (mass_electron**2 + data.p_1**2) )
data["dE_2_from_p"] = np.sqrt( data.p_2**2 * data.dp_2**2 /
                               (mass_electron**2 + data.p_2**2) )

### INVARIANT MASS
data["inv_mass"] = np.sqrt(2 * mass_electron**2 + \
                           2 * (data.E_1_from_p * data.E_2_from_p - \
                                np.abs(data.p_1 * data.p_2) * data.cos))


data["inv_mass_me0"] = np.sqrt(2 * np.abs(data.p_1 * data.p_2) * (1 - data.cos))

data["dinv_mass_me0"] = np.sqrt((1 - data.cos) * np.abs(data.p_2) * data.dp_1**2 /
                                (2 * np.abs(data.p_1)) +
                                (1 - data.cos) * np.abs(data.p_1) * data.dp_2**2 /
                                (2 * np.abs(data.p_2)))

data["dinv_mass"] = data.dinv_mass_me0


### LATEX
from scripts.tools import round

# Invariant masses
out = data[data.event != 24]


pt1 = out[["event", "p_T_1", "dp_T_1", "eta_1", "phi0_1", "p_1", "dp_1"]]
pt2 = out[["event", "p_T_2", "dp_T_2", "eta_2", "phi0_2", "p_2", "dp_2"]]
pt3 = out[["event", "cos", "inv_mass", "dinv_mass", "inv_mass_me0", "dinv_mass_me0"]]


pt1.columns = ["{Ereignis}",
               "{$p_{\mathrm{T}^{(1)}}$ / \si{\GeV}}",
               "{$\sigma_{p_\mathrm{T}^{(1)}}$ / \si{\GeV}}",
               "{$\eta^{(1)}$}",
               "{$\phi^{(1)}$ / \si{\degree}}",
               "{$p^{(1)}$ / \si{\GeV}}",
               "{$\sigma_{p^{(1)}}$ / \si{\GeV}}"]

pt2.columns = ["{Ereignis}",
               "{$p_{\mathrm{T}^{(2)}}$ / \si{\GeV}}",
               "{$\sigma_{p_\mathrm{T}^{(2)}}$ / \si{\GeV}}",
               "{$\eta^{(2)}$}",
               "{$\phi^{(2)}$ / \si{\degree}}",
               "{$p^{(2)}$ / \si{\GeV}}",
               "{$\sigma_{p^{(2)}}$ / \si{\GeV}}"]

pt3.columns = ["{Ereignis}",
               "{$\cos\theta$}",
               "{$m_\mathrm{inv}^\mathrm{exakt}$ / \si{\GeV}}",
               "{$\sigma_{m_\mathrm{inv}^\mathrm{exakt}}$ / \si{\GeV}}",
               "{$m_\mathrm{inv}$ / \si{\GeV}}",
               "{$\sigma_{m_\mathrm{inv}}$ / \si{\GeV}}"]


pt1.to_latex("tables/zee_inv_mass_pt1.tex", index=False,
             formatters=[str, round(2), round(2), round(2), round(1), round(1), round(1)],
             column_format="rSSSSSS", escape=False)

pt2.to_latex("tables/zee_inv_mass_pt2.tex", index=False,
             formatters=[str, round(2), round(2), round(2), round(1), round(1), round(1)],
             column_format="rSSSSSS", escape=False)

pt3.to_latex("tables/zee_inv_mass_pt3.tex", index=False,
             formatters=[str, round(3), round(1), round(1), round(1), round(1)],
             column_format="rSSSSS", escape=False)

