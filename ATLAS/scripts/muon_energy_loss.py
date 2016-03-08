import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
data["p_loss"] = data.p_in - data.p_out

# Error: dp_loss = sqrt(dp_in^2 + dp_out^2)
data["dp_loss"] = np.sqrt(data.dp_in**2 + data.dp_out**2)


### ENERGY
muon_mass = 0.105658  # GeV
data["E_in"] = np.sqrt(muon_mass**2 + data.p_in**2)
data["E_out"] = np.sqrt(muon_mass*2 + data.p_out**2)

# Error: dE = sqrt(p^2 dp^2 / (m^2 + p^2))
data["dE_in"] = np.sqrt(data.p_in**2 * data.dp_in**2 / (muon_mass**2 + data.p_in**2))
data["dE_out"] = np.sqrt(data.p_out**2 * data.dp_out**2 / (muon_mass**2 + data.p_out**2))

### ENERGY LOSS
data["E_loss"] = data.E_in - data.E_out
data["dE_loss"] = np.sqrt(data.dE_in**2 + data.dE_out**2)


### LATEX

### PLOTS
