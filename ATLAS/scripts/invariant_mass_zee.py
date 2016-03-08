import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# raw data
raw_data = pd.read_csv("data/atlantis/invariant_mass_zee.csv", comment="#")

# data for analysis
data = raw_data.drop(["p_1", "p_2"], axis=1)


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
data["theta_2"] = 2* np.arctan(np.exp(-data.eta_2))


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


### INVARIANT MASS
mass_electron = 0.000511  # GeV
data["inv_mass"] = np.sqrt(2 * mass_electron**2 + \
                           2 * (data.E_1 * data.E_2 - \
                                np.abs(data.p_1 * data.p_2) * data.cos))

data["inv_mass_me0"] = np.sqrt(2 * np.abs(data.p_1 * data.p_2) * (1 - data.cos))
