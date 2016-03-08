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
