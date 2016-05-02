import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from lmfit.models import LorentzianModel
from mcerp import Binomial

np.random.seed(12856712)

### Monte-Carlo Data

# After precut
precut = np.array([93802, 94381, 79214, 98563], dtype=np.float64)

# After |cos|<0.9 cut (only electron and muon)
cos_cut = np.array([37191, 81847, 79214, 98563], dtype=np.float64)

# After s-channel cut -0.9 < cos < 0.5 (only electron)
s_cut = np.array([20531, 81847, 79214, 98563], dtype=np.float64)


# Final cuts
e_cut = np.array([20229, 1, 117, 0], dtype=np.float64)
mu_cut = np.array([0, 78537, 549, 0], dtype=np.float64)
tau_cut = np.array([692, 4089, 76651, 934], dtype=np.float64)
hadr_cut = np.array([6, 0, 503, 97533], dtype=np.float64)


### Calculate s-channel cut efficiency
a, b = (-0.9, 0.5)
korr_s_cut = 8 / 3 * (-a - a**3 / 3 + b + b**3 / 3)**-1

### Total number of simulated s-channel events
precut_factor = 100000 / precut[0]
true_s_events = e_cut[0] * precut_factor * korr_s_cut


# Calculate efficiencies
e_cut[0] /= true_s_events
mu_cut[0] /= true_s_events
tau_cut[0] /= true_s_events
hadr_cut[0] /= true_s_events

e_cut[1:] /= 100000
mu_cut[1:] /= 100000
tau_cut[1:] /= 100000
hadr_cut[1:] /= 100000

# Efficiency matrix
efficiency = np.array([e_cut, mu_cut, tau_cut, hadr_cut])
inv_efficiency = np.linalg.inv(efficiency)

# Calculate errors
efficiency_err = np.zeros((4,4))

for i in range(4):
    efficiency_err[i,0] = np.sqrt((efficiency[i,0]-efficiency[i,0]**2) / true_s_events)

for i in range(1,4):
    for j in range(1,4):
        efficiency_err[i,j] = np.sqrt((efficiency[i,j]-efficiency[i,j]**2) / 100000)


# Monte Carlo Error Propagation
matrices = []
for i in range(10000):
    dE = np.random.randn(16)
    dE.resize((4,4))
    dE *= efficiency_err
    E = efficiency + dE
    matrices.append(np.linalg.inv(E).flatten())

matrices = np.array(matrices)
matrices = matrices.transpose()

inv_efficiency_err = np.array([m.std() for m in matrices])
inv_efficiency_err = inv_efficiency_err.reshape((4,4))


### Cross section measurement
data = pd.read_csv("data/cross_sec_data.csv", dtype=np.float64)


# apply efficiency correction
pvecs = data[["e", "mu", "tau", "hadr"]].get_values()
pvecs = np.array([np.dot(inv_efficiency, pvec) for pvec in pvecs])
corr = pd.DataFrame(pvecs, columns=["e_corr", "mu_corr", "tau_corr", "hadr_corr"])
data = pd.concat([data, corr], axis=1)

data["sig_e"] = data.e_corr / data.lumi + data.corr_lep
data["sig_mu"] = data.mu_corr / data.lumi + data.corr_lep
data["sig_tau"] = data.tau_corr / data.lumi + data.corr_lep
data["sig_hadr"] = data.hadr_corr / data.lumi + data.corr_hadr


# BreitWignerFits
model = LorentzianModel()

E_cm = data.E_cm.get_values()
e_fit = model.fit(data.sig_e.get_values(), x=E_cm,
                  amplitude=8.0, center= 91.2,  sigma=1.0)
mu_fit = model.fit(data.sig_mu.get_values(), x=E_cm,
                   amplitude=8.0, center= 91.2,  sigma=1.0)
tau_fit = model.fit(data.sig_tau.get_values(), x=E_cm,
                   amplitude=8.0, center= 91.2,  sigma=1.0)
hadr_fit = model.fit(data.sig_hadr.get_values(), x=E_cm,
                   amplitude=160.0, center= 91.2,  sigma=1.0)


print(efficiency)
print(e_fit.fit_report())
print(mu_fit.fit_report())
print(tau_fit.fit_report())
print(hadr_fit.fit_report())
