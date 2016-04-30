import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sp

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
E = np.array([e_cut, mu_cut, tau_cut, hadr_cut])
