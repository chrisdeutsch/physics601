import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.odr import Model, RealData, ODR, polynomial
from scipy.linalg import block_diag
from scripts.tools import set_mpl_comma

plt.style.use("publication")

# raw data
raw_data = pd.read_csv("data/root/wmass/w_mass.csv")

### MC GAUGE CURVE
mc = raw_data[raw_data.set == "MC"]

### ODR FIT TO MC-GAUGE
fit_f = lambda B, x: B[0] * x + B[1]
linear_model = Model(fit_f)
fit_data = RealData(mc.half_height, mc.mc_w_mass, sx=mc.dhalf_height, sy=0.1)
odr = ODR(fit_data, linear_model, beta0=[2.0, -3.0])
odr_res = odr.run()

# Fitted parameters
m, b = odr_res.beta[0], odr_res.beta[1]

# Covariance matrix
cov_mat = odr_res.res_var * odr_res.cov_beta

### CALCULATION OF W MASS
w_meas = raw_data[raw_data.set == "ATLAS_W"]
HH = w_meas.half_height.values[0]
dHH = w_meas.dhalf_height.values[0]

# Error propagation
# covariance matrix for 3 parameters m, b, HH
cov_gauge = block_diag(cov_mat, dHH**2)
cov_gauge = np.matrix(cov_gauge)

# Jacobian J = (HH, 1.0, m) for gauge curve
jac = np.matrix([[HH, 1.0, m]])

# W mass
w_mass = fit_f(odr_res.beta, HH)
dw_mass = np.sqrt(jac * cov_gauge * jac.transpose())[0,0]

print("\nW mass: " + str(w_mass) + " +- " + str(dw_mass))

### PLOTTING
def make_gauge_plot():
    set_mpl_comma()
    plt.xlim(40.8, 42.2)
    plt.ylim(78.5, 81.5)
    plt.xlabel(r"Position des Halbhöhepunktes~$h$ / $\mathrm{GeV}$")
    plt.ylabel(r"simulierte Masse~$M_\mathrm{MC}$ / $\mathrm{GeV}$")
    
    # Data
    plt.errorbar(mc.half_height, mc.mc_w_mass, xerr=mc.dhalf_height, fmt="o",
                 zorder=2, label="MC-Simulation")
    
    # Fit
    x_fit = np.linspace(*plt.xlim(), 200)
    y_fit = fit_f(odr_res.beta, x_fit)
    plt.plot(x_fit, y_fit, "-", zorder=1, label="Anpassung")
    
    plt.legend(loc="best")
    
    plt.savefig("figures/wmass/gauge.pdf")
    plt.close()