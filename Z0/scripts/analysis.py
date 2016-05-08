import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from lmfit.models import LorentzianModel, Model
from uncertainties import ufloat, umath, correlated_values, unumpy
from scipy.stats import chi2

plt.style.use("publication")
np.random.seed(12856712)


# Breit Wigner
def breit_wigner(s, peak, fwhm, mass):
    return peak * s * fwhm**2 / ( (s - mass**2)**2 + s**2 * fwhm**2 / mass**2)



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
for i in range(100000):
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

# error on counts (binomial)
which = ["e", "mu", "tau", "hadr"]

for s in which:
    data["{}_err".format(s)] = np.sqrt(data.events * (data[s] / data.events - (data[s] / data.events)**2))


# error prop for efficiency
for s, i in zip(which, range(4)):
    sum = 0.0
    for c, k in zip(which, range(4)):
        sum += data[c]**2 * inv_efficiency_err[i,k]**2 +\
               inv_efficiency[i,k]**2 * data["{}_err".format(c)]**2
    
    data["{}_corr_err".format(s)] = np.sqrt(sum)


# calculate cross section
data["sig_e"] = data.e_corr / data.lumi + data.corr_lep
data["sig_mu"] = data.mu_corr / data.lumi + data.corr_lep
data["sig_tau"] = data.tau_corr / data.lumi + data.corr_lep
data["sig_hadr"] = data.hadr_corr / data.lumi + data.corr_hadr

data["sig_e_err"] = np.sqrt(data.e_corr_err**2 / data.lumi**2 + data.e_corr**2 / data.lumi**4 * data.lumi_tot**2)
data["sig_mu_err"] = np.sqrt(data.mu_corr_err**2 / data.lumi**2 + data.mu_corr**2 / data.lumi**4 * data.lumi_tot**2)
data["sig_tau_err"] = np.sqrt(data.tau_corr_err**2 / data.lumi**2 + data.tau_corr**2 / data.lumi**4 * data.lumi_tot**2)
data["sig_hadr_err"] = np.sqrt(data.hadr_corr_err**2 / data.lumi**2 + data.hadr_corr**2 / data.lumi**4 * data.lumi_tot**2)

# Forward backward asymmetry
data["afb_events"] = data.fwd + data.bwd
data["fwd_err"] = np.sqrt(data.afb_events * (data.fwd / data.afb_events - (data.fwd / data.afb_events)**2))
data["bwd_err"] = np.sqrt(data.afb_events * (data.bwd / data.afb_events - (data.bwd / data.afb_events)**2))


data["afb"] = (data.fwd - data.bwd) / (data.fwd + data.bwd) + data.corr_afb
data["afb_err"] = 2.0 / (data.fwd + data.bwd)**2 * np.sqrt(data.bwd**2 * data.fwd_err**2 + data.fwd**2 * data.bwd_err**2)

afb_data = data[["afb", "afb_err"]].get_values()
afb = ufloat(afb_data[3,0], afb_data[3,1])

weinberg = (1.0 - umath.sqrt(afb / 3.0)) / 4.0


# BreitWignerFits
model = LorentzianModel()
E_cm = data.E_cm.get_values()

e_fit = model.fit(data.sig_e.get_values(), x=E_cm, weights=data.sig_e_err**-1,
                  amplitude=8.0, center= 91.2, sigma=1.0)
mu_fit = model.fit(data.sig_mu.get_values(), x=E_cm, weights=data.sig_mu_err**-1,
                   amplitude=8.0, center= 91.2,  sigma=1.0)
tau_fit = model.fit(data.sig_tau.get_values(), x=E_cm, weights=data.sig_tau_err**-1,
                   amplitude=8.0, center= 91.2,  sigma=1.0)
hadr_fit = model.fit(data.sig_hadr.get_values(), x=E_cm, weights=data.sig_hadr_err**-1,
                   amplitude=160.0, center= 91.2,  sigma=1.0)


# chi square test
ndf = 4
e_pval = 1 - chi2.cdf(e_fit.chisqr, ndf)
mu_pval = 1 - chi2.cdf(mu_fit.chisqr, ndf)
tau_pval = 1 - chi2.cdf(tau_fit.chisqr, ndf)
hadr_pval = 1 - chi2.cdf(hadr_fit.chisqr, ndf)


# Calculate peak cross sections
(e_sigma, e_mu, e_A) = correlated_values([e_fit.values[name] for name in ["sigma", "center", "amplitude"]], e_fit.covar)
(mu_sigma, mu_mu, mu_A) = correlated_values([mu_fit.values[name] for name in ["sigma", "center", "amplitude"]], mu_fit.covar)
(tau_sigma, tau_mu, tau_A) = correlated_values([tau_fit.values[name] for name in ["sigma", "center", "amplitude"]], tau_fit.covar)
(hadr_sigma, hadr_mu, hadr_A) = correlated_values([hadr_fit.values[name] for name in ["sigma", "center", "amplitude"]], hadr_fit.covar)

e_peak = e_A / (np.pi * e_sigma)
mu_peak = mu_A / (np.pi * mu_sigma)
tau_peak = tau_A / (np.pi * tau_sigma)
hadr_peak = hadr_A / (np.pi * hadr_sigma)

nb_gev = 2.5681897e-6 # GeV^-2 per nbarn

# error weighted mean total decay width
Z_fwhm = np.array([2.0 * e_sigma, 2.0 * mu_sigma, 2.0 * tau_sigma, 2.0 * hadr_sigma])
weights = unumpy.std_devs(Z_fwhm)**-2
var = weights.sum()**-1
mean = unumpy.nominal_values(Z_fwhm) * weights
mean = mean.sum()
mean *= var

Z_fwhm = ufloat(mean, np.sqrt(var))

#error weighted mean mass
Z_mass = np.array([e_mu, mu_mu, tau_mu, hadr_mu])
weights = unumpy.std_devs(Z_mass)**-2
var = weights.sum()**-1
mean = unumpy.nominal_values(Z_mass) * weights
mean = mean.sum()
mean *= var

Z_mass = ufloat(mean, np.sqrt(var))

# partial decay width
e_part_fwhm = umath.sqrt(Z_fwhm**2 * e_peak * nb_gev * Z_mass**2 /(12.0 * np.pi))
mu_part_fwhm = Z_fwhm**2 * mu_peak * nb_gev * Z_mass**2 / (12.0 * np.pi * e_part_fwhm)
tau_part_fwhm = Z_fwhm**2 * tau_peak * nb_gev * Z_mass**2 / (12.0 * np.pi * e_part_fwhm)
hadr_part_fwhm = Z_fwhm**2 * hadr_peak * nb_gev * Z_mass**2 / (12.0 * np.pi * e_part_fwhm)

# Lepton universality
sig_mu_e = mu_peak / e_peak
sig_tau_e = tau_peak / e_peak

r_mu_e = mu_part_fwhm / e_part_fwhm
r_tau_e = tau_part_fwhm / e_part_fwhm

# ratio of hadronic to lepton
r_hadr_lepton = hadr_part_fwhm / (e_part_fwhm + mu_part_fwhm + tau_part_fwhm)





# number of neutrino generations
neutrino_fwhm = 0.1678 # GeV
N_nu = (Z_fwhm - e_part_fwhm - mu_part_fwhm - tau_part_fwhm - hadr_part_fwhm) / neutrino_fwhm


print(e_fit.fit_report())
print(mu_fit.fit_report())
print(tau_fit.fit_report())
print(hadr_fit.fit_report())


### Plots

# Cross sections
fig = plt.gcf()
fig.set_size_inches(5.4, 4.8)


x = np.linspace(88.0, 94.0, 1000)

ax1 = plt.subplot(221)
plt.errorbar(data.E_cm.get_values(), data.sig_e.get_values(), yerr=data.sig_e_err.get_values(), fmt="o", zorder=2, markersize=2.5)
plt.plot(x, e_fit.eval(x=x), "-", zorder=1)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.ylim((0.0, 2.0))
plt.title("electrons", fontsize=11)

ax2 = plt.subplot(222, sharex=ax1)
plt.errorbar(data.E_cm.get_values(), data.sig_mu.get_values(), yerr=data.sig_mu_err.get_values(), fmt="o", zorder=2, markersize=2.5)
plt.plot(x, mu_fit.eval(x=x), "-", zorder=1)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.ylim((0.0, 2.0))
plt.title("muons", fontsize=11)


ax3 = plt.subplot(223)
plt.errorbar(data.E_cm.get_values(), data.sig_tau.get_values(), yerr=data.sig_tau_err.get_values(), fmt="o", zorder=2, markersize=2.5)
plt.plot(x, tau_fit.eval(x=x), "-", zorder=1)
plt.ylim((0.0, 2.2))
plt.title("tauons", fontsize=11)

ax4 = plt.subplot(224, sharex=ax2)
plt.errorbar(data.E_cm.get_values(), data.sig_hadr.get_values(), yerr=data.sig_hadr_err.get_values(), fmt="o", zorder=2, markersize=2.5)
plt.plot(x, hadr_fit.eval(x=x), "-", zorder=1)
plt.ylim((0.0, 50.0))
plt.title("hadrons", fontsize=11)

axes = [ax1, ax2, ax3, ax4]
yaxes = [ax1, ax3]
xaxes = [ax3, ax4]
for ax in xaxes:
    ax.set_xlabel(r"$\sqrt{s}$ / GeV", fontsize=11)

for ax in yaxes:
    ax.set_ylabel(r"$\sigma$ / nb", fontsize=11)

for ax in axes:
    ax.tick_params(axis="both", labelsize=9)


plt.tight_layout(pad=0.5)
plt.savefig("figures/cross_sections.pdf")
plt.close()


# AFB plot
fig = plt.gcf()
fig.set_size_inches(0.8 * fig.get_size_inches())

plt.xlabel(r"$\sqrt{s}$ / GeV")
plt.ylabel(r"$A_\mathrm{FB}$")
plt.errorbar(data.E_cm.get_values(), data.afb.get_values(), yerr=data.afb_err.get_values(), fmt="o")

plt.tight_layout(pad=0.2)
plt.savefig("figures/afb.pdf")
plt.close()

