import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit.models import LorentzianModel, PolynomialModel, ConstantModel, GaussianModel
from uncertainties import ufloat, umath

plt.style.use("publication")

# Spectrum Fit Model
bg = PolynomialModel(4, prefix="bg_")
L1 = LorentzianModel(prefix="L1_")
L2 = LorentzianModel(prefix="L2_")
L3 = LorentzianModel(prefix="L3_")

spec_model = L1 + L2 + L3 + bg

# Cooling Peak Model
bg = ConstantModel()
G = GaussianModel()

cool_model = G + bg


### SIGNAL:

# Spektroskopie
sig_spec = pd.read_csv("data/detuning/cooling/signal/F0095CH1.csv",
                  usecols=[3, 4], header=None, names=["t", "I"])
# Cut data
sig_spec = sig_spec[sig_spec.t < 22.3]

# Fluoreszenz
sig_fluo = pd.read_csv("data/detuning/cooling/signal/F0095CH2.csv",
                  usecols=[3, 4], header=None, names=["t", "I"])

# Cut data
sig_fluo = sig_fluo[sig_fluo.t < 22.3]


### Fit to Spectrum
t = sig_spec.t.get_values()
I = sig_spec.I.get_values()

sig_spec_fit = spec_model.fit(I, x=t, L1_amplitude=1.0, L1_sigma=1.0, L1_center=5.8,
                              L2_amplitude=1.0, L2_sigma=1.0, L2_center=9.9,
                              L3_amplitude=0.5, L3_sigma=1.0, L3_center=18.2,
                              bg_c0=9.0, bg_c1=0.0, bg_c2=0.0, bg_c3=0.0, bg_c4=0.0)



### BACKGROUND:

# Spektroskopie
bg_spec = pd.read_csv("data/detuning/cooling/background/F0096CH1.csv",
                  usecols=[3, 4], header=None, names=["t", "I"])
# Cut data
bg_spec = bg_spec[bg_spec.t > 3.2]

# Fluoreszenz
bg_fluo = pd.read_csv("data/detuning/cooling/background/F0096CH2.csv",
                  usecols=[3, 4], header=None, names=["t", "I"])

# Cut data
bg_fluo = bg_fluo[bg_fluo.t > 3.2]


# Fit to spectrum
t = bg_spec.t.get_values()
I = bg_spec.I.get_values()

bg_spec_fit = spec_model.fit(I, x=t, L1_amplitude=1.0, L1_sigma=1.0, L1_center=9.4,
                             L2_amplitude=1.0, L2_sigma=1.0, L2_center=13.4,
                             L3_amplitude=0.5, L3_sigma=1.0, L3_center=22.0,
                             bg_c0=9.0, bg_c1=0.0, bg_c2=0.0, bg_c3=0.0, bg_c4=0.0)


# Match sig & bg spectra
sig_param = sig_spec_fit.params
bg_param = bg_spec_fit.params

param_names = ["L1_center", "L2_center", "L3_center"]
diffs = []
for p in param_names:
    diffs.append(sig_param[p].value - bg_param[p].value)

mean_diff = np.array(diffs).mean()

bg_spec.t += mean_diff
bg_fluo.t += mean_diff


# Find matching range
low = max(bg_spec.t.min(), sig_spec.t.min())
high = min(bg_spec.t.max(), sig_spec.t.max())

# Cut to range
sig_spec = sig_spec[(sig_spec.t > low) & (sig_spec.t < high)].reset_index()
sig_fluo = sig_fluo[(sig_fluo.t > low) & (sig_fluo.t < high)].reset_index()
bg_spec = bg_spec[(bg_spec.t > low) & (bg_spec.t < high)].reset_index()
bg_fluo = bg_fluo[(bg_fluo.t > low) & (bg_fluo.t < high)].reset_index()


# Refit Spectra
t = sig_spec.t.get_values()
I = sig_spec.I.get_values()

sig_spec_fit = spec_model.fit(I, x=t, L1_amplitude=1.0, L1_sigma=1.0, L1_center=5.8,
                              L2_amplitude=1.0, L2_sigma=1.0, L2_center=9.9,
                              L3_amplitude=0.5, L3_sigma=1.0, L3_center=18.2,
                              bg_c0=9.0, bg_c1=0.0, bg_c2=0.0, bg_c3=0.0, bg_c4=0.0)

t = bg_spec.t.get_values()
I = bg_spec.I.get_values()

bg_spec_fit = spec_model.fit(I, x=t, L1_amplitude=1.0, L1_sigma=1.0, L1_center=5.7,
                             L2_amplitude=1.0, L2_sigma=1.0, L2_center=9.7,
                             L3_amplitude=0.5, L3_sigma=1.0, L3_center=18.3,
                             bg_c0=9.0, bg_c1=0.0, bg_c2=0.0, bg_c3=0.0, bg_c4=0.0)


# Remove background
fluo = pd.DataFrame()
fluo["I"] = sig_fluo.I - bg_fluo.I
fluo["t"] = sig_fluo.t


# Fit cooling peak
t = fluo.t.get_values()
I = fluo.I.get_values()
cool_fit = cool_model.fit(fluo.I.get_values(), amplitude=0.06, sigma=1.0, center=16.8, c=0.0, x=t)

# Frequency scaling
spec_params = sig_spec_fit.params
cool_params = cool_fit.params

# For calibration: F->F' = 3->2,4
left_peak = ufloat(spec_params["L1_center"].value, spec_params["L1_center"].stderr)

# Cooling Transition F->F' = 3->4
cooling_transition = ufloat(spec_params["L3_center"].value, spec_params["L3_center"].stderr)
cooling_transition_fwhm = ufloat(spec_params["L3_fwhm"].value, spec_params["L3_fwhm"].stderr)

# Frequency band between both peaks
freq = 92.0 # MHz

# Proportionality constant
freq_prop = freq / (cooling_transition - left_peak)

def detuning(t):
    ret = t - cooling_transition.n
    ret *= freq_prop.n
    return ret

def detuning_inverse(f):
    t = f
    t = t / freq_prop.n
    t += cooling_transition.n
    return t

# This describes the detuning
cooling_freq = ufloat(cool_params["center"].value, cool_params["center"].stderr)

# Calculate detuning
detuning = cooling_freq - cooling_transition
detuning *= freq_prop

# Calculate FWHM
cooling_transition_fwhm *= freq_prop


# Plot
plt.xlim(0.0, 21.2)
plt.ylim(-0.05, 0.45)

plt.xlabel(r"Zeit $t$ / \si{s}")
plt.ylabel(r"Signal / beliebige Einheiten")


plt.tick_params(axis="y", which="both", left="off", right="off", labelleft="off")


# Daten
plt.plot(sig_spec.t, sig_spec.I / 10.0 - 0.6, "+", label="Rb-Spektrum", markersize=2.0, markeredgewidth=0.5)
plt.plot(sig_fluo.t, sig_fluo.I + 0.1, "+", label="Fluoresz. MOT\n(mit Untergrund)", markersize=2.0, markeredgewidth=0.5)
plt.plot(bg_fluo.t, bg_fluo.I + 0.1, "+", label="Fluoresz.\nUntergrund", markersize=2.0, markeredgewidth=0.5)
plt.plot(fluo.t, fluo.I, "+", label="Fluoresz. MOT\n(ohne Untergrund)", markersize=2.0, markeredgewidth=0.5)

# Fits
plt.plot(fluo.t, cool_fit.best_fit, "-", label="Fluoresz. MOT\n(Anpassung)")
plt.plot(sig_spec.t, sig_spec_fit.best_fit / 10.0 - 0.6, "-", label="Rb-Spektrum\n(Anpassung)")


# Vertical line at transition
plt.axvline(x=cooling_transition.n, c="k", zorder=1)
plt.axvline(x=cooling_freq.n, c="k", zorder=1, ls=":")

plt.legend(loc="lower left", framealpha=0.85, fontsize=7, markerscale=2.0)


# Scaled axis
ax1 = plt.gca()
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())

tic_loc = detuning_inverse(np.arange(-130, 20, 10))
ax2.set_xticks(tic_loc)
ax2.set_xticklabels(np.arange(-130, 20, 10), fontsize=7)
ax2.set_xlabel(r"Verstimmung $\delta$ / \si{MHz}")


plt.tight_layout(pad=0.2)
plt.savefig("figures/detuning_cooling.pdf")
plt.close()
