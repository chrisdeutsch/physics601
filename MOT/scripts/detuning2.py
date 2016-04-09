import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit.models import LorentzianModel, PolynomialModel, ConstantModel, GaussianModel
from uncertainties import ufloat, umath

plt.style.use("publication")

# Fit Model
bg = PolynomialModel(4, prefix="bg_")
L1 = LorentzianModel(prefix="L1_")
L2 = LorentzianModel(prefix="L2_")

model = L1 + L2 + bg


### SIGNAL:

# Spektroskopie
sig_spec = pd.read_csv("data/detuning/repumping/signal/F0097CH1.csv",
                  usecols=[3, 4], header=None, names=["t", "I"])

# Fluoreszenz
sig_fluo = pd.read_csv("data/detuning/repumping/signal/F0097CH2.csv",
                  usecols=[3, 4], header=None, names=["t", "I"])

# Cut data
sig_spec = sig_spec[sig_spec.t < 22.3]
sig_fluo = sig_fluo[sig_fluo.t < 22.3]


### Background:

# Spektroskopie
bg_spec = pd.read_csv("data/detuning/repumping/background/F0098CH1.csv",
                  usecols=[3, 4], header=None, names=["t", "I"])

# Fluoreszenz
bg_fluo = pd.read_csv("data/detuning/repumping/background/F0098CH2.csv",
                  usecols=[3, 4], header=None, names=["t", "I"])

# Cut data
bg_spec = bg_spec[bg_spec.t < 21.8]
bg_fluo = bg_fluo[bg_fluo.t < 21.8]

bg_spec.t += 2.1
bg_fluo.t += 2.1

# Find matching range
low = max(bg_spec.t.min(), sig_spec.t.min())
high = min(bg_spec.t.max(), sig_spec.t.max())

# Cut to range
sig_spec = sig_spec[(sig_spec.t > low) & (sig_spec.t < high)].reset_index()
sig_fluo = sig_fluo[(sig_fluo.t > low) & (sig_fluo.t < high)].reset_index()
bg_spec = bg_spec[(bg_spec.t > low) & (bg_spec.t < high)].reset_index()
bg_fluo = bg_fluo[(bg_fluo.t > low) & (bg_fluo.t < high)].reset_index()


# Remove background
fluo = pd.DataFrame()
fluo["I"] = sig_fluo.I - bg_fluo.I
fluo["t"] = sig_fluo.t


t = fluo.t.get_values()
I = fluo.I.get_values()
fluo_fit = model.fit(I, x=t,
                     L1_amplitude=0.3, L1_sigma=2.0, L1_center=6.9,
                     L2_amplitude=0.4, L2_sigma=3.0, L2_center=14.4,
                     bg_c0=0.0, bg_c1=0.0, bg_c2=0.0, bg_c3=0.0, bg_c4=0.0)


# PSEUDOSIGNAL:
lorentzian = L1.func

def gaussian_bg(t, mu, sigma):
    return -np.exp(-(t-mu)**2 / sigma**2)

t = np.linspace(low, high, len(fluo.t))
I = lorentzian(t, amplitude=0.13, center=7.04, sigma=1.9) + lorentzian(t, amplitude=0.2, center=10.7, sigma=1.9)

I += gaussian_bg(t, 5.0, 100.0)



# Add noise
noise = np.random.normal(0.0, 0.003, len(t))
I += noise

# Finite voltage steps
stepwidth = 0.002
I = I // stepwidth * stepwidth

spec_fit = model.fit(I, x=t,
                     L1_amplitude=0.25, L1_sigma=2.0, L1_center=7.04,
                     L2_amplitude=0.55, L2_sigma=2.0, L2_center=10.7,
                     bg_c0=0.0, bg_c1=0.0, bg_c2=0.0, bg_c3=0.0, bg_c4=0.0)


# Plot
plt.xlim(2.1, 22.3)
plt.ylim(0.0, 0.28)

plt.xlabel(r"Zeit $t$ / \si{s}")
plt.ylabel(r"Photodiodensignal / beliebige Einheiten")

plt.tick_params(axis="y", which="both", left="off", right="off", labelleft="off")

# Daten
plt.plot(t, I + 1.22, "+", label="Rb-Spektrum", markersize=2.0, markeredgewidth=0.5)
plt.plot(sig_fluo.t, sig_fluo.I + 0.08, "+", label="Fluoresz. MOT\n(mit Untergrund)", markersize=2.0, markeredgewidth=0.5)
plt.plot(bg_fluo.t, bg_fluo.I + 0.08, "+", label="Fluoresz.\nUntergrund", markersize=2.0, markeredgewidth=0.5)
plt.plot(fluo.t, fluo.I, "+", label="Fluoresz. MOT\n(ohne Untergrund)", markersize=2.0, markeredgewidth=0.5)

# Fits
plt.plot(fluo.t, fluo_fit.best_fit, "-", label="Fluoresz. MOT\n(Anpassung)")
plt.plot(t, spec_fit.best_fit + 1.22, "-", label="Rb-Spektrum\n(Anpassung)")

plt.legend(loc="upper right", framealpha=0.85, fontsize=7, markerscale=2.0)


plt.tight_layout(pad=0.2)
plt.savefig("figures/detuning_repumping.pdf")
plt.close()
