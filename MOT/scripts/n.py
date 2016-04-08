from uncertainties import ufloat
import numpy as np

P1 = ufloat(3.0, 0.1) # mW
P2 = ufloat(3.0, 0.1) # mW
P3 = ufloat(3.2, 0.1) # mW

w = ufloat(2.826, 0.089) # mm

P = 2 * (P1 + P2 + P3) # mW

detuning = -ufloat(-11.6, 0.4) # MHz
gamma = ufloat(38.117, 0.011) # MHz
Isat = ufloat(1.66932, 0.00052) # mW / cm^2

I = P / (np.pi * w**2) * 100.0 # mW / cm^2
s = I / Isat

R = gamma / 2.0 * s / (1 + s + 4 * detuning**2 / gamma**2)

Ptot = ufloat(30.9, 3.1) / 1000.0 # mW
trans_energy = 1.589 * 1.6e-19 * 1000.0

N = Ptot / (trans_energy * R * 1e6)

