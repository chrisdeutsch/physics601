import numpy as np
from uncertainties import ufloat
import uncertainties.umath as umath


kb = 1.38064852e-23 # J / K
m_rb = 84.911794 * 1.660539040e-27 # kg

rb_radius = 248.0e-12 # m
geom_cross_sec = np.pi * rb_radius**2


tau = ufloat(99.6, 4.9) / 1000.0 # s
p = ufloat(1.1e-7, 0.2e-7) * 100.0 # Pa
T = ufloat(24.6, 0.5) + 273.15 # K


n = p / (kb * T)
v = umath.sqrt(8.0 * kb * T / (np.pi * m_rb))

cross_sec = 1.0 / (n * v * tau)