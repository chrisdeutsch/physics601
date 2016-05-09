import numpy as np
import matplotlib.pyplot as plt

plt.style.use("publication")

costheta = np.linspace(-1.0, 1.0, 1000)

plt.plot(costheta, 1.0 + costheta**2, "-", label="$s$-channel")
plt.plot(costheta, (1.0 - costheta)**-2, "-", label="$t$-channel")


plt.ylim((0, 5))
plt.setp(plt.gca().get_yticklabels(), visible=False)
plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega}$ / a.u.")
plt.xlabel(r"$\cos\theta_\mathrm{W}$")
plt.legend(loc=0)

plt.tight_layout(pad=0.5)
plt.savefig("figures/s_t_channel.pdf")
plt.close()