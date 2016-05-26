import numpy as np
import matplotlib.pyplot as plt

plt.style.use("publication")

fig = plt.gcf()
fig.set_size_inches(4.6, 2.8)

costheta = np.linspace(-1.0, 1.0, 1000)
s_channel = 1.0 + costheta**2
t_channel = (1.0 - costheta)**-2

plt.plot(costheta, s_channel, "-", label="$s$-channel")
plt.plot(costheta, t_channel, "-", label="$t$-channel")
plt.plot(costheta, s_channel + t_channel, "-", label="total")

plt.xlabel(r"$\cos\theta$")
plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega}$ / a.u.")
plt.ylim((0, 10))
plt.legend(loc=0)

plt.tight_layout(pad=0.5)
plt.savefig("talkfigs/pdf/s_t_channel.pdf")
plt.close()