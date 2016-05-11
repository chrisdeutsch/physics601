import numpy as np
import matplotlib.pyplot as plt

plt.style.use("publication")

fig = plt.gcf()
fig.set_size_inches(4.6, 5.0)

costheta = np.linspace(-1.0, 1.0, 1000)
s_channel = 1.0 + costheta**2
t_channel = (1.0 - costheta)**-2


ax1 = plt.subplot(211)
plt.plot(costheta, s_channel, "-", label="$s$-channel (total)")

plt.title(r"$\mathrm{e}^+ \mathrm{e}^- \rightarrow \mathrm{\mu}^+ \mathrm{\mu}^-$")
plt.ylim((0, 3))
plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega}$ / a.u.")
plt.setp(ax1.get_xticklabels(), visible=False)
plt.legend(loc=0)


ax2 = plt.subplot(212, sharex=ax1)
plt.plot(costheta, s_channel, "-", label="$s$-channel")
plt.plot(costheta, t_channel, "-", label="$t$-channel")
plt.plot(costheta, s_channel + t_channel, "-", label="total")

plt.title(r"$\mathrm{e}^+ \mathrm{e}^- \rightarrow \mathrm{e}^+ \mathrm{e}^-$")
plt.xlabel(r"$\cos\theta$")
plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega}$ / a.u.")
plt.ylim((0, 10))
plt.legend(loc=0)

plt.tight_layout(pad=0.5)
plt.savefig("figures/s_t_channel.pdf")
plt.close()