import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Model
from scripts.tools import set_mpl_comma

plt.style.use("publication")

def charge(t, t0, amp, tau, offset):
    charging = t > t0
    ret = np.zeros(len(t))
    ret[charging] = amp * (1.0 - np.exp(- (t[charging] - t0) / tau))
    return ret + offset

files = glob.glob("data/loading/*/*.csv")

fit_params = {"t0": [], "dt0": [], "amp": [], "damp": [], "tau": [], "dtau": [],
              "offset": [], "doffset": []}

for filename, i in zip(files, range(len(files))):
    data = pd.read_csv(filename, usecols=[3, 4], header=None, names=["t", "I"])
    data.I *= 1000
    
    # Cut data
    data = data[data.I > 2.8]
    cut = data.t > 0.5 * data.t.max()
    cut &= data.I < 0.5 * (data.I.max() + data.I.min())
    cut = ~cut
    data = data[cut]
    
    # Fit
    mod = Model(charge)
    fit_res = mod.fit(data.I.get_values(), t=data.t.get_values(), 
                      t0=0.2, amp=30.0, tau=0.2, offset=30.0)
    
    # Save fit params
    params = fit_res.params
    for pname in params:
        p = params[pname]
        fit_params[pname].append(p.value)
        fit_params["d" + pname].append(p.stderr)
    
    # Plots
    set_mpl_comma()
    plt.xlim(0.0, 1.0)
    plt.xlabel(r"Zeit $t$ / \si{\second}")
    plt.ylabel(r"Photodiodensignal $U$ / \si{mV}")
    
    x = np.linspace(0.0, 1.0, 1000)
    plt.plot(data.t, data.I, "-", label="Messwerte")
    plt.plot(x, fit_res.eval(t=x), "-", label="Anpassung")
    
    plt.legend(loc="lower right")
    
    plt.tight_layout()
    plt.savefig("figures/loading/loading{}.pdf".format(i))
    plt.close()

fit_params = pd.DataFrame(fit_params)
fit_params = fit_params[fit_params.index > 0]

# Weighted mean of tau
fit_params["w"] = fit_params.dtau**-2

var_tau = 1.0 / fit_params.w.sum()
mean_tau = (fit_params.tau * fit_params.w).sum() * var_tau
stderr_tau = np.sqrt(var_tau)

print("tau = {} +- {}".format(mean_tau, stderr_tau))

# Latex Tabelle
out = fit_params[["amp", "damp", "tau", "dtau", "t0", "dt0", "offset", "doffset"]]

# to ms
out[["tau", "dtau", "t0", "dt0"]] *= 1000.0

out.columns = [r"{$U_0$ / \si{mV}}", r"{$\Delta U_0$ / \si{mV}}",
               r"{$\tau$ / \si{ms}}", r"{$\Delta \tau$ / \si{ms}}",
               r"{$t_0$ / \si{ms}}", r"{$\Delta t_0$ / \si{ms}}",
               r"{$\mathrm{BG}$ / \si{mV}}", r"{$\Delta \mathrm{BG}$ / \si{mV}}"]

from scripts.tools import round
out.to_latex("tables/loading.tex", index=False,
             formatters=[round(2), round(2),
                         round(2), round(2),
                         round(2), round(2),
                         round(2), round(2)], 
             column_format="S[table-format=2.2]S[table-format=1.2]S[table-format=3.2]S[table-format=1.2]S[table-format=3.2]S[table-format=1.2]S[table-format=2.2]S[table-format=1.2]",
             escape=False)