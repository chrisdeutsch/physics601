import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("publication")

files = ["ELECTRONS", "MUONS", "TAUS", "HADRONS"]
data = dict((f, pd.read_csv("./data/tag1/" + f + ".csv")) for f in files)

variables = [("Ctrk(Sump)", 10), ("Ctrk(N)", 10), ("Ecal(SumE)", 10), ("Hcal(SumE)", 10)]

for var, b in variables:
    fig = plt.figure(figsize=(4.05, 3.4))
    sets = [data[key][var] for key in data]
    plt.hist(sets, bins=b, label=[key.lower() for key in data])
    plt.xlabel(var)
    plt.legend(loc=0, framealpha=0.5)
    plt.minorticks_on()
    plt.tight_layout(pad=0.5)
    plt.savefig('./figures/event_display/'+str(var)+'.pdf')
