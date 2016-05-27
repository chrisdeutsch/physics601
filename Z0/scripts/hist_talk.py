import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("publication")

files = ["ELECTRONS", "MUONS", "TAUS", "HADRONS"]
data = dict((f, pd.read_csv("./data/tag1/" + f + ".csv")) for f in files)

variables = [("Ctrk(Sump)", "Total momentum of charged tracks / GeV", 10), ("Ctrk(N)", "Number of charged tracks", 15), ("Ecal(SumE)", "Energy in e.m. calorimeter / GeV", 10), ("Hcal(SumE)", "Energy in hadron calorimeter / GeV", 10)]

for var, detail, b in variables:
    fig = plt.figure(figsize=(4.05, 2.7))
    sets = [data[key][var] for key in data]
    plt.hist(sets, bins=b, label=[key.lower() for key in data], edgecolor=None, linewidth=0.5)
    plt.xlabel(format(detail))
    plt.legend(loc=0, framealpha=0.5)
    plt.minorticks_on()
    plt.tight_layout(pad=0.5)
    plt.savefig('./talkfigs/pdf/'+str(var)+'.pdf')
    plt.close()