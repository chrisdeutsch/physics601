import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("publication")


# Natrium links
data = pd.read_csv("data/spectrum/Na_links.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)

plt.errorbar(data.index, data.cnt, yerr=data.dcnt, fmt="|")

plt.annotate(r"$1275 \, \si{\keV}$", xy=(5913, 133), xytext=(5913, 300),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$511 \, \si{\keV}$", xy=(2250, 960), xytext=(1400, 1100),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.xlim(0, 7000)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/na_links.pdf")
plt.close()


# Natrium rechts
data = pd.read_csv("data/spectrum/Na_rechts.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)

plt.errorbar(data.index, data.cnt, yerr=data.dcnt, fmt="|")

plt.annotate(r"$1275 \, \si{\keV}$", xy=(6730, 120), xytext=(6730, 290),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$511 \, \si{\keV}$", xy=(2420, 820), xytext=(1400, 900),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.xlim(0, 8000)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/na_rechts.pdf")
plt.close()


# LYSO links
data = pd.read_csv("data/spectrum/LYSO_links.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)

plt.errorbar(data.index, data.cnt, yerr=data.dcnt, fmt="|")

plt.annotate(r"$597 \, \si{\keV}$", xy=(3000, 900), xytext=(4000, 930),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$290 \, \si{\keV}$", xy=(1400, 380), xytext=(1400, 600),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$395 \, \si{\keV}$", xy=(1900, 340), xytext=(1900, 500),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$88 \, \si{\keV}$", xy=(480, 80), xytext=(480, 200),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$202 \, \si{\keV}$", xy=(922, 134), xytext=(922, 300),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.xlim(0, 6000)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/LYSO_links.pdf")
plt.close()


### LYSO rechts
data = pd.read_csv("data/spectrum/LYSO_rechts.txt", sep="\t", index_col=0,
                   names=["ch", "cnt"])
data["dcnt"] = np.sqrt(data.cnt)

plt.errorbar(data.index, data.cnt, yerr=data.dcnt, fmt="|")

plt.annotate(r"$597 \, \si{\keV}$", xy=(3200, 780), xytext=(4400, 800),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$290 \, \si{\keV}$", xy=(1420, 340), xytext=(1420, 500),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$395 \, \si{\keV}$", xy=(1980, 300), xytext=(1980, 420),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$88 \, \si{\keV}$", xy=(420, 75), xytext=(420, 220),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.annotate(r"$202 \, \si{\keV}$", xy=(886, 110), xytext=(886, 300),
             arrowprops=dict(arrowstyle="->", facecolor="black"), ha="center")

plt.xlim(0, 7000)
plt.xlabel("Kanal")
plt.ylabel("Ereignisse~$N$")

plt.savefig("figures/LYSO_rechts.pdf")
plt.close()
