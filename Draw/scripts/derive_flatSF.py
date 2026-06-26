import uproot
import numpy as np
import argparse
import matplotlib.pyplot as plt
import mplhep as hep
import os

plt.style.use(hep.style.CMS)
plt.rcParams.update({"font.size": 16})


parser = argparse.ArgumentParser()
parser.add_argument("--file", type=str, required=True, help="Input ROOT file")
parser.add_argument("--label", type=str, required=True, help="label")
args = parser.parse_args()

node = "mt_cp_inclusive_mTLt65"
subtract = ["ZL", "W", "VVJ", "QCD"]
add = ["ZTT", "TTT", "VVT"]

f = uproot.open(args.file)

vals, bins = f[f"{node}/data_obs"].to_numpy()
centers = 0.5 * (bins[:-1] + bins[1:])

bkg_vals = sum(f[f"{node}/{s}"].values() for s in subtract)
mc_vals = sum(f[f"{node}/{s}"].values() for s in add)
subtracted = vals - bkg_vals

mask = (centers >= 40) & (centers <= 80)
print(f"Data - non genuine : {subtracted.sum():.2f}")
print(f"Genuine: {mc_vals.sum():.2f}")
print(f"Ratio (40-80) : {subtracted[mask].sum()/mc_vals[mask].sum():.4f}")

ratio = np.divide(subtracted, mc_vals, out=np.ones_like(subtracted), where=mc_vals != 0)
ratio_err = np.divide(np.sqrt(vals), mc_vals, out=np.zeros_like(subtracted), where=mc_vals != 0)

fig, (ax, rax) = plt.subplots(2, 1, figsize=(6, 6), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
fig.subplots_adjust(hspace=0.05)

ax.errorbar(centers, subtracted, yerr=np.sqrt(vals), fmt="ko", label="Data - Fake")
ax.stairs(mc_vals, bins, color="orange", label="Total Genuine")
ax.set_ylabel("Events")
ax.set_xlim(0,200)
ax.legend()
ratio_40_80 = subtracted[mask].sum() / mc_vals[mask].sum()
ax.text(0.95, 0.8, f"Correction (40-80 GeV): {ratio_40_80:.4f}", transform=ax.transAxes,
        ha="right", va="top", fontweight="bold", fontsize=12)

ax.text(0.95, 0.9, args.label, transform=ax.transAxes,
        ha="right", va="top", fontweight="bold", fontsize=14)

rax.errorbar(centers, ratio, yerr=ratio_err, fmt="ko")
rax.axhline(1, color="orange", linestyle="--")
rax.set_xlabel(r"$m_{vis}$")
rax.set_ylabel("Ratio")
rax.set_ylim(0.5, 1.5)

lumi = 109.08 if "2024" in args.file else 110.58
hep.cms.label(ax=ax, text="Work in progress", data=True, lumi=lumi, com=13.6, fontsize=16)

fig.savefig(args.file.replace('.root', "_correction.pdf"), bbox_inches="tight")
plt.close(fig)
