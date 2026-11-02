import ROOT
import numpy as np
from pprint import pprint

PATH = "Draw/plots/a1-investigation/pre-meeting/Run3_2024/control/mt/datacard_pt_2_mTLt65_sel_mua1_mt_Run3_2024.root"
NAME = "mt_sel_mua1_mTLt65"
VAR = "pt_2(40,0,200)"


def load_histograms(path, name):
    file = ROOT.TFile(path)
    hist = file.Get(name)
    obs = hist.Get("data_obs")
    exp = hist.Get("total_bkg")
    dy = hist.Get("ZTT")
    return file, obs, exp, dy


def get_bin_edges(var):
    n, low, high = var.split("(")[1][:-1].split(",")
    n = int(n)
    low = float(low)
    high = float(high)
    return np.linspace(low, high, n + 1)


file, obs, exp, dy = load_histograms(PATH, NAME)
bin_edges = get_bin_edges(VAR)
corrections = []
for i in range(len(bin_edges) - 1):
    obs_bin = obs.GetBinContent(i + 1)
    exp_bin = exp.GetBinContent(i + 1)
    dy_bin = dy.GetBinContent(i + 1)
    if dy_bin > 0:
        corrections.append(1 + (obs_bin - exp_bin) / dy_bin)
    else:
        corrections.append(1)

compressed_corrections = {
    "20-25": corrections[4],
    "25-30": corrections[5],
    "30-35": corrections[6],
    "35-40": corrections[7],
    "40-45": corrections[8],
    "45-50": corrections[9],
}

print("======================================================================")
print("Detailed corrections:")
for i in range(len(bin_edges) - 1):
    print(f"{int(bin_edges[i])}—{int(bin_edges[i + 1])} GeV: {corrections[i]}")
print("======================================================================")
print("Compressed corrections:")
pprint(compressed_corrections)
