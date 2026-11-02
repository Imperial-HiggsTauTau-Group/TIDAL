import ROOT
import numpy as np
import matplotlib.pyplot as plt

BASE = "Draw/plots/a1-investigation/post-post-meeting"
NAME = "tt_cp_problem_region"
VAR = "m_vis(50,0,400)"
SUBDIRS = ["3-PNet-Medium", "3-PNet-Tight", "3-PNet-VTight",
           "4-DeepTau-Medium", "4-DeepTau-Tight", "4-DeepTau-VTight"]


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


def get_data_minus_fakes(base, subdir, name, var):
    path = f"{base}/{subdir}/Run3_2024/control/tt/datacard_m_vis_cp_problem_region_tt_Run3_2024.root"
    file, obs, exp, dy = load_histograms(path, name)
    bin_edges = get_bin_edges(var)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    obs_minus_fakes_over_dy = []
    for i in range(len(bin_edges) - 1):
        obs_bin = obs.GetBinContent(i + 1)
        exp_bin = exp.GetBinContent(i + 1)
        dy_bin = dy.GetBinContent(i + 1)
        if dy_bin > 0:
            obs_minus_fakes_over_dy.append(1 + (obs_bin - exp_bin) / dy_bin)
        else:
            obs_minus_fakes_over_dy.append(1)
    return bin_centers, obs_minus_fakes_over_dy


plt.figure(figsize=(9, 5))
for i, subdir in enumerate(SUBDIRS):
    bin_centers, obs_minus_fakes_over_dy = get_data_minus_fakes(BASE, subdir, NAME, VAR)
    plt.subplot(2, 3, i + 1)
    plt.scatter(bin_centers, obs_minus_fakes_over_dy, label=subdir, marker='+')
    plt.xlabel("m_vis [GeV]")
    plt.ylabel("(Data - Fakes) / ZTT")
    plt.title(subdir[2:])
    plt.xlim(50, 125)
    plt.ylim(0.5, 1.5)
plt.tight_layout()
plt.savefig("Draw/a1_investigation/check_data_minus_fakes.png")

