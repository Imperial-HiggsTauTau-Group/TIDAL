import json
import ROOT
import numpy as np

BASE = "Draw/plots/a1-investigation/era-investigations/derive_residual_corr"
ERAS = ["CDE", "F", "G", "HI", "inclusive"]
CATS = ["sel_mupi", "sel_murho", "sel_mua11pr", "sel_mua1"]
VAR = "pt_2[0.0,35.0,40.0,50.0,200.0]"


def get_path_and_name(base, era, cat):
    path = f"{base}/{era}/Run3_2024/control/mt/datacard_pt_2_mTLt65_{cat}_mt_Run3_2024.root"
    name = f"mt_{cat}_mTLt65"
    return path, name


def load_histograms(path, name):
    file = ROOT.TFile(path)
    hist = file.Get(name)
    obs = hist.Get("data_obs")
    exp = hist.Get("total_bkg")
    dy = hist.Get("ZTT")
    return file, obs, exp, dy


def get_bin_edges(var):
    if "(" in var:
        n, low, high = var.split("(")[1][:-1].split(",")
        n = int(n)
        low = float(low)
        high = float(high)
        return np.linspace(low, high, n + 1)
    elif "[" in var:
        edges = var.split("[")[1][:-1].split(",")
        return [float(edge) for edge in edges]


def derive_corrections(base, eras, cats, var):
    corrections = {}
    for era in eras:
        corrections[era] = {}
        for cat in cats:
            corrections_by_bin = {}
            path, name = get_path_and_name(base, era, cat)
            file, obs, exp, dy = load_histograms(path, name)
            bin_edges = get_bin_edges(var)
            for i in range(len(bin_edges) - 1):
                obs_bin = obs.GetBinContent(i + 1)
                exp_bin = exp.GetBinContent(i + 1)
                dy_bin = dy.GetBinContent(i + 1)
                if dy_bin > 0:
                    correction = 1.0 + (obs_bin - exp_bin) / dy_bin
                else:
                    correction = 1.0
                corrections_by_bin[f"{int(bin_edges[i])}-{int(bin_edges[i + 1])}"] = correction
            corrections[era][cat] = corrections_by_bin
            print(f"\033[32mDerived corrections for era '{era}', category {cat}\033[0m")
    return corrections


def main():
    corrections = derive_corrections(BASE, ERAS, CATS, VAR)
    with open("Draw/residual_correcting/corrections.json", "w") as f:
        json.dump(corrections, f, indent=4)


if __name__ == "__main__":
    main()

