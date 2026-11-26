import ROOT
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import mplhep as hep
from Draw.python.era_utilities import LUMI
from Draw.compare_eras.utils import *
from Draw.compare_eras.compare_data import get_data_subnodes
from Draw.compare_eras.compare_data import get_counts as get_data_counts
from Draw.compare_eras.compare_fakes import get_counts as get_fake_counts

plt.style.use(hep.style.CMS)


def get_args():
    parser = argparse.ArgumentParser(description="Compare data from different eras")
    parser.add_argument(
        "--config",
        type=str,
        default="Draw/compare_eras/config.json",
        help="Path to the configuration .json file",
    )
    return parser.parse_args()


def save_fig(config, hist_dict, lumi_dict, var):
    reference_era = config["reference_era"]
    year = config["year"]
    channel = config["channel"]
    output_path = config["output_path"]

    if reference_era not in hist_dict:
        raise ValueError(f"Reference era {reference_era} not in {list(hist_dict)}")

    var_name = var.split("(")[0]
    edges = next(iter(hist_dict.values()))[0]

    total_lumi_fb = sum(lumi_dict.values()) / 1000.0

    # Lumi-normalised reference and its relative stat. uncertainty
    ref_data_counts = hist_dict[reference_era][1]
    ref_fake_counts = hist_dict[reference_era][2]
    ref_counts = ref_data_counts - ref_fake_counts
    ref_lumi_fb = lumi_dict[reference_era] / 1000.0
    ref = ref_counts / ref_lumi_fb
    with np.errstate(divide="ignore", invalid="ignore"):
        ref_err = np.sqrt(ref_data_counts + ref_fake_counts)
        ref_rel_err = np.where(ref_counts > 0, ref_err / ref_counts, 0.0)

    fig, (ax, rax) = plt.subplots(
        2, 1, figsize=(10, 10), sharex=True,
        gridspec_kw={"height_ratios": (3, 1), "hspace": 0.05},
    )

    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, (era, (bin_edges, data_counts, fake_counts)) in enumerate(hist_dict.items()):
        colour = colours[i % len(colours)]
        lumi_fb = lumi_dict[era] / 1000.0
        counts = data_counts - fake_counts
        norm = counts / lumi_fb
        err = np.sqrt(data_counts + fake_counts) / lumi_fb
        label = f"{era} ({lumi_fb:.1f} fb$^{{-1}}$)"
        if era == reference_era:
            label += " [ref]"
        hep.histplot(norm, bin_edges, yerr=err, ax=ax, histtype="step", color=colour, label=label)

        if era == reference_era:
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(ref > 0, norm / ref, np.nan)
            ratio_err = np.where(ref > 0, err / ref, np.nan)
        hep.histplot(ratio, bin_edges, yerr=ratio_err, ax=rax, histtype="errorbar", marker='_', markersize=4, color=colour)

    # Reference era's stat. uncertainty shown as a grey band around 1
    rax.stairs(1 + ref_rel_err, edges, baseline=1 - ref_rel_err, fill=True,
               color="grey", alpha=0.3, linewidth=0)
    rax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    rax.set_ylabel(f"Era / {reference_era}", fontsize=20)
    rax.set_ylim(*RATIO_RANGE)
    # Drop the top and bottom tick labels (avoids clash with upper panel, keeps symmetry)
    rax.yaxis.set_major_locator(MaxNLocator(nbins=5, prune="both"))
    rax.set_xlabel(VAR_LABELS.get(var_name, var_name))
    rax.set_xlim(edges[0], edges[-1])

    hep.cms.label("Preliminary", data=True, ax=ax, com=13.6, year=year.split("_")[-1],
                  lumi=f"{total_lumi_fb:.1f}")
    ax.text(0.04, 0.90, f"{CHANNEL_LABELS.get(channel, channel)}",
            transform=ax.transAxes, fontsize=20)

    ax.set_ylabel(rf"Data $-$ Jet Fakes Est. / fb$^{{-1}}$")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
    ax.legend(ncol=1, fontsize=16, title_fontsize=16, loc="upper right")

    for ext in ("png", "pdf"):
        fig.savefig(output_path.rsplit(".", 1)[0] + f".{ext}", bbox_inches="tight")
    plt.close(fig)


def main():
    args = get_args()
    with open(args.config, "r") as f:
        config = json.load(f)

    base_path = config["input_path"]
    file, datacard_name, datasets = get_data_subnodes(base_path)
    var = config["variable_definitions"][config["variable"]]
    hist_dict = {}

    for era, runs in config["era_definitions"].items():
        bin_edges, data_counts = get_data_counts(file, datacard_name, datasets, runs, var)
        bin_edges, fake_counts = get_fake_counts(base_path, era, var)
        hist_dict[era] = (bin_edges, data_counts, fake_counts)

    lumi_by_run = LUMI[config["year"]]
    lumi_dict = build_lumi_dict(config["era_definitions"], lumi_by_run)
    save_fig(config, hist_dict, lumi_dict, var)


if __name__ == "__main__":
    main()

