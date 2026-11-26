import os
import ROOT
import numpy as np
from Draw.compare_eras.utils import build_lumi_dict
from Draw.python.era_utilities import LUMI
from prettytable import PrettyTable, TableStyle

BASE = "/vols/cms/dmw25/HiggsDNA/output/WiktorFFs"
YEAR = "Run3_2024"
CHANNEL = "tt"
ERAS = {
    "2024CDE": ["Run2024C", "Run2024D", "Run2024E"],
    "2024F": ["Run2024F"],
    "2024G": ["Run2024G"],
    "2024HI": ["Run2024H", "Run2024I"],
}
WPS = {"Loose": 4, "Medium": 5, "Tight": 6, "VTight": 7}


def get_datasets(base, year, channel, eras):
    datasets = {}
    samples = os.listdir(os.path.join(base, year, channel))

    for era, runs in eras.items():
        datasets[era] = [
            os.path.join(
                base, year, channel, sample
            ) for sample in samples if any(run in sample for run in runs)
        ]

    return datasets


def get_events_passing_cuts(dataset, cut_strings):
    print(f"\033[92mProcessing dataset:\033[0m {dataset}")

    file_path = os.path.join(dataset, "nominal", "merged.root")
    file = ROOT.TFile(file_path)
    df = ROOT.RDataFrame("ntuple", file)

    counts = []
    for cut_string in cut_strings:
        print(f"Applying cut: {cut_string}")
        filtered_df = df.Filter(cut_string)
        counts.append(filtered_df.Count().GetValue())
    return np.array(counts)


def make_table(title, eras, values, errors, fmt):
    table = PrettyTable(["Era", *WPS])
    table.title = title
    table.set_style(TableStyle.SINGLE_BORDER)
    table.align = "r"
    table.align["Era"] = "l"
    for era, vals, errs in zip(eras, values, errors):
        table.add_row([era, *(f"{v:{fmt}} ± {e:{fmt}}" for v, e in zip(vals, errs))])
    return table


def main():
    datasets = get_datasets(BASE, YEAR, CHANNEL, ERAS)
    lumi_by_run = LUMI[YEAR]
    lumi_dict = build_lumi_dict(ERAS, lumi_by_run)

    cut_strings = [
        f"(trg_doubletau || trg_doubletauandjet) && (idPNetVSjet_1 >= {wp} && idPNetVSjet_2 >= {wp})"
        for wp in WPS.values()
    ]

    eras = list(datasets)
    raw_counts = np.zeros((len(eras), len(cut_strings)))
    for i, era in enumerate(eras):
        for dataset in datasets[era]:
            raw_counts[i] += get_events_passing_cuts(dataset, cut_strings)

    # Events per unit lumi, with Poisson uncertainty sqrt(N) / lumi
    lumis = np.array([lumi_dict[era] for era in eras])[:, None]
    rate = raw_counts / lumis
    rate_err = np.sqrt(raw_counts) / lumis

    # Ratio to the first era; uncorrelated Poisson errors on numerator and
    # denominator, except the reference row, which is 1 by construction
    ratio = rate / rate[0]
    ratio_err = ratio * np.sqrt(1 / raw_counts + 1 / raw_counts[0])
    ratio_err[0] = 0.0

    print()
    print(make_table(
        "Events / lumi [pb]: trigger + both taus pass PNet VSjet WP",
        eras, rate, rate_err, ".2f",
    ))
    print(make_table(
        f"Ratio to {eras[0]}",
        eras, ratio, ratio_err, ".3f",
    ))


if __name__ == "__main__":
    main()
