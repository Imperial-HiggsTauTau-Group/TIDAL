import ROOT
import argparse
from prettytable import PrettyTable

overview = [
    "Higgs_flat_htt125",
    "ZTT",
    "JetFakes",
]

signal = [
    "ggH_flat_prod_sm_htt125",
    "qqH_flat_htt125",
]


def describe_datacard(obj: ROOT.TH1):
    """
    Returns the root mean square of the errors the 4 BDT score windows.
    """
    N_bins = obj.GetNbinsX()
    errors = []
    for i in range(1, N_bins + 1):
        errors.append(
            obj.GetBinError(i) / obj.GetBinContent(i) if obj.GetBinContent(i) != 0 else 0
        )
    N_per_window = N_bins // 4
    RMS_per_window = []
    for i in range(4):
        start_bin = i * N_per_window + 1
        end_bin = (i + 1) * N_per_window
        window_errors = errors[start_bin - 1 : end_bin]
        RMS_per_window.append(
            (sum(e**2 for e in window_errors) / N_per_window)**0.5
        )
    return RMS_per_window


def summarise(directory: ROOT.TDirectory,
              of_interest: list[str]
) -> dict[str, list[float]]:
    """
    Summarises the statistics of the variables in the given TDirectory.
    """
    errors = {}
    for key in directory.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        if name in of_interest:
            errors[name] = describe_datacard(obj)
    return errors


def print_error_table(args, of_interest):
    # Open the ROOT file
    file = ROOT.TFile(args.file)

    # List all the TDirectories in the file
    all_TDirectories = [key.GetName() for key in file.GetListOfKeys() if key.IsFolder()]
    #print("TDirectories in the file:", all_TDirectories)

    # Remove directories we are not interested in
    TDirectories = all_TDirectories.copy()
    for dir_name in all_TDirectories:
        if dir_name.startswith("tt_mva_") or dir_name.endswith("_aiso"):
            TDirectories.remove(dir_name)
    #print("TDirectories after filtering:", TDirectories)

    # Create summary table
    table = PrettyTable()
    table.field_names = [
        "TDirectory", "MC Sample", "BDT Window 1",
        "BDT Window 2", "BDT Window 3", "BDT Window 4"
    ]
    table.title = f"Root mean square of errors for each BDT score window ({args.era})"
    for dir_name in TDirectories:
        directory = file.Get(dir_name)
        errors = summarise(directory, of_interest)
        i = 0
        for sample in of_interest:
            table.add_row([
                dir_name if i == 0 else "",
                sample,
                f"{errors[sample][0]:.4f}",
                f"{errors[sample][1]:.4f}",
                f"{errors[sample][2]:.4f}",
                f"{errors[sample][3]:.4f}",
            ])
            i += 1
        table.add_divider()
    print(table)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check statistics of variables in a ROOT file.")
    parser.add_argument("--file", type=str, required=True, help="Path to the ROOT file")
    parser.add_argument("--era", type=str, required=True, help="Era of the data (e.g., 2016, 2017, 2018)")
    parser.add_argument("--overview", action="store_true", help="Print statistics for overview samples")
    parser.add_argument("--signal", action="store_true", help="Print statistics for signal samples")
    args = parser.parse_args()

    if args.overview:
        print_error_table(args, overview)
    if args.signal:
        print_error_table(args, signal)

