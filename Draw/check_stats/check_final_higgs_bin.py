import ROOT
import argparse
from prettytable import PrettyTable

of_interest = [
    "ZTT_0J",
    "ZTT_1J",
    "ZTT_2J"
]


def get_final_bin_error(hist):
    # Get the number of bins in the histogram
    n_bins = hist.GetNbinsX()
    
    # Get the content and error of the last bin
    last_bin_content = hist.GetBinContent(n_bins)
    last_bin_error = hist.GetBinError(n_bins)
    
    return last_bin_error / last_bin_content, last_bin_content


def main(args):
    table = PrettyTable()
    table.title = args.title
    table.field_names = ["Sample", "Final Bin Error", "Final Bin Content"]

    file = ROOT.TFile(args.file, "READ")
    directory = file.Get('tt_mva_higgs')
    
    for sample in of_interest:
        hist = directory.Get(sample)
        if hist:
            error, content = get_final_bin_error(hist)
            table.add_row([sample, f"{error:.4f}", f"{content:.4f}"])
        else:
            print(f"Histogram for {sample} not found in the file.")
    
    print(table)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="Path to the input ROOT file")
    parser.add_argument("--title", type=str, default="Final Bin Errors", help="Title for the output table")
    args = parser.parse_args()
    main(args)