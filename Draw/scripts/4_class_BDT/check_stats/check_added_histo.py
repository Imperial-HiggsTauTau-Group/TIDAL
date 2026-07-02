import ROOT
import argparse

samples_of_interest = [
    "ggH_sm_prod_sm_htt125",
    "qqH_sm_htt125",
]

directories_of_interest = [
    "tt_ggH_rhorho",
    "tt_VBF_rhorho",
]


def get_final_five_bins(hist):
    # Get the number of bins in the histogram
    n_bins = hist.GetNbinsX()
    
    # Get the content and error of the last five bins
    bin_contents = [hist.GetBinContent(i) for i in range(n_bins - 4, n_bins + 1)]
    bin_abs_errors = [hist.GetBinError(i) for i in range(n_bins - 4, n_bins + 1)]
    bin_rel_errors = [bin_abs_errors[i] / bin_contents[i] if bin_contents[i] != 0 else 0 for i in range(5)]
   
    return bin_rel_errors, bin_contents


def main(args):

    for j, histogram_name in enumerate(directories_of_interest):
        sample = samples_of_interest[j]
        file = ROOT.TFile(args.file, "READ")
        histogram = file.Get(histogram_name)
        hist = histogram.Get(sample)
        if hist:
            error, content = get_final_five_bins(hist)
            print(f"{histogram_name} - {sample}:")
            for i, (rel_error, cont) in enumerate(zip(error, content)):
                print(f"  Bin {i}: Relative Error = {rel_error:.2%}, Content = {cont:.2f}")
        else:
            print(f"Histogram for {sample} not found in the file.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="Path to the input ROOT file")
    parser.add_argument("--title", type=str, default="Final Bin Errors", help="Title for the output table")
    args = parser.parse_args()
    main(args)