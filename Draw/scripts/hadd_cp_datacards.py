import ROOT
import argparse
import os
import yaml
from tqdm import tqdm
from Draw.python.PlotHistograms import HTT_Histogram
from Draw.scripts.makeDatacards import create_bins


def find_variable(dir_name, channel, config):
    category = dir_name[3:] # remove the channel prefix
    if dir_name.endswith('_aiso'):
        category = category[:-5] # remove the _aiso suffix
    
    plotting_variable = None
    for entry in config['cpdecay'][channel]:
        if entry['category'][0] == category:
            plotting_variable = entry['plotting_variable'][0]
            break
    
    try:
        var = config['cpdecay']['variable_definitions'][plotting_variable]
        return create_bins(var)
    except KeyError:
        return 'Bin number'


def find_config(input_files):
    channel_dir = os.path.dirname(input_files[0])
    scheme_dir = os.path.dirname(channel_dir)
    year_dir = os.path.dirname(scheme_dir)
    output_dir = os.path.dirname(year_dir)
    ls = os.listdir(output_dir)
    for f in ls:
        if f.endswith('.yaml') and f.startswith('cpdecay_datacards'):
            return os.path.join(output_dir, f)
    return 'Draw/scripts/cpdecay_datacards.yaml'  # default if not found


def hadd_root_files(
        input_files,
        output_file,
        dir_combinations,
        channel,
        config,
        exp_num=None,
):
    """
    Combine ROOT files and merge specific directories by adding their histograms.
    
    Args:
        input_files (list of str): List of input ROOT files to combine.
        output_file (str): Name of the output ROOT file.
        dir_combinations (dict): Dictionary where keys are new directory names and values are lists
                                 of directories to combine.
                                 E.g., {'higgs_pirho': ['tt_higgs_pirho', 'tt_higgs_rhopi']}
    """
    # Create an empty output ROOT file
    output = ROOT.TFile(output_file, "RECREATE")
    print("Created output file:", output_file)
    # First, create a map of (hist_name, dir_name, file_name) to histogram
    histogram_dir_file_map = {}

    for file_name in tqdm(input_files):

        input_root = ROOT.TFile.Open(file_name, "READ")
        
        # Loop over all the keys in the input file
        for key in input_root.GetListOfKeys():

            obj = key.ReadObj()

            # Check if the object is a directory
            if isinstance(obj, ROOT.TDirectoryFile):
                dir_name = obj.GetName()

                # Loop over all histograms in the directory
                input_root.cd(dir_name)
                for hist_key in ROOT.gDirectory.GetListOfKeys():
                    hist_name = hist_key.GetName()
                    hist = hist_key.ReadObj()
                        
                    # Check if it's a histogram
                    if isinstance(hist, ROOT.TH1):
                        out_key = (hist_name,dir_name,file_name)
                        histogram_dir_file_map[out_key] = hist.Clone()
                        histogram_dir_file_map[out_key].SetDirectory(0)
    print("Mapped histograms from input files")

    # now we convert this map into a histogram_dir_map with format (hist_name, dir_name) : hist
    # we check that a histogram exists for all files, if it is a systematic histogram identifed by ending in 'Up', 'Up_2D', 'Down', or 'Down_2D' then it might not exist in all files
    # in this case we will add the nominal histogram from that file instead

    # first create empty histogram_dir_map by looping over the keys in histogram_dir_file_map
    # we also get a set of the nominal histograms names which are these ones not ending in 'Up', 'Up_2D', 'Down', or 'Down_2D'
    histogram_dir_map = {}
    nominal_hist_names = set()

    for (hist_name, dir_name, file_name) in histogram_dir_file_map.keys():
        out_key = (hist_name, dir_name)
        if out_key not in histogram_dir_map:
            histogram_dir_map[out_key] = None
        if not (hist_name.endswith('Up') or hist_name.endswith('Down') or hist_name.endswith('Up_2D') or hist_name.endswith('Down_2D')):
            nominal_hist_names.add(hist_name)

    print("Getting Keys")
    for key in tqdm(histogram_dir_map.keys()):
        hist_name, dir_name = key
        hist = None
        hist_count = 0
        for file_name in input_files:
            file_key = (hist_name, dir_name, file_name)
        
            if file_key in histogram_dir_file_map:
                if hist is None:
                    hist = histogram_dir_file_map[file_key].Clone()
                    hist.SetDirectory(0)
                    hist_count += 1
                else:
                    hist.Add(histogram_dir_file_map[file_key])
                    hist_count += 1
            else:
                # check if the histogram is a systematic histogram
                if hist_name.endswith('Up') or hist_name.endswith('Down') or hist_name.endswith('Up_2D') or hist_name.endswith('Down_2D'):
                    # get the nominal histogram name, loop over the nominal_hist_names and find the one that matches
                    nominal_hist_name = None
                    for nh in nominal_hist_names:
                        if hist_name.startswith(nh.replace('_2D','')) and (hist_name.endswith('_2D') == nh.endswith('_2D')):
                            nominal_hist_name = nh
                            break

                    if nominal_hist_name is not None:
                        # add the nominal histogram from this file instead

                        nominal_file_key = (nominal_hist_name, dir_name, file_name)
                        if nominal_file_key in histogram_dir_file_map:
                            if hist is None:
                                hist = histogram_dir_file_map[nominal_file_key].Clone()
                                hist.SetDirectory(0)
                                hist_count += 1
                            else:
                                hist.Add(histogram_dir_file_map[nominal_file_key])
                                hist_count += 1
                    
        if hist is None:
            print(f"Warning: Histogram {hist_name} in directory {dir_name} not found in any input files. Not adding to output.")
        elif exp_num is None or (hist_count == exp_num):
            histogram_dir_map[key] = hist
        elif hist_count != exp_num:
            print(f"Warning: Histogram {hist_name} in directory {dir_name} found in {hist_count} out of expected {exp_num} input files. Not adding to output.")

    print("Processed Keys")

    # now we take care of combining directories together if needed based on dir_combinations
    if dir_combinations:
        combined_histogram_dir_map = {}
        for (hist_name, dir_name), hist in histogram_dir_map.items():
            if hist is None: continue
            combined = False
            for combined_dir, dirs_to_combine in dir_combinations.items():
                if dir_name in dirs_to_combine:
                    combined = True
                    out_key = (hist_name, combined_dir)
                    if out_key not in combined_histogram_dir_map:
                        combined_histogram_dir_map[out_key] = hist.Clone()
                        combined_histogram_dir_map[out_key].SetDirectory(0)
                    else:
                        combined_histogram_dir_map[out_key].Add(hist)
            if not combined:
                out_key = (hist_name, dir_name)
                if out_key not in combined_histogram_dir_map:
                    combined_histogram_dir_map[out_key] = hist.Clone()
                    combined_histogram_dir_map[out_key].SetDirectory(0)
                else:
                    combined_histogram_dir_map[out_key].Add(hist)
        histogram_dir_map = combined_histogram_dir_map

    dir_names = []
    print("Combined directories")

    # now we write the histograms to the output file
    for (hist_name, dir_name), hist in histogram_dir_map.items():
        if hist is None: continue
        # Create the directory in the output file if it doesn't exist
        if not output.GetDirectory(dir_name):
            output.mkdir(dir_name)
            dir_names.append(dir_name)
        output.cd(dir_name)
        hist.SetName(hist_name)
        hist.Write(hist_name)

    # Uncomment the below if you want to manually create 'JetFakes' histograms by summing MC jet backgrounds
    if channel in ['mt', 'et']: # manually create jet fakes while without FFs
        for key in output.GetListOfKeys():
            if key.GetClassName() != "TDirectoryFile":
                continue
            dirname = key.GetName()
            print(f"Processing directory: {dirname}")
            dir_obj = output.GetDirectory(dirname)
            hists = [dir_obj.Get(hn) for hn in ["TTJ", "VVJ", "W", "QCD", "ZJ"]]
            hists_to_sum = [h for h in hists if h is not None]
            if not hists_to_sum:
                continue
            # Clone first histo and add others
            h_sum = hists_to_sum[0].Clone("JetFakes")
            h_sum.Reset() # clear
            for h in hists_to_sum:
                h_sum.Add(h)
            # Write to directory
            dir_obj.cd()
            h_sum.Write("JetFakes")


    # Close the output file
    output.Close()

    if len(config['eras']) == 1:
        era = config['eras'][0]
    elif set(config['eras']) == {'Run3_2022', 'Run3_2022EE', 'Run3_2023', 'Run3_2023BPix'}:
        era = 'earlyrun3'
    elif set(config['eras']) == {'Run3_2024', 'Run3_2025'}:
        era = 'laterun3'
    else:
        era = '...'

    print("Written output file:", output_file)

    aco_categories = []
    for entry in config['cpdecay'][channel]:
        if 'aco' in entry['plotting_variable'][0]:
            aco_categories.append(ch + '_' + entry['category'][0])
            aco_categories.append(ch + '_' + entry['category'][0] + '_aiso')
    aco_categories = set(aco_categories)

    for dir_name in dir_names:

        print('\033[1;34m=============================================\033[0m')
        print(f"\033[1;34mPlotting directory: {dir_name}\033[0m")

        blind = True
        is2Dunrolled = False
        var_name = "Bin number"
        
        if 'mva_fake' in dir_name or 'mva_tau' in dir_name or 'aiso' in dir_name or '_ss' in dir_name:
            blind = False 
        if 'mva_fake' in dir_name or 'mva_tau' in dir_name or 'mva_higgs' in dir_name:
            var_name = "BDT score"
        else:
            for aco_category in aco_categories:
                if aco_category in dir_name:
                    is2Dunrolled = True
                    var_name = find_variable(aco_category, channel, config)
                    break
            
        method = 6 # method that plots jetfakes
        # make a plot of the combined histograms
        Histo_Plotter = HTT_Histogram(
            output_file,
            dir_name,
            channel,
            era,
            var_name,
            method,
            blind=blind,
            log_y=False,
            is2Dunrolled=is2Dunrolled,
            save_name=output_file.replace('.root', f'_{dir_name}',)
        )
        Histo_Plotter.plot_1D_histo()



if __name__ == "__main__":
    # Argument parser for command-line arguments
    parser = argparse.ArgumentParser()
    
    # Add input and output files arguments
    parser.add_argument('-i', '--input', nargs='+', required=True, help="List of input ROOT files")
    parser.add_argument('-o', '--output', required=True, help="Name of the output ROOT file")
    parser.add_argument('-e', '--expected', type=int, default=None, help="Expected number of input files each histogram should appear in i.e the number of eras")
    parser.add_argument('-c', '--channel', default='tt', help="Channel to process (default: 'tt')")
    parser.add_argument('--config', type=str, default=None, help="Path to the configuration file (default: 'Draw/scripts/cpdecay_datacards.yaml')")

    # Parse arguments
    args = parser.parse_args()

    # List of input ROOT files from command-line
    input_files = args.input

    # Output ROOT file from command-line
    output_file = args.output

    # channel to precess
    ch = args.channel

    # expected number of input files each histogram should appear in i.e the number of eras
    exp_num = args.expected

    # Load the configuration file
    if args.config is None:
        args.config = find_config(args.input)
        print(
f"""\033[1;91mWARNING: No config file provided, using default: \
{args.config} — luminosity and BDT labels may not be \
correct!\033[0m"""
        )
        
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # Define which directories to combine (customize this based on your case)
    if ch == 'tt':
        dir_combinations = {
            'tt_higgs_pirho': ['tt_higgs_pirho', 'tt_higgs_rhopi'],
            'tt_higgs_rhoa1': ['tt_higgs_a1rho', 'tt_higgs_rhoa1'],
            'tt_higgs_pia1': ['tt_higgs_a1pi', 'tt_higgs_pia1'],
            'tt_higgs_pia11pr': ['tt_higgs_pia11pr', 'tt_higgs_a11prpi'],
            'tt_higgs_a11pra1': ['tt_higgs_a11pra1', 'tt_higgs_a1a11pr'],

            'tt_fake_pirho': ['tt_fake_pirho', 'tt_fake_rhopi'],
            'tt_fake_rhoa1': ['tt_fake_a1rho', 'tt_fake_rhoa1'],
            'tt_fake_pia1': ['tt_fake_a1pi', 'tt_fake_pia1'],
            'tt_fake_pia11pr': ['tt_fake_pia11pr', 'tt_fake_a11prpi'],
            'tt_fake_a11pra1': ['tt_fake_a11pra1', 'tt_fake_a1a11pr'],

            'tt_tau_pirho': ['tt_tau_pirho', 'tt_tau_rhopi'],
            'tt_tau_rhoa1': ['tt_tau_a1rho', 'tt_tau_rhoa1'],
            'tt_tau_pia1': ['tt_tau_a1pi', 'tt_tau_pia1'],
            'tt_tau_pia11pr': ['tt_tau_pia11pr', 'tt_tau_a11prpi'],
            'tt_tau_a11pra1': ['tt_tau_a11pra1', 'tt_tau_a1a11pr'],
        }
        # add aiso directories to the combinations
        dir_combinations_extra = {}
        for dir_name in dir_combinations.keys():
            dir_combinations_extra[dir_name + '_aiso'] = [d + '_aiso' for d in dir_combinations[dir_name]]
            dir_combinations_extra[dir_name + '_aco_aiso'] = [d + '_aco_aiso' for d in dir_combinations[dir_name]]
            dir_combinations_extra[dir_name + '_BDT_score_aiso'] = [d + '_BDT_score_aiso' for d in dir_combinations[dir_name]]
        dir_combinations.update(dir_combinations_extra)

    elif ch in ['et', 'mt']:
        dir_combinations = None
        # {
        #     'mt_higgs_mupi': ['mt_higgs_mupi'],
        #     'mt_higgs_murho': ['mt_higgs_murho'],
        #     'mt_higgs_mua11pr': ['mt_higgs_mua11pr'],
        #     'mt_higgs_mua1': ['mt_higgs_mua1'],
        # }


    # Call the hadd function with the provided arguments
    hadd_root_files(input_files, output_file, dir_combinations, ch, config, exp_num=exp_num)
   
