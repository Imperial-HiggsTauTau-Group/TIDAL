import os
import shutil
import shlex
import yaml
import argparse
from Draw.python.PlotHistograms import HTT_Histogram

def find_sh_file(root_file):
    # Find the shell file with parameters used to make the plot
    job_dir = os.path.join(os.path.dirname(root_file), 'logs')
    log_file_name = root_file.split('/')[-1].split('datacard_')[1].split(f'_{channel}_Run3')[0] + ".sh"
    return os.path.join(job_dir, log_file_name)

def combine_histos(root_files, output_file):
    # Combine multiple ROOT files into one using hadd
    command = f"hadd -v 1 -f {output_file} " + " ".join(root_files)
    print("Combining histograms")
    os.system(command)

def plot_combined(cmb_file, tree_name, channel, era, variable, blind=False):
    # Plot the combined root files
    print(f"\n>> PLOTTING FILE: {cmb_file}")
    Histo_Plotter = HTT_Histogram(
        cmb_file,
        tree_name,
        channel,
        era,
        variable,
        blind=blind,
        log_y=False,
        is2Dunrolled=False,
    )
    Histo_Plotter.plot_1D_histo()

def combine_eras(eras_to_combine, channel, scheme, directory, f_name):
    # example args: [Run3_2022, Run3_2022EE], et, control, /path/to/dir, datacard_m_vis_cp_inclusive_tt_Run3_2022.root
    print(f"\nERAS TO COMBINE: {eras_to_combine}")

    # Get the files
    files_to_combine = []
    old_args = {}
    for era in eras_to_combine:
        src_file = os.path.join(directory, era, scheme, channel, f_name.replace("$ERA", era))
        files_to_combine.append(src_file)
        print(f">> FILE: {src_file}")
        log_file = find_sh_file(src_file)
        with open(log_file, 'r') as f:
            tokens = shlex.split(f.read())
        # Build a dictionary with all arguments
        args = {
            tokens[i].split("\n--")[1]: tokens[i + 1]
            for i in range(len(tokens) - 1)
            if tokens[i].startswith("\n--")
        }
        # Check that all parameters match in the arguments!
        if old_args:
            if args['channel'] != old_args['channel'] or args['category'] != old_args['category'] or args['var'] != old_args['var']:
                raise ValueError("Different parameters found in the files to combine!")
        else:
            old_args = args

    # Combine histograms
    out_prefix = f"datacard_{args['var'].split('[')[0]}_{args['category']}_{args['channel']}"
    if ("Run3_2022" in eras_to_combine) and ("Run3_2022EE" in eras_to_combine) and ("Run3_2023" in eras_to_combine) and ("Run3_2023BPix" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_earlyRun3")
        cmb_file = os.path.join(cmb_directory, out_prefix + "_full2223.root")
        era = "earlyrun3"
    elif ("Run3_2022" in eras_to_combine) and ("Run3_2022EE" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_2022")
        cmb_file = os.path.join(cmb_directory, out_prefix + "_full22.root")
        era = "full22"
    elif ("Run3_2023" in eras_to_combine) and ("Run3_2023BPix" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_2023")
        cmb_file = os.path.join(cmb_directory, out_prefix + "_full23.root")
        era = "full23"
    else:
        raise ValueError("Combination of eras not supported!")
    # make output directory
    os.makedirs(cmb_directory, exist_ok=True)

    # Check for extra name in node
    if 'nodename' in args.keys():
        tree_name = args['channel'] + '_' + args['category'] + args['nodename']
    else:
        tree_name = args['channel'] + '_' + args['category']

    # combine the histograms
    combine_histos(files_to_combine, cmb_file)

    # plot the combined histograms
    plot_combined(cmb_file, tree_name, channel, era, args['var'])


def get_args():
    parser = argparse.ArgumentParser(description="Combine datacards from different eras and plot them")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration YAML file')
    parser.add_argument('--all', action='store_true', help='Combine all eras')
    parser.add_argument('--only22', action='store_true', help='Combine 2022 eras')
    parser.add_argument('--only23', action='store_true', help='Combine 2023 eras')
    input_args = parser.parse_args()
    return input_args

if __name__ == "__main__":

    input_args = get_args()
    config = input_args.config
    if input_args.all:
        eras_to_combine = ['Run3_2022', 'Run3_2022EE', 'Run3_2023', 'Run3_2023BPix']
    elif input_args.only22:
        eras_to_combine = ['Run3_2022', 'Run3_2022EE']
    elif input_args.only23:
        eras_to_combine = ['Run3_2023', 'Run3_2023BPix']
    else:
        raise ValueError("Please specify which eras to combine using --all, --only22, or --only23")

    with open(config) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    directory = config['output_path']
    for scheme in config['schemes']:
        scheme_config = config[scheme]
        for channel in config['channels']:
            plot_config = scheme_config[channel]
            for plot in plot_config: # can have several lines in the config
                for category in plot['category']:
                    for variable in plot['plotting_variable']:
                        if channel in ['et', 'mt']:
                            f_name = f"datacard_{variable}_mTLt65_{category}_{channel}_$ERA.root"
                        else:  # tt channel
                            f_name = f"datacard_{variable}_{category}_{channel}_$ERA.root"
                        print("\n" + "*"*100)
                        print(f"-> Datacards to merge: {f_name}")
                        combine_eras(eras_to_combine, channel, scheme, directory, f_name)

