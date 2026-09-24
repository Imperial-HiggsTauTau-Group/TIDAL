import os
import shutil
import shlex
import yaml
import argparse
from Draw.python.PlotHistograms import HTT_Histogram
from Draw.scripts.makeDatacards import create_bins, format_first_selection

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

def plot_combined(cmb_file, tree_name, channel, method, era, variable, blind, unrolled):
    # Plot the combined root files
    print(f"\n>> PLOTTING FILE: {cmb_file}")
    Histo_Plotter = HTT_Histogram(
        cmb_file,
        tree_name,
        channel,
        era,
        variable,
        method,
        blind=blind,
        log_y=False,
        is2Dunrolled=unrolled,
    )
    Histo_Plotter.plot_1D_histo()

#: names kept for the combinations people already refer to by name; any other
#: set of eras is named after the eras themselves
NAMED_COMBINATIONS = {
    ("Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"): ("combined_earlyRun3", "full2223", "earlyrun3"),
    ("Run3_2022", "Run3_2022EE"): ("combined_2022", "full22", "full22"),
    ("Run3_2023", "Run3_2023BPix"): ("combined_2023", "full23", "full23"),
}


def combination_names(eras_to_combine):
    """Directory, file postfix and era label for a set of eras.

    The label is what the plotter uses to find the luminosity: a known name
    for the combinations that have one, otherwise the eras joined with "+",
    which the plotter sums.
    """
    known = ("Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix", "Run3_2024", "Run3_2025")
    # chronological for the eras we know, then anything else alphabetically
    key = tuple(era for era in known if era in eras_to_combine)
    key += tuple(sorted(era for era in eras_to_combine if era not in known))
    if key in NAMED_COMBINATIONS:
        return NAMED_COMBINATIONS[key]
    joined = "_".join(key)
    return f"combined_{joined}", joined, "+".join(key)


def get_files_to_combine(eras_to_combine, channel, scheme, directory, f_name):
    print(f"\nERAS TO COMBINE: {eras_to_combine}")
    if len(eras_to_combine) < 2:
        raise ValueError("at least two eras are needed to combine")

    # Get the files
    files_to_combine = []
    missing = []
    for era in eras_to_combine:
        src_file = os.path.join(directory, era, scheme, channel, f_name.replace("$ERA", era))
        files_to_combine.append(src_file)
        print(f">> FILE: {src_file}")
        if not os.path.isfile(src_file):
            missing.append(src_file)
    if missing:
        raise FileNotFoundError("cannot combine, these inputs do not exist: " + ", ".join(missing))
    cmb_dirname, cmb_postfix, era = combination_names(eras_to_combine)
    cmb_directory = os.path.join(directory, cmb_dirname, channel)
    # output file name
    cmb_file = os.path.join(cmb_directory, f_name.replace("$ERA", "").replace(".root", "") + cmb_postfix + ".root")
    # make output directory
    os.makedirs(cmb_directory, exist_ok=True)

    return files_to_combine, cmb_file, era


def main(args, eras):

    # NEW
    with open(args.config) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    combined = []  # what was written, for the record

    input_folder = config["input_folder"]
    output_path = config["output_path"]
    channels = config["channels"]
    # the eras named on the command line, else all the eras the configuration lists
    eras = eras or config["eras"]
    parameter_path = config["parameter_path"]
    schemes = config["schemes"]
    run_systematics = config["run_systematics"]

    for channel in channels:
        for scheme in schemes:
            settings = config[scheme]
            if "Aliases" in settings:
                available_aliases = settings["Aliases"]
            else:
                available_aliases = {}
            for setting in settings[channel]:
                method = setting.get("method", "1")
                category = setting.get("category", "[inclusive]")
                variables = setting.get("plotting_variable", "[m_vis]")
                blind = setting.get("blind", False)
                auto_rebin = setting.get("auto_rebin", False)
                additional_selections = setting.get("additional_selections", ["(1)"])
                if isinstance(additional_selections, str):
                    additional_selections = [additional_selections]
                set_alias = setting.get("set_alias", "")
                additional_weight = setting.get("additional_weight", "(1)")
                extra_identifier = setting.get("extra_identifier", "")
                aiso = setting.get("aiso", False)
                same_sign = setting.get("same_sign", False)
                unroll = setting.get("unroll", False)
                rename_procs = setting.get("rename_procs", False)
                if set_alias and set_alias in available_aliases:
                    set_alias = available_aliases[set_alias]
                for cat in category:
                    if set_alias:
                        alias = set_alias.replace("category", cat)
                    else:
                        alias = None
                    for additional_selection in additional_selections:
                        for variable in variables:
                            variable = settings["variable_definitions"][variable]
                            variable = create_bins(variable)
                            variable_name = variable.split("[")[0]
                            variable_name = variable_name.replace(",", "_vs_")
                            nodename = ""
                            if additional_selection and additional_selection != "" and additional_selection != "(1)":
                                variable_name = variable_name + '_' + format_first_selection(additional_selection)
                                nodename = format_first_selection(additional_selection)
                            else:
                                variable_name = variable_name
                            nodename += setting.get("nodename", "")
                            if extra_identifier:
                                variable_name = variable_name + "_" + extra_identifier
                            if aiso:
                                variable_name = variable_name + "_aiso"
                                nodename = nodename + "_aiso"
                            if same_sign:
                                variable_name = variable_name + "_ss"
                                nodename = nodename + "_ss"
                            if nodename != "":
                                if nodename[0] != "_":
                                    nodename = "_" + nodename
                            filename = f'{variable_name}_{cat}'
                            f_name = f"datacard_{filename}_{channel}_$ERA.root"
                            if same_sign:
                                f_name = f_name.replace('$ERA', '$ERA_ss')
                            print(f"Identified files to merge: {f_name}")
                            # Get files to combine and target file
                            files_to_combine, cmb_file, era_name = get_files_to_combine(eras, channel, scheme, output_path, f_name)
                            # Actually combine the histos
                            tree_name = f"{channel}_{cat}{nodename}"
                            combine_histos(files_to_combine, cmb_file)
                            combined.append({"inputs": files_to_combine, "output": cmb_file,
                                             "channel": channel, "scheme": scheme, "era": era_name})
                            # plot the combined histograms
                            if variable_name.count(",") > 1 and not unroll: # check if unrolled (for CP)
                                print("2D variable cannot be plotted directly, skipping plotting step.")
                            else:
                                plot_combined(cmb_file, tree_name, channel, method, era_name, variable, blind, unroll)

    if args.record:
        # A record of what was actually combined, for whatever drives this
        # script; a wrapper cannot know from the outside which files came out.
        import json
        os.makedirs(os.path.dirname(os.path.abspath(args.record)), exist_ok=True)
        with open(args.record, "w") as handle:
            json.dump({"eras": list(eras), "config": os.path.abspath(args.config),
                       "combined": combined}, handle, indent=2)
        print(f"Record written to {args.record}")


def get_args():
    parser = argparse.ArgumentParser(description="Combine datacards from different eras and plot them")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration YAML file')
    parser.add_argument('--eras', nargs='+', default=None,
                        help='Eras to combine; default: every era the configuration lists')
    parser.add_argument('--record', type=str, default=None,
                        help='Write a JSON record of the combined files to this path')
    parser.add_argument('--all', action='store_true', help='Combine Run3_2022, Run3_2022EE, Run3_2023 and Run3_2023BPix')
    parser.add_argument('--only22', action='store_true', help='Combine the 2022 eras')
    parser.add_argument('--only23', action='store_true', help='Combine the 2023 eras')
    input_args = parser.parse_args()
    return input_args

if __name__ == "__main__":

    input_args = get_args()
    if input_args.all:
        eras_to_combine = ['Run3_2022', 'Run3_2022EE', 'Run3_2023', 'Run3_2023BPix']
    elif input_args.only22:
        eras_to_combine = ['Run3_2022', 'Run3_2022EE']
    elif input_args.only23:
        eras_to_combine = ['Run3_2023', 'Run3_2023BPix']
    else:
        eras_to_combine = input_args.eras  # None means: whatever the configuration lists

    main(input_args, eras_to_combine)
