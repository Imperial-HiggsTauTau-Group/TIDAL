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

def plot_combined(cmb_file, tree_name, channel, era, variable, blind, unrolled):
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
        is2Dunrolled=unrolled,
    )
    Histo_Plotter.plot_1D_histo()

def get_files_to_combine(eras_to_combine, channel, scheme, directory, f_name):
    # example args: [Run3_2022, Run3_2022EE], et, control, /path/to/dir, datacard_m_vis_cp_inclusive_tt_Run3_2022.root
    print(f"\nERAS TO COMBINE: {eras_to_combine}")

    # Get the files
    files_to_combine = []
    old_args = {}
    for era in eras_to_combine:
        src_file = os.path.join(directory, era, scheme, channel, f_name.replace("$ERA", era))
        files_to_combine.append(src_file)
        print(f">> FILE: {src_file}")
    # Combine histograms
    if ("Run3_2022" in eras_to_combine) and ("Run3_2022EE" in eras_to_combine) and ("Run3_2023" in eras_to_combine) and ("Run3_2023BPix" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_earlyRun3", channel)
        cmb_postfix = "full2223.root"
        era = "earlyrun3"
    elif ("Run3_2022" in eras_to_combine) and ("Run3_2022EE" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_2022", channel)
        cmb_postfix = "full22.root"
        era = "full22"
    elif ("Run3_2023" in eras_to_combine) and ("Run3_2023BPix" in eras_to_combine):
        cmb_directory = os.path.join(directory, "combined_2023", channel)
        cmb_postfix = "full23.root"
        era = "full23"
    else:
        raise ValueError("Combination of eras not supported!")
    # output file name
    cmb_file = os.path.join(cmb_directory, f_name.replace("$ERA", "").replace(".root", "") + cmb_postfix)
    # make output directory
    os.makedirs(cmb_directory, exist_ok=True)

    return files_to_combine, cmb_file, era

def main(args, eras):

    # NEW
    with open(args.config) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    input_folder = config["input_folder"]
    output_path = config["output_path"]
    channels = config["channels"]
    eras = config["eras"]
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
                            print(f"Identified files to merge: {f_name}")
                            # Get files to combine and target file
                            files_to_combine, cmb_file, era_name = get_files_to_combine(eras, channel, scheme, output_path, f_name)
                            # Actually combine the histos
                            tree_name = f"{channel}_{cat}{nodename}"
                            combine_histos(files_to_combine, cmb_file)
                            # plot the combined histograms
                            if variable_name.count(",") > 1 and not unroll: # check if unrolled (for CP)
                                print("2D variable cannot be plotted directly, skipping plotting step.")
                            else:
                                plot_combined(cmb_file, tree_name, channel, era_name, variable, blind, unroll)


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
    if input_args.all:
        eras_to_combine = ['Run3_2022', 'Run3_2022EE', 'Run3_2023', 'Run3_2023BPix']
    elif input_args.only22:
        eras_to_combine = ['Run3_2022', 'Run3_2022EE']
    elif input_args.only23:
        eras_to_combine = ['Run3_2023', 'Run3_2023BPix']
    else:
        raise ValueError("Please specify which eras to combine using --all, --only22, or --only23")

    main(input_args, eras_to_combine)

