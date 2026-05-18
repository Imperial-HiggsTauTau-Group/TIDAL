import argparse
import yaml
import numpy
import subprocess
import os
import re
import shutil

def create_bins(variable: str) -> str:
    if "(" in variable:
        variable_name = variable.split("(")[0]

        matches = re.findall(r'\((.*?)\)', variable)
        full_binning = ""
        for index,match in enumerate(matches):
            number_of_bins = match.split(",")[0]
            lower_bound = match.split(",")[1].strip()
            upper_bound = match.split(",")[2].strip()

            # Function to handle pi expressions
            def evaluate_bound(bound: str) -> float:
                if "pi" in bound:
                    bound = bound.replace("pi", "*numpy.pi")
                    if bound.startswith("*"):
                        bound = bound[1:]  # Remove leading '*' for valid syntax
                    if bound.startswith("-*"):
                        bound = "-" + bound[2:]  # Fix '-*' syntax
                return eval(bound)  # Safely evaluate

            # Evaluate bounds
            lower_bound = str(evaluate_bound(lower_bound))
            upper_bound = str(evaluate_bound(upper_bound))

            # Create the binning array
            binning = numpy.linspace(
                float(lower_bound), float(upper_bound), int(number_of_bins) + 1
            )
            if index == 0:
                full_binning = f'{full_binning}[{",".join([str(i) for i in binning])}]'
            else:
                full_binning = f'{full_binning},[{",".join([str(i) for i in binning])}]'

        variable = f"{variable_name}{full_binning}"

    return variable


def create_condor_submit_file(
    logs_path: str, variable_name: str, submit_file: str, script_path: str, era: str = "", channel: str = "", run_systematics: bool = False
):
    condor_template = f"""
executable = {script_path}
priority = 15
output = {logs_path}/condor_{variable_name}.out
error = {logs_path}/condor_{variable_name}.err
log = {logs_path}/condor_{variable_name}.log
request_memory = 8000
request_cpus = 1
getenv = True
+MaxRuntime = 10500
queue
"""
    if run_systematics and channel in ["et", "mt"]:
        print("ASSIGNING EXTRA RUNTIME AND MEMORY (et/mt with systematics)")
        condor_template = condor_template.replace("request_cpus = 1", "request_cpus = 3")
        condor_template = condor_template.replace("request_memory = 8000", "request_memory = 12000")
        if era in ['Run3_2022EE']:
            condor_template = condor_template.replace("+MaxRuntime = 10500", "+MaxRuntime = 100000")
        else:
            condor_template = condor_template.replace("+MaxRuntime = 10500", "+MaxRuntime = 35800")
    elif run_systematics and channel in ['tt'] and era in ["Run3_2022EE", "Run3_2023"]:
        print("ASSIGNING EXTRA RUNTIME AND MEMORY (tt with systematics)")
        condor_template = condor_template.replace("+MaxRuntime = 10500", "+MaxRuntime = 35800") # one core should still be enough
    elif not run_systematics: # limit runtime to 7000 s to get in the shorter queues
        condor_template = condor_template.replace("+MaxRuntime = 10500", "+MaxRuntime = 7000")
    with open(submit_file, "w") as f:
        f.write(condor_template)
    os.system(f"chmod +x {submit_file}")


def create_shell_script(
    input_folder,
    output_folder,
    parameter_file,
    channel,
    era,
    method,
    category,
    variable,
    additional_selection,
    additional_weight,
    datacard_name,
    script_path,
    run_systematics=False,
    systematics_to_run=[],
    blind=False,
    auto_rebin=False,
    set_alias="",
    aiso=False,
    same_sign=False,
    unroll=False,
    rename_procs=False,
    dy_LO=False,
    dy_NLO=False,
    use_filtered_DY=False,
    nodename="",
):
    shell_script = f"""
#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh
python3 Draw/scripts/HiggsTauTauPlot.py \\
--parameter_file {parameter_file} \\
--input_folder {input_folder} \\
--output_folder {output_folder} \\
--channel {channel} \\
--era {era} \\
--method {method} \\
--category {category} \\
--var {variable} \\
--sel '{additional_selection}' \\
--add_weight '{additional_weight}' \\
--datacard_name {datacard_name}"""

    if blind:
        shell_script += " \\\n--blind"
    if auto_rebin:
        shell_script += " \\\n--auto_rebin"
    if run_systematics:
        if systematics_to_run:
            shell_script += " \\\n--run_systematics"
            for syst in systematics_to_run:
                shell_script += f" \\\n--{syst}"
    if set_alias:
        shell_script += f" \\\n--set_alias='{set_alias}'"
    if aiso:
        shell_script += " \\\n--do_aiso"
    if same_sign:
        shell_script += " \\\n--do_ss"
    if unroll:
        shell_script += " \\\n--do_unrolling"
    if rename_procs:
        shell_script += " \\\n--rename_procs"
    if dy_LO:
        shell_script += " \\\n--LO_DY"
    if dy_NLO:
        shell_script += " \\\n--NLO_DY"
    if use_filtered_DY:
        shell_script += " \\\n--use_filtered_DY"
    if nodename != "":
        shell_script += f" \\\n--nodename {nodename}"

    with open(script_path, "w") as script_file:
        print(shell_script)
        script_file.write(shell_script)
    os.system(f"chmod +x {script_path}")


def format_first_selection(selection):
    # Extract the first condition using regex
    match = re.match(r'([^&|]+)', selection.strip())
    if match:
        condition = match.group(1).strip()
        # Replace operators with readable words
        formatted = condition.replace(">=", "_GTE_").replace("<=", "_LTE_")\
                              .replace(">", "_GT_").replace("<", "_LT_")\
                              .replace("==", "_EQ_").replace("!=", "_NEQ_")
        # Replace spaces with underscores
        readable_format = formatted.replace(" ", "").removeprefix("(").removesuffix(")")
        format_dict = {
            'mt_1_LT_65': 'mTLt65',
            'mt_1_GT_70': 'mTGt70',
        }
        if readable_format in format_dict:
            return format_dict[readable_format]

        return readable_format
    return None

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--config", type=str, help="Configuration file")
    parser.add_argument("--batch", action="store_true", help="Run in batch mode")

    args = parser.parse_args()

    config_file = args.config

    # read config file which is a yaml file
    with open(config_file) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    input_folder = config["input_folder"]
    output_path = config["output_path"]

    # save config
    os.makedirs(output_path, exist_ok=True)
    shutil.copy(config_file, output_path)

    channels = config["channels"]
    eras = config["eras"]
    parameter_path = config["parameter_path"]
    schemes = config["schemes"]
    run_systematics = config["run_systematics"]

    available_channels = ["mm", "ee", "mt", "tt", "et"]
    for channel in channels:
        if channel not in available_channels:
            raise ValueError(
                f"Channel {channel} is not a valid channel. Please choose from {available_channels}"
            )

    available_eras = ["Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"]
    for era in eras:
        if era not in available_eras:
            raise ValueError(
                f"Era {era} is not a valid era. Please choose from {available_eras}"
            )

    available_schemes = ["sf_calculation", "control", "cpdecay", "cp_acoplanarity", "cpdecay_fakefactors_control", "cp_alpha_angles", "cp_chi2_tests"]

    for scheme in schemes:
        if scheme not in available_schemes:
            raise ValueError(
                f"Scheme {scheme} is not a valid scheme. Please choose from {available_schemes}"
            )

    for era in eras:
        parameter_file = f"{parameter_path}/{era}/params.yaml"
        for channel in channels:
            for scheme in schemes:
                settings = config[scheme]
                output_folder = f"{output_path}/{era}/{scheme}/{channel}"
                subprocess.run(["mkdir", "-p", output_folder])

                # check if settings Aliases exists
                if "Aliases" in settings:
                    available_aliases = settings["Aliases"]
                else:
                    available_aliases = {}

                systematics_to_run = []
                if run_systematics and "systematics" in config[scheme]:
                    if "per_era" in config[scheme]["systematics"]:
                        if (
                            era not in config[scheme]["systematics"]["per_era"]
                            and "default" in config[scheme]["systematics"]["per_era"]
                        ):
                            if (
                                channel
                                in config[scheme]["systematics"]["per_era"]["default"]
                            ):
                                channel_systematics = config[scheme]["systematics"][
                                    "per_era"
                                ]["default"][channel]
                                for syst in channel_systematics:
                                    systematics_to_run.append(f"{syst[0]}={syst[1]}")
                            if (
                                "common"
                                in config[scheme]["systematics"]["per_era"]["default"]
                            ):
                                common_systematics = config[scheme]["systematics"][
                                    "per_era"
                                ]["default"]["common"]
                                for syst in common_systematics:
                                    systematics_to_run.append(f"{syst[0]}={syst[1]}")

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
                    dy_LO= setting.get("dy_LO", False)
                    dy_NLO = setting.get("dy_NLO", False)
                    use_filtered_DY = setting.get('use_filtered_DY', False)
                    if use_filtered_DY:
                        print("WARNING: To use filtered DY, different parameters will be sourced (automatically)")
                        parameter_file = f"{parameter_path}/{era}/params_DYfiltered.yaml"
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
                                if dy_LO:
                                    variable_name = variable_name + "_dy_LO"
                                if dy_NLO:
                                    variable_name = variable_name + "_dy_NLO"

                                if nodename != "":
                                    if nodename[0] != "_":
                                        nodename = "_" + nodename

                                filename = f'{variable_name}_{cat}'

                                logs = f"{output_folder}/logs"
                                subprocess.run(["mkdir", "-p", logs])
                                script_path = os.path.join(
                                    logs, f"{filename}.sh"
                                )
                                create_shell_script(
                                    input_folder,
                                    output_folder,
                                    parameter_file,
                                    channel,
                                    era,
                                    method,
                                    cat,
                                    variable,
                                    additional_selection,
                                    additional_weight,
                                    variable_name,
                                    script_path,
                                    run_systematics=run_systematics,
                                    systematics_to_run=systematics_to_run,
                                    blind=blind,
                                    auto_rebin=auto_rebin,
                                    set_alias=alias,
                                    aiso=aiso,
                                    same_sign=same_sign,
                                    unroll=unroll,
                                    rename_procs=rename_procs,
                                    dy_LO=dy_LO,
                                    dy_NLO=dy_NLO,
                                    use_filtered_DY=use_filtered_DY,
                                    nodename=nodename,
                                )

                                submit_file = os.path.join(
                                logs, f"submit_{filename}.sub"
                                )
                                create_condor_submit_file(
                                logs, filename, submit_file, script_path, era, channel, run_systematics
                                )
                                if args.batch:
                                    subprocess.run(["condor_submit", submit_file])
                                else:
                                    os.chmod(script_path, 0o755)
                                    subprocess.run(["/bin/bash", script_path])