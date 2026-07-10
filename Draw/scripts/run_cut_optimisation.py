from glob import glob
import subprocess
import argparse
import yaml
import os

config_files = {
    "control": "Draw/scripts/config_plot_LateRun3.yaml",
    "cpdecay": "Draw/scripts/cpdecay_datacards_LateRun3.yaml",
    "hadd": None
}


def get_args():
    parser = argparse.ArgumentParser(description="Run cut optimisation for HiggsTauTauPlot")
    parser.add_argument("--output", "-o", type=str, default="Draw/plots/cuts_for_LateRun3")
    parser.add_argument("--step", type=str, help="Choose from: 'control', 'cpdecay', 'hadd'", required=True)
    parser.add_argument("--channels", type=str, help="e.g. et,mt,tt", required=True)
    args = parser.parse_args()
    args.channels = args.channels.split(",")
    return args


def to_string(value: float) -> str:
    return str(value).replace(".", "p")


def create_temp_config(config, output, step, channels, IPsig, Esplit):
    with open(config, "r") as f:
        config_dict = yaml.safe_load(f)

    temp_config_dir = "Draw/scripts/temp_configs"
    os.makedirs(temp_config_dir, exist_ok=True)
    temp_config = os.path.join(
        temp_config_dir,
        f"{step}_IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}.yaml"
    )

    config_dict["output_path"] = f"{output}/IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}"
    config_dict["channels"] = channels
    config_dict["IPsig"] = IPsig
    config_dict["Esplit"] = Esplit

    with open(temp_config, "w") as f:
        yaml.dump(config_dict, f, sort_keys=False)

    return temp_config


def main(args):
    config = config_files[args.step]
    IPsig_values = [1.25]
    Esplit_values = [0.20]

    for IPsig in IPsig_values:
        for Esplit in Esplit_values:
            if args.step in ["control", "cpdecay"]:
                temp_config = create_temp_config(config, args.output, args.step, args.channels, IPsig, Esplit)
                subprocess.run([
                    "python",
                    "Draw/scripts/makeDatacards.py",
                    "--config",
                    temp_config,
                    "--batch"
                ])
            
            elif args.step == "hadd":
                os.makedirs(f"{args.output}/IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}/Combined", exist_ok=True)
                for channel in args.channels:
                    hadd_command = (
                        ["python", "Draw/scripts/hadd_cp_datacards.py", "-i"]
                        + glob(f"{args.output}/IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}/Run3_*/cpdecay/{channel}/*.root")
                        + ["-o", f"{args.output}/IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}/Combined/added_histo_{channel}.root"]
                        + ["-c", channel]
                    )
                    subprocess.run(hadd_command, check=True)

    if args.step in ["control", "cpdecay"]:
        # Clean up temporary config files
        temp_config_dir = "Draw/scripts/temp_configs"
        for temp_file in glob(os.path.join(temp_config_dir, "*.yaml")):
            os.remove(temp_file)


if __name__ == "__main__":
    args = get_args()
    main(args)