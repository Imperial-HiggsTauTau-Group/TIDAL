from glob import glob
import subprocess
import argparse
import yaml
import os

config_files = {
    "control": "Draw/scripts/config_plot_LateRun3.yaml",
    "cpdecay": "Draw/scripts/cpdecay_datacards_LateRun3.yaml",
}


def get_args():
    parser = argparse.ArgumentParser(description="Run cut optimisation for HiggsTauTauPlot")
    parser.add_argument("--output", "-o", type=str, default="Draw/plots/cuts_for_LateRun3")
    parser.add_argument("--step", type=str, help="Choose from: 'control', 'cpdecay', 'hadd'")
    return parser.parse_args()


def to_string(value: float) -> str:
    return str(value).replace(".", "p")


def create_temp_config(config, output, IPsig, Esplit):
    with open(config, "r") as f:
        config_dict = yaml.safe_load(f)

    temp_config_dir = "Draw/scripts/temp_configs"
    os.makedirs(temp_config_dir, exist_ok=True)
    temp_config = os.path.join(
        temp_config_dir,
        f"IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}.yaml"
    )

    config_dict["output_path"] = f"{output}/IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}"
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
                temp_config = create_temp_config(config, args.output, IPsig, Esplit)
                subprocess.run([
                    "python",
                    "Draw/scripts/makeDatacards.py",
                    "--config",
                    temp_config,
                    "--batch"
                ])
            
            elif args.step == "hadd":
                os.makedirs(f"{args.output}/IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}/Combined", exist_ok=True)
                hadd_command = (
                    ["python", "Draw/scripts/hadd_cp_datacards.py", "-i"]
                    + glob(f"{args.output}/IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}/Run3_*/cpdecay/*/*.root")
                    + ["-o", f"{args.output}/IPsig_{to_string(IPsig)}_Esplit_{to_string(Esplit)}/Combined/added_histo.root"]
                )
                subprocess.run(hadd_command)

    if args.step in ["control", "cpdecay"]:
        # Clean up temporary config files
        temp_config_dir = "Draw/scripts/temp_configs"
        for temp_file in glob(os.path.join(temp_config_dir, "*.yaml")):
            os.remove(temp_file)


if __name__ == "__main__":
    args = get_args()
    main(args)