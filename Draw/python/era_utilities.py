import yaml

LUMI2024 = {
    "Run2024B":   130,
    "Run2024C":  7260,
    "Run2024D":  7980,
    "Run2024E": 11420,
    "Run2024F": 28040,
    "Run2024G": 38070,
    "Run2024H":  5490,
    "Run2024I": 11560,
}

LUMI2025 = {
    "Run2025B":   250,
    "Run2025C": 21560,
    "Run2025D": 25820,
    "Run2025E": 14050,
    "Run2025F": 26690,
    "Run2025G": 22250,
}

LUMI = {
    "Run3_2024": LUMI2024,
    "Run3_2025": LUMI2025,
}


def build_temp_config(parameter_file, era, specific_eras):
    with open(parameter_file, "r") as f:
        config = yaml.safe_load(f)
    lumi = sum(LUMI[era][e] for e in specific_eras)
    config["lumi"] = lumi
    return config

