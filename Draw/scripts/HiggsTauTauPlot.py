# HiggsTauTauPlot.py

# Methods:
# 1. Basic method using MC samples & SS method for QCD estimation
# 2. Basic method using MC samples & SS method for QCD estimation && WJets shape method (High mT control region)
# 3. Basic method using MC samples & Flat FF method for QCD estimation
# 4. Basic method using MC samples && FF method for QCD estimation
# 5. Basic method using MC samples && SS method for QCD estimation but, relaxed isolation

import argparse
from collections import OrderedDict
from prettytable import PrettyTable
import copy
import numpy as np
import ROOT
import re
from Draw.python import Analysis
from Draw.python import Plotting
from Draw.python.nodes import (
    BuildCutString,
    GenerateZTT,
    GenerateZLL,
    GenerateTop,
    GenerateVV,
    GenerateW,
    GenerateQCD,
    GenerateFakes,
    GenerateReweightedCPSignal,
    GenerateEWKZ
)
from Draw.python.HiggsTauTauPlot_utilities import (
    PrintSummary,
    GetTotals,
    FixBins,
    FindRebinning,
    RebinHist,
    UnrollHist2D,
    RenameDatacards,
    Total_Uncertainty
)
from Draw.scripts.systematics.systematics import generate_systematics_dict
from Draw.python.PlotHistograms import HTT_Histogram

ROOT.TH1.SetDefaultSumw2(True)

parser = argparse.ArgumentParser()
# ------------------------------------------------------------------------------------------------------------------------
# Main Options:
# parameter_file: yaml file containing the scaling parameters (e.g. cross-sections, luminosity, etc.)
# input/output folders, channel, era, method, category, variable, run_systematics, systematics_file
parser.add_argument("--bypass_plotter", action="store_true", help="Bypass plotter")
parser.add_argument("--parameter_file", type=str, help="Parameter file")
parser.add_argument("--input_folder", type=str, help="Input folder")
parser.add_argument("--output_folder", default="output", help="Output folder")
parser.add_argument("--channel", default="mm", help="Channel to run on")
parser.add_argument("--era", default="2016", help="Era to run on")
parser.add_argument("--method", default=1, help="Method to run on")
parser.add_argument("--category", default="inclusive", help="Category to run on")
parser.add_argument("--var", type=str, help="Variable to plot")
parser.add_argument("--run_systematics", action="store_true", help="Run systematics")

# Available Systematic Options:
# ------------------------------------------------------------------------------------------------------------------------
systematic_options = [
    ["Muon_ID_iso", "Muon ID and isolation systematics"],
    ["Electron_ID_Reco", "Electron ID and Reconstruction systematic"],
    ["Tau_ID", "Tau ID systematic"],
    ["Tau_FakeRate_e","Tau Fake Rate systematic for genuine electrons misidentified as taus"],
    ["Tau_FakeRate_mu","Tau Fake Rate systematic for genuine muons misidentified as taus"],
    ["Tau_EnergyScale_PNet_TSCALE", "Tau Energy Scale systematic for genuine taus"],
    ["Tau_EnergyScale_PNet_ESCALE", "Tau Energy Scale systematic for genuine electrons misidentified as taus"],
    ["Tau_EnergyScale_PNet_MUSCALE", "Tau Energy Scale systematic for genuine muons misidentified as taus"],
    ['MET_Recoil', 'MET Recoil systematic'],
    ['Trigger', 'Trigger efficiency systematic'],
    ["Jet_EnergyScale_Total", "Jet Energy Scale Total systematic"],
    ["Jet_EnergyScale_Regrouped", "Jet Energy Scale Regrouped systematics"],
    ["Jet_EnergyResolution", "Jet Energy Resolution systematic"],
    ["Electron_Scale_Smearing", "Electron Scale and Smearing systematics"],
    ["QCD_Background", "QCD Background systematic"],
    ["DY_Shape", "DY Shape systematic from ZpT reweighting"],
    ["DY_Shape_Imperial", "DY Shape systematic from ZpT reweighting (Imperial)"],
    ["TTbar_Shape", "TTbar Shape systematic from top pT reweighting"],
    ["Fake_Flat_Uncertainty", "flat fake uncertainty"],
    ["Tau_ID_PNet", "Tau ID systematic"],
    ["Fake_Factors", "fake factor related uncertainties"],
    ["Pileup", "pileup systematic"],
    ["Signal_Theory", "theoretical uncertainties on the signal"],
    ["IP_Calibration", "Uncertainty on the IP calibration"],
    ["SV_Resolution", "Uncertainty on the SV resolution"],
    ['IP_Significance', 'Uncertainty on the IP significance cut SFs']
]


# The nargs="?" option allows the user to provide an optional value for the systematic e.g. overriding the name assigned to the histogram
for systematic, description in systematic_options:
    parser.add_argument(
        f"--systematic_{systematic}",
        nargs="?",  # Accepts an optional value
        const=systematic,  # Default to the systematic name if no value is provided
        help=description,  # Use the description from the list
    )
# ------------------------------------------------------------------------------------------------------------------------

# Additional Options:
parser.add_argument("--LO_DY", action="store_true", help="Use LO")
parser.add_argument("--NLO_DY", action="store_true", help="Use NLO DY")
parser.add_argument("--use_filtered_DY", action="store_true", help="Include filtered DY")
parser.add_argument("--sel", type=str, help="Additional Selection to apply", default="")
parser.add_argument(
    "--set_alias",
    action="append",
    dest="set_alias",
    type=str,
    default=None,
    help="Overwrite alias selection using this options. Specify with the form --set_alias=nameofaliastoreset:newselection",
)
parser.add_argument("--add_weight", default="", help="Additional weight to apply")
parser.add_argument("--do_aiso", action="store_true", help="Do Anti-Isolated")
parser.add_argument("--do_ss", action="store_true", help="Do SS")
parser.add_argument("--blind", action="store_true", help="Blind the plot (remove data)")
parser.add_argument(
    "--masses", default="125", help="Mass points to process, seperated by commas"
)
parser.add_argument("--do_unrolling", action="store_true", help="Unroll 2D histograms")
parser.add_argument(
    "--rename_procs", action="store_true", help="Rename processes in the datacards"
)
parser.add_argument("--datacard_name", help="Override the datacard name")
parser.add_argument("--nodename", help="Override the nodename")
parser.add_argument(
    "--auto_rebin", action="store_true", help="Automatically rebin histograms"
)

# ------------------------------------------------------------------------------------------------------------------------
args = parser.parse_args()

masses = args.masses.split(",")

available_channels = ["ee", "mm", "em", "mt", "et", "tt"]

if args.channel not in available_channels:
    raise ValueError(
        "Invalid channel. Please choose from: {}".format(available_channels)
    )
available_methods = ["1", "2", "3", "4", "5", "6"]
if args.method not in available_methods:
    raise ValueError("Invalid method. Please choose from: {}".format(available_methods))

table = PrettyTable()
table.field_names = ["Details", "Choices"]
table.add_row(["Parameter File", args.parameter_file])
table.add_row(["Input Folder", args.input_folder])
table.add_row(["Output Folder", args.output_folder])
table.add_row(["Channel", args.channel])
table.add_row(["Era", args.era])
table.add_row(["Method", args.method])
table.add_row(["Category", args.category])
table.add_row(["Variable", args.var])
table.add_row(["Run Systematics", args.run_systematics])
table.add_row(["Selection", args.sel])
table.add_row(["Additional Weight", args.add_weight])
table.add_row(["Do SS", args.do_ss])
table.add_row(["Blind", args.blind])
table.add_row(["Masses", args.masses])
table.add_row(["Datacard Name", args.datacard_name])
table.add_row(["Auto Rebin", args.auto_rebin])

method = int(args.method)

# ------------------------------------------------------------------------------------------------------------------------
# Define baseline selections and different categories
# TODO: add option to change triggers
categories = {}
if args.era in ["Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"]:
    if args.channel == "ee":
        categories["baseline"] = (
            "(iso_1<0.15 && iso_2<0.15 && (trg_singleelectron && pt_1 > 31 && abs(eta_1) < 2.1))"
        )
    if args.channel == "mm":
        categories["baseline"] = (
            "(m_vis > 50 && iso_1<0.15 && iso_2<0.15 && (trg_singlemuon && pt_1 > 26 && abs(eta_1) < 2.4))"
        )
    if args.channel == "mt":
        mt_cross_only = "(trg_mt_cross && pt_1 > 21 && pt_1 <= 26 && abs(eta_1) < 2.1 && pt_2 > 32 && abs(eta_2) < 2.1)"
        single_muon_only = "(trg_singlemuon && pt_1 > 26  && abs(eta_1) < 2.4)"
        trg_full = single_muon_only #"(%s || %s)" % (mt_cross_only, single_muon_only)
        categories["baseline"] = (
            "(m_vis>40 && iso_1 < 0.15 && idDeepTau2018v2p5VSjet_2 >= 7 && idDeepTau2018v2p5VSe_2 >= 2 && idDeepTau2018v2p5VSmu_2 >= 4 && %s)"
            % trg_full
        )
        if args.do_aiso:
            categories["baseline"] = re.sub(
                "iso_1\s*<\s*0.15", "iso_1 > 0.05 && iso_1 < 0.2", categories["baseline"] # aiso cut based on what is used for the FF DR
            )
            
        categories["lt_ff_AR"] = categories["baseline"].replace(
            "idDeepTau2018v2p5VSjet_2 >= 7",
            "idDeepTau2018v2p5VSjet_2 < 7 && idDeepTau2018v2p5VSjet_2 >= 1",
        )
        
    if args.channel == "et":
        # et_cross_only = "(trg_et_cross && pt_1 > 25 && pt_1 < 31 && abs(eta_1) < 2.1 && pt_2 > 35 && abs(eta_2) < 2.1)"
        single_electron_only = "(trg_singleelectron && pt_1 >= 32 && abs(eta_1) < 2.1 )"
        trg_full = single_electron_only # remove et cross trigger until Nanoprod v3
        categories["baseline"] = ( # Tight VSe for et
            "(m_vis>40 && iso_1 < 0.15 && idDeepTau2018v2p5VSjet_2 >= 7 && idDeepTau2018v2p5VSe_2 >= 6 && idDeepTau2018v2p5VSmu_2 >= 4 && %s)"
            % trg_full
        )
        if args.do_aiso:
            categories["baseline"] = re.sub(
                "iso_1\s*<\s*0.15", "iso_1 > 0.05 && iso_1 < 0.2", categories["baseline"]  # aiso cut based on what is used for the FF DR
            )

        categories["lt_ff_AR"] = categories["baseline"].replace(
            "idDeepTau2018v2p5VSjet_2 >= 7",
            "idDeepTau2018v2p5VSjet_2 < 7 && idDeepTau2018v2p5VSjet_2 >= 1",
        )

    if args.channel == "tt":
        doubletau_only_trg = "(trg_doubletau && pt_1 > 40 && pt_2 > 40)"
        doubletaujet_only_trg = "(trg_doubletauandjet && pt_1 > 35 && pt_2 > 35 && jpt_1 > 60)"  # might need to revise jet cut later on
        trg_full = "(%s || %s)" % (doubletau_only_trg, doubletaujet_only_trg)
        categories["baseline"] = (
            "(m_vis > 40 && idDeepTau2018v2p5VSjet_1 >= 7 && idDeepTau2018v2p5VSjet_2 >= 7 && idDeepTau2018v2p5VSe_1 >= 2 && idDeepTau2018v2p5VSe_2 >= 2 && idDeepTau2018v2p5VSmu_1 >= 4 && idDeepTau2018v2p5VSmu_2 >= 4 && %s)"
            % trg_full
        )

        if args.do_aiso:
            categories["baseline"] = categories["baseline"].replace(
                "idDeepTau2018v2p5VSjet_2 >= 7",
                "idDeepTau2018v2p5VSjet_2 < 7 && idDeepTau2018v2p5VSjet_2 >= 1",
            )

        categories["tt_qcd_norm"] = categories["baseline"].replace(
            "idDeepTau2018v2p5VSjet_1 >= 7",
            "idDeepTau2018v2p5VSjet_1 < 7 && idDeepTau2018v2p5VSjet_1 >= 1",
        )
        categories["tt_ff_AR"] = categories["baseline"].replace(
            "idDeepTau2018v2p5VSjet_1 >= 7",
            "idDeepTau2018v2p5VSjet_1 < 7 && idDeepTau2018v2p5VSjet_1 >= 1",
        )
        categories["subleadfake"] = (
            categories["baseline"] + "&& genPartFlav_1 != 0 && genPartFlav_2 == 0"
        )

categories["inclusive"] = "(1)"
categories["nobtag"] = "(n_bjets==0)"
categories["btag"] = "(n_bjets>=1)"
categories["w_sdb"] = "mt_1>70."
categories["w_shape"] = ""
categories["qcd_loose_shape"] = re.sub(
    "iso_1\s*<\s*0.15", "iso_1 < 0.3", categories["baseline"]
)

categories["ip_control"] = "(m_vis > 40 && m_vis < 90 && mt_1 < 40)"

categories["xt_dM0"] = "(decayMode_2 == 0)"
categories["xt_dM1"] = "(decayMode_2 == 1)"
categories["xt_dM10"] = "(decayMode_2 == 10)"
categories["xt_dM11"] = "(decayMode_2 == 11)"

if args.channel == "tt":

    sel_pi = "decayModePNet_X==0 && ip_LengthSig_X>=1.25"
    sel_rho = "decayMode_X==1 && decayModePNet_X==1 && pion_E_split_X>0.2"
    sel_a11pr = "decayMode_X==1 && decayModePNet_X==2 && pion_E_split_X>0.2"
    sel_a1 = "decayModePNet_X==10 && hasRefitSV_X"

    sel_pi_1 = sel_pi.replace("X", "1")
    sel_pi_2 = sel_pi.replace("X", "2")
    sel_rho_1 = sel_rho.replace("X", "1")
    sel_rho_2 = sel_rho.replace("X", "2")
    sel_a1_1 = sel_a1.replace("X", "1")
    sel_a1_2 = sel_a1.replace("X", "2")
    sel_a11pr_1 = sel_a11pr.replace("X", "1")
    sel_a11pr_2 = sel_a11pr.replace("X", "2")

    categories["cp_inclusive"] = (
        f"( ({sel_pi_1}) || ({sel_rho_1}) || ({sel_a1_1}) || ({sel_a11pr_1}) ) && ( ({sel_pi_2}) || ({sel_rho_2}) || ({sel_a1_2}) || ({sel_a11pr_2}) )"
    )

    categories["inclusive_PNet_rhorho"] = "(%(sel_rho_1)s && %(sel_rho_2)s)" % vars()
    categories["inclusive_PNet_pipi"] = "(%(sel_pi_1)s && %(sel_pi_2)s)" % vars()
    categories["inclusive_PNet_a1a1"] = "(%(sel_a1_1)s && %(sel_a1_2)s)" % vars()
    categories["inclusive_PNet_rhoa11pr"] = (
        "((%(sel_rho_1)s && %(sel_a11pr_2)s) || (%(sel_a11pr_1)s && %(sel_rho_2)s) || (%(sel_a11pr_1)s && %(sel_a11pr_2)s))"
        % vars()
    )

    categories["inclusive_PNet_pirho"] = "(%(sel_pi_1)s && %(sel_rho_2)s)" % vars()
    categories["inclusive_PNet_a1rho"] = "(%(sel_a1_1)s && %(sel_rho_2)s)" % vars()
    categories["inclusive_PNet_a1pi"] = "(%(sel_a1_1)s && %(sel_pi_2)s)" % vars()
    categories["inclusive_PNet_pia11pr"] = "(%(sel_pi_1)s && %(sel_a11pr_2)s)" % vars()
    categories["inclusive_PNet_a1a11pr"] = "(%(sel_a1_1)s && %(sel_a11pr_2)s)" % vars()

    categories["inclusive_PNet_rhopi"] = "(%(sel_rho_1)s && %(sel_pi_2)s)" % vars()
    categories["inclusive_PNet_rhoa1"] = "(%(sel_rho_1)s && %(sel_a1_2)s)" % vars()
    categories["inclusive_PNet_pia1"] = "(%(sel_pi_1)s && %(sel_a1_2)s)" % vars()
    categories["inclusive_PNet_a11prpi"] = "(%(sel_a11pr_1)s && %(sel_pi_2)s)" % vars()
    categories["inclusive_PNet_a11pra1"] = "(%(sel_a11pr_1)s && %(sel_a1_2)s)" % vars()

    categories["mva_higgs"] = f"(BDT_pred_class==1) && ( ({sel_pi_1}) || ({sel_rho_1}) || ({sel_a1_1}) || ({sel_a11pr_1}) ) && ( ({sel_pi_2}) || ({sel_rho_2}) || ({sel_a1_2}) || ({sel_a11pr_2}) )"
    categories["mva_fake"] = f"(BDT_pred_class==2) && ( ({sel_pi_1}) || ({sel_rho_1}) || ({sel_a1_1}) || ({sel_a11pr_1}) ) && ( ({sel_pi_2}) || ({sel_rho_2}) || ({sel_a1_2}) || ({sel_a11pr_2}) )"
    categories["mva_tau"] = f"(BDT_pred_class==0) && ( ({sel_pi_1}) || ({sel_rho_1}) || ({sel_a1_1}) || ({sel_a11pr_1}) ) && ( ({sel_pi_2}) || ({sel_rho_2}) || ({sel_a1_2}) || ({sel_a11pr_2}) )"

    tt_channels = [
        "rhorho",
        "pirho",
        "rhopi",
        "a1rho",
        "rhoa1",
        "a1pi",
        "pia1",
        "a1a1",
        "pipi",
        "pia11pr",
        "a11prpi",
        "rhoa11pr",
        "a1a11pr",
        "a11pra1",
    ]

    for c in tt_channels:
        categories[f"higgs_{c}"] = f"((BDT_pred_class==1) && {categories[f'inclusive_PNet_{c}']})"
        categories[f"tau_{c}"] = f"((BDT_pred_class==0) && {categories[f'inclusive_PNet_{c}']})"
        categories[f"fake_{c}"] = f"((BDT_pred_class==2) && {categories[f'inclusive_PNet_{c}']})"


    # Define MVA categories by channels instead
    # categories['mva_tau'] = f"((BDT_pred_class==0) && ({categories['inclusive_PNet_rhorho']} || {categories['inclusive_PNet_pirho']} || {categories['inclusive_PNet_rhopi']} || {categories['inclusive_PNet_a1rho']} || {categories['inclusive_PNet_rhoa1']} || {categories['inclusive_PNet_a1pi']} || {categories['inclusive_PNet_pia1']} || {categories['inclusive_PNet_a1a1']} || {categories['inclusive_PNet_pipi']} || {categories['inclusive_PNet_pia11pr']} || {categories['inclusive_PNet_a11prpi']} || {categories['inclusive_PNet_rhoa11pr']} || {categories['inclusive_PNet_a1a11pr']} || {categories['inclusive_PNet_a11pra1']}) )"
    # categories['mva_fake'] = f"((BDT_pred_class==2) && ({categories['inclusive_PNet_rhorho']} || {categories['inclusive_PNet_pirho']} || {categories['inclusive_PNet_rhopi']} || {categories['inclusive_PNet_a1rho']} || {categories['inclusive_PNet_rhoa1']} || {categories['inclusive_PNet_a1pi']} || {categories['inclusive_PNet_pia1']} || {categories['inclusive_PNet_a1a1']} || {categories['inclusive_PNet_pipi']} || {categories['inclusive_PNet_pia11pr']} || {categories['inclusive_PNet_a11prpi']} || {categories['inclusive_PNet_rhoa11pr']} || {categories['inclusive_PNet_a1a11pr']} || {categories['inclusive_PNet_a11pra1']}) )"
    # categories['mva_higgs'] = f"((BDT_pred_class==1) && ({categories['inclusive_PNet_rhorho']} || {categories['inclusive_PNet_pirho']} || {categories['inclusive_PNet_rhopi']} || {categories['inclusive_PNet_a1rho']} || {categories['inclusive_PNet_rhoa1']} || {categories['inclusive_PNet_a1pi']} || {categories['inclusive_PNet_pia1']} || {categories['inclusive_PNet_a1a1']} || {categories['inclusive_PNet_pipi']} || {categories['inclusive_PNet_pia11pr']} || {categories['inclusive_PNet_a11prpi']} || {categories['inclusive_PNet_rhoa11pr']} || {categories['inclusive_PNet_a1a11pr']} || {categories['inclusive_PNet_a11pra1']}) )"


    # alpha minus validation angles
    for opt, opt_cut in zip(['high', 'low'], ['>', '<']):
        for obj in [1, 2]:
            categories[f'aminus_rho_rho_{obj}_{opt}'] = f"(alphaAngle_rho_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_rhorho']}))"
            categories[f'aminus_pi_rho_{obj}_{opt}'] = f"(alphaAngle_pi_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_pirho']}))"
            categories[f'aminus_rho_pi_{obj}_{opt}'] = f"(alphaAngle_rho_pi_{obj} {opt_cut} {np.pi/4} && ({categories['tau_rhopi']}))"
            categories[f'aminus_a1_rho_{obj}_{opt}'] = f"(alphaAngle_a1_rho_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['tau_a1rho']}))"
            categories[f'aminus_rho_a1_{obj}_{opt}'] = f"(alphaAngle_rho_a1_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['tau_rhoa1']}))"
            categories[f'aminus_a1_pi_{obj}_{opt}'] = f"(alphaAngle_a1_pi_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['tau_a1pi']}))"
            categories[f'aminus_pi_a1_{obj}_{opt}'] = f"(alphaAngle_pi_a1_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['tau_pia1']}))"
            categories[f'aminus_a1_a1_{obj}_{opt}'] = f"(alphaAngle_a1_a1_{obj} {opt_cut} {np.pi/4} && ({categories['tau_a1a1']}))"
            categories[f'aminus_pi_pi_{obj}_{opt}'] = f"(alphaAngle_pi_pi_{obj} {opt_cut} {np.pi/4} && ({categories['tau_pipi']}))"
            categories[f'aminus_pi_a11pr_{obj}_{opt}'] = f"(alphaAngle_pi_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_pia11pr']}))"
            categories[f'aminus_a11pr_pi_{obj}_{opt}'] = f"(alphaAngle_rho_pi_{obj} {opt_cut} {np.pi/4} && ({categories['tau_a11prpi']}))"
            categories[f'aminus_rho_a11pr_{obj}_{opt}'] = f"(alphaAngle_rho_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_rhoa11pr']}))"
            categories[f'aminus_a11pr_rho_{obj}_{opt}'] = f"(alphaAngle_rho_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_a1a11pr']}))"
            categories[f'aminus_a1_a11pr_{obj}_{opt}'] = f"(alphaAngle_a1_rho_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['tau_a11pra1']}))"
            categories[f'aminus_a11pr_a1_{obj}_{opt}'] = f"(alphaAngle_rho_a1_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['tau_a11pra1']}))"



elif args.channel == "mt":
    # PNet DM separated categories for SFs (with CP selections)
    categories["DM0_tau_cp"] = "decayModePNet_2 == 0 && ip_LengthSig_2 >= 1.25"
    categories["DM1_tau_cp"] = "decayMode_2==1 && decayModePNet_2 == 1 && pion_E_split_2 > 0.2"
    categories["DM2_tau_cp"] = "decayMode_2==1 && decayModePNet_2 == 2 && pion_E_split_2 > 0.2"
    categories["DM10_tau_cp"] = "decayModePNet_2 == 10 && hasRefitSV_2"
    categories["DM11_tau_cp"] = "decayModePNet_2 == 11 && hasRefitSV_2"

    # HPS DM separate categories (inclusive)
    categories["DM0_tau"] = "decayMode_2 == 0"
    categories["DM1_tau"] = "decayMode_2==1"
    categories["DM10_tau"] = "decayMode_2==10"
    categories["DM11_tau"] = "decayMode_2==11"

    # CP measurement categories
    categories["sel_mupi"] = "ip_LengthSig_1>1.0 && decayModePNet_2 == 0 && ip_LengthSig_2 >= 1.25"
    categories["sel_murho"] = "ip_LengthSig_1>1.0 && decayMode_2==1 && decayModePNet_2 == 1 && pion_E_split_2 > 0.2"
    categories["sel_mua11pr"] = "ip_LengthSig_1>1.0 && decayMode_2==1 && decayModePNet_2 == 2 && pion_E_split_2 > 0.2"
    categories["sel_mua1"] = "ip_LengthSig_1>1.0 && decayModePNet_2 == 10 && hasRefitSV_2"
    # inclusive category including an mT cut to suppress W+jets
    categories['cp_inclusive'] = f"(({categories['sel_mupi']}) || ({categories['sel_murho']}) || ({categories['sel_mua11pr']}) || ({categories['sel_mua1']}))"

    categories["mva_higgs"] = f"(BDT_pred_class==1) && ({categories['cp_inclusive']})"
    categories["mva_fake"] = f"(BDT_pred_class==2) && ({categories['cp_inclusive']})"
    categories["mva_tau"] = f"(BDT_pred_class==0) && ({categories['cp_inclusive']})"

    lt_channels = [
        "mupi",
        "murho",
        "mua11pr",
        "mua1"
    ]
    for c in lt_channels:
        categories[f"higgs_{c}"] = f"({categories['mva_higgs']} && {categories[f'sel_{c}']})"
        categories[f"tau_{c}"] = f"({categories['mva_tau']} && {categories[f'sel_{c}']})"
        categories[f"fake_{c}"] = f"({categories['mva_fake']} && {categories[f'sel_{c}']})"

    # alpha minus validation angles
    for opt, opt_cut in zip(['high', 'low'], ['>', '<']):
        for obj in [1, 2]:
            categories[f'aminus_mu_pi_{obj}_{opt}'] = f"(alphaAngle_mu_pi_{obj} {opt_cut} {np.pi/4} && ({categories['tau_mupi']}))"
            categories[f'aminus_mu_rho_{obj}_{opt}'] = f"(alphaAngle_mu_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_murho']}))"
            categories[f'aminus_mu_a11pr_{obj}_{opt}'] = f"(alphaAngle_mu_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_mua11pr']}))"
            categories[f'aminus_mu_a1_{obj}_{opt}'] = f"(alphaAngle_mu_a1_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['tau_mua1']}))"
            categories[f'BDTinc_aminus_mu_pi_{obj}_{opt}'] = f"(alphaAngle_mu_pi_{obj} {opt_cut} {np.pi/4} && ({categories['sel_mupi']}))"
            categories[f'BDTinc_aminus_mu_rho_{obj}_{opt}'] = f"(alphaAngle_mu_rho_{obj} {opt_cut} {np.pi/4} && ({categories['sel_murho']}))"
            categories[f'BDTinc_aminus_mu_a11pr_{obj}_{opt}'] = f"(alphaAngle_mu_rho_{obj} {opt_cut} {np.pi/4} && ({categories['sel_mua11pr']}))"
            categories[f'BDTinc_aminus_mu_a1_{obj}_{opt}'] = f"(alphaAngle_mu_a1_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['sel_mua1']}))"


elif args.channel == "et":
    # PNet DM separated categories for SFs (with CP selections)
    categories["DM0_tau_cp"] = "decayModePNet_2 == 0 && ip_LengthSig_2 >= 1.25"
    categories["DM1_tau_cp"] = "decayMode_2==1 && decayModePNet_2 == 1 && pion_E_split_2 > 0.2"
    categories["DM2_tau_cp"] = "decayMode_2==1 && decayModePNet_2 == 2 && pion_E_split_2 > 0.2"
    categories["DM10_tau_cp"] = "decayModePNet_2 == 10 && hasRefitSV_2"
    categories["DM11_tau_cp"] = "decayModePNet_2 == 11 && hasRefitSV_2"

    # HPS DM separate categories (inclusive)
    categories["DM0_tau"] = "decayMode_2 == 0"
    categories["DM1_tau"] = "decayMode_2==1"
    categories["DM10_tau"] = "decayMode_2==10"
    categories["DM11_tau"] = "decayMode_2==11"

    # CP measurement categories
    categories["sel_epi"] = "ip_LengthSig_1>1.0 && decayModePNet_2 == 0 && ip_LengthSig_2 >= 1.25"
    categories["sel_erho"] = "ip_LengthSig_1>1.0 && decayMode_2==1 && decayModePNet_2 == 1 && pion_E_split_2 > 0.2"
    categories["sel_ea11pr"] = "ip_LengthSig_1>1.0 && decayMode_2==1 && decayModePNet_2 == 2 && pion_E_split_2 > 0.2"
    categories["sel_ea1"] = "ip_LengthSig_1>1.0 && decayModePNet_2 == 10 && hasRefitSV_2"
    # inclusive category including an mT cut to suppress W+jets
    categories['cp_inclusive'] = f"(({categories['sel_epi']}) || ({categories['sel_erho']}) || ({categories['sel_ea11pr']}) || ({categories['sel_ea1']}))"

    categories["mva_higgs"] = f"(BDT_pred_class==1) && ({categories['cp_inclusive']})"
    categories["mva_fake"] = f"(BDT_pred_class==2) && ({categories['cp_inclusive']})"
    categories["mva_tau"] = f"(BDT_pred_class==0) && ({categories['cp_inclusive']})"

    lt_channels = [
        "epi",
        "erho",
        "ea11pr",
        "ea1"
    ]
    for c in lt_channels:
        categories[f"higgs_{c}"] = f"({categories['mva_higgs']} && {categories[f'sel_{c}']})"
        categories[f"tau_{c}"] = f"({categories['mva_tau']} && {categories[f'sel_{c}']})"
        categories[f"fake_{c}"] = f"({categories['mva_fake']} && {categories[f'sel_{c}']})"

    # alpha minus validation angles
    for opt, opt_cut in zip(['high', 'low'], ['>', '<']):
        for obj in [1, 2]:
            categories[f'aminus_e_pi_{obj}_{opt}'] = f"(alphaAngle_e_pi_{obj} {opt_cut} {np.pi/4} && ({categories['tau_epi']}))"
            categories[f'aminus_e_rho_{obj}_{opt}'] = f"(alphaAngle_e_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_erho']}))"
            categories[f'aminus_e_a11pr_{obj}_{opt}'] = f"(alphaAngle_e_rho_{obj} {opt_cut} {np.pi/4} && ({categories['tau_ea11pr']}))"
            categories[f'aminus_e_a1_{obj}_{opt}'] = f"(alphaAngle_e_a1_FASTMTT_MassConstraint_{obj} {opt_cut} {np.pi/4} && ({categories['tau_ea1']}))"


# if args.set_alias is not None then overwrite the categories with the selection provided

if args.set_alias is not None:
    for i in args.set_alias:
        cat_to_overwrite = i.split(":")[0]
        cat_to_overwrite = cat_to_overwrite.replace('"', "")
        overwrite_with = i.split(":")[1]
        overwrite_with = overwrite_with.replace('"', "")
        start_index = overwrite_with.find("{")
        end_index = overwrite_with.find("}")
        while start_index > 0:
            replace_with = overwrite_with[start_index : end_index + 1]
            replace_with = replace_with.replace("{", "")
            replace_with = replace_with.replace("}", "")
            replace_string = categories[replace_with]
            overwrite_with = (
                overwrite_with[0:start_index]
                + replace_string
                + overwrite_with[end_index + 1 :]
            )
            start_index = overwrite_with.find("{")
            end_index = overwrite_with.find("}")

        print(
            'Overwriting alias: "'
            + cat_to_overwrite
            + '" with selection: "'
            + overwrite_with
            + '"'
        )
        if cat_to_overwrite == "sel":
            args.sel = overwrite_with
        else:
            categories[cat_to_overwrite] = overwrite_with

# ------------------------------------------------------------------------------------------------------------------------
# Define the samples (Data and MC (Background & Signal))
if args.era in ["Run3_2022", "Run3_2022EE", "Run3_2023", "Run3_2023BPix"]:
    samples_dict = {}
    # Data Samples
    if args.era in ["Run3_2022"]:
        if args.channel in ["ee", "et"]:
            data_samples = ["EGamma_Run2022C", "EGamma_Run2022D"]
        elif args.channel in ["mm", "mt"]:
            data_samples = ["SingleMuon_Run2022C", "Muon_Run2022C", "Muon_Run2022D"]
        elif args.channel == "tt":
            data_samples = ["Tau_Run2022C", "Tau_Run2022D"]
    elif args.era in ["Run3_2022EE"]:
        if args.channel in ["ee", "et"]:
            data_samples = ["EGamma_Run2022E", "EGamma_Run2022F", "EGamma_Run2022G"]
        elif args.channel in ["mm", "mt"]:
            data_samples = ["Muon_Run2022E", "Muon_Run2022F", "Muon_Run2022G"]
        elif args.channel == "tt":
            data_samples = ["Tau_Run2022E", "Tau_Run2022F", "Tau_Run2022G"]
    elif args.era in ["Run3_2023"]:
        if args.channel in ["ee", "et"]:
            data_samples = [
                "EGamma0_Run2023C_v1",
                "EGamma0_Run2023C_v2",
                "EGamma0_Run2023C_v3",
                "EGamma0_Run2023C_v4",
                "EGamma1_Run2023C_v1",
                "EGamma1_Run2023C_v2",
                "EGamma1_Run2023C_v3",
                "EGamma1_Run2023C_v4",
            ]
        elif args.channel in ["mm", "mt"]:
            data_samples = [
                "Muon0_Run2023C_v1",
                "Muon0_Run2023C_v2",
                "Muon0_Run2023C_v3",
                "Muon0_Run2023C_v4",
                "Muon1_Run2023C_v1",
                "Muon1_Run2023C_v2",
                "Muon1_Run2023C_v3",
                "Muon1_Run2023C_v4",
            ]
        elif args.channel == "tt":
            data_samples = [
                "Tau_Run2023C_v1",
                "Tau_Run2023C_v2",
                "Tau_Run2023C_v3",
                "Tau_Run2023C_v4",
            ]
    elif args.era in ["Run3_2023BPix"]:
        if args.channel in ["ee", "et"]:
            data_samples = [
                "EGamma0_Run2023D_v1",
                "EGamma0_Run2023D_v2",
                "EGamma1_Run2023D_v1",
                "EGamma1_Run2023D_v2",
            ]
        elif args.channel in ["mm", "mt"]:
            data_samples = [
                "Muon0_Run2023D_v1",
                "Muon0_Run2023D_v2",
                "Muon1_Run2023D_v1",
                "Muon1_Run2023D_v2",
            ]
        elif args.channel == "tt":
            data_samples = ["Tau_Run2023D_v1", "Tau_Run2023D_v2"]

    samples_dict["data_samples"] = data_samples

    # MC Samples
    if args.LO_DY:
        print("WARNING: Using LO DY samples")
        ztt_samples = [
            "DYto2L_M_50_madgraphMLM",
            "DYto2L_M_50_madgraphMLM_ext1",
            "DYto2L_M_50_1J_madgraphMLM",
            "DYto2L_M_50_2J_madgraphMLM",
            "DYto2L_M_50_3J_madgraphMLM",
            "DYto2L_M_50_4J_madgraphMLM",
        ]
        if args.era in ["Run3_2023", "Run3_2023BPix"]:
            ztt_samples.remove("DYto2L_M_50_madgraphMLM_ext1")

    elif args.NLO_DY:
        print("WARNING: Using NLO DY samples")
        ztt_samples = [
            "DYto2L_M_50_amcatnloFXFX",
            "DYto2L_M_50_amcatnloFXFX_ext1",
            "DYto2L_M_50_0J_amcatnloFXFX",
            "DYto2L_M_50_1J_amcatnloFXFX",
            "DYto2L_M_50_2J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_40to100_1J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_100to200_1J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_200to400_1J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_400to600_1J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_600_1J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_40to100_2J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_100to200_2J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_200to400_2J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_400to600_2J_amcatnloFXFX",
            "DYto2L_M_50_PTLL_600_2J_amcatnloFXFX",
        ]  # use NLO samples
        if args.era in ["Run3_2023", "Run3_2023BPix"]:
            ztt_samples.remove("DYto2L_M_50_amcatnloFXFX_ext1")

    else:
        print("Using New DY samples")
        ztt_samples = [
            "DYto2Tau_MLL_50_0J_amcatnloFXFX",
            "DYto2Tau_MLL_50_1J_amcatnloFXFX",
            "DYto2Tau_MLL_50_2J_amcatnloFXFX"
        ]
        if args.use_filtered_DY:
            print(f"WARNING: Will use filtered DY, and read effective events from alternate file")
            ztt_samples += [
                "DYto2Tau_MLL_50_0J_Filtered_amcatnloFXFX",
                "DYto2Tau_MLL_50_1J_Filtered_amcatnloFXFX",
                "DYto2Tau_MLL_50_2J_Filtered_amcatnloFXFX"
            ]
        zll_samples = [
            "DYto2L_M_50_amcatnloFXFX",
            "DYto2L_M_50_amcatnloFXFX_ext1",
            "DYto2L_M_50_0J_amcatnloFXFX",
            "DYto2L_M_50_1J_amcatnloFXFX",
            "DYto2L_M_50_2J_amcatnloFXFX"
        ]
        if args.era in ["Run3_2023", "Run3_2023BPix"]:
            zll_samples.remove("DYto2L_M_50_amcatnloFXFX_ext1")

    top_samples = [
        "TTto2L2Nu",
        "TTto2L2Nu_ext1",
        "TTtoLNu2Q",
        "TTtoLNu2Q_ext1"
    ]
    vv_samples = [
        "WW",
        "WZ",
        "ZZ",
        "ST_t_channel_top_4f_InclusiveDecays",
        "ST_t_channel_antitop_4f_InclusiveDecays",
        "ST_tW_top_2L2Nu",
        "ST_tW_top_2L2Nu_ext1",
        "ST_tW_antitop_2L2Nu",
        "ST_tW_antitop_2L2Nu_ext1",
        "ST_tW_top_LNu2Q",
        "ST_tW_top_LNu2Q_ext1",
        "ST_tW_antitop_LNu2Q",
        "ST_tW_antitop_LNu2Q_ext1",
    ]
    ewkz_samples = [
        "EWKZ_MLL-50_TuneCP5_13p6TeV_madgraph-pythia8"
    ]
    wjets_samples = [
        "WtoLNu_madgraphMLM",
        "WtoLNu_madgraphMLM_ext1",
        "WtoLNu_1J_madgraphMLM",
        "WtoLNu_2J_madgraphMLM",
        "WtoLNu_3J_madgraphMLM",
        "WtoLNu_4J_madgraphMLM",
    ]

    if args.era in ["Run3_2023", "Run3_2023BPix"]:
        top_samples.remove("TTto2L2Nu_ext1")
        top_samples.remove("TTtoLNu2Q_ext1")
        vv_samples.remove("ST_tW_top_2L2Nu_ext1")
        vv_samples.remove("ST_tW_antitop_2L2Nu_ext1")
        vv_samples.remove("ST_tW_top_LNu2Q_ext1")
        vv_samples.remove("ST_tW_antitop_LNu2Q_ext1")
        wjets_samples.remove("WtoLNu_madgraphMLM_ext1")

    if args.channel in ["et", "mt", "tt"]:
        signal_samples = {
            # Unfiltered samples
            # "qqH_sm_unfiltered_htt*": "VBFHToTauTau_UncorrelatedDecay_UnFiltered",
            # "WH_sm_unfiltered_htt*": [
            #     "WplusHToTauTau_UncorrelatedDecay_UnFiltered",
            #     "WminusHToTauTau_UncorrelatedDecay_UnFiltered",
            # ],
            # "ZH_sm_unfiltered_htt*": "ZHToTauTau_UncorrelatedDecay_UnFiltered",

            # VBF samples
            "qqH_sm_htt*": "VBFHToTauTau_UncorrelatedDecay_Filtered",
            "qqH_ps_htt*": "VBFHToTauTau_UncorrelatedDecay_Filtered",
            "qqH_mm_htt*": "VBFHToTauTau_UncorrelatedDecay_Filtered",

            # ggH samples (SM production)
            "ggH_sm_prod_sm_htt*": "GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay",
            "ggH_ps_prod_sm_htt*": "GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay",
            "ggH_mm_prod_sm_htt*": "GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay",
            # ggH samples (PS production)
            "ggH_sm_prod_ps_htt*": "GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay",
            "ggH_ps_prod_ps_htt*": "GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay",
            "ggH_mm_prod_ps_htt*": "GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay",
            # ggH samples (MM production)
            "ggH_sm_prod_mm_htt*": "GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay",
            "ggH_ps_prod_mm_htt*": "GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay",
            "ggH_mm_prod_mm_htt*": "GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay",

            # WH samples
            "WH_sm_htt*": [
                "WplusHToTauTau_UncorrelatedDecay_Filtered",
                "WminusHToTauTau_UncorrelatedDecay_Filtered",
            ],
            "WH_ps_htt*": [
                "WplusHToTauTau_UncorrelatedDecay_Filtered",
                "WminusHToTauTau_UncorrelatedDecay_Filtered",
            ],
            "WH_mm_htt*": [
                "WplusHToTauTau_UncorrelatedDecay_Filtered",
                "WminusHToTauTau_UncorrelatedDecay_Filtered",
            ],

            # ZH samples
            "ZH_sm_htt*": "ZHToTauTau_UncorrelatedDecay_Filtered",
            "ZH_ps_htt*": "ZHToTauTau_UncorrelatedDecay_Filtered",
            "ZH_mm_htt*": "ZHToTauTau_UncorrelatedDecay_Filtered",

            # Higgs Flat samples (uncorrelated decay)
            "Higgs_flat_htt*": [
                "VBFHToTauTau_UncorrelatedDecay_Filtered",
                "GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay",
                "GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay",
                "GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay",
                "ZHToTauTau_UncorrelatedDecay_Filtered",
                "WplusHToTauTau_UncorrelatedDecay_Filtered",
                "WminusHToTauTau_UncorrelatedDecay_Filtered"
            ],
            "qqH_flat_htt*": "VBFHToTauTau_UncorrelatedDecay_Filtered",
            "ggH_flat_prod_sm_htt*": "GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay",
            "ggH_flat_prod_ps_htt*": "GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay",
            "ggH_flat_prod_mm_htt*": "GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay",
            "WH_flat_htt*": [
                "WplusHToTauTau_UncorrelatedDecay_Filtered",
                "WminusHToTauTau_UncorrelatedDecay_Filtered",
            ],
            "ZH_flat_htt*": "ZHToTauTau_UncorrelatedDecay_Filtered",

            # combined ggH (events from all production mechanisms)
            # "ggH_sm_prod_sm_reweight_htt*": [
            #     'GluGluHTo2Tau_UncorrelatedDecay_CPodd_Filtered_ProdAndDecay',
            #     'GluGluHTo2Tau_UncorrelatedDecay_MM_Filtered_ProdAndDecay',
            #     'GluGluHTo2Tau_UncorrelatedDecay_SM_Filtered_ProdAndDecay'
            # ],
        }


    else:
        signal_samples = {}

    # # TODO: REMOVE THIS IS TEMPORARY FOR TAU ID SFs
    # signal_samples = {}

    samples_dict["ztt_samples"] = ztt_samples
    samples_dict["zll_samples"] = zll_samples
    samples_dict["top_samples"] = top_samples
    samples_dict["vv_samples"] = vv_samples
    samples_dict["wjets_samples"] = wjets_samples
    samples_dict["ewkz_samples"] = ewkz_samples
    samples_dict["signal_samples"] = signal_samples
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Define gen selections to separate MC samples by Gen Flags

# Hint:
# Muon Gen Matching:
# 1 = prompt muon
# 15 = muon from prompt tau
# 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched

# Electron Gen Matching:
# 1 = prompt electron (including gamma*->mu mu)
# 15 = electron from prompt tau
# 22 = prompt photon (likely conversion)
# 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched

# Tau Gen Matching:
# 1 = prompt electron
# 2 = prompt muon
# 3 = tau->e decay
# 4 = tau->mu decay
# 5 = hadronic tau decay
# 0 = unknown or unmatched

gen_sels = {}
gen_sels_dict = {}
if args.channel in ["ee", "mm"]:
    gen_sels["ll_sel"] = "(genPartFlav_1==1 & genPartFlav_2==1)"
    gen_sels["tt_sel"] = "(genPartFlav_1==15 & genPartFlav_2==15)"
    gen_sels["j_sel"] = (
        "(!(" + gen_sels["ll_sel"] + ") && !(" + gen_sels["tt_sel"] + "))"
    )

    z_sels = {}
    z_sels["ztt_sel"] = gen_sels["tt_sel"]
    z_sels["zl_sel"] = gen_sels["ll_sel"]
    z_sels["zj_sel"] = gen_sels["j_sel"]

    top_sels = {}
    top_sels["ttt_sel"] = gen_sels["ll_sel"]
    top_sels["ttj_sel"] = "!(" + gen_sels["ll_sel"] + ")"

    vv_sels = {}
    vv_sels["vvt_sel"] = gen_sels["ll_sel"]
    vv_sels["vvj_sel"] = "!(" + gen_sels["ll_sel"] + ")"

if args.channel in ["mt", "et"]:
    gen_sels["tt_sel"] = "(genPartFlav_1==15 && genPartFlav_2==5)"
    gen_sels["ll_sel"] = f"(genPartFlav_2!=0 && !{gen_sels['tt_sel']})"
    gen_sels["j_sel"] = "(genPartFlav_2==0)"

    z_sels = {}
    z_sels["ztt_sel"] = gen_sels["tt_sel"]
    z_sels["zl_sel"] = gen_sels["ll_sel"]
    z_sels["zj_sel"] = gen_sels["j_sel"]

    top_sels = {}
    top_sels["ttt_sel"] = "!(" + gen_sels["j_sel"] + ")"
    top_sels["ttj_sel"] = gen_sels["j_sel"]

    vv_sels = {}
    vv_sels["vvt_sel"] = "!(" + gen_sels["j_sel"] + ")"
    vv_sels["vvj_sel"] = gen_sels["j_sel"]

if args.channel in ["tt"]:
    gen_sels["tt_sel"] = "(genPartFlav_1==5 && genPartFlav_2==5)"
    gen_sels["ll_sel"] = (
        f"(!(genPartFlav_1==0 || genPartFlav_2==0) && !{gen_sels['tt_sel']})"
    )
    gen_sels["j_sel"] = "(genPartFlav_1==0 || genPartFlav_2==0)"

    z_sels = {}
    z_sels["ztt_sel"] = gen_sels["tt_sel"]
    z_sels["zl_sel"] = gen_sels["ll_sel"]
    z_sels["zj_sel"] = gen_sels["j_sel"]

    top_sels = {}
    top_sels["ttt_sel"] = "!(" + gen_sels["j_sel"] + ")"
    top_sels["ttj_sel"] = gen_sels["j_sel"]

    vv_sels = {}
    vv_sels["vvt_sel"] = "!(" + gen_sels["j_sel"] + ")"
    vv_sels["vvj_sel"] = gen_sels["j_sel"]

gen_sels_dict["z_sels"] = z_sels
gen_sels_dict["top_sels"] = top_sels
gen_sels_dict["vv_sels"] = vv_sels
# ------------------------------------------------------------------------------------------------------------------------

# RunPlotting handles how each process is added to the analysis


def RunPlotting(
    ana,
    nodename,
    samples_dict,
    gen_sels_dict,
    systematic="",
    cat_name="",
    categories={},
    categories_unmodified={},
    sel="",
    add_name="",
    wt="wt",
    do_data=True,
    qcd_factor=1.0,
    method=1,
    nodes_to_skip=[],
):
    """
    RunPlotting handles how each process is added to the analysis
    ana: Analysis object
    nodename: name of the node
    samples_dict: dictionary containing the samples
    gen_sels_dict: dictionary containing the gen selections
    systematic: systematic variation
    cat_name: the name of the category to be used for the data selection
    categories: dictionary containing the categories
    categories_unmodified: dictionary containing the categories - unmodified by any systematic variations - to be applied for data only
    sel: additional selection
    add_name: additional name to be added to the process (e.g. ZTT + xxx)
    wt: weight to be applied
    do_data: boolean to decide if data should be added
    qcd_factor: factor to be applied to QCD
    """

    cat = categories["cat"]
    cat_data = categories_unmodified["cat"]

    doZL = True if "ZL" not in nodes_to_skip else False
    doZJ = True if "ZJ" not in nodes_to_skip else False
    doTTT = True if "TTT" not in nodes_to_skip else False
    doTTJ = True if "TTJ" not in nodes_to_skip else False
    doVVT = True if "VVT" not in nodes_to_skip else False
    doVVJ = True if "VVJ" not in nodes_to_skip else False


    if method in [3,4,6]: # jet fake estimate so don't include other MC jet fakes:
        doTTJ = False
        doZJ = False
        doVVJ = False

    if do_data:
        if args.do_ss:
            OSSS = "!os"
        else:
            OSSS = "os"
        weight = "weight"
        full_selection = BuildCutString(weight, sel, cat_data, OSSS)
        ana.nodes[nodename].AddNode(
            ana.SummedFactory("data_obs", data_samples, plot_unmodified, full_selection)
        )

    if "ZTT" not in nodes_to_skip:
        GenerateZTT(
            ana,
            nodename,
            add_name,
            samples_dict["ztt_samples"],
            plot,
            wt,
            sel,
            cat,
            gen_sels_dict["z_sels"],
            not args.do_ss,
        )
    if "ZLL" not in nodes_to_skip:
        GenerateZLL(
            ana,
            nodename,
            add_name,
            samples_dict["ztt_samples"] + samples_dict["zll_samples"],
            plot,
            wt,
            sel,
            cat,
            gen_sels_dict["z_sels"],
            not args.do_ss,
            doZL,
            doZJ,
        )
    if "TT" not in nodes_to_skip:
        GenerateTop(
            ana,
            nodename,
            add_name,
            samples_dict["top_samples"],
            plot,
            wt,
            sel,
            cat,
            gen_sels_dict["top_sels"],
            not args.do_ss,
            doTTT,
            doTTJ,
        )
    if "VV" not in nodes_to_skip:
        GenerateVV(
            ana,
            nodename,
            add_name,
            samples_dict["vv_samples"],
            plot,
            wt,
            sel,
            cat,
            gen_sels_dict["vv_sels"],
            not args.do_ss,
            doVVT,
            doVVJ,
        )
    if "EWKZ" not in nodes_to_skip:
        GenerateEWKZ( # just genuine for now
            ana,
            nodename,
            add_name,
            samples_dict["ewkz_samples"],
            plot,
            wt,
            sel,
            cat,
            gen_sels_dict["z_sels"],
            not args.do_ss,
        )
    if "W" not in nodes_to_skip:
        if method in [1, 2, 5]: # only generate W if no jetfakes
            GenerateW(
                ana,
                nodename,
                add_name,
                samples_dict,
                gen_sels_dict,
                plot,
                plot_unmodified,
                wt,
                sel,
                cat_name,
                categories,
                categories_unmodified=categories_unmodified,
                method=method,
                qcd_factor=qcd_factor,
                get_os=not args.do_ss,
            )
    if "QCD" not in nodes_to_skip and method in [1, 2, 5]:  # QCD estimate
        GenerateQCD(
            ana,
            nodename,
            add_name,
            samples_dict,
            gen_sels_dict,
            systematic,
            plot,
            plot_unmodified,
            wt,
            sel,
            cat_name,
            categories=categories,
            categories_unmodified=categories_unmodified,
            method=method,
            qcd_factor=qcd_factor,
            get_os=not args.do_ss,
        )
    elif "JetFakes" not in nodes_to_skip and method in [3,4,6]:  # Jet Fakes
        if method == 6 and "BDT_pred_score,aco" in args.var:
            print("WARNING: For CP datacards with variable names matching: `BDT_pred_score,aco*`, jet fake fractions are computed per BDT bin")
            flatten_y = True
        else:
            flatten_y = False
        GenerateFakes(
            ana,
            nodename,
            add_name,
            samples_dict,
            gen_sels_dict,
            systematic,
            plot,
            plot_unmodified,
            wt,
            sel,
            cat_name,
            categories=categories,
            categories_unmodified=categories_unmodified,
            method=method,
            qcd_factor=qcd_factor,
            get_os=not args.do_ss,
            flatten_y = flatten_y
            )

    if "signal" not in nodes_to_skip:
        # generate correct signal
        # TODO: add scheme or similar flat to determine which ones to use
        GenerateReweightedCPSignal(
            ana,
            nodename,
            add_name,
            samples_dict["signal_samples"],
            masses,
            plot,
            wt,
            sel,
            cat,
            not args.do_ss,
        )

# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Defining name of the output file
is_2d = False
is_3d = False
var_name = args.var.split("[")[0]
var_name = var_name.split("(")[0]
var_name = var_name.replace("/", "_over_")
if var_name.count(",") == 1:
    is_2d = True
    var_name = var_name.split(",")[0] + "_vs_" + var_name.split(",")[1]
if var_name.count(",") == 2:
    is_3d = True
    var_name = (
        var_name.split(",")[0]
        + "_vs_"
        + var_name.split(",")[1]
        + "_vs_"
        + var_name.split(",")[2]
    )

category_name = args.category
if args.datacard_name:
    output_name = f"{args.output_folder}/datacard_{args.datacard_name}_{category_name}_{args.channel}_{args.era}.root"
else:
    output_name = f"{args.output_folder}/datacard_{var_name}_{category_name}_{args.channel}_{args.era}.root"

if args.nodename:
    nodename = args.channel + "_" + args.category + args.nodename
else:
    nodename = args.channel + "_" + args.category

if args.do_ss:
    output_name = output_name.replace(".root", "_ss.root")

if args.bypass_plotter:
    outfile = ROOT.TFile(output_name, "UPDATE")
else:
    outfile = ROOT.TFile(output_name, "RECREATE")
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# Define qcd factor, systematics
if args.channel in ["ee", "mm"]:
    qcd_factor = 1.07
elif args.channel == "mt":
    qcd_factor = 1.12
elif args.channel == "et":
    qcd_factor = 1.13
else:
    qcd_factor = 1.0

weight = "(weight)"
if args.add_weight:
    weight += "*" + args.add_weight
# weight += "/(w_Tau_e_FakeRate*w_Tau_mu_FakeRate)"
# set systematics:
# - 1st index sets folder name contaning systematic samples
# - 2nd index sets string to be appended to output histograms
# - 3rd index specifies the weight to be applied
# - 4th lists nodes that should be skipped
# - 5th specifies if this is a FF systematic (string specifying which one)
# - 6th specifies if variable to plot needs to be changed
systematics = OrderedDict()
if args.channel == "mt":
    systematics["nominal"] = ("nominal", "", f"({weight})", [], None, None)
elif args.channel == "et":
    systematics["nominal"] = ("nominal", "", f"({weight})", [], None, None)
elif args.channel in ["ee", "mm"]:
    systematics["nominal"] = ("nominal", "", f"({weight})", [], None, None)
elif args.channel == "tt":
    systematics["nominal"] = ("nominal", "", f"({weight})", [], None, None)

if args.run_systematics and not args.do_aiso: # we dont run systematics for anti-iso region
    enabled_systematics = {
        systematic: getattr(args, "systematic_" + systematic)
        for systematic, _ in systematic_options
        if getattr(args, "systematic_" + systematic) is not None
    }

    for syst in enabled_systematics.keys():
        if syst == enabled_systematics[syst]:
            specific_systematic_name = ""
        else:
            specific_systematic_name = enabled_systematics[syst]

        systematics_dict = generate_systematics_dict(
            specific_era=args.era,
            specific_channel=args.channel,
            specific_systematic=syst,
            specific_name=specific_systematic_name,
            variable_to_plot=args.var,
        )

        for available_systematic in systematics_dict.keys():
            systematics[available_systematic] = systematics_dict[available_systematic]

    # loop over systematics and replace weight_to_replace with weight
    for syst in systematics.keys():
        systematics[syst] = (
            systematics[syst][0],
            systematics[syst][1],
            systematics[syst][2].replace("weight_to_replace", weight),
            systematics[syst][3],
            systematics[syst][4],
            systematics[syst][5],
        )
# ------------------------------------------------------------------------------------------------------------------------

categories["cat"] = (
    "(" + categories[args.category] + ")*(" + categories["baseline"] + ")"
)

# ------------------------------------------------------------------------------------------------------------------------

# Loop over systematics & run plotting, etc
systematic_suffixes = []
max_systematics_per_pass = 10

if not args.bypass_plotter:
    while len(systematics) > 0:
        analysis = Analysis.Analysis()
        analysis.nodes.AddNode(Analysis.ListNode(nodename))
        analysis.remaps = {}

        if args.channel in ["mm", "mt"]:
            analysis.remaps["Muon"] = "data_obs"
        if args.channel == "tt":
            analysis.remaps["Tau"] = "data_obs"

        previous_systematic_variation = None
        for index, systematic in enumerate(
            list(systematics.keys())[:max_systematics_per_pass]
        ):
            if (
                previous_systematic_variation is not None
                and systematics[systematic][0] != previous_systematic_variation
            ):
                continue
            previous_systematic_variation = systematics[systematic][0]
            print("Processing:", systematic)
            print("")

            sel = args.sel
            plot = args.var if systematics[systematic][5] is None else systematics[systematic][5]
            # use plot_unmodified and categories_unmodified in cases where the data and MC get different selections due to a systematic variation
            plot_unmodified = plot
            categories_unmodified = copy.deepcopy(categories)
            systematic_folder_name = systematics[systematic][0]
            systematic_suffix = systematics[systematic][1]
            weight = systematics[systematic][2]
            nodes_to_skip = systematics[systematic][3]
            ff_syst = systematics[systematic][4]
            if not isinstance(ff_syst, str): ff_syst = None

            systematic_suffixes.append(systematic_suffix)

            for sample_name in data_samples:
                analysis.AddSamples(
                    f"{args.input_folder}/{args.era}/{args.channel}/{sample_name}/nominal/merged.root",
                    "ntuple",
                    None,
                    sample_name,
                )

            for sample_name in ztt_samples + zll_samples + top_samples + vv_samples + wjets_samples + ewkz_samples:
                analysis.AddSamples(
                    f"{args.input_folder}/{args.era}/{args.channel}/{sample_name}/{systematic_folder_name}/merged.root",
                    "ntuple",
                    None,
                    sample_name,
                )

            for key, value in signal_samples.items():
                if not isinstance(value, (list,)):
                    value = [value]
                for samp in value:
                    for mass in masses:
                        sample_name = samp.replace("*", mass)
                        analysis.AddSamples(
                            f"{args.input_folder}/{args.era}/{args.channel}/{sample_name}/{systematic_folder_name}/merged.root",
                            "ntuple",
                            None,
                            sample_name,
                        )

            analysis.AddInfo(args.parameter_file, scaleTo="data_obs")

            if systematic == "nominal":
                do_data = True
            else:
                do_data = False
            RunPlotting(
                analysis,
                nodename,
                samples_dict,
                gen_sels_dict,
                ff_syst if ff_syst else systematic,
                args.category,
                categories,
                categories_unmodified,
                sel,
                systematic_suffix,
                weight,
                do_data,
                qcd_factor,
                method,
                nodes_to_skip,
            )

            del systematics[systematic]

        analysis.Run()
        analysis.nodes.Output(outfile)

        FixBins(analysis, nodename, outfile)
        for suffix in systematic_suffixes:
            GetTotals(analysis, nodename, suffix, samples_dict, outfile)
        PrintSummary(
            analysis,
            nodename,
            ["data_obs"],
            add_names=systematic_suffixes,
            channel=args.channel,
            samples_dict=samples_dict,
        )
# ------------------------------------------------------------------------------------------------------------------------

# unroll 2D histograms into 1D histograms but store both versions

if is_2d and args.do_unrolling:
    x_lines = []
    y_labels = []
    first_hist = True
    # loop over all TH2Ds and for each one unroll to produce TH1D and add to datacard
    directory = outfile.Get(nodename)
    outfile.cd(nodename)
    hists_to_add = []
    for key in directory.GetListOfKeys():
        hist_name = key.GetName()
        hist = directory.Get(hist_name).Clone()

        if not isinstance(hist, ROOT.TDirectory):
            include_of = False

            h1d = UnrollHist2D(hist, include_of)
            hists_to_add.append(h1d)
            if first_hist:
                first_hist = False
                Nxbins = hist.GetNbinsX()
                for i in range(1, hist.GetNbinsY() + 1):
                    x_lines.append(Nxbins * i)
                for j in range(1, hist.GetNbinsY() + 1):
                    y_labels.append(
                        [
                            hist.GetYaxis().GetBinLowEdge(j),
                            hist.GetYaxis().GetBinLowEdge(j + 1),
                        ]
                    )
                if include_of:
                    y_labels.append(
                        [hist.GetYaxis().GetBinLowEdge(hist.GetNbinsY() + 1), -1]
                    )
    for hist in hists_to_add:
        hist_2d = directory.Get(hist.GetName()).Clone()
        hist_2d.SetName(hist.GetName() + "_2D")
        hist_2d.Write("", ROOT.TObject.kOverwrite)
        hist.Write("", ROOT.TObject.kOverwrite)

# --------------------
# Full Uncertainty Band

directory = outfile.Get(nodename)
keys = [key.GetName() for key in directory.GetListOfKeys()]

h0 = directory.Get('total_bkg')
hists=[]

# first process systematics affecting normalisation
normalisation_systematics = {}
normalisation_systematics['lumi'] = ((0.014), ["ZTT", "ZL", "TTT", "VVT"])
normalisation_systematics['dy_xs'] = ((0.016, 0.013), ["ZTT", "ZL"])
normalisation_systematics['top_xs'] = ((0.05), ["TTT"])

for norm_syst, (value, processes) in normalisation_systematics.items():
    if isinstance(value, tuple) and len(value) == 2:
        down_shift = value[0]
        up_shift = value[1]
    else:
        down_shift = value
        up_shift = value

    h1 = h0.Clone()
    h2 = h0.Clone()
    h1.SetName(h0.GetName() + "_" + norm_syst + "Up")
    h2.SetName(h0.GetName() + "_" + norm_syst + "Down")
    for proc in processes:
        if proc not in keys:
            print(f"Process {proc} not found in the directory.")
            if proc in ["TTT", "VVT"]:
                proc = proc.replace("TTT", "TTL").replace("VVT", "VVL")
                if proc not in keys:
                    print(f"Process {proc} not found in the directory.")
                    continue
            else:
                continue

        hup = directory.Get(proc).Clone()
        hdown = directory.Get(proc).Clone()

        h1.Add(hup, -1.0)
        h2.Add(hdown, -1.0)

        hup.Scale(1+up_shift)
        hdown.Scale(1-down_shift)

        h1.Add(hup)
        h2.Add(hdown)

    hists.append(h1.Clone())
    hists.append(h2.Clone())

# now process systematics affecting shape
for hist in directory.GetListOfKeys():
    if ".subnodes" in hist.GetName():
        continue

    processes = ["ZTT", "ZL", "ZJ", "TTT", "TTJ","VVT","VVJ", "W", "QCD", "JetFakes", "JetFakesSublead", "ggH_sm_prod_sm_htt125","qqH_sm_htt125","WH_sm_htt125","ZH_sm_htt125"]
    if hist.GetName().endswith("Up") or hist.GetName().endswith("Down"):
        for proc in processes:
            if hist.GetName().startswith(proc + '_'):
                if "signal" in hist.GetName() and proc in ["WH_sm_htt125", "ZH_sm_htt125"]:
                    print("->Skipping signal systematic for VH (not valid):", hist.GetName())
                    continue
                else:
                    print(f"Adding {hist.GetName()} to total uncertainty")
                    no_syst_name = proc
                    temp_hist = h0.Clone()
                    temp_hist.Add(directory.Get(no_syst_name),-1)
                    temp_hist.Add(directory.Get(hist.GetName()))
                    temp_hist.SetName(hist.GetName())
                    hists.append(temp_hist)

(uncert, up, down) = Total_Uncertainty(h0, hists)
outfile.cd(nodename)
uncert.Write()
up.Write()
down.Write()

# --------------------

outfile.Close()
plot_file = ROOT.TFile(output_name, "READ")

if args.auto_rebin:
    outfile_rebin = ROOT.TFile(
        output_name.replace(".root", "_rebinned.root"), "RECREATE"
    )
    outfile_rebin.mkdir(nodename)
    outfile_rebin.cd(nodename)
    total_bkghist = plot_file.Get(nodename + "/total_bkg").Clone()
    binning = FindRebinning(total_bkghist, BinThreshold=100, BinUncertFraction=0.5)

    print("New binning:", binning)
    hists_done = []
    for i in plot_file.Get(nodename).GetListOfKeys():
        if i.GetName() not in hists_done:
            if ".subnodes" not in i.GetName():
                RebinHist(
                    plot_file.Get(nodename + "/" + i.GetName()).Clone(), binning
                ).Write()
                hists_done.append(i.GetName())

    outfile_rebin.Close()
    plot_file = ROOT.TFile(output_name.replace(".root", "_rebinned.root"))

titles = Plotting.SetAxisTitles(args.var, args.channel)
x_title = titles[0]
y_title = titles[1]


if args.rename_procs:
    outfile = ROOT.TFile(output_name, "UPDATE")
    if args.channel in ["mm", "ee"]:
        RenameDatacards(outfile, nodename)
    outfile.Close()

# new plotting available for 1D histograms (NB only 2D unrolled histograms are supported)
if (not is_2d) or (is_2d and args.do_unrolling):
    Histo_Plotter = HTT_Histogram(
        output_name,
        nodename,
        args.channel,
        args.era,
        args.var,
        method,
        blind=args.blind,
        log_y=False,
        is2Dunrolled=is_2d,
    )

    Histo_Plotter.plot_1D_histo()

