RATIO_RANGE = (0.6, 1.4)

VAR_LABELS = {
    "m_vis": r"$m_{\mathrm{vis}}$ [GeV]",
    "pt_1": r"$p_{\mathrm{T}}^{\tau_1}$ [GeV]",
    "pt_2": r"$p_{\mathrm{T}}^{\tau_2}$ [GeV]",
    "eta_1": r"$\eta^{\tau_1}$",
    "eta_2": r"$\eta^{\tau_2}$",
    "phi_1": r"$\phi^{\tau_1}$",
    "phi_2": r"$\phi^{\tau_2}$",
    "dR": r"$\Delta R(\tau_1, \tau_2)$",
    "n_jets": r"$N_\mathrm{jets}$",
    "met_pt": r"$p_{\mathrm{T}}^{\mathrm{miss}}$ [GeV]",
    "pt_tt": r"$p_{\mathrm{T}}^{\tau\tau}$ [GeV]",
    "l1pt_1": r"Leading $\tau_h$ L1 $p_{\mathrm{T}}$ [GeV]",
    "l1pt_2": r"Sub-leading $\tau_h$ L1 $p_{\mathrm{T}}$ [GeV]",
    "decayModePNet_1": r"$\tau_1$ decay mode (PNet)",
    "decayModePNet_2": r"$\tau_2$ decay mode (PNet)",
}

CHANNEL_LABELS = {
    "tt": r"$\tau_h\tau_h$",
    "mt": r"$\mu\tau_h$",
    "et": r"$e\tau_h$",
}


def build_lumi_dict(era_definitions, lumi_by_run):
    """Return a dictionary of integrated luminosities for each era."""
    lumi_dict = {}
    for era, runs in era_definitions.items():
        lumi_dict[era] = sum(lumi_by_run[run] for run in runs)
    return lumi_dict

