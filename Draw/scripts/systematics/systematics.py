# script to generate systematics dictionary
from collections import OrderedDict


def generate_systematics_dict(specific_era='Run3_2022', specific_channel='mt', specific_systematic='MuonID', specific_name='', variable_to_plot='m_vis'):
    systematics = OrderedDict()

    # we add eras and channels to names when the name contains *year or *channel
    # note could do the same thing for other binning e.g dm
    specific_name = specific_name.replace("*year", specific_era.split("Run3_")[1]) # keep the year part only (CMS/combine convention)
    specific_name = specific_name.replace("*channel", specific_channel)

    # Muon ID/Isolation systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Muon_ID_iso':

        nodes_to_skip = [
             'JetFakes'
        ]
        for kind in ['ID', 'Isolation']:
            up_var = f'w_Muon_{kind}_Up'
            down_var = f'w_Muon_{kind}_Down'

            for updown in [up_var,down_var]:
                formula_leading = (
                    f"(({updown}) * (genPartFlav_1 == 1 || genPartFlav_1 == 15) + "
                    f"(!(genPartFlav_1 == 1 || genPartFlav_1 == 15)))"
                )

                formula_subleading = (
                    f"(({updown}) * (genPartFlav_2 == 1 || genPartFlav_2 == 15) + "
                    f"(!(genPartFlav_2 == 1 || genPartFlav_2 == 15)))"
                )

                if kind == "ID":
                    extension = 'id'
                elif kind == "Isolation":
                    extension = 'iso'

                systematic_name = 'Muon_' + kind + '_' + updown.split('_')[-1]
                if specific_name == '':
                    histogram_name = 'syst_muon_' + extension + updown.split('_')[-1]
                else:
                    histogram_name = '_' + specific_name.replace("*kind", extension) + updown.split('_')[-1]

                if specific_channel == 'mm':
                    systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading}) * ({formula_subleading})", nodes_to_skip, None, None)
                elif specific_channel == 'mt':
                    systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading})", nodes_to_skip, None, None)

            del up_var, down_var
    # ----------------------------------------------------------------------------------------------------

    # Electron ID systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Electron_ID_Reco':

        nodes_to_skip = [
             'JetFakes'
        ]

        for kind in ['ID', 'Reco']:
            up_var = f'w_Electron_{kind}_Up'
            down_var = f'w_Electron_{kind}_Down'

            for updown in [up_var,down_var]:
                formula_leading = (
                    f"(({updown}) * (genPartFlav_1 == 1 || genPartFlav_1 == 15) + "
                    f"(!(genPartFlav_1 == 1 || genPartFlav_1 == 15)))"
                )

                formula_subleading = (
                    f"(({updown}) * (genPartFlav_2 == 1 || genPartFlav_2 == 15) + "
                    f"(!(genPartFlav_2 == 1 || genPartFlav_2 == 15)))"
                )

                systematic_name = 'Electron_ID_' + kind + updown.split('_')[-1]
                if specific_name == '':
                    histogram_name = 'syst_electron_id' + updown.split('_')[-1]
                else:
                    histogram_name = '_' + specific_name.replace('*kind', kind.lower()) + updown.split('_')[-1]

                if specific_channel == 'ee':
                    systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading}) * ({formula_subleading})", nodes_to_skip, None, None)
                elif specific_channel == 'et':
                    systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({formula_leading})", nodes_to_skip, None, None)

            del up_var, down_var
    # ----------------------------------------------------------------------------------------------------

    # Tau ID systematics
    # ----------------------------------------------------------------------------------------------------

    # First Kind: stat1, stat2
    # should be uncorrelated across DMs and eras
    if specific_systematic == 'Tau_ID_PNet':

        nodes_to_skip = ['JetFakes', 'QCD']
        kinds = ['stat1','stat2']
        decay_modes = ["0", "1", "2", "10"]
        era = specific_era # no longer run all eras (can deal with this in hadding)

        for kind in kinds:
            for dm in decay_modes:
                up_weights = []
                down_weights = []
                for obj_index, obj_type in enumerate(specific_channel):
                    if obj_type == 't':
                        up_var = f'w_Tau_ID_PNet_{obj_index+1}_{kind}_Up'
                        down_var = f'w_Tau_ID_PNet_{obj_index+1}_{kind}_Down'

                        formula = (
                            f"((variation_to_replace) * (decayModePNet_{obj_index+1} == {dm}) + "
                            f"(!(decayModePNet_{obj_index+1} == {dm})))"
                        )

                        up_weights.append(formula.replace('variation_to_replace', up_var))
                        down_weights.append(formula.replace('variation_to_replace', down_var))
                        del up_var, down_var

                systematic_name = f'Tau_ID_PNet_{kind.replace("_era_", "_")}_DM{dm}_{era}'
                if specific_name == '':
                    histogram_name =  f'syst_tau_id_pnet_{kind.replace("_era_", "_")}_DM{dm}_{era}'
                else:
                    histogram_name = specific_name.replace("*group", f"{kind}_DM{dm}PNet_{specific_era.split('Run3_')[1]}")

                if specific_channel in ["et","mt","tt"]:
                    systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), nodes_to_skip, None, None)
                    systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), nodes_to_skip, None, None)

        del up_weights, down_weights

        # Second Kind: syst_era
        # should be correlated across DMs but uncorrelated across eras

        up_weights = []
        down_weights = []
        for obj_index, obj_type in enumerate(specific_channel):
            if obj_type == 't':
                up_var = f'w_Tau_ID_PNet_{obj_index+1}_syst_era_Up'
                down_var = f'w_Tau_ID_PNet_{obj_index+1}_syst_era_Down'

                up_weights.append(f"({up_var})")
                down_weights.append(f"({down_var})")

                del up_var, down_var

        systematic_name = f'Tau_ID_PNet_syst_era_{era}'
        if specific_name == '':
            histogram_name =  f'syst_tau_id_pnet_{era}'
        else:
            histogram_name = specific_name.replace("*group", f"syst_{specific_era.split('Run3_')[1]}")

        if specific_channel in ["et","mt","tt"]:
            systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), nodes_to_skip, None, None)
            systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), nodes_to_skip, None, None)

        del up_weights, down_weights

        # Third Kind: syst_all_eras
        # should be correlated across DMs and eras

        up_weights = []
        down_weights = []
        for obj_index, obj_type in enumerate(specific_channel):
            if obj_type == 't':
                up_var = f'w_Tau_ID_PNet_{obj_index+1}_syst_all_eras_Up'
                down_var = f'w_Tau_ID_PNet_{obj_index+1}_syst_all_eras_Down'

                up_weights.append(f"({up_var})")
                down_weights.append(f"({down_var})")

                del up_var, down_var

        systematic_name = 'Tau_ID_PNet_syst_all_eras'
        if specific_name == '':
            histogram_name = 'syst_tau_id_pnet_all_eras'
        else:
            histogram_name = specific_name.replace("*group", f"syst_alleras")

        if specific_channel in ["et","mt","tt"]:
            systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), nodes_to_skip, None, None)
            systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), nodes_to_skip, None, None)

        del up_weights, down_weights
    # ----------------------------------------------------------------------------------------------------

    # Tau Fakerate (e) systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Tau_FakeRate_e':
        nodes_to_skip =  ['ZTT','ZJ','VVT','VVJ','TTT','TTJ','QCD','JetFakes','signal','W']  # only ZLL
        eta_bins = ["0.0", "1.5", "2.5"]
        era = specific_era

        for i, eta in enumerate(eta_bins[:-1]):
            up_weights = []
            down_weights = []
            for obj_index, obj_type in enumerate(specific_channel):
                if obj_type == 't':
                    up_var = f"w_Tau_e_FakeRate_{obj_index+1}_Up"
                    down_var = f"w_Tau_e_FakeRate_{obj_index+1}_Down"

                    formula = (
                        f"((variation_to_replace * ((fabs(eta_{obj_index+1}) >= {eta_bins[i]}) && (fabs(eta_{obj_index+1}) < {eta_bins[i+1]}))) + "
                        f"(!((fabs(eta_{obj_index+1}) >= {eta_bins[i]}) && (fabs(eta_{obj_index+1}) < {eta_bins[i+1]}))))"
                    )

                    up_weights.append(formula.replace('variation_to_replace', up_var))
                    down_weights.append(formula.replace('variation_to_replace', down_var))

                    del up_var, down_var

            systematic_name = f'Tau_e_FakeRate_{era}_eta_{eta}'
            if specific_name == '':
                histogram_name = f'syst_etau_fakerate_{era}_eta_{eta}'
            else:
                histogram_name = '_' + specific_name.replace("*year", era).replace('*eta', eta.replace('.','p'))

            if specific_channel in ["et","mt","tt"]:
                systematics[systematic_name + '_up'] = ('nominal', histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), nodes_to_skip, None, None)
                systematics[systematic_name + '_down'] = ('nominal', histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), nodes_to_skip, None, None)

        del up_weights, down_weights
    # ----------------------------------------------------------------------------------------------------

    # Tau Fakerate (mu) systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Tau_FakeRate_mu':
        nodes_to_skip =  ['ZTT','ZJ','VVT','VVJ','TTT','TTJ','QCD','JetFakes','signal','W']  #only ZLL
        eta_bins = ["0.0", "0.4", "0.8", "1.2", "1.7", "2.4"]
        era = specific_era

        for i, eta in enumerate(eta_bins[:-1]):
            up_weights = []
            down_weights = []
            for obj_index, obj_type in enumerate(specific_channel):
                if obj_type == 't':
                    up_var = f"w_Tau_mu_FakeRate_{obj_index+1}_Up"
                    down_var = f"w_Tau_mu_FakeRate_{obj_index+1}_Down"

                    formula = (
                        f"((variation_to_replace * ((fabs(eta_{obj_index+1}) >= {eta_bins[i]}) && (fabs(eta_{obj_index+1}) < {eta_bins[i+1]}))) + "
                        f"(!((fabs(eta_{obj_index+1}) >= {eta_bins[i]}) && (fabs(eta_{obj_index+1}) < {eta_bins[i+1]}))))"
                    )

                    up_weights.append(formula.replace('variation_to_replace', up_var))
                    down_weights.append(formula.replace('variation_to_replace', down_var))

                    del up_var, down_var

            systematic_name = f'Tau_mu_FakeRate_{era}_eta_{eta}'
            if specific_name == '':
                histogram_name = f'syst_mutau_fakerate_{era}_eta_{eta}'
            else:
                histogram_name = specific_name.replace("*year", era).replace('*eta', eta.replace('.','p'))

            if specific_channel in ["et","mt","tt"]:
                systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), nodes_to_skip, None, None)
                systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), nodes_to_skip, None, None)

        del up_weights, down_weights
    # ----------------------------------------------------------------------------------------------------

    # Tau Energy Scale systematics (This are recommended to be uncorrelated across eras but we will leave it for now)
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic in ['Tau_EnergyScale_PNet_TSCALE', 'Tau_EnergyScale_PNet_ESCALE', 'Tau_EnergyScale_PNet_MUSCALE']:
        kinds = {
            'DM0PNet': '1PRONG',
            'DM1PNet': '1PRONG_1PI0',
            'DM2PNet': '1PRONG_2PI0',
            'DM10PNet': '3PRONG'
        }

        # Genuine Taus
        if specific_systematic == 'Tau_EnergyScale_PNet_TSCALE':
            prefixes = ['Tau_EnergyScale_PNet_TSCALE_']
            nodes_to_skip = ["JetFakes", "QCD", "ZLL"]
        # Genuine electrons misidentified as taus
        elif specific_systematic == 'Tau_EnergyScale_PNet_ESCALE':
            prefixes = ['Tau_EnergyScale_PNet_ESCALE_']
            nodes_to_skip = ['ZTT','ZJ','VVT','VVJ','TTT','TTJ','QCD','JetFakes','signal','W']
        # Genuine muons misidentified as taus
        elif specific_systematic == 'Tau_EnergyScale_PNet_MUSCALE':
            prefixes = ['Tau_EnergyScale_PNet_MUSCALE_']
            nodes_to_skip = ['ZTT','ZJ','VVT','VVJ','TTT','TTJ','QCD','JetFakes','signal','W']

        for name, folder_suffix in kinds.items():
            for prefix in prefixes:
                for updown in ['up', 'down']:
                    scale_type = prefix.removesuffix('_').split('_')[-1]
                    systematic_name = 'syst_tau_escale_' + scale_type + '_' + name + '_' + updown
                    folder_name = prefix + folder_suffix + '_' + updown
                    if specific_name == '':
                        histogram_name = 'syst_tau_escale_' + scale_type + '_' + name + updown.capitalize()
                    else:
                        histogram_name = '_' + specific_name.replace("*group", name) + updown.capitalize()

                    if specific_channel in ["et","mt","tt"]:
                        systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', nodes_to_skip, None, None)

    # ----------------------------------------------------------------------------------------------------

    # Trigger systematics
    # ----------------------------------------------------------------------------------------------------

    # should be uncorrelated across DMs and eras
    if specific_systematic == 'Trigger':

        nodes_to_skip = ['JetFakes', 'QCD']
        decay_modes = ["0", "1", "2", "10"]
        era = specific_era

        if specific_channel in ['tt']:

            # DOUBLE TAU - TAU VARIATION
            for dm in decay_modes:
                up_weights = []
                down_weights = []
                for tau_number in [1, 2]:  # two taus in tt channel
                    up_var = 'w_Trigger_doubletau_tauUp'
                    down_var = 'w_Trigger_doubletau_tauDown'
                    formula = (
                        f"((variation_to_replace) * (decayModePNet_{tau_number} == {dm}) + "
                        f"(!(decayModePNet_{tau_number} == {dm})))"
                    )
                    up_weights.append(formula.replace('variation_to_replace', up_var))
                    down_weights.append(formula.replace('variation_to_replace', down_var))
                    del up_var, down_var

                systematic_name = f'syst_Tau_Trigger_doubletau_tau_DM{dm}_{era}'
                if specific_name == '':
                    histogram_name =  f'syst_Tau_Trigger_doubletau_tau_DM{dm}_{era}'
                else:
                    histogram_name = specific_name.replace("*obj", 't').replace("*trigger", "ditau").replace("*group", f'VTight_DM{dm}PNet_{specific_era.split("Run3_")[1]}')

                systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), nodes_to_skip, None, None)
                systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), nodes_to_skip, None, None)


            # DOUBLE TAU JET - TAU VARIATION
            for dm in decay_modes:
                up_weights = []
                down_weights = []
                for tau_number in [1, 2]:  # two taus in tt channel
                    up_var = 'w_Trigger_doubletaujet_tauUp'
                    down_var = 'w_Trigger_doubletaujet_tauDown'
                    formula = (
                        f"((variation_to_replace) * (decayModePNet_{tau_number} == {dm}) + "
                        f"(!(decayModePNet_{tau_number} == {dm})))"
                    )
                    up_weights.append(formula.replace('variation_to_replace', up_var))
                    down_weights.append(formula.replace('variation_to_replace', down_var))
                    del up_var, down_var

                systematic_name = f'syst_Tau_Trigger_doubletaujet_tau_DM{dm}_{era}'
                if specific_name == '':
                    histogram_name =  f'syst_Tau_Trigger_doubletaujet_tau_DM{dm}_{era}'
                else:
                    histogram_name = specific_name.replace("*obj", 't').replace("*trigger", "ditaujet").replace("*group", f'VTight_DM{dm}PNet_{specific_era.split("Run3_")[1]}')

                systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), nodes_to_skip, None, None)
                systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), nodes_to_skip, None, None)

            # DOUBLE TAU JET - JET VARIATION
            up_weight = '(w_Trigger_doubletaujet_jetUp)'
            down_weight = '(w_Trigger_doubletaujet_jetDown)'

            systematic_name = f'syst_Tau_Trigger_doubletaujet_jet_DM{dm}_{era}'
            if specific_name == '':
                histogram_name =  f'syst_Tau_Trigger_doubletaujet_jet_DM{dm}_{era}'
            else:
                histogram_name = specific_name.replace("*obj", 'j').replace("*trigger", "ditaujet").replace("*group", specific_era.split("Run3_")[1])

            systematics[systematic_name + '_up'] = ('nominal', '_' + histogram_name + 'Up', 'weight_to_replace*' + '*'.join(up_weights), nodes_to_skip, None, None)
            systematics[systematic_name + '_down'] = ('nominal', '_' + histogram_name + 'Down', 'weight_to_replace*' + '*'.join(down_weights), nodes_to_skip, None, None)

            del up_weights, down_weights

        else:  # Trigger systematics are simpler in the et and mt cases

            up_var = '(w_TriggerUp)'
            down_var = '(w_TriggerDown)'

            for updown in ["up", "down"]:

                systematic_name = 'syst_Trigger_' + specific_channel + updown
                if specific_name == '':
                    histogram_name = 'syst_trigger_' + extension + updown.capitalize()
                else:
                    histogram_name = '_' + specific_name + updown.capitalize()

                weight_updown = up_var if updown == "up" else down_var
                systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({weight_updown})", nodes_to_skip, None, None)

            del up_var, down_var

    if specific_systematic == "IP_Significance":
        nodes_to_skip = ['JetFakes', 'QCD']

        for type in ['prompt', 'tauDecay']:
            for eta in ['Lt1p0', '1p0to1p6', 'Gt1p6']:
                for updown in ['Up', 'Down']:
                    # get alternative weight
                    weight_updown = f'w_IPSig_correction_{type}_eta{eta}_stat{updown}'

                    systematic_name = f'syst_ip_significance_{type}_{eta}_{updown}'
                    if specific_name == '':
                        histogram_name = f'_syst_ip_significance_{type}_{eta}{updown.capitalize()}'
                    else:
                        histogram_name = f'_{specific_name.replace("*obj", type).replace("*eta", eta)}{updown.capitalize()}'
                    systematics[systematic_name] = ('nominal', histogram_name, f'weight_to_replace * ({weight_updown})', nodes_to_skip, None, None)

    # ----------------------------------------------------------------------------------------------------

    # Jet Energy Scale systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Jet_EnergyScale_Total':

        nodes_to_skip = ['JetFakes', 'QCD']

        prefix = 'jec_syst_Total'
        name = 'syst_jet_scale_Total'
        for updown in ['up', 'down']:
            systematic_name = name + '_' + updown
            folder_name = prefix + '_' + updown
            if specific_name == '':
                histogram_name = '_' + prefix + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()
            systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', nodes_to_skip, None, None)

    # ----------------------------------------------------------------------------------------------------

    # Jet Energy Resolution systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Jet_EnergyResolution':

        nodes_to_skip = ['JetFakes', 'QCD']

        prefix = 'jer_syst'
        name = 'syst_jet_resolution'
        for updown in ['up', 'down']:
            systematic_name = name + '_' + updown
            folder_name = prefix + '_' + updown
            if specific_name == '':
                histogram_name = '_' + prefix + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()
            systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', nodes_to_skip, None, None)

    # ----------------------------------------------------------------------------------------------------

    # Electron Energy Scale + Smearing systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'Electron_Scale_Smearing':

        nodes_to_skip = [
             'JetFakes'
        ]

        for kind in ['Scale', 'Smearing']:
            for updown in ['up', 'down']:
                systematic_name = 'syst_electron_' + kind.lower() + '_' + updown
                folder_name = 'Electron_' + kind + '_' + updown

                if kind == "Scale":
                    extension = 'scale'
                elif kind == "Smearing":
                    extension = 'res'

                if specific_name == '':
                    histogram_name = '_' +'syst_electron_' + kind.lower() + updown.capitalize()
                else:
                    histogram_name = '_' + specific_name.replace('*kind', extension) + updown.capitalize()

                if specific_channel in ["ee","et"]:
                    systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', nodes_to_skip, None, None)

    # ----------------------------------------------------------------------------------------------------

    # TTbar pT reweighting systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'TTbar_Shape':
        nodes_to_skip = [
            "ZTT", "ZLL", "ZL", "ZJ",
            "VV", "VVT", "VVJ",
            "W","signal", 'QCD', 'JetFakes'
        ]
        up_var = 'w_Top_pt_Reweighting'
        down_var = '(1/w_Top_pt_Reweighting)'

        for updown in ["up", "down"]:
            systematic_name = 'syst_ttbar_shape_' + updown
            if specific_name == '':
                histogram_name = '_syst_ttbar_shape' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            weight_updown = up_var if updown == "up" else down_var
            systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({weight_updown})", nodes_to_skip, None, None)

        del up_var, down_var
    # ----------------------------------------------------------------------------------------------------

    # DY pT reweighting systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == 'DY_Shape' or specific_systematic == 'DY_Shape_Imperial':
        nodes_to_skip = [
            "TT", "TTT", "TTJ",
            "VV", "VVT", "VVJ",
            "W","signal", 'QCD', 'JetFakes'
        ]

        if specific_systematic == 'DY_Shape':
            up_var = 'w_Zpt_Reweighting'
            down_var = '(1/w_Zpt_Reweighting)'
        # else:
        #     up_var = 'w_Zpt_Reweighting_Imperial'
        #     down_var = '(1/w_Zpt_Reweighting_Imperial)'

        for updown in ["up", "down"]:
            systematic_name = 'syst_dy_shape_' + updown
            if specific_name == '':
                histogram_name = '_syst_dy_shape' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            weight_updown = up_var if updown == "up" else down_var
            systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({weight_updown})", nodes_to_skip, None, None)

        del up_var, down_var

    if specific_systematic == 'Fake_Flat_Uncertainty':
        nodes_to_skip = [
            "TT", "TTT", "TTJ",
            "ZTT", "ZLL", "ZL", "ZJ",
            "VV", "VVT", "VVJ",
            "W","signal"
        ]

        for updown in ["up", "down"]:
            systematic_name = 'flat_fake_sub_' + updown
            if specific_name == '':
                histogram_name = '_flat_fake_sub' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            systematics[systematic_name] = ('nominal', histogram_name, "weight_to_replace", nodes_to_skip, None, None)

    if specific_systematic == 'Fake_Factors':
        # note the way this is setup at the moment won't work of the box if doing it for mt and et
        nodes_to_skip = [
            "TT", "TTT", "TTJ",
            "ZTT", "ZLL", "ZL", "ZJ",
            "VV", "VVT", "VVJ",
            "W","signal"
        ]


        if specific_channel == 'tt':

            if specific_name == '':
                name_prefix = 'ff'
            else:
                name_prefix = specific_name

            for updown in ["up", "down"]:

                # we split the stat error by the dm of the leading tau as this is how the FFs are measured
                for dm in [0,1,2,10]:
                    for njets in [0,1,2]:
                        systematic_name = f'ff_stat_dm{dm}_{njets}j_{updown}'
                        histogram_name = f'_{name_prefix}_stat_dm{dm}_{njets}j{updown.capitalize()}'

                        if njets in [0,1]:
                            jet_selection = f'n_jets=={njets}'
                            anti_jet_selection = f'n_jets!={njets}'
                        elif njets == 2: # greater or equal
                            jet_selection = f'n_jets>=2'
                            anti_jet_selection = f'n_jets<2'

                        systematics[systematic_name] = ('nominal', histogram_name, "weight_to_replace", nodes_to_skip, f"FF_syst:(*ff_nom*((decayModePNet_1 != {dm}) || ({anti_jet_selection})) + *ff_nom_{updown}*(decayModePNet_1 == {dm} && {jet_selection}))", None)

        elif specific_channel in ['mt','et']:
            # three fake factors to vary in this case:
            if specific_name == '':
                name_prefix = 'ff_' + specific_channel
            else:
                name_prefix = specific_name

            # Statistical uncertainties, uncorrelated across DM and Njets
            for updown in ["up", "down"]:
                for ff_type in ['wj','qcd','mc_top']:
                    nominal_weight = f"FF_{ff_type}_noipcut_nom" # nominal weight in jet fake calculation
                    alternative_weight = f"FF_{ff_type}_noipcut_{updown}" # weight to replace in jet fake calculation
                    for dm in [0,1,2,10]:
                        for njets in [0,1,2]:
                            if njets in [0,1]:
                                jet_selection = f'n_jets=={njets}'
                                anti_jet_selection = f'n_jets!={njets}'
                            elif njets == 2: # greater or equal
                                jet_selection = f'n_jets>=2'
                                anti_jet_selection = f'n_jets<2'

                            systematic_name = f'ff_{ff_type}_stat_dm{dm}_{njets}j_{updown}'
                            histogram_name = f'_{name_prefix}_{ff_type}_stat_dm{dm}_{njets}j{updown.capitalize()}'

                            ff_string = f"FF_uct_{ff_type}_stat:({nominal_weight}*((decayModePNet_2 != {dm}) || ({anti_jet_selection})) + {alternative_weight}*(decayModePNet_2 == {dm} && {jet_selection}))"
                            systematics[systematic_name] = ('nominal', histogram_name, "weight_to_replace", nodes_to_skip, ff_string, None)

            # Systematic uncertainties
            for updown in ["up", "down"]:
                for ff_type in ['wj','qcd','mc_top']:
                    if updown == 'up': # placeholder systematic uncertainties
                        alternative_weight = f"1.1*FF_{ff_type}_noipcut_nom" # weight to replace in jet fake calculation
                    elif updown == 'down': # placeholder systematic uncertainties
                        alternative_weight = f"0.9*FF_{ff_type}_noipcut_nom" # weight to replace in jet fake calculation
                    systematic_name = f'ff_{ff_type}_syst_{updown}'
                    histogram_name = f'_{name_prefix}_{ff_type}_syst{updown.capitalize()}'
                    ff_string = f"FF_uct_{ff_type}_syst:{alternative_weight}"
                    systematics[systematic_name] = ('nominal', histogram_name, "weight_to_replace", nodes_to_skip, ff_string, None)

        # subtraction uncertainty
        systematic_name = 'ff_sub_syst_' + updown
        histogram_name = f'_{name_prefix}_sub_syst' + updown.capitalize()
        systematics[systematic_name] = ('nominal', histogram_name, "weight_to_replace", nodes_to_skip, None, None)

    # ----------------------------------------------------------------------------------------------------

    # QCD Background systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == "QCD_Background":
        nodes_to_skip = [
            "ZTT", "ZLL", "ZL", "ZJ",
            "TT", "TTT", "TTJ",
            "VV", "VVT", "VVJ",
            "W","signal", 'JetFakes'
        ]
        for updown in ['up','down']:
            systematic_name = 'qcd_sub_' + updown
            if specific_name == '':
                histogram_name = '_qcd_sub' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            systematics[systematic_name] = ('nominal', histogram_name, 'weight_to_replace', nodes_to_skip, None, None)

    if specific_systematic == "Pileup":
        nodes_to_skip = ['QCD', 'JetFakes'
        ]

        up_var = 'w_Pileup_Up'
        down_var = 'w_Pileup_Up'

        for updown in ["up", "down"]:
            systematic_name = 'syst_pileup_' + updown
            if specific_name == '':
                histogram_name = '_syst_pileup' + updown.capitalize()
            else:
                histogram_name = '_' + specific_name + updown.capitalize()

            weight_updown = up_var if updown == "up" else down_var
            systematics[systematic_name] = ('nominal', histogram_name, f"weight_to_replace * ({weight_updown})", nodes_to_skip, None, None)

        del up_var, down_var


    # ----------------------------------------------------------------------------------------------------

    # MET Recoil systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == "MET_Recoil":
        # exists for W, DY and Signal Samples
        nodes_to_skip = [
            "TT", "TTT", "TTJ",
            "VV", "VVT", "VVJ",
            "QCD", 'JetFakes'
        ]
        variations = [{'dir': 'MET_Recoil_Resolution', 'name': 'syst_met_recoil_resolution', 'type': 'res'},
                      {'dir': 'MET_Recoil_Response', 'name': 'syst_met_recoil_response', 'type': 'scale'}]

        for updown in ['up','down']:
            for syst_info in variations:
                systematic_name = syst_info['name'] + '_' + updown
                folder_name = syst_info['dir'] + '_' + updown
                if specific_name == '':
                    histogram_name = '_' + syst_info['name'] + updown.capitalize()
                else:
                    histogram_name = '_' + specific_name.replace('*type', syst_info['type']) + updown.capitalize()

                systematics[systematic_name] = (folder_name, histogram_name, 'weight_to_replace', nodes_to_skip, None, None)

    # ----------------------------------------------------------------------------------------------------

    # Signal Theory systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == "Signal_Theory":
        nodes_to_skip = [
            "ZTT", "ZLL", "ZL", "ZJ",
            "TT", "TTT", "TTJ",
            "VV", "VVT", "VVJ",
            "W", "QCD", 'JetFakes'
        ]

        variations = ["Scale_muR", "Scale_muF",
                   "PS_ISR", "PS_FSR",
        ]

        for updown in ['up','down']:
            for var in variations:
                weight = "w_Signal_Theory_" + var + updown.capitalize()
                systematic_name = f'signal_theory_{var}_{updown}'
                if specific_name == '':
                    histogram_name = f'_signal_theory_{var}{updown.capitalize()}'
                else:
                    if var == "Scale_muR":
                        histogram_name = '_'+ specific_name.replace("*type", "QCDscale_ren") + '_ACCEPT' + updown.capitalize()
                    elif var == "Scale_muF":
                        histogram_name = '_'+ specific_name.replace("*type", "QCDscale_fac") + '_ACCEPT' + updown.capitalize()
                    elif var == "PS_ISR":
                        histogram_name = '_'+ specific_name.replace("*type", "ps_isr")+ updown.capitalize()
                    elif var == "PS_FSR":
                        histogram_name = '_'+ specific_name.replace("*type", "ps_fsr")+ updown.capitalize()

                systematics[systematic_name] = ('nominal', histogram_name, f'weight_to_replace * ({weight})', nodes_to_skip, None, None)
    # ----------------------------------------------------------------------------------------------------


    # IP Calibration Systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == "IP_Calibration":
        nodes_to_skip = [
            "ZJ", "TTJ", "VVJ",
            "W","QCD", 'JetFakes'
        ]

        variables_to_consider = [
            "aco_pi_pi", "aco_pi_rho", "aco_pi_a1_FASTMTT_MassConstraint",
            "aco_rho_pi", "aco_a1_pi_FASTMTT_MassConstraint",
            "aco_mu_pi", "aco_mu_rho", "aco_mu_a1_FASTMTT_MassConstraint",
            "aco_e_pi", "aco_e_rho", "aco_e_a1_FASTMTT_MassConstraint"
        ]

        for var in variables_to_consider:
            if var in variable_to_plot:
                for updown in ['up','down']:
                    systematic_name = f'IP_Calibration_{updown}'
                    if specific_name == '':
                        histogram_name = f'_IP_Calibration{updown.capitalize()}'
                    else:
                        histogram_name = '_' + specific_name + f'{updown.capitalize()}'

                    new_variable_to_plot = variable_to_plot.replace(var, f'{var}_IP_syst_{updown}')

                    systematics[systematic_name] = ('nominal', histogram_name, 'weight_to_replace', nodes_to_skip, None, new_variable_to_plot)
    # ----------------------------------------------------------------------------------------------------

    # Secondary Vertex Systematics
    # ----------------------------------------------------------------------------------------------------
    if specific_systematic == "SV_Resolution":
        nodes_to_skip = [
            "ZJ", "TTJ", "VVJ",
            "W","QCD", 'JetFakes'
        ]

        variables_to_consider = [
            "aco_pi_a1_FASTMTT_MassConstraint", "aco_a1_pi_FASTMTT_MassConstraint",
            "aco_a1_rho_FASTMTT_MassConstraint", "aco_rho_a1_FASTMTT_MassConstraint",
            "aco_a1_a1",
            "aco_mu_a1_FASTMTT_MassConstraint", "aco_e_a1_FASTMTT_MassConstraint"
        ]

        for var in variables_to_consider:
            if var in variable_to_plot:
                for updown in ['up','down']:
                    systematic_name = f'SV_Resolution_{updown}'
                    if specific_name == '':
                        histogram_name = f'_SV_Resolution{updown.capitalize()}'
                    else:
                        histogram_name = '_' + specific_name + f'{updown.capitalize()}'

                    new_variable_to_plot = variable_to_plot.replace(var, f'{var}_SV_syst_{updown}')

                    systematics[systematic_name] = ('nominal', histogram_name, 'weight_to_replace', nodes_to_skip, None, new_variable_to_plot)
    # ----------------------------------------------------------------------------------------------------

    return systematics
