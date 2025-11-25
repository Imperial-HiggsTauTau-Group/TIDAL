import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
import uproot
import os

plt.style.use(hep.style.CMS)
plt.rcParams.update({"font.size": 16})


class HTT_Histogram:

    def __init__(self, file, category, channel, era, variable,blind=False, log_y=False, is2Dunrolled=False, save_name=None, for_combine=False):
        self.file = uproot.open(file)
        self.file_name = file
        self.directory = os.path.dirname(file)
        self.category = category
        self.channel = channel
        self.era = era
        self.blind = blind
        self.variable = variable.split("[")[0]
        self.log_y = log_y
        self.is2Dunrolled = is2Dunrolled
        self.save_name = save_name
        self.for_combine = for_combine
        self.initialize_plotting()
        self.initialize_nodes()
        self.get_labels()

        if self.is2Dunrolled:
            self.var_dim_1 = np.fromstring(variable.split('[')[1][:-2], sep=',', dtype=float)
            self.var_dim_2 = np.fromstring(variable.split('[')[2][:-1], sep=',', dtype=float)
        # get yields
        self.get_backgrounds()
        self.get_data()

    def initialize_plotting(self):
        # get bin information
        self.bin_edges = self.file[f"{self.category}/ZTT"].axis().edges()
        self.bin_centers = self.file[f"{self.category}/ZTT"].axis().centers()
        self.bin_widths = np.diff(self.bin_edges)
        self.step_edges = np.append(self.bin_edges,2*self.bin_edges[-1]-self.bin_edges[-2]) # for outline
        # track height of stacked backgrounds
        self.stacked_block = np.zeros(len(self.bin_centers))
        self.stacked_step = np.zeros(len(self.bin_centers)+2)
        # define colour scheme
        self.colors = {
                "yellow": "#ffa90e",
                "red": "#e42536",
                "pink": "#f6b5fd", #facaff",
                "lightblue": "#3f90da",
                "violet": "#882ebd",
                "purple": "#964a8b",
                "brown": "#a96b59",
                "green": "#a0c172", # #b1cf86
                "grey": "#94a4a2",
                "ash": "#717581",
            }
        # channel label height:
        self.ch_label_height = 0.915

        # define channel label (inclusive)
        label_map = {
                "tt": "$\\tau_h\\tau_h$",
                "mt": "$\\mu\\tau_h$",
                "et": "$e\\tau_h$",
                "em": "$e\\mu$",
                "ee": "$ee$",
                "mm": "$\\mu\\mu$",
            }
        self.channel_label = label_map[self.channel]

        if self.channel == 'tt':
            # define CP channel label
            if "pipi" in self.category or "pi_pi" in self.category:
                self.channel_label = r"$\pi\pi$"
            elif "pirho" in self.category or "pi_rho" in self.category:
                self.channel_label = r"$\pi\rho$"
            elif "rhopi" in self.category or "rho_pi" in self.category:
                self.channel_label = r"$\rho\pi$"
            elif "rhorho" in self.category or "rho_rho" in self.category:
                self.channel_label = r"$\rho\rho$"
            elif "pia11pr" in self.category or "pi_a11pr" in self.category:
                self.channel_label = r"$\pi a_1^{1pr}$"
            elif "pia1" in self.category or "pi_a1" in self.category:
                self.channel_label = r"$\pi a_1^{3pr}$"
            elif "a11prpi" in self.category or "a11pr_pi" in self.category:
                self.channel_label = r"$a_1^{1pr}\pi$"
            elif "a1pi" in self.category or "a1_pi" in self.category:
                self.channel_label = r"$a_1^{3pr} \pi$"
            elif "rhoa11pr" in self.category or "rho_a11pr" in self.category:
                self.channel_label = r"""$\rho a_1^{1pr}$ / $a_1^{1pr}\rho$ / $a_1^{1pr}a_1^{1pr}$"""
                # self.ch_label_height = 0.825
            elif "rhoa1" in self.category or "rho_a1" in self.category:
                self.channel_label = r"$\rho a_1^{3pr}$"
            elif "a1rho" in self.category or "a1_rho" in self.category:
                self.channel_label = r"$a_1^{3pr}\rho$"
            elif "a11pra1" in self.category or "a11pr_a1" in self.category:
                self.channel_label = r"$a_1^{1pr}a_1^{3pr}$"
            elif "a1a11pr" in self.category or "a1_a11pr" in self.category:
                self.channel_label = r"$a_1^{3pr}a_1^{1pr}$"
            elif "a1a1" in self.category or "a1_a1" in self.category:
                self.channel_label = r"$a_1^{3pr}a_1^{3pr}$"
        elif self.channel == 'mt':
            if "DM0_" in self.category:
                self.channel_label = r"$\mu\pi$"
            elif "DM1_" in self.category:
                self.channel_label = r"$\mu\rho$"
            elif "DM2_" in self.category:
                self.channel_label = r"$\mu a_1^{1pr}$"
            elif "DM10_" in self.category:
                self.channel_label = r"$\mu a_1^{3pr}$"
            # define CP channel label
            # if ("higgs_" in self.category):
            if "mupi" in self.category or "mu_pi" in self.category:
                self.channel_label = r"$\mu\pi$"
            elif "murho" in self.category or "mu_rho" in self.category:
                self.channel_label = r"$\mu\rho$"
            elif "mua11pr" in self.category or "mu_a11pr" in self.category:
                self.channel_label = r"$\mu a_1^{1pr}$"
            elif "mua1" in self.category or "mu_a1" in self.category:
                self.channel_label = r"$\mu a_1^{3pr}$"
        elif self.channel == 'et':
            if "DM0_" in self.category:
                self.channel_label = r"$e\pi$"
            elif "DM1_" in self.category:
                self.channel_label = r"$e\rho$"
            elif "DM2_" in self.category:
                self.channel_label = r"$e a_1^{1pr}$"
            elif "DM10_" in self.category:
                self.channel_label = r"$e a_1^{3pr}$"
            # define CP channel label
            if "epi" in self.category or "e_pi" in self.category:
                self.channel_label = r"$e\pi$"
            elif "erho" in self.category or "e_rho" in self.category:
                self.channel_label = r"$e\rho$"
            elif "ea11pr" in self.category or "e_a11pr" in self.category:
                self.channel_label = r"$e a_1^{1pr}$"
            elif "ea1" in self.category or "e_a1" in self.category:
                self.channel_label = r"$e a_1^{3pr}$"



    def initialize_nodes(self):
        # get total background nodes
        self.total_background_nodes = {
            "nominal": "total_bkg",
            "up": "total_bkg_full_uncerts_up",
            "down": "total_bkg_full_uncerts_down",
        }
        if self.channel == "tt":
            self.backgrounds = {
                                "Jet$\\to\\tau_h$": {"nodes": ["JetFakes", "JetFakesSublead"], "color": "green"},
                                "Z$\\to\\ell\\ell$": {"nodes": ["ZL"], "color": "lightblue"},
                                "Genuine $\\tau$": {"nodes": ["ZTT", "TTT", "VVT", "qqH_sm_htt125","ggH_sm_prod_sm_htt125","WH_sm_htt125","ZH_sm_htt125"], "color": "yellow"},
                            }
            self.lep1 = "\\tau_1"
            self.lep2 = "\\tau_2"
        elif self.channel == "mt":
            self.backgrounds = {
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVJ", "W", "VVT"], "color": "red"},
                                "$t\\bar{t}$": {"nodes": ["TTJ", "TTT"], "color": "violet"},
                                "Z$\\to\\ell\\ell$": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            if self.for_combine: # Jet Fake Cat test
                self.backgrounds = {
                                    "Jet$\\to\\tau_h$": {"nodes": ["JetFakes",], "color": "green"},
                                    "Z$\\to\\ell\\ell$": {"nodes": ["ZL"], "color": "lightblue"},
                                    "Genuine $\\tau$": {"nodes": ["ZTT", "TTT", "VVT"], "color": "yellow"},
                                }
            self.lep1 = "\\mu"
            self.lep2 = "\\tau"

        elif self.channel == "et":
            self.backgrounds = {
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVJ", "W", "VVT"], "color": "red"},
                                "$t\\bar{t}$": {"nodes": ["TTJ", "TTT"], "color": "violet"},
                                "Z$\\to\\ell\\ell$": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            if self.for_combine: # Jet Fake Cat test
                self.backgrounds = {
                                    "Jet$\\to\\tau_h$": {"nodes": ["JetFakes",], "color": "green"},
                                    "Z$\\to\\ell\\ell$": {"nodes": ["ZL"], "color": "lightblue"},
                                    "Genuine $\\tau$": {"nodes": ["ZTT", "TTT", "VVT"], "color": "yellow"},
                                }
            self.lep1 = "e"
            self.lep2 = "\\tau"
        elif self.channel == "em":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to\\ell\\ell": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "e"
            self.lep2 = "\\mu"
        elif self.channel == "ee":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to$ ee": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "e_1"
            self.lep2 = "e_2"
        elif self.channel == "mm":
            self.backgrounds = {
                                "$t\\bar{t}$": {"nodes": ["TTT", "TTJ"], "color": "violet"},
                                "QCD": {"nodes": ["QCD"], "color": "pink"},
                                "Electroweak": {"nodes": ["VVT", "VVJ", "W"], "color": "red"},
                                "Z$\\to\\mu\\mu$": {"nodes": ["ZL", "ZJ"], "color": "lightblue"},
                                "Z$\\to\\tau\\tau$": {"nodes": ["ZTT"], "color": "yellow"},
                            }
            self.lep1 = "\\mu_1"
            self.lep2 = "\\mu_2"
        # get luminosity
        if self.era == "Run3_2022":
            self.lumi = 7.98
        elif self.era == "Run3_2022EE":
            self.lumi = 26.67
        elif self.era == "Run3_2023":
            self.lumi = 18.06
        elif self.era == "Run3_2023BPix":
            self.lumi = 9.69
        elif self.era == 'full22':
            self.lumi = 34.65
        elif self.era == 'full23':
            self.lumi = 27.76
        elif self.era == 'earlyrun3':
            self.lumi = 62.41
        else: 
            self.lumi = 0.0
        # get color for each background
        for bkg, info in self.backgrounds.items():
            info['color'] = self.colors[info["color"]]  # replace color name with hex code

    def get_labels(self):
        label_map = {
            "m_vis": r"m$_{vis}$ (GeV)",
            "eta_1": rf"$\eta^{{{self.lep1}}}$",
            "eta_2": rf"$\eta^{{{self.lep2}}}$",
            "phi_1": rf"$\phi^{{{self.lep1}}}$",
            "phi_2": rf"$\phi^{{{self.lep2}}}$",
            "pt_1": rf"$p_T^{{{self.lep1}}}$ (GeV)",
            "pt_2": rf"$p_T^{{{self.lep2}}}$ (GeV)",
            "met_pt": r"$p_T^{miss}$ (GeV)",
            "met_phi": r"$\phi^{miss}$",
            "met_dphi_1": rf"$\Delta\phi({self.lep1},$" + r"$p_T^{miss})$",
            "met_dphi_2": rf"$\Delta\phi({self.lep2},$" + r"$p_T^{miss})$",
            "dR": rf"$\Delta R({self.lep1}, {self.lep2})$",
            "dphi": rf"$\Delta\phi({self.lep1}, {self.lep2})$",
            "pt_tt": r"$p_T$" + rf"$^{{{self.lep1}}}$" +rf"$^{{{self.lep2}}}$ (GeV)",
            "pt_vis": r"$p_T^{vis}$ (GeV)",
            "FastMTT_mass": r"$m_{FastMTT}$ (GeV)",
            "FastMTT_mass_constraint": r"$m_{FastMTT}^{constr.}$ (GeV)",
            "mt_1": rf"$m_T({self.lep1},$" + r"$p_T^{miss})$",
            "mt_2": rf"$m_T({self.lep2},$" + r"$p_T^{miss})$",
            "mt_lep": rf"$m_T({self.lep1}, {self.lep2})$",
            "mt_tot": r"$m_T^{tot}$ (GeV)",
            "jpt_1": r"$p_T^{j_1}$ (GeV)",
            "jpt_2": r"$p_T^{j_2}$ (GeV)",
            "jeta_1": r"$\eta^{j_1}$",
            "jeta_2": r"$\eta^{j_2}$",
            "mjj": r"$m_{jj}$ (GeV)",
            "jdeta": r"$\Delta\eta_{jj}$",
            "dijetpt": r"$p_T^{jj}$ (GeV)",
            "n_jets": r"$N_{jets}$",
            "n_bjets": r"$N_{bjets}$",

        }
        self.variable_label = label_map.get(self.variable, self.variable)
        if self.is2Dunrolled:
            if "BDT_pred_score,aco" in self.variable:
                self.variable_label = "Acoplanarity Bin Number"


    def get_counts_errors(self, node):
        # get count and error for a given node
        path = f"{self.category}/{node}"
        histo = self.file[path]
        # We want these to be differential (divide by bin width)
        counts = histo.values()/self.bin_widths
        errors = histo.errors()/self.bin_widths
        return counts, errors

    def get_backgrounds(self):
        # total background template
        bkg_counts, _ = self.get_counts_errors(self.total_background_nodes["nominal"])
        bkg_counts_up, _ = self.get_counts_errors(self.total_background_nodes["up"])
        bkg_counts_down, _ = self.get_counts_errors(self.total_background_nodes["down"])

        self.total_bkg_error_up = np.abs(bkg_counts_up - bkg_counts)
        self.total_bkg_error_down = np.abs(bkg_counts_down - bkg_counts)
        self.total_bkg_error_pce_up = np.divide(self.total_bkg_error_up, bkg_counts, where=bkg_counts > 0, out=np.zeros_like(bkg_counts))
        self.total_bkg_error_pce_down = np.divide(self.total_bkg_error_down, bkg_counts, where=bkg_counts > 0, out=np.zeros_like(bkg_counts))

        self.total_background = {"counts": bkg_counts, "error_up": self.total_bkg_error_up, "error_down": self.total_bkg_error_down, "pc_error_up": self.total_bkg_error_pce_up, "pc_error_down": self.total_bkg_error_pce_down}

        for bkg, info in self.backgrounds.items():
            bkg_counts = np.zeros(len(self.bin_centers))
            bkg_sq_errors = np.zeros(len(self.bin_centers))
            for contribution in info["nodes"]:
                counts, errors = self.get_counts_errors(contribution)
                bkg_counts += counts
                bkg_sq_errors += errors**2
            # store counts and yields for each contribution
            self.backgrounds[bkg]["counts"] = bkg_counts
            self.backgrounds[bkg]["errors"] = np.sqrt(bkg_sq_errors)

        return True

    def get_data(self):
        # get counts for the data
        data_counts, data_errors = self.get_counts_errors("data_obs")
        pce_data = np.divide(data_errors, data_counts, where=data_counts > 0, out=np.zeros_like(data_counts))
        self.data = {"counts": data_counts, "errors": data_errors, "pc_error": pce_data}
        # get data/MC ratio
        count_ratio =  np.divide(data_counts, self.total_background["counts"], where=self.total_background["counts"] > 0, out=np.zeros_like(data_counts))
        self.data_MC_ratio = {"counts": count_ratio, "error_data": self.data["pc_error"], "error_MC_up": self.total_background["pc_error_up"], "error_MC_down": self.total_background["pc_error_down"]}
        return True

    def stack_background(self, ax, name):
        # add background to stacked histo
        counts = self.backgrounds[name]["counts"]
        steps = np.append(np.insert(counts,0,0.0),0.0)
        # plot block and outline
        self.ax.bar(self.bin_centers, counts, width = self.bin_widths, bottom = self.stacked_block,
                color = self.backgrounds[name]["color"], label = rf"{name}")
        self.ax.step(self.step_edges, steps + self.stacked_step, color='black', linewidth = 0.5)
        # Update the bottom for the next process
        self.stacked_block += counts
        self.stacked_step += steps
        return True

    # def get_signal(self):
    #     # total signal
    #     for sig, info in self.signal.items():
    #         sig_counts = np.zeros(len(self.bin_centers))
    #         sig_sq_errors = np.zeros(len(self.bin_centers))
    #         for contribution in info["nodes"]:
    #             counts, errors = self.get_counts_errors(contribution)
    #             sig_counts += counts
    #             sig_sq_errors += errors**2
    #         # store counts and yields for each contribution
    #         self.signal[sig]["counts"] = sig_counts
    #         self.signal[sig]["errors"] = np.sqrt(sig_sq_errors)

    #     return True


    def plot_1D_histo(self, ratio_min=0.5, ratio_max=1.5):
        print("Plotting 1D histogram")
        # plot 1D histogram
        self.fig, (self.ax, self.ax_ratio) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]}, sharex=True, figsize=(9, 8))

        ## MAIN PLOT
        # add background contributions
        for bkg in self.backgrounds.keys():
            self.stack_background(self.ax, bkg)
        # add uncertainty on total MC
        self.stacked_block = np.insert(self.stacked_block, len(self.stacked_block), 0)
        self.ax.fill_between(self.bin_edges,
                    self.stacked_block-np.insert(self.total_background["error_down"], len(self.total_background["error_down"]), 0),
                    self.stacked_block+np.insert(self.total_background["error_up"], len(self.total_background["error_up"]), 0),
                    step="post", facecolor='none', hatch='////////', edgecolor='grey', linewidth=0, label = "Bkg. Uncert.")
        if not self.blind:
            # add data
            self.ax.errorbar(self.bin_centers, self.data['counts'], label='Observation', yerr=self.data['errors'], fmt='o', color = 'black', markersize=3, linewidth=0.6)
            self.ax.errorbar(self.bin_centers, self.data['counts'],xerr=self.bin_widths/2, fmt='o', color = 'black', markersize=3, linewidth=0.6) # add width marker


        # TODO: MAKE OPTION to add signal

        ## RATIO PLOT
        self.ax_ratio.axhline(1, color='black', linestyle=':')
        for l in [0.6, 0.8, 1.2, 1.4]: # add horizontal lines
            self.ax_ratio.axhline(l, color='darkgray', linestyle='dotted')
        self.fig.subplots_adjust(hspace=0.05)
        # add data/MC ratio as points
        if not self.blind:
            self.ax_ratio.errorbar(self.bin_centers, self.data_MC_ratio['counts'], yerr=self.data_MC_ratio['error_data'], xerr=self.bin_widths/2, fmt='o', color = 'black', markersize=3, linewidth=0.6)
        # add MC uncertainty as shaded band
        self.ax_ratio.fill_between(self.bin_edges,
            1 - np.insert(self.data_MC_ratio['error_MC_down'], len(self.data_MC_ratio['error_MC_down']), 0),
            1 + np.insert(self.data_MC_ratio['error_MC_up'], len(self.data_MC_ratio['error_MC_up']), 0),
            step="post", facecolor='none', hatch='////////', edgecolor='grey', linewidth=0)

        # legends and labels
        hep.cms.label(ax=self.ax, label="Preliminary", data=True, lumi=self.lumi, com=13.6, fontsize=16)
        self.ax.text(0.035, self.ch_label_height, self.channel_label, fontsize=18, fontweight="bold", transform=self.ax.transAxes)
        handles, labels = self.ax.get_legend_handles_labels()

        # add vertical lines if 2D unrolled:
        if self.is2Dunrolled:
            n_bins = len(self.bin_centers)
            nrows = len(self.var_dim_1)-1
            ncols = len(self.var_dim_2)-1
            # draw boundaries
            boundaries = np.arange(0, (nrows + 1) * ncols, ncols)
            for b in boundaries:
                self.ax_ratio.axvline(b, color='black', linestyle='--', linewidth=2)
                self.ax.axvline(b, color='black', linestyle='--', linewidth=2)
            if "BDT_pred_score," in self.variable:
                # add text for binning of variable 1
                label_loc = (np.arange(nrows)/nrows) + 0.03
                for i, l in zip(range(nrows), label_loc):
                    # print(self.var_dim_1[i], self.var_dim_1[i+1])
                    self.ax.text(l, 0.84, f"BDT ({self.var_dim_1[i]}, {self.var_dim_1[i+1]})", fontsize=16, transform=self.ax.transAxes)
            # place legend outside of plot
            plt.legend(handles[::-1], labels[::-1], loc='upper left', frameon=1, framealpha=1, bbox_to_anchor=(1.005, 5.2))
        else:
            if len(handles) <= 5:
                ncols=1
            else:
                ncols=2
            self.ax.legend(handles[::-1], labels[::-1], loc='upper right', frameon=1, framealpha=1, bbox_to_anchor=(0.98, 0.98), ncol=ncols)

        # main plot
        if "(GeV)" in self.variable_label:
            self.ax.set_ylabel(f"Events / {round(self.bin_widths[0],2)} GeV")
        else:
            self.ax.set_ylabel(f"Events / {round(self.bin_widths[0],2)}")
        if self.log_y:
            self.ax.set_yscale('log')
            self.ax.set_ylim(0.1, 10*np.max(self.stacked_block))
        else:
            self.ax.set_ylim(0, 1.9*np.max(self.stacked_block))
        self.ax.set_xlim(self.bin_edges[0], self.bin_edges[-1])
        # ratio plot
        self.ax_ratio.set_ylabel("Obs/Exp")
        self.ax_ratio.set_xlabel(self.variable_label)
        self.ax_ratio.set_ylim(ratio_min, ratio_max)

        # Save to pdf and png
        if self.save_name:
            save_path_png = self.save_name+".png"
            save_path_pdf = self.save_name+".pdf"
        else:
            save_path_png = self.file_name.split(os.sep)
            save_path_png.insert(-1, 'pngs')
            os.makedirs(os.sep.join(save_path_png[:-1]), exist_ok=True) # make dir if not exists
            save_path_png = os.sep.join(save_path_png).replace(".root", ".png")
            save_path_pdf = self.file_name.split(os.sep)
            save_path_pdf.insert(-1, 'pdfs')
            os.makedirs(os.sep.join(save_path_pdf[:-1]), exist_ok=True) # make dir if not exists
            save_path_pdf = os.sep.join(save_path_pdf).replace(".root", ".pdf")
        plt.savefig(save_path_png, bbox_inches='tight')
        print(f'Saved histogram to {save_path_png}')

        plt.savefig(save_path_pdf, bbox_inches='tight')
        print(f'Saved histogram to {save_path_pdf}')


if __name__ == "__main__":
    # histo = HTT_Histogram("/vols/cms/lcr119/offline/HiggsCP/TIDAL/Draw/ZpT_Recoil_v4_TEST/Run3_2022/control/mm/datacard_pt_tt_inclusive_mm_Run3_2022.root", "mm_inclusive", "mm", "Run3_2022", "pt_tt", log_y=False, blind=False)
    # histo.plot_1D_histo()



    histo = HTT_Histogram("/vols/cms/lcr119/offline/HiggsCP/TIDAL/Draw/test_datacards_newBDT/Run3_2022/cpdecay/tt/datacard_BDT_pred_score_vs_aco_rho_rho_higgs_rhorho_tt_Run3_2022.root", "tt_higgs_rhorho", "tt", "Run3_2022", "BDT_pred_score,aco_rho_rho[0.,0.7,0.8,0.9,1.0],[0.0,0.6283185307179586,1.2566370614359172,1.8849555921538759,2.5132741228718345,3.141592653589793,3.7699111843077517,4.39822971502571,5.026548245743669,5.654866776461628,6.283185307179586]", log_y=False, blind=True, is2Dunrolled=True)
    histo.plot_1D_histo()

    # histo = HTT_Histogram("/vols/cms/lcr119/offline/HiggsCP/TIDAL/Draw/DeriveIDSFs_May25/AntiIso_0p3/Run3_2022/sf_calculation/mt/datacard_m_vis_mTLt65_aiso_inclusive_mt_Run3_2022.root", "mt_inclusive_mTLt65_aiso", "mt", "Run3_2022", "m_vis", log_y=False, blind=False)

    # histo = HTT_Histogram("/vols/cms/lcr119/offline/HiggsCP/TIDAL/Draw/3104/BinnedSFs_DoubleTauJet/Run3_2023/control/tt/datacard_m_vis_1_cp_inclusive_tt_Run3_2023.root", "tt_higgs_rhorho", "tt", "Run3_2022", "BDT_pred_score,aco_rho_rho[0.,0.7,0.8,0.9,1.0],[0.0,0.6283185307179586,1.2566370614359172,1.8849555921538759,2.5132741228718345,3.141592653589793,3.7699111843077517,4.39822971502571,5.026548245743669,5.654866776461628,6.283185307179586]", blind=True, is2Dunrolled=True)
    # histo.plot_1D_histo()


