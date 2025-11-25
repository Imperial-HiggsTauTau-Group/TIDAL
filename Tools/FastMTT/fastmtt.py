import math
import vector
import numpy as np
import numba as nb
import os
import awkward as ak
import argparse
from alive_progress import alive_bar
import pyarrow.parquet as pq


def fastmtt(input_file: str, channel_: str, chunk_number: int, start_idx: int, end_idx: int, delta=1/1.15, reg_order=6, constrain=True, constrain_setting="Window", constrain_window=np.array([123, 127])):

    if constrain and constrain_setting == "Window":
        output_dir_name = f"fastmtt_{constrain_setting}_{constrain_window[0]}_{constrain_window[1]}"
    elif constrain and constrain_setting == "BreitWigner":
        output_dir_name = f"fastmtt_{constrain_setting}"
    else:
        output_dir_name = "fastmtt"

    if channel_ == "mt" or channel_ == "et":
        columns = [
            "met_pt", "met_phi", "met_covXX", "met_covXY", "met_covYY",
            "pt_1", "eta_1", "phi_1", "mass_1",
            "pt_2", "eta_2", "phi_2", "mass_2", "decayMode_2",
            "event", "run", "lumi"
        ]
    elif channel_ == "tt":
        columns = [
            "met_pt", "met_phi", "met_covXX", "met_covXY", "met_covYY",
            "pt_1", "eta_1", "phi_1", "mass_1", "decayMode_1",
            "pt_2", "eta_2", "phi_2", "mass_2", "decayMode_2",
            "event", "run", "lumi"
        ]

    events = ak.from_parquet(input_file, columns=columns)[start_idx:end_idx]
    events = ak.fill_none(events, -9999)
    nevents = len(events)

    if channel_ == "et":
        decay_type_1 = np.zeros(nevents, dtype=np.uint8)
        decay_type_2 = np.zeros(nevents, dtype=np.uint8)
    elif channel_ == "mt":
        decay_type_1 = np.ones(nevents, dtype=np.uint8)
        decay_type_2 = np.zeros(nevents, dtype=np.uint8)
    elif channel_ == "tt":
        decay_type_1 = np.ones(nevents, dtype=np.uint8) * 2
        decay_type_2 = np.ones(nevents, dtype=np.uint8) * 2

    # Create a directory for fastmtt outputs if it doesn't exist
    output_dir = os.path.join(os.path.dirname(input_file), output_dir_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # initialize global parameters
    m_ele = 0.51100e-3
    m_muon = 0.10566
    m_tau = 1.77685
    m_pion = 0.13957

    N = nevents

    pt_1 = (events.pt_1).to_numpy()
    eta_1 = (events.eta_1).to_numpy()
    phi_1 = (events.phi_1).to_numpy()
    mass1 = (events.mass_1).to_numpy()
    pt_2 = (events.pt_2).to_numpy()
    eta_2 = (events.eta_2).to_numpy()
    phi_2 = (events.phi_2).to_numpy()
    mass2 = (events.mass_2).to_numpy()
    met_pt = (events.met_pt).to_numpy()
    met_phi = (events.met_phi).to_numpy()
    metcov_xx = (events.met_covXX).to_numpy()
    metcov_xy = (events.met_covXY).to_numpy()
    metcov_yx = (events.met_covXY).to_numpy()
    metcov_yy = (events.met_covYY).to_numpy()

    met_x = met_pt * np.cos(met_phi)
    met_y = met_pt * np.sin(met_phi)

    fastmttMass_values, fastmttPt_values, fastmttPt1_values, fastmttPt2_values = compute_fastmtt(
        N,
        pt_1,
        eta_1,
        phi_1,
        mass1,
        pt_2,
        eta_2,
        phi_2,
        mass2,
        met_x,
        met_y,
        metcov_xx,
        metcov_xy,
        metcov_yx,
        metcov_yy,
        decay_type_1,
        decay_type_2,
        m_ele,
        m_muon,
        m_tau,
        m_pion,
        delta,
        reg_order,
        constrain,
        constrain_setting,
        constrain_window
    )

    # save the results to a parquet file
    columns_to_store = ["event", "run", "lumi"]
    chunk_events = ak.from_parquet(input_file, columns=columns_to_store)[start_idx:end_idx]
    chunk_events = ak.with_field(chunk_events, fastmttMass_values, "FastMTT_mass")
    chunk_events = ak.with_field(chunk_events, fastmttPt_values, "FastMTT_pt")
    chunk_events = ak.with_field(chunk_events, fastmttPt1_values, "FastMTT_pt_1")
    chunk_events = ak.with_field(chunk_events, fastmttPt2_values, "FastMTT_pt_2")

    output_file = os.path.join(output_dir, f"fastmtt_chunk{chunk_number}.parquet")
    ak.to_parquet(chunk_events, output_file)
    print("Job is done")


@nb.jit(nopython=True, parallel=False)
def compute_fastmtt(
    N,
    pt_1,
    eta_1,
    phi_1,
    mass1,
    pt_2,
    eta_2,
    phi_2,
    mass2,
    met_x,
    met_y,
    metcov_xx,
    metcov_xy,
    metcov_yx,
    metcov_yy,
    decay_type_1,
    decay_type_2,
    m_ele,
    m_muon,
    m_tau,
    m_pion,
    delta,
    reg_order,
    constrain,
    constrain_setting,
    constrain_window
):
    fastmttMass_values = np.zeros(N, dtype=np.float32)
    fastmttPt_values = np.zeros(N, dtype=np.float32)
    fastmttPt1_values = np.zeros(N, dtype=np.float32)
    fastmttPt2_values = np.zeros(N, dtype=np.float32)

    mass_dict = {0: m_ele, 1: m_muon, 2: m_tau}

    for i in range(N):

        # grab the correct masses based on tau decay type
        # tau decay_type: 0 ==> leptonic to electron,
        #                 1 ==> leptonic to muon,
        #                 2 ==> leptonic to hadronic
        if (decay_type_1[i] != 2):
            m1 = mass_dict[decay_type_1[i]]
        else:
            m1 = mass1[i]
        if (decay_type_2[i] != 2):
            m2 = mass_dict[decay_type_2[i]]
        else:
            m2 = mass2[i]

        # store visible masses
        m_vis_1 = m1
        m_vis_2 = m2

        # determine minimum and maximum possible masses
        m_vis_min_1, m_vis_max_1 = 0, 0
        m_vis_min_2, m_vis_max_2 = 0, 0
        if (decay_type_1[i] == 0):
            m_vis_min_1, m_vis_max_1 = m_ele, m_ele
        if (decay_type_1[i] == 1):
            m_vis_min_1, m_vis_max_1 = m_muon, m_muon
        if (decay_type_1[i] == 2):
            m_vis_min_1, m_vis_max_1 = m_pion, 1.5
        if (decay_type_2[i] == 0):
            m_vis_min_2, m_vis_max_2 = m_ele, m_ele
        if (decay_type_2[i] == 1):
            m_vis_min_2, m_vis_max_2 = m_muon, m_muon
        if (decay_type_2[i] == 2):
            m_vis_min_2, m_vis_max_2 = m_pion, 1.5
        if (m_vis_1 < m_vis_min_1):
            m_vis_1 = m_vis_min_1
        if (m_vis_1 > m_vis_max_1):
            m_vis_1 = m_vis_max_1
        if (m_vis_2 < m_vis_min_2):
            m_vis_2 = m_vis_min_2
        if (m_vis_2 > m_vis_max_2):
            m_vis_2 = m_vis_max_2

        # store both tau candidate four vectors
        leg1 = vector.obj(pt=pt_1[i], eta=eta_1[i], phi=phi_1[i], mass=m_vis_1)
        leg2 = vector.obj(pt=pt_2[i], eta=eta_2[i], phi=phi_2[i], mass=m_vis_2)

        # store visible mass of ditau pair
        m_vis = math.sqrt(2*leg1.pt*leg2.pt*(math.cosh(leg1.eta - leg2.eta) -
                                             math.cos(leg1.phi - leg2.phi)))

        # correct initial visible masses
        if (decay_type_1[i] == 2 and m_vis_1 > 1.5):
            m_vis_1 = 0.3
        if (decay_type_2[i] == 2 and m_vis_2 > 1.5):
            m_vis_2 = 0.3

        # invert met covariance matrix, calculate determinant
        metcovinv_xx, metcovinv_yy = metcov_yy[i], metcov_xx[i]
        metcovinv_xy, metcovinv_yx = -metcov_xy[i], -metcov_yx[i]
        metcovinv_det = (metcovinv_xx*metcovinv_yy -
                         metcovinv_yx*metcovinv_xy)
        if (metcovinv_det<1e-10):
                print("Warning! Ill-conditioned MET covariance at event index", i)
                continue

        # perform likelihood scan
        # see http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_032_v3.pdf
        met_const = 1/(2*math.pi*math.sqrt(metcovinv_det))
        min_likelihood, x1_opt, x2_opt = 999, 0.01, 0.01
        mass_likelihood, met_transfer = 0, 0

        initialise = True

        # scan over weights for each ditau four-vector
        for x1 in np.arange(0.01, 1, 0.01):
            for x2 in np.arange(0.01, 1, 0.01):
                x1_min = min(1, math.pow((m_vis_1/m_tau), 2))
                x2_min = min(1, math.pow((m_vis_2/m_tau),2))
                if ((x1 < x1_min) or (x2 < x2_min)):
                    continue

                # test weighted four-vectors
                leg1_x1, leg2_x2 = leg1*(1/x1), leg2*(1/x2)
                ditau_test = vector.obj(px=leg1_x1.px+leg2_x2.px,
                                     py=leg1_x1.py+leg2_x2.py,
                                     pz=leg1_x1.pz+leg2_x2.pz,
                                     E=leg1_x1.E+leg2_x2.E)
                nu_test = vector.obj(px=ditau_test.px-leg1.px-leg2.px,
                                  py=ditau_test.py-leg1.py-leg2.py,
                                  pz=ditau_test.pz-leg1.pz-leg2.pz,
                                  E=ditau_test.E-leg1.E-leg2.E)
                test_mass = ditau_test.mass

                if constrain_setting == "Window":
                    if (((test_mass < constrain_window[0]) or
                        (test_mass > constrain_window[1])) and
                        constrain):
                        continue

                # calculate mass likelihood integral
                m_shift = test_mass * delta
                if (m_shift < m_vis):
                    continue
                x1_min = min(1.0, math.pow((m_vis_1/m_tau),2))
                x2_min = max(math.pow((m_vis_2/m_tau),2),
                             math.pow((m_vis/m_shift),2))
                x2_max = min(1.0, math.pow((m_vis/m_shift),2)/x1_min)
                if (x2_max < x2_min):
                    continue
                J = 2*math.pow(m_vis,2) * math.pow(m_shift, -reg_order)
                I_x2 = math.log(x2_max) - math.log(x2_min)
                I_tot = I_x2
                if (decay_type_1[i] != 2):
                    I_m_nunu_1 = math.pow((m_vis/m_shift),2) * (math.pow(x2_max,-1) - math.pow(x2_min,-1))
                    I_tot += I_m_nunu_1
                if (decay_type_2[i] != 2):
                    I_m_nunu_2 = math.pow((m_vis/m_shift),2) * I_x2 - (x2_max - x2_min)
                    I_tot += I_m_nunu_2
                mass_likelihood = 1e9 * J * I_tot

                # calculate MET transfer function
                residual_x = met_x[i] - nu_test.x
                residual_y = met_y[i] - nu_test.y
                pull2 = (residual_x*(metcovinv_xx*residual_x +
                                     metcovinv_xy*residual_y) +
                         residual_y*(metcovinv_yx*residual_x +
                                     metcovinv_yy*residual_y))
                pull2 /= metcovinv_det
                met_transfer = met_const*math.exp(-0.5*pull2)

                # calculate final likelihood, store if minimum
                likelihood = -met_transfer * mass_likelihood

                if constrain and constrain_setting == "BreitWigner":
                    mH = 125.0
                    GammaH = 0.004
                    deltaM = test_mass*test_mass - mH*mH
                    mG = test_mass*GammaH
                    BreitWigner_likelihood = 1/(deltaM*deltaM + mG*mG)
                    likelihood = likelihood*BreitWigner_likelihood

                if initialise:
                    min_likelihood = likelihood
                    x1_opt, x2_opt = x1, x2
                    initialise = False
                else:
                    if (likelihood < min_likelihood):
                        min_likelihood = likelihood
                        x1_opt, x2_opt = x1, x2

        leg1_x1, leg2_x2 = leg1*(1/x1_opt), leg2*(1/x2_opt)
        p4_ditau_opt = vector.obj(px=leg1_x1.px+leg2_x2.px,
                               py=leg1_x1.py+leg2_x2.py,
                               pz=leg1_x1.pz+leg2_x2.pz,
                               E=leg1_x1.E+leg2_x2.E)

        mass_opt = p4_ditau_opt.mass
        pt_opt = p4_ditau_opt.pt
        pt1_opt = pt_1[i]/x1_opt
        pt2_opt = pt_2[i]/x2_opt

        fastmttMass_values[i] = mass_opt
        fastmttPt_values[i] = pt_opt
        fastmttPt1_values[i] = pt1_opt
        fastmttPt2_values[i] = pt2_opt

    return fastmttMass_values, fastmttPt_values, fastmttPt1_values, fastmttPt2_values


def run_processes(base_dir, constrain, constrain_setting, constrain_window, use_condor=False, chunk_size=10000):
    variations = [
"Electron_Smearing_down", "Tau_EnergyScale_PNet_ESCALE_3PRONG_1PI0_down", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_1PI0_up",
"Electron_Smearing_up", "Tau_EnergyScale_PNet_ESCALE_3PRONG_1PI0_up", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_down",
"jec_syst_Total_down", "Tau_EnergyScale_PNet_ESCALE_3PRONG_down", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_up",
"jec_syst_Total_up", "Tau_EnergyScale_PNet_ESCALE_3PRONG_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_1PI0_down",
"jer_syst_down",
"Tau_EnergyScale_PNet_JSCALE_down", "Tau_EnergyScale_PNet_TSCALE_1PRONG_1PI0_up",
"jer_syst_up",
"Tau_EnergyScale_PNet_JSCALE_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_2PI0_down",
"nominal", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_1PI0_down",  "Tau_EnergyScale_PNet_TSCALE_1PRONG_2PI0_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_1PI0_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_1PI0_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_1PI0_up", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_2PI0_down",   "Tau_EnergyScale_PNet_TSCALE_1PRONG_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_2PI0_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_2PI0_up", "Tau_EnergyScale_PNet_TSCALE_3PRONG_1PI0_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_2PI0_up", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_down", "Tau_EnergyScale_PNet_TSCALE_3PRONG_1PI0_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_up", "Tau_EnergyScale_PNet_TSCALE_3PRONG_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_up", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_1PI0_down", "Tau_EnergyScale_PNet_TSCALE_3PRONG_up"
    ]
    for channel in os.listdir(base_dir):
        if channel not in ["mt","et","tt"]:
           continue
        channel_dir = os.path.join(base_dir, channel)
        if os.path.isdir(channel_dir):
            for process in os.listdir(channel_dir):
                process_dir = os.path.join(channel_dir, process)
                if os.path.isdir(process_dir):
                    for variation in variations:
                        variation_dir = os.path.join(process_dir, variation)
                        print(variation_dir)
                        if os.path.isdir(variation_dir):
                            parquet_file = os.path.join(variation_dir, 'merged.parquet')
                            if os.path.exists(parquet_file):
                                run_process(parquet_file, constrain, constrain_setting, constrain_window, channel, use_condor, chunk_size)


def run_process(parquet_file, constrain, constrain_setting, constrain_window, channel, use_condor, chunk_size):
    events = ak.from_parquet(parquet_file, columns=["event"])  # Load only the event column to determine size
    nevents = len(events)
    del events

    for chunk_number, start_idx in enumerate(range(0, nevents, chunk_size)):
        end_idx = min(start_idx + chunk_size, nevents)

        if not use_condor:
            fastmtt(parquet_file, channel, chunk_number, start_idx, end_idx, constrain=constrain, constrain_setting=constrain_setting, constrain_window=constrain_window)
        else:
            shell_script_path = create_shell_script(parquet_file, constrain, constrain_setting, constrain_window, channel, chunk_number, start_idx, end_idx)
            submission_file_path = create_condor_submission_file(shell_script_path, chunk_number)
            submit_condor_job(submission_file_path)


def create_condor_submission_file(shell_script_path, chunk_number):
    directory = os.path.dirname(shell_script_path)
    submission_file_content = f"""
executable = {shell_script_path}
output = {directory}/fastmtt_chunk{chunk_number}_$(CLUSTER).out
error = {directory}/fastmtt_chunk{chunk_number}_$(CLUSTER).err
log = {directory}/fastmtt_chunk{chunk_number}_$(CLUSTER).log
request_memory = 8000
getenv = True
+MaxRuntime = 10800
queue
"""
    submission_file_path = os.path.join(directory, f'condor_submit_{shell_script_path.split("_")[-1].replace(".sh", "")}.sub')
    with open(submission_file_path, 'w') as f:
        f.write(submission_file_content)

    return submission_file_path


def create_shell_script(input_file: str, constrain, constrain_setting, constrain_window, channel: str, chunk_number: int, start_idx: int, end_idx: int):
    if constrain and constrain_setting == "Window":
        output_dir_name = f"fastmtt_{constrain_setting}_{constrain_window[0]}_{constrain_window[1]}"
    elif constrain and constrain_setting == "BreitWigner":
        output_dir_name = f"fastmtt_{constrain_setting}"
    else:
        output_dir_name = "fastmtt"
    shell_script = f"""
#!/bin/bash
python3 Tools/FastMTT/fastmtt.py \\
--input_file {input_file} \\
--channel {channel} \\
--chunk_number {chunk_number} \\
--start_idx {start_idx} \\
--end_idx {end_idx} \\
--running_dir """

    if constrain:
        shell_script += f"--constrain --constrain_setting {constrain_setting} --constrain_window {constrain_window[0]},{constrain_window[1]}"

    shell_script_path = os.path.join(os.path.dirname(input_file), output_dir_name, 'condor', f'fast_chunk{chunk_number}.sh')
    if not os.path.exists(os.path.dirname(shell_script_path)):
        os.makedirs(os.path.dirname(shell_script_path))
    with open(shell_script_path, 'w') as f:
        f.write(shell_script)
    os.system(f"chmod +x {shell_script_path}")
    return shell_script_path


def submit_condor_job(submission_file_path):
    print(f"Submitting job {submission_file_path}")
    os.system(f"condor_submit {submission_file_path}")


def merge_fastmtt_chunks(dir, output_dir_name, channels):
    variations = [
"Electron_Smearing_down", "Tau_EnergyScale_PNet_ESCALE_3PRONG_1PI0_down", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_1PI0_up",
"Electron_Smearing_up", "Tau_EnergyScale_PNet_ESCALE_3PRONG_1PI0_up", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_down",
"jec_syst_Total_down", "Tau_EnergyScale_PNet_ESCALE_3PRONG_down", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_up",
"jec_syst_Total_up", "Tau_EnergyScale_PNet_ESCALE_3PRONG_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_1PI0_down",
"jer_syst_down",
"Tau_EnergyScale_PNet_JSCALE_down", "Tau_EnergyScale_PNet_TSCALE_1PRONG_1PI0_up",
"jer_syst_up",
"Tau_EnergyScale_PNet_JSCALE_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_2PI0_down",
"nominal", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_1PI0_down",  "Tau_EnergyScale_PNet_TSCALE_1PRONG_2PI0_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_1PI0_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_1PI0_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_1PI0_up", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_2PI0_down",   "Tau_EnergyScale_PNet_TSCALE_1PRONG_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_2PI0_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_2PI0_up", "Tau_EnergyScale_PNet_TSCALE_3PRONG_1PI0_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_2PI0_up", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_down", "Tau_EnergyScale_PNet_TSCALE_3PRONG_1PI0_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_up", "Tau_EnergyScale_PNet_TSCALE_3PRONG_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_up", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_1PI0_down", "Tau_EnergyScale_PNet_TSCALE_3PRONG_up"
    ]
    for channel in channels.split(","):
        print("-"*50)
        print("Processing channel:", channel)
        directory = os.path.join(dir, channel)
        for process in os.listdir(directory):
            print(f"Merging for process: {process}")
            process_dir = os.path.join(directory, process)
            if os.path.isdir(process_dir):
                for variation in variations:
                    variation_dir = os.path.join(process_dir, variation)
                    if os.path.isdir(variation_dir):
                        fastmtt_dir = os.path.join(variation_dir, output_dir_name)
                        parquet_files = [os.path.join(fastmtt_dir, file) for file in os.listdir(fastmtt_dir) if file.endswith(".parquet")]
                        if len(parquet_files) == 1:
                            print("-"*50)
                            print(f"Copying single file for {process} and renaming to {output_dir_name}.parquet")
                            os.system(f"cp {parquet_files[0]} {variation_dir}/{output_dir_name}.parquet")
                        if len(parquet_files) > 1:
                            print("-"*50)
                            print(f"Merging fastmtt chunks for {process}")
                            schema = None
                            writer = None
                            with alive_bar(len(parquet_files)) as bar:
                                for f in parquet_files:
                                    dataset = pq.ParquetDataset(f)
                                    table = dataset.read()
                                    if schema is None:
                                        schema = table.schema
                                        writer = pq.ParquetWriter(os.path.join(variation_dir, f"{output_dir_name}.parquet"), schema)
                                    else:
                                        table = table.cast(schema)
                                    writer.write_table(table)
                                    bar()
                            if writer is not None:
                                writer.close()


def check_logs(dir, output_dir_name, channels):
    variations = [
"Electron_Smearing_down", "Tau_EnergyScale_PNet_ESCALE_3PRONG_1PI0_down", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_1PI0_up",
"Electron_Smearing_up", "Tau_EnergyScale_PNet_ESCALE_3PRONG_1PI0_up", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_down",
"jec_syst_Total_down", "Tau_EnergyScale_PNet_ESCALE_3PRONG_down", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_up",
"jec_syst_Total_up", "Tau_EnergyScale_PNet_ESCALE_3PRONG_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_1PI0_down",
"jer_syst_down",
"Tau_EnergyScale_PNet_JSCALE_down", "Tau_EnergyScale_PNet_TSCALE_1PRONG_1PI0_up",
"jer_syst_up",
"Tau_EnergyScale_PNet_JSCALE_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_2PI0_down",
"nominal", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_1PI0_down",  "Tau_EnergyScale_PNet_TSCALE_1PRONG_2PI0_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_1PI0_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_1PI0_up", "Tau_EnergyScale_PNet_TSCALE_1PRONG_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_1PI0_up", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_2PI0_down",   "Tau_EnergyScale_PNet_TSCALE_1PRONG_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_2PI0_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_2PI0_up", "Tau_EnergyScale_PNet_TSCALE_3PRONG_1PI0_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_2PI0_up", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_down", "Tau_EnergyScale_PNet_TSCALE_3PRONG_1PI0_up",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_down", "Tau_EnergyScale_PNet_MUSCALE_1PRONG_up", "Tau_EnergyScale_PNet_TSCALE_3PRONG_down",
"Tau_EnergyScale_PNet_ESCALE_1PRONG_up", "Tau_EnergyScale_PNet_MUSCALE_3PRONG_1PI0_down", "Tau_EnergyScale_PNet_TSCALE_3PRONG_up"
    ]
    for channel in channels.split(","):
        print("-"*50)
        print("Processing channel:", channel)

        count_failed_jobs = 0
        count_running_jobs = 0
        count_done = 0
        count_sub_files = 0

        directory = os.path.join(dir, channel)
        for process in os.listdir(directory):
            print(f"Checking process: {process}")
            process_dir = os.path.join(directory, process)
            if os.path.isdir(process_dir):
                for variation in variations:
                    variation_dir = os.path.join(process_dir, variation)
                    if os.path.isdir(variation_dir):
                        log_dir = os.path.join(variation_dir, output_dir_name, "condor")
                        # remove log_dir
                        #if os.path.exists(log_dir):
                        #    os.system(f"rm -r {log_dir}")
                        #continue
                        if not os.path.exists(os.path.join(variation_dir, output_dir_name, "resubmit")):
                            os.makedirs(os.path.join(variation_dir, output_dir_name, "resubmit"))
                        for log_file in os.listdir(log_dir):
                            job_status_running = "Total for query: 1 jobs; 0 completed, 0 removed, 0 idle, 1 running, 0 held, 0 suspended"
                            job_status_idle = "Total for query: 1 jobs; 0 completed, 0 removed, 1 idle, 0 running, 0 held, 0 suspended"
                            job_status_held = "Total for query: 1 jobs; 0 completed, 0 removed, 0 idle, 0 running, 1 held, 0 suspended"
                            if log_file.endswith(".sub"):
                                count_sub_files += 1
                            if log_file.endswith(".log"):
                                chunk_number = log_file.split("_")[-2].replace("chunk", "")
                                sub_file = f"condor_submit_chunk{chunk_number}.sub"
                                out_file = log_file.replace(".log", ".out")
                                clusterID = log_file.split("_")[-1].replace(".log", "")
                                condor_q = os.popen(f"condor_q {clusterID}").read()

                                if os.path.exists(os.path.join(log_dir, out_file)):
                                    with open(os.path.join(log_dir, out_file), 'r') as f:
                                        lines = f.readlines()
                                        if len(lines) > 0:
                                            last_line = lines[-1]

                                            if "Job is done" not in last_line:
                                                if job_status_running not in condor_q and job_status_idle not in condor_q:
                                                    count_failed_jobs += 1

                                                    if job_status_held not in condor_q:
                                                        os.system(f"mv {os.path.join(log_dir, log_file)} {os.path.join(variation_dir, output_dir_name, 'resubmit')}")
                                                        os.system(f"mv {os.path.join(log_dir, log_file.replace('.log', '.err'))} {os.path.join(variation_dir, output_dir_name, 'resubmit')}")
                                                        os.system(f"mv {os.path.join(log_dir, out_file)} {os.path.join(variation_dir, output_dir_name, 'resubmit')}")
                                                        os.system(f"condor_submit {os.path.join(log_dir, sub_file)}")

                                                elif job_status_running in condor_q:
                                                    count_running_jobs += 1

                                            elif "Job is done" in last_line:
                                                if job_status_running in condor_q:
                                                    os.popen(f"condor_rm {clusterID}")
                                                    print(f"Job {clusterID} is done but still running. Removing the job")
                                                count_done += 1
                                        else:
                                            if job_status_running in condor_q or job_status_idle in condor_q:
                                                count_running_jobs += 1
                                            else:
                                                count_failed_jobs += 1
                                                if job_status_held not in condor_q:
                                                    os.system(f"mv {os.path.join(log_dir, log_file)} {os.path.join(variation_dir, output_dir_name, 'resubmit')}")
                                                    os.system(f"mv {os.path.join(log_dir, log_file.replace('.log', '.err'))} {os.path.join(variation_dir, output_dir_name, 'resubmit')}")
                                                    os.system(f"mv {os.path.join(log_dir, out_file)} {os.path.join(variation_dir, output_dir_name, 'resubmit')}")
                                                    os.system(f"condor_submit {os.path.join(log_dir, sub_file)}")
                                else:
                                    if job_status_running in condor_q or job_status_idle in condor_q:
                                        count_running_jobs += 1
                                    else:
                                        count_failed_jobs += 1
                                        if job_status_held not in condor_q:
                                            os.system(f"mv {os.path.join(log_dir, log_file)} {os.path.join(variation_dir, output_dir_name, 'resubmit')}")
                                            os.system(f"condor_submit {os.path.join(log_dir, sub_file)}")

        print(f"Number of jobs currently running: {count_running_jobs}")
        print(f"Number of finished jobs: {count_done}")
        print(f"Total number of jobs: {count_sub_files}")
        print(f"Number of failed jobs: {count_failed_jobs}")


def main():
    parser = argparse.ArgumentParser(description="Run SVFit")
    parser.add_argument("--source_dir", type=str, help="Source directory with parquet files e.g. /vols/cms/ks1021/offline/HiggsDNA/IC/outputs/production_v3110/Run3_2022/")
    parser.add_argument("--output_dir_name", type=str, help="Output directory name")
    parser.add_argument("--constrain", action="store_true", help="Constrain the mass")
    parser.add_argument("--constrain_setting", type=str, default="Window", help="Constrain setting (Window or BreitWigner)")
    parser.add_argument("--constrain_window", type=str, default="123,127", help="Constrain window (used by condor submission)")

    parser.add_argument("--chunk_size", type=int, default=500000, help="Number of events per chunk")
    parser.add_argument("--use_condor", action="store_true", help="Use condor for parallel processing")
    parser.add_argument("--input_file", type=str, help="Input file (used by condor submission)")
    parser.add_argument("--channel", type=str, help="Channel (used by condor submission)")
    parser.add_argument("--chunk_number", type=int, help="Chunk number (used by condor submission)")
    parser.add_argument("--start_idx", type=int, help="Start index for the chunk (used by condor submission)")
    parser.add_argument("--end_idx", type=int, help="End index for the chunk (used by condor submission)")
    parser.add_argument("--running_dir", action="store_true", help="Running directory (used by condor submission)")

    parser.add_argument("--check_logs", action="store_true", help="Check logs for failed jobs")
    parser.add_argument("--merge_chunks", action="store_true", help="Merge svfit chunks")
    parser.add_argument("--channels", type=str, help="Channels to check/merge e.g. mt,et,tt")

    args = parser.parse_args()

    constrain_window = np.array([float(x) for x in args.constrain_window.split(",")])

    if args.check_logs:
        check_logs(args.source_dir, args.output_dir_name, args.channels)
    elif args.merge_chunks:
        merge_fastmtt_chunks(args.source_dir, args.output_dir_name, args.channels)
    else:
        if not args.running_dir:
            run_processes(args.source_dir, args.constrain, args.constrain_setting, constrain_window, args.use_condor, args.chunk_size)
        else:
            fastmtt(args.input_file, args.channel, args.chunk_number, args.start_idx, args.end_idx, constrain=args.constrain, constrain_setting=args.constrain_setting, constrain_window=constrain_window)

if __name__ == "__main__":
    main()


