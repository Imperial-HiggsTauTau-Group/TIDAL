import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from TauAnalysis.ClassicSVfit.wrapper.pybind_wrapper import ClassicSVfit, TauTauHistogramAdapter, MeasuredTauLepton
import ROOT as root
import vector
import awkward as ak
import numpy as np
import os
from alive_progress import alive_bar
import argparse
import pyarrow.parquet as pq


def calculate_p4(pt: ak.Array, eta: ak.Array, phi: ak.Array, mass: ak.Array) -> ak.Array:
    return vector.awk(
        ak.zip(
            {
                "pt": pt,
                "eta": eta,
                "phi": phi,
                "mass": mass
            },
            with_name="Momentum4D"
        )
    )


def is_valid_p4(p4: ak.Array) -> bool:
    is_nan_or_inf = np.isnan(p4.pt) | np.isnan(p4.eta) | np.isnan(p4.phi) | np.isnan(p4.mass) | np.isinf(p4.pt) | np.isinf(p4.eta) | np.isinf(p4.phi) | np.isinf(p4.mass)
    return not is_nan_or_inf


def process_event(lep1, lep2, mode, kappa, metx, mety, covMatrix,
                  lep1_MeasuredTauLepton, lep2_MeasuredTauLepton,
                  do_massConstraint, mass_constraint, verbosity):
    if is_valid_p4(lep1) and is_valid_p4(lep2):
        svFitAlgo = ClassicSVfit(verbosity)
        svFitAlgo.setHistogramAdapter(TauTauHistogramAdapter())
        svFitAlgo.addLogM_fixed(True, kappa)
        if do_massConstraint:
            svFitAlgo.setDiTauMassConstraint(mass_constraint)

        measuredTauLeptons = [lep1_MeasuredTauLepton, lep2_MeasuredTauLepton]

        svFitAlgo.integrate(measuredTauLeptons, metx, mety, covMatrix)

        isValidSolution = svFitAlgo.isValidSolution()

        tau1_P4 = svFitAlgo.getHistogramAdapter().GetFittedTau1LV()
        tau2_P4 = svFitAlgo.getHistogramAdapter().GetFittedTau2LV()

        if isValidSolution:
            return (svFitAlgo.getHistogramAdapter().getMass(), svFitAlgo.getHistogramAdapter().getMassErr(), tau1_P4.Pt(), tau2_P4.Pt())
        else:
            return (-1, -1, -1, -1)
    else:
        return (-1, -1, -1, -1)


def run_svfit(input_file: str, channel_: str, chunk_number: int, start_idx: int, end_idx: int):
    # Load the specific chunk of data from the parquet file

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

    # Create a directory for svfit outputs if it doesn't exist
    output_dir = os.path.join(os.path.dirname(input_file), "svfit")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    channel = ak.Array([channel_] * nevents)
    dummy = ak.Array([0] * nevents)

    mode = ak.where(
        channel == "tt",
        ak.ones_like(dummy),
        ak.where(
            channel == "et",
            ak.ones_like(dummy) * 2,
            ak.where(
                channel == "mt",
                ak.ones_like(dummy) * 3,
                ak.zeros_like(dummy)
            )
        )
    )

    kappa = ak.where(
        channel == "em",
        ak.ones_like(dummy) * 3,
        ak.where(
            channel == "mt",
            ak.ones_like(dummy) * 4,
            ak.where(
                channel == "et",
                ak.ones_like(dummy) * 4,
                ak.where(
                    channel == "tt",
                    ak.ones_like(dummy) * 5,
                    ak.ones_like(dummy) * 3
                )
            )
        )
    )

    metx = events.met_pt * np.cos(events.met_phi)
    mety = events.met_pt * np.sin(events.met_phi)

    covXX = events.met_covXX
    covXY = events.met_covXY
    covYY = events.met_covYY

    covMET = np.zeros((len(covXX), 2, 2))
    covMET[:, 0, 0] = covXX
    covMET[:, 1, 0] = covXY
    covMET[:, 0, 1] = covXY
    covMET[:, 1, 1] = covYY

    rootMETMatrices = []

    for cov in covMET:
        matrix = root.TMatrixD(2, 2)
        matrix[0][0] = cov[0][0]
        matrix[1][0] = cov[1][0]
        matrix[0][1] = cov[0][1]
        matrix[1][1] = cov[1][1]
        rootMETMatrices.append(matrix)

    lep_1 = calculate_p4(events.pt_1, events.eta_1, events.phi_1, events.mass_1)
    lep_2 = calculate_p4(events.pt_2, events.eta_2, events.phi_2, events.mass_2)

    if channel_ == "tt":
        leps1_MeasuredTauLepton = np.array(
            [MeasuredTauLepton(*args) for args in zip(mode, lep_1.pt, lep_1.eta, lep_1.phi, lep_1.mass, events.decayMode_1)]
        )
    elif channel_ == "et" or channel_ == "mt":
        dummy_decayMode = ak.Array([-1] * nevents)
        leps1_MeasuredTauLepton = np.array(
            [MeasuredTauLepton(*args) for args in zip(mode, lep_1.pt, lep_1.eta, lep_1.phi, lep_1.mass, dummy_decayMode)]
        )

    leps2_MeasuredTauLepton = np.array(
        [MeasuredTauLepton(*args) for args in zip(mode, lep_2.pt, lep_2.eta, lep_2.phi, lep_2.mass, events.decayMode_2)]
    )

    del events

    verbosity = 0
    do_massConstraint = False
    mass_constraint = -1

    svfitMass_values = [0] * nevents
    svfitMass_err_values = [0] * nevents
    pt_1_values = [0] * nevents
    pt_2_values = [0] * nevents

    with alive_bar(nevents) as bar:
        for index in range(nevents):
            mass, mass_err, pt_1, pt_2 = process_event(
                lep_1[index], lep_2[index], mode[index], kappa[index],
                metx[index], mety[index], rootMETMatrices[index],
                leps1_MeasuredTauLepton[index], leps2_MeasuredTauLepton[index],
                do_massConstraint, mass_constraint, verbosity
            )
            svfitMass_values[index] = mass
            svfitMass_err_values[index] = mass_err
            pt_1_values[index] = pt_1
            pt_2_values[index] = pt_2

            bar()

    columns_to_store = ["event", "run", "lumi"]
    chunk_events = ak.from_parquet(input_file, columns=columns_to_store)[start_idx:end_idx]
    chunk_events = ak.with_field(chunk_events, svfitMass_values, "svfitMass")
    chunk_events = ak.with_field(chunk_events, svfitMass_err_values, "svfitMass_err")
    chunk_events = ak.with_field(chunk_events, pt_1_values, "tau1_pt_svfit")
    chunk_events = ak.with_field(chunk_events, pt_2_values, "tau2_pt_svfit")

    output_file = os.path.join(output_dir, f"svfit_chunk{chunk_number}.parquet")
    ak.to_parquet(chunk_events, output_file)
    print("Job is done")
    sys.exit(0)


def run_processes(base_dir, use_condor=False, chunk_size=10000):
    for channel in os.listdir(base_dir):
        channel_dir = os.path.join(base_dir, channel)
        if os.path.isdir(channel_dir):
            for process in os.listdir(channel_dir):
                process_dir = os.path.join(channel_dir, process)
                if os.path.isdir(process_dir):
                    variation_dir = os.path.join(process_dir, "nominal")
                    if os.path.isdir(variation_dir):
                        parquet_file = os.path.join(variation_dir, 'merged.parquet')
                        if os.path.exists(parquet_file):
                            run_process(parquet_file, channel, use_condor, chunk_size)


def run_process(parquet_file, channel, use_condor, chunk_size):
    events = ak.from_parquet(parquet_file, columns=["event"])  # Load only the event column to determine size
    nevents = len(events)
    del events

    for chunk_number, start_idx in enumerate(range(0, nevents, chunk_size)):
        end_idx = min(start_idx + chunk_size, nevents)

        if not use_condor:
            print(f"Running chunk {chunk_number} for {channel} from {start_idx} to {end_idx}")
            run_svfit(parquet_file, channel, chunk_number, start_idx, end_idx)
        else:
            shell_script_path = create_shell_script(parquet_file, channel, chunk_number, start_idx, end_idx)
            submission_file_path = create_condor_submission_file(shell_script_path, chunk_number)
            submit_condor_job(submission_file_path)

def create_shell_script(parquet_file, channel, chunk_number, start_idx, end_idx):
    script_content = f"""
#!/bin/bash
LD_LIBRARY_PATH=/usr/local/lib
CPPYY_BACKEND_LIBRARY=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/libcppyy_backend3_9.so
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/bin/thisroot.sh
export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/

python3 scripts/run_svfit.py --input_file {parquet_file} --channel {channel} --chunk_number {chunk_number} --start_idx {start_idx} --end_idx {end_idx} --running_dir
"""
    shell_script_path = os.path.join(os.path.dirname(parquet_file), 'svfit', 'condor', f'svfit_chunk{chunk_number}.sh')
    if not os.path.exists(os.path.dirname(shell_script_path)):
        os.makedirs(os.path.dirname(shell_script_path))
    with open(shell_script_path, 'w') as f:
        f.write(script_content)
    #print(f"Shell script created at {shell_script_path}")
    os.system(f"chmod +x {shell_script_path}")
    return shell_script_path


def create_condor_submission_file(shell_script_path, chunk_number):
    directory = os.path.dirname(shell_script_path)
    submission_file_content = f"""
executable = {shell_script_path}
output = {directory}/svfit_chunk{chunk_number}_$(CLUSTER).out
error = {directory}/svfit_chunk{chunk_number}_$(CLUSTER).err
log = {directory}/svfit_chunk{chunk_number}_$(CLUSTER).log
request_memory = 8G
getenv = True
+MaxRuntime = 10800
queue
"""
    submission_file_path = os.path.join(directory, f'condor_submit_{shell_script_path.split("_")[-1].replace(".sh", "")}.sub')
    with open(submission_file_path, 'w') as f:
        f.write(submission_file_content)

    return submission_file_path


def submit_condor_job(submission_file_path):
    print(f"Submitting job {submission_file_path}")
    os.system(f"condor_submit {submission_file_path}")


def check_logs(directory, channel):
    count_failed_jobs = 0
    count_running_jobs = 0
    count_done = 0
    count_sub_files = 0

    directory = os.path.join(directory, channel)
    for process in os.listdir(directory):
        process_dir = os.path.join(directory, process)
        if os.path.isdir(process_dir):
            variation_dir = os.path.join(process_dir, "nominal")
            if os.path.isdir(variation_dir):
                log_dir = os.path.join(variation_dir, "svfit", "condor")
                if not os.path.exists(os.path.join(variation_dir, "svfit", "resubmit")):
                    os.makedirs(os.path.join(variation_dir, "svfit", "resubmit"))
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
                                                os.system(f"mv {os.path.join(log_dir, log_file)} {os.path.join(variation_dir, 'svfit', 'resubmit')}")
                                                os.system(f"mv {os.path.join(log_dir, log_file.replace('.log', '.err'))} {os.path.join(variation_dir, 'svfit', 'resubmit')}")
                                                os.system(f"mv {os.path.join(log_dir, out_file)} {os.path.join(variation_dir, 'svfit', 'resubmit')}")
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
                                            os.system(f"mv {os.path.join(log_dir, log_file)} {os.path.join(variation_dir, 'svfit', 'resubmit')}")
                                            os.system(f"mv {os.path.join(log_dir, log_file.replace('.log', '.err'))} {os.path.join(variation_dir, 'svfit', 'resubmit')}")
                                            os.system(f"mv {os.path.join(log_dir, out_file)} {os.path.join(variation_dir, 'svfit', 'resubmit')}")
                                            os.system(f"condor_submit {os.path.join(log_dir, sub_file)}")
                        else:
                            if job_status_running in condor_q or job_status_idle in condor_q:
                                count_running_jobs += 1
                            else:
                                count_failed_jobs += 1
                                if job_status_held not in condor_q:
                                    os.system(f"mv {os.path.join(log_dir, log_file)} {os.path.join(variation_dir, 'svfit', 'resubmit')}")
                                    os.system(f"condor_submit {os.path.join(log_dir, sub_file)}")

    print(f"Number of jobs currently running: {count_running_jobs}")
    print(f"Number of finished jobs: {count_done}")
    print(f"Total number of jobs: {count_sub_files}")
    print(f"Number of failed jobs: {count_failed_jobs}")


def merge_svfit_chunks(directory, channel):
    directory = os.path.join(directory, channel)
    for process in os.listdir(directory):
        process_dir = os.path.join(directory, process)
        if os.path.isdir(process_dir):
            variation_dir = os.path.join(process_dir, "nominal")
            if os.path.isdir(variation_dir):
                svfit_dir = os.path.join(variation_dir, "svfit")
                parquet_files = [os.path.join(svfit_dir, file) for file in os.listdir(svfit_dir) if file.endswith(".parquet")]
                if len(parquet_files) == 1:
                    print(f"Copying single file for {process} and renaming to svfit.parquet")
                    os.system(f"cp {parquet_files[0]} {variation_dir}/svfit.parquet")
                if len(parquet_files) > 1:
                    print(f"Merging svfit chunks for {process}")
                    schema = None
                    writer = None
                    with alive_bar(len(parquet_files)) as bar:
                        for f in parquet_files:
                            dataset = pq.ParquetDataset(f)
                            table = dataset.read()
                            if schema is None:
                                schema = table.schema
                                writer = pq.ParquetWriter(os.path.join(variation_dir, "svfit.parquet"), schema)
                            else:
                                table = table.cast(schema)
                            writer.write_table(table)
                            bar()
                    if writer is not None:
                        writer.close()


def main():
    parser = argparse.ArgumentParser(description="Run SVFit")
    parser.add_argument("--source_dir", type=str, help="Source directory with parquet files e.g. /vols/cms/ks1021/offline/HiggsDNA/IC/output/test/Run3_2022/")
    parser.add_argument("--use_condor", action="store_true", help="Use condor submission")
    parser.add_argument("--input_file", type=str, help="Input file (used by condor submission)")
    parser.add_argument("--channel", type=str, help="Channel (used by condor submission)")
    parser.add_argument("--chunk_number", type=int, help="Chunk number (used by condor submission)")
    parser.add_argument("--start_idx", type=int, help="Start index for the chunk (used by condor submission)")
    parser.add_argument("--end_idx", type=int, help="End index for the chunk (used by condor submission)")
    parser.add_argument("--running_dir", action="store_true", help="Running directory (used by condor submission)")
    parser.add_argument("--chunk_size", type=int, default=10000, help="Number of events per chunk")
    parser.add_argument("--check_logs", action="store_true", help="Check logs for failed jobs")
    parser.add_argument("--merge_chunks", action="store_true", help="Merge svfit chunks")
    args = parser.parse_args()

    if args.check_logs:
        check_logs(args.source_dir, args.channel)
    elif args.merge_chunks:
        merge_svfit_chunks(args.source_dir, args.channel)
    else:
        if not args.running_dir:
            run_processes(args.source_dir, args.use_condor, args.chunk_size)
        else:
            run_svfit(args.input_file, args.channel, args.chunk_number, args.start_idx, args.end_idx)

if __name__ == "__main__":
    main()
