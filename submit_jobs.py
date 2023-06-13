from JobSubmitter import JobSubmitter
import glob
import os


def submit_sample(script_dir, base_dir, name, files, debug):

    workdir = "{}/workdir_{}".format(base_dir, name)
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    submitter = JobSubmitter(workdir)

    json_dir = "{}/mass_parser_jsons".format(workdir)
    if not os.path.exists(json_dir):
        os.makedirs(json_dir)

    common_wrapper_lines = [
        # "cd {}".format(script_dir),
        # "source setup_env.sh",
        # "cd {}/workdir".format(base_dir),
        "#!/bin/bash",
        "TMP_DIR=$(mktemp -d -t ${USER}.XXXX)",
        "cd $TMP_DIR",
        'echo "created tmp workdir in ${TMP_DIR}"',
        "pwd",
        "source /cvmfs/cms.cern.ch/cmsset_default.sh",
        "cmsrel CMSSW_11_3_4 && cd CMSSW_11_3_4/src && eval `scramv1 runtime -sh`",
        "pwd",
    ]

    wrapper_lines = [] + common_wrapper_lines
    if name == "signal":
        wrapper_lines += [
            "ls -lah",
            "cp {}/dummy_nanoAOD_plotter.py .".format(script_dir),
            "cp {}/hvt_mass_hypotheses_util.py .".format(script_dir),
            'echo "processing signal tree"',
            "python3 dummy_nanoAOD_plotter.py $1 $2",
            "ls -lah",
            'echo "copying outputfiles"',
            "cp *.root {}".format(workdir),
            "cp mass_parser_jsons/*.json {}".format(json_dir),
        ]
    else:
        wrapper_lines += [
            "ls -lah",
            "cp {}/dummy_nanoAOD_plotter_bkg.py .".format(script_dir),
            'echo "processing background tree"',
            "python3 dummy_nanoAOD_plotter_bkg.py $1 $2",
            "ls -lah",
            'echo "copying outputfiles"',
            "cp *.root {}".format(workdir),
        ]

    submitter.add_job_array(
        name="{}_TreeMaker".format(name),
        subworkdir=".",
        wrapper_lines=wrapper_lines,
        arguments=["$(InputFile) $(MyJobIndex)"],
        queue="MyJobIndex InputFile from ({}\n)".format("\n".join(" ".join(map(str, job)) for job in enumerate(files))),
    )
    submitter.submit_jobs(debug=debug)


if __name__ == "__main__":
    import argparse

    script_dir = "/afs/desy.de/user/s/schmelch/GenLevelStudies"
    base_dir = "/nfs/dust/cms/user/schmelch/customTrees"

    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", help="Do not submit any jobs but create all necessary files.")
    parser.add_argument("--hadd", action="store_true", help="Hadd all trees.")
    parser.add_argument("--submit", action="store_true", help="Submit jobs to write trees.")

    args = parser.parse_args()

    signal_samples = glob.glob(
        "/pnfs/desy.de/cms/tier2/store/user/anmehta/resSrch/XToYYprime_YToQQ_YprimeToQQ_narrow_TuneCP5_PSWeights_13TeV"
        "-madgraph-pythia8_RunIIAutumn18/*.root"
    )

    background_samples = glob.glob(
        "/pnfs/desy.de/cms/tier2/store/user/anmehta/resSrch/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_"
        "RunIIAutumn18/*.root"
    )
    if args.submit:
        submit_sample(script_dir, base_dir, "signal", signal_samples, debug=args.debug)
        submit_sample(script_dir, base_dir, "background", background_samples, debug=args.debug)
    if args.hadd:
        os.system("hadd -f {}/signal.root {}/workdir_signal/*.root".format(base_dir, base_dir))
        os.system("hadd -f {}/background.root {}/workdir_background/*.root".format(base_dir, base_dir))
