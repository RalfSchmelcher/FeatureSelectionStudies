from JobSubmitter import JobSubmitter
import glob
import os

def submit_sample(name, files):
    script_dir = "/afs/desy.de/user/s/schmelch/GenLevelStudies"
    base_dir = "/nfs/dust/cms/user/schmelch/customTrees"

    workdir = "{}/workdir_{}".format(base_dir,name)
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
            "echo \"created tmp workdir in ${TMP_DIR}\"",
            "pwd",
            "source /cvmfs/cms.cern.ch/cmsset_default.sh",
            "cmsrel CMSSW_11_3_4 && cd CMSSW_11_3_4/src && eval `scramv1 runtime -sh`",
            "pwd",
    ]
    

    for job_index, fname in enumerate(files):
        wrapper_lines = common_wrapper_lines
        if name == "signal":
            wrapper_lines += [            
                "cp {}/dummy_nanoAOD_plotter.py .".format(script_dir),
                "cp {}/hvt_mass_hypotheses_util.py .".format(script_dir),
                "python3 dummy_nanoAOD_plotter.py {} {}".format(fname, job_index),
                "cp *.root {}".format(workdir),
                "cp mass_parser_jsons/*.json {}".format(json_dir)
            ]
        else:
            wrapper_lines += [
                "cp {}/dummy_nanoAOD_plotter_bkg.py .".format(script_dir),
                "python3 dummy_nanoAOD_plotter_bkg.py {} {}".format(fname, job_index),
                "cp *.root {}".format(workdir),

            ]

        submitter.add_job("{}_TreeMaker_{}".format(name, job_index), ".", wrapper_lines)

    submitter.submit_jobs(debug=False)

submit_sample("signal", glob.glob("/pnfs/desy.de/cms/tier2/store/user/anmehta/resSrch/XToYYprime_YToQQ_YprimeToQQ_narrow_TuneCP5_PSWeights_13TeV-madgraph-pythia8_RunIIAutumn18/*.root"))
submit_sample("background", glob.glob("/pnfs/desy.de/cms/tier2/store/user/anmehta/resSrch/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18/*.root"))


#hadd RunIIAutumn18_all.root workdir/RunIIAutumn18_*.root