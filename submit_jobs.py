#!/usr/bin/env python3
from __future__ import print_function
from JobSubmitter import JobSubmitter
import ROOT
import glob
import os
import logging
import json
logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def test_root_file(fname, tree_name="tree"):
    try:
        f = ROOT.TFile.Open(fname)
        tree = f.Get(tree_name)
        n_entries = tree.GetEntries() # noqa
        return True
    except BaseException as e:
        logger.warn(e)
        return False


class SampleManager(object):
    def __init__(self, name, file_pattern, base_dir, script_dir, signal=True):
        self._output_file_prefix = "RunIIAutumn18"
        self._name = name
        self._file_pattern = file_pattern
        self._base_dir = base_dir
        self._script_dir = script_dir
        self._signal = signal
        self.workdir = "{}/workdir_{}".format(self._base_dir, self._name)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        self._resubmit_counter = 0
        self._merged = False
        self.status_lut = {0: "unsubmitted", 1: "submitted", 2: "finished"}

        self._manager_info_path = "{}/manager_info.json".format(self.workdir)
        if os.path.isfile(self._manager_info_path):
            self.load_info()
        else:
            self.files = glob.glob(self._file_pattern)
            self.njobs = len(self.files)
            self.jobids = range(self.njobs)
            self.status = [0] * self.njobs
        self.dump_info()

        self._wrapper_lines = [
            "#!/bin/bash",
            "TMP_DIR=$(mktemp -d -t ${USER}.XXXX)",
            "cd $TMP_DIR",
            'echo "created tmp workdir in ${TMP_DIR}"',
            "pwd",
            "source /cvmfs/cms.cern.ch/cmsset_default.sh",
            "cmsrel CMSSW_11_3_4 && cd CMSSW_11_3_4/src && eval `scramv1 runtime -sh`",
            "pwd",
        ]

        if self._signal:
            json_dir = "{}/mass_parser_jsons".format(self.workdir)
            if not os.path.exists(json_dir):
                os.makedirs(json_dir)
            self._wrapper_lines += [
                "ls -lah",
                "cp {}/dummy_nanoAOD_plotter.py .".format(self._script_dir),
                "cp {}/hvt_mass_hypotheses_util.py .".format(self._script_dir),
                'echo "processing signal tree"',
                "python3 dummy_nanoAOD_plotter.py $1 $2",
                "ls -lah",
                'echo "copying outputfiles"',
                "cp *.root {}".format(self.workdir),
                "cp mass_parser_jsons/*.json {}".format(json_dir),
            ]
        else:
            self._wrapper_lines += [
                "ls -lah",
                "cp {}/dummy_nanoAOD_plotter_bkg.py .".format(self._script_dir),
                'echo "processing background tree"',
                "python3 dummy_nanoAOD_plotter_bkg.py $1 $2",
                "ls -lah",
                'echo "copying outputfiles"',
                "cp *.root {}".format(self.workdir),
            ]

    def print_status(self):
        self.check_jobs()
        status_blocks = [
            "({}/{}) {}".format(self.status.count(status_id), self.njobs, status_text)
            for status_id, status_text in self.status_lut.items()
            if self.status.count(status_id) > 0
        ]
        if self._merged:
            status_blocks.append(" MERGED")
        print(self._name, ", ".join(status_blocks))

    def dump_info(self):
        with open(self._manager_info_path, "w") as manager_info_file:
            manager_info_file.write(json.dumps(
                {
                    "files": self.files,
                    "ids": self.jobids,
                    "status": self.status,
                    "resubmit_counter": self._resubmit_counter,
                    "merged": self._merged,
                },
                sort_keys=False,
                indent=2
            ))

    def load_info(self):
        manager_info = json.load(open(self._manager_info_path, "r"))
        self.files = manager_info["files"]
        self.njobs = len(self.files)
        self.jobids = manager_info["ids"]
        self.status = manager_info["status"]
        self._resubmit_counter = manager_info["resubmit_counter"]
        self._merged = manager_info["merged"]

    def check_jobs(self):
        for jobid in self.jobids:
            outputfile_path = "{}/{}_{}.root".format(self.workdir, self._output_file_prefix, jobid)
            if os.path.isfile(outputfile_path):
                if test_root_file(outputfile_path):
                    self.status[jobid] = 2
            final_file = "{}/{}.root".format(self.workdir, self._name)
            if os.path.isfile(final_file):
                if test_root_file(final_file):
                    self._merged = True
        self.dump_info()

    def submit_jobs(self, debug, jobids=[], batchname_suffix=""):
        submitter = JobSubmitter(self.workdir)

        if len(jobids) == 0:
            # submitting all jobs
            jobs = zip(self.jobids, self.files)
        else:
            # only submitting some jobs
            jobs = [(jobid, self.files[jobid]) for jobid in jobids]

        if not debug:
            for jobid in jobids:
                self.status[jobid] = 1

        submitter.add_job_array(
            name="{}_TreeMaker{}".format(self._name, batchname_suffix),
            subworkdir=".",
            wrapper_lines=self._wrapper_lines,
            arguments=["$(InputFile) $(MyJobIndex)"],
            queue="MyJobIndex InputFile from ({}\n)".format(
                "\n".join(" ".join(map(str, job)) for job in jobs)
            )
        )
        submitter.submit_jobs(debug=debug)
        self.dump_info()

    def resubmit_jobs(self, debug):
        self.check_jobs()
        unfinished_jobs = [jobid for jobid in self.jobids if self.status[jobid] != 2]
        if len(unfinished_jobs) == 0:
            return
        self.submit_jobs(debug, unfinished_jobs, "_resubmit_{}".format(self._resubmit_counter))
        if not debug:
            self._resubmit_counter += 1

    def merge_jobs(self):
        if sum(self.status) == self.njobs*2:
            print("merging files for {}".format(self._name))
            os.system(
                "hadd -f {}/{}.root {}".format(
                    self.workdir,
                    self._name,
                    " ".join(
                        "{}/{}_{}.root".format(self.workdir, self._output_file_prefix, jobid) for jobid in self.jobids
                    ),
                )
            )


if __name__ == "__main__":
    import argparse

    script_dir = "/afs/desy.de/user/s/schmelch/GenLevelStudies"
    base_dir = "/nfs/dust/cms/user/schmelch/customTrees"

    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", help="Do not submit any jobs but create all necessary files.")
    parser.add_argument("--hadd", action="store_true", help="Hadd all trees.")
    parser.add_argument("--submit", action="store_true", help="Submit jobs to write trees.")
    parser.add_argument("--resubmit", action="store_true", help="Resubmit jobs that failed previously.")

    args = parser.parse_args()
    signal_manager = SampleManager(
        name="signal",
        file_pattern=(
            "/pnfs/desy.de/cms/tier2/store/user/anmehta/resSrch/XToYYprime_YToQQ_YprimeToQQ_narrow_TuneCP5"
            "_PSWeights_13TeV-madgraph-pythia8_RunIIAutumn18/*.root"
        ),
        base_dir=base_dir,
        script_dir=script_dir,
        signal=True,
    )

    background_manager = SampleManager(
        name="background",
        file_pattern=(
            "/pnfs/desy.de/cms/tier2/store/user/anmehta/resSrch/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_"
            "RunIIAutumn18/*.root"
        ),
        base_dir=base_dir,
        script_dir=script_dir,
        signal=False,
    )
    if args.submit:
        signal_manager.submit_jobs(debug=args.debug)
        background_manager.submit_jobs(debug=args.debug)
    if args.resubmit:
        signal_manager.resubmit_jobs(debug=args.debug)
        background_manager.resubmit_jobs(debug=args.debug)

    if args.hadd:
        signal_manager.merge_jobs()
        background_manager.merge_jobs()

    signal_manager.print_status()
    background_manager.print_status()
