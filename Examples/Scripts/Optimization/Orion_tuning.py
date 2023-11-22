#!/usr/bin/env python3
import sys

import sys
import os
import yaml
import pprint
import time
import datetime
import warnings

import logging
import uproot

import pathlib
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import random
import subprocess
import multiprocessing
import numpy as np
import json
import array
import sys
import argparse
import pandas as pd

from typing import Optional, Union
from pathlib import Path

from orion.client import build_experiment

srcDir = Path(__file__).resolve().parent


def run_ckf(params, names, outDir):

    if len(params) != len(names):
        raise Exception("Length of Params must equal names")

    ckf_script = srcDir / "ckf.py"
    nevts = "--nEvents=1"
    indir = "--indir=" + str(srcDir)
    outdir = "--output=" + str(outDir)

    ret = ["python"]
    ret.append(ckf_script)
    ret.append(nevts)
    ret.append(indir)
    ret.append(outdir)

    i = 0
    for param in params:
        arg = "--sf_" + names[i] + "=" + str(param)
        ret.append(arg)
        i += 1

    # Run CKF for the given parameters
    subprocess.call(ret)


class Objective:
    def __init__(self, k_dup, k_time):
        self.res = {
            "eff": [],
            "fakerate": [],
            "duplicaterate": [],
            "runtime": [],
        }

        self.k_dup = k_dup
        self.k_time = k_time

    def __call__(
        self,
        maxSeedsPerSpM,
        cotThetaMax,
        sigmaScattering,
        radLengthPerSeed,
        impactMax,
        maxPtScattering,
        deltaRMin,
        deltaRMax,
        ckf_perf=True,
    ):

        params = [
            maxSeedsPerSpM,
            cotThetaMax,
            sigmaScattering,
            radLengthPerSeed,
            impactMax,
            maxPtScattering,
            deltaRMin,
            deltaRMax,
        ]
        keys = [
            "maxSeedsPerSpM",
            "cotThetaMax",
            "sigmaScattering",
            "radLengthPerSeed",
            "impactMax",
            "maxPtScattering",
            "deltaRMin",
            "deltaRMax",
        ]

        if ckf_perf:
            outDirName = "Output_CKF"
            outputfile = srcDir / outDirName / "performance_ckf.root"
            effContName = "particles"
            contName = "tracks"
        else:
            outDirName = "Output_Seeding"
            outputfile = srcDir / outDirName / "performance_seeding.root"
            effContName = "seeds"
            contName = "seeds"

        outputDir = Path(srcDir / outDirName)
        outputDir.mkdir(exist_ok=True)
        run_ckf(params, keys, outputDir)
        rootFile = uproot.open(outputfile)
        self.res["eff"].append(rootFile["eff_" + effContName].member("fElements")[0])
        self.res["fakerate"].append(
            rootFile["fakerate_" + contName].member("fElements")[0]
        )
        self.res["duplicaterate"].append(
            rootFile["duplicaterate_" + contName].member("fElements")[0]
        )

        timingfile = srcDir / outDirName / "timing.tsv"
        timing = pd.read_csv(timingfile, sep="\t")

        if ckf_perf:
            time_ckf = float(
                timing[
                    timing["identifier"].str.match("Algorithm:TrackFindingAlgorithm")
                ]["time_perevent_s"]
            )
        else:
            time_ckf = 0

        time_seeding = float(
            timing[timing["identifier"].str.match("Algorithm:SeedingAlgorithm")][
                "time_perevent_s"
            ]
        )
        self.res["runtime"].append(time_ckf + time_seeding)

        efficiency = self.res["eff"][-1]
        penalty = (
            self.res["fakerate"][-1]
            + self.res["duplicaterate"][-1] / self.k_dup
            + self.res["runtime"][-1] / self.k_time
        )

        return [
            {"name": "objective", "type": "objective", "value": -(efficiency - penalty)}
        ]


def main():

    k_dup = 5
    k_time = 5

    # Initializing the objective (score) function
    objective = Objective(k_dup, k_time)

    # Defining the parameter space
    space = {
        "maxSeedsPerSpM": "uniform(0,10,discrete=True)",
        "cotThetaMax": "uniform(5.0,10.0)",
        "sigmaScattering": "uniform(0.2,50.0)",
        "radLengthPerSeed": "uniform(.001,0.1)",
        "impactMax": "uniform(0.1,25.0)",
        "maxPtScattering": "uniform(1.0, 50.0)",
        "deltaRMin": "uniform(0.25, 30.0)",
        "deltaRMax": "uniform(50.0,300.0)",
    }

    # Remove storage file if already exists (conflicts with present run if not deleted)
    if os.path.exists("./db.pkl"):
        os.remove("./db.pkl")

    # location to store metadata
    storage = {
        "type": "legacy",
        "database": {
            "type": "pickleddb",
            "host": "./db.pkl",
        },
    }

    # Build new orion experiment
    experiment = build_experiment(
        "orion_new",
        space=space,
        storage=storage,
    )

    # Start Optimization
    experiment.workon(objective, max_trials=3)

    outputDir = Path("OrionResults")
    outputDir.mkdir(exist_ok=True)

    # fetching trials in a dataframe
    df = experiment.to_pandas()
    df.to_csv(outputDir / "results.txt")

    # Getting the best parameters
    df_imp = df[
        [
            "objective",
            "maxSeedsPerSpM",
            "cotThetaMax",
            "sigmaScattering",
            "radLengthPerSeed",
            "impactMax",
            "maxPtScattering",
            "deltaRMin",
            "deltaRMax",
        ]
    ]
    df_obj = df["objective"]
    min_obj = df_obj.min()
    df_final = df_imp[df_imp["objective"] == min_obj]
    print("Best Score = %s" % (df_final["objective"]))
    print("maxSeedsPerSpM = %s" % (df_final["maxSeedsPerSpM"]))
    print("cotThetaMax = %s" % (df_final["cotThetaMax"]))
    print("sigmaScattering = %s" % (df_final["sigmaScattering"]))
    print("radLengthPerSeed = %s" % (df_final["radLengthPerSeed"]))
    print("impactMax = %s" % (df_final["impactMax"]))
    print("maxPtScattering = %s" % (df_final["maxPtScattering"]))
    print("deltaRMin = %s" % (df_final["deltaRMin"]))
    print("deltaRMax = %s" % (df_final["deltaRMax"]))


if __name__ == "__main__":
    main()
