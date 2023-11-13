#!/usr/bin/env python3
import sys

import sys
import os
import yaml
import pprint
import time
import datetime
import warnings

import optuna
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

from optuna.visualization import plot_contour
from optuna.visualization import plot_edf
from optuna.visualization import plot_intermediate_values
from optuna.visualization import plot_optimization_history
from optuna.visualization import plot_parallel_coordinate
from optuna.visualization import plot_param_importances
from optuna.visualization import plot_slice

srcDir = Path(__file__).resolve().parent


def run_ckf(params, names, outDir):

    if len(params) != len(names):
        raise Exception("Length of Params must equal names")

    ckf_script = srcDir / "ckf.py"
    nevts = "--nEvents=10"
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

    def __call__(self, trial):
        params = []

        deltaPhiMax = trial.suggest_float("deltaPhiMax1", 0.030, 0.050)
        params.append(deltaPhiMax)
        cotThetaMax = trial.suggest_float("cotThetaMax", 20.0, 30.0)
        params.append(cotThetaMax)
        collisionRegionMin = trial.suggest_float("collisionRegionMin", -280., -220.)
        params.append(collisionRegionMin)
        collisionRegionMax = trial.suggest_float("collisionRegionMax", 220., 280.)
        params.append(collisionRegionMax)
        deltaRMinTopSP = trial.suggest_float("deltaRMinTopSP", 20, 28.)
        params.append(deltaRMinTopSP)
        deltaRMaxTopSP = trial.suggest_float("deltaRMaxTopSP", 100., 200.)
        params.append(deltaRMaxTopSP)
        deltaRMinBottomSP = trial.suggest_float("deltaRMinBottomSP", 20, 30.)
        params.append(deltaRMinBottomSP)
        deltaRMaxBottomSP = trial.suggest_float("deltaRMaxBottomSP", 70., 120.)
        params.append(deltaRMaxBottomSP)
        keys = [
            "deltaPhiMax1",
            "cotThetaMax",
            "collisionRegionMin",
            "collisionRegionMax",
            "deltaRMinTopSP",
            "deltaRMaxTopSP",
            "deltaRMinBottomSP",
            "deltaRMaxBottomSP",
        ]

        outputDir = Path(srcDir / "Output_CKF")
        outputfile = srcDir / "Output_CKF/performance_ckf.root"
        outputDir.mkdir(exist_ok=True)
        run_ckf(params, keys, outputDir)
        rootFile = uproot.open(outputfile)
        self.res["eff"].append(rootFile["eff_particles"].member("fElements")[0])
        self.res["fakerate"].append(rootFile["fakerate_tracks"].member("fElements")[0])
        self.res["duplicaterate"].append(
            rootFile["duplicaterate_tracks"].member("fElements")[0]
        )

        timingfile = srcDir / "Output_CKF/timing.tsv"
        timing = pd.read_csv(timingfile, sep="\t")
        time_ckf = float(
            timing[timing["identifier"].str.match("Algorithm:TrackFindingAlgorithm")][
                "time_perevent_s"
            ]
        )
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
            + 2 * self.res["runtime"][-1] / self.k_time
        )

        return efficiency - penalty


def main():

    k_dup = 5
    k_time = 5

    # Initializing the objective (score) function
    objective = Objective(k_dup, k_time)


    start_values = {
        "deltaPhiMax1": 0.025,
        "cotThetaMax": 27.2899,
        "collisionRegionMin": -200,
        "collisionRegionMax": 200,
				"deltaRMinTopSP": 6,
        "deltaRMaxTopSP": 280,
        "deltaRMinBottomSP": 6,
        "deltaRMaxBottomSP": 150
    }


    # Optuna logger
    optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))
    study_name = "test_study"
    storage_name = "sqlite:///{}.db".format(study_name)

    # creating a new optuna study
    study = optuna.create_study(
        study_name=study_name,
        storage="sqlite:///{}.db".format(study_name),
        direction="maximize",
        load_if_exists=True,
    )

    #study.enqueue_trial(start_values)
    # Start Optimization
    study.optimize(objective, n_trials=30)

    # Printout the best trial values
    print("Best Trial until now", flush=True)
    for key, value in study.best_trial.params.items():
        print(f"    {key}: {value}", flush=True)

    outputDir = Path("OptunaResults")
    outputDir.mkdir(exist_ok=True)

    with open(outputDir / "results.json", "w") as fp:
        json.dump(study.best_params, fp)


if __name__ == "__main__":
    main()
                                                   
