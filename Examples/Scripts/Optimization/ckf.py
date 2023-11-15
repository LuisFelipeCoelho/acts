#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union
from collections import namedtuple
import argparse
import sys
import os

from acts.examples import Sequencer, GenericDetector, RootParticleReader

import acts

from acts import UnitConstants as u

import pathlib, acts, acts.examples, acts.examples.itk
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    TruthSeedRanges,
    addCKFTracks,
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
)

ttbar_pu200 = True

def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(description="Command line arguments for CKF")
    parser.add_argument(
        "-i",
        "--indir",
        dest="indir",
        help="Directory with input root files",
        default="./",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="outdir",
        help="Output directory for new ntuples",
        default="./",
    )
    parser.add_argument(
        "-n", "--nEvents", dest="nEvts", help="Number of events to run over", default=1
    )
    parser.add_argument(
        "--sf_maxSeedsPerSpM",
        dest="sf_maxSeedsPerSpM",
        help="Number of compatible seeds considered for middle seed",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--sf_cotThetaMax",
        dest="sf_cotThetaMax",
        help="cot of maximum theta angle",
        type=float,
        default=7.40627,
    )
    parser.add_argument(
        "--sf_sigmaScattering",
        dest="sf_sigmaScattering",
        help="How many sigmas of scattering to include in seeds",
        type=float,
        default=5,
    )
    parser.add_argument(
        "--sf_radLengthPerSeed",
        dest="sf_radLengthPerSeed",
        help="Average Radiation Length",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--sf_impactMax",
        dest="sf_impactMax",
        help="max impact parameter in mm",
        type=float,
        default=3.0,
    )
    parser.add_argument(
        "--sf_maxPtScattering",
        dest="sf_maxPtScattering",
        help="maximum Pt for scattering cut in GeV",
        type=float,
        default=10.0,
    )
    parser.add_argument(
        "--sf_deltaRMinTopSP",
        dest="sf_deltaRMinTopSP",
        help="minimum value for deltaR separation in mm",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--sf_deltaRMaxTopSP",
        dest="sf_deltaRMaxTopSP",
        help="maximum value for deltaR separation in mm",
        type=float,
        default=60.0,
    )
    parser.add_argument(
        "--sf_deltaRMinBottomSP",
        dest="sf_deltaRMinBottomSP",
        help="minimum value for deltaR separation in mm",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--sf_deltaRMaxBottomSP",
        dest="sf_deltaRMaxBottomSP",
        help="maximum value for deltaR separation in mm",
        type=float,
        default=60.0,
    )
    parser.add_argument(
        "--sf_deltaPhiMax1",
        dest="sf_deltaPhiMax1",
        help="maximum value for deltaPhiMax separation for the Orthogonal Seeding",
        type=float,
        default=0.025,
    )
    parser.add_argument(
        "--sf_collisionRegionMin",
        dest="sf_collisionRegionMin",
        help="",
        type=float,
        default=-300,
    )
    parser.add_argument(
        "--sf_collisionRegionMax",
        dest="sf_collisionRegionMax",
        help="",
        type=float,
        default=300,
    )
    parser.add_argument(
        "--sf_deltaZMax",
        dest="sf_deltaZMax",
        help="",
        type=float,
        default=300,
    )

    return parser


def runCKFTracks(
    trackingGeometry,
    decorators,
    geometrySelection: Path,
    digiConfigFile: Path,
    field,
    outputDir: Path,
    NumEvents=1,
    truthSmearedSeeded=False,
    truthEstimatedSeeded=False,
    outputCsv=True,
    inputParticlePath: Optional[Path] = None,
    s=None,
    MaxSeedsPerSpM=1,
    CotThetaMax=7.40627,
    SigmaScattering=5,
    RadLengthPerSeed=0.1,
    ImpactMax=3.0,
    MaxPtScattering=10.0,
    DeltaPhiMax=0.025,
    CollisionRegionMin=-200,
    CollisionRegionMax=200,
    DeltaRMinTopSP=1.0,
    DeltaRMaxTopSP=60.0,
    DeltaRMinBottomSP=1.0,
    DeltaRMaxBottomSP=60.0,
    DeltaZMax=1500,
):
    from acts.examples.simulation import (
        addParticleGun,
        EtaConfig,
        PhiConfig,
        ParticleConfig,
        addFatras,
        addDigitization,
    )

    from acts.examples.reconstruction import (
        addSeeding,
        TruthSeedRanges,
        ParticleSmearingSigmas,
        SeedFinderConfigArg,
        SeedFinderOptionsArg,
        SeedingAlgorithm,
        TruthEstimatedSeedingAlgorithmConfigArg,
        addCKFTracks,
    )

    s = s or acts.examples.Sequencer(
        events=int(NumEvents),
        numThreads=-1,
        logLevel=acts.logging.INFO,
        outputDir=outputDir,
    )
    for d in decorators:
        s.addContextDecorator(d)
    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    if inputParticlePath is None:
        if not ttbar_pu200:
          addParticleGun(
            s,
            EtaConfig(-2.0, 2.0),
            ParticleConfig(10, acts.PdgParticle.eMuon, True),
            PhiConfig(0.0, 360.0 * u.degree),
            multiplicity=2,
            rnd=rnd,
          )
        else:
          addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=200,
            vtxGen=acts.examples.GaussianVertexGenerator(
              stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
              mean=acts.Vector4(0, 0, 0, 0),
            ),
            rnd=rnd,
            outputDirRoot=outputDir,
          )
    else:
        acts.logging.getLogger("CKFExample").info(
            "Reading particles from %s", inputParticlePath.resolve()
        )
        assert inputParticlePath.exists()
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection="particles_input",
                orderedEvents=False,
            )
        )

    addFatras(
        s,
        trackingGeometry,
        field,
        preSelectParticles=ParticleSelectorConfig(
          rho=(0.0 * u.mm, 28.0 * u.mm),
          absZ=(0.0 * u.mm, 1.0 * u.m),
          eta=(-4.0, 4.0),
          pt=(150 * u.MeV, None),
          removeNeutral=True,
        ),
        rnd=rnd,
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None)),
        ParticleSmearingSigmas(pRel=0.01),  # only used by SeedingAlgorithm.TruthSmeared
        *acts.examples.itk.itkSeedingAlgConfig(
          acts.examples.itk.InputSpacePointsType.PixelSpacePoints,
          highOccupancyConfig=False,
          #deltaPhiMax=DeltaPhiMax,
          #cotThetaMax=CotThetaMax,
          #collisionRegionMin=CollisionRegionMin,
          #collisionRegionMax=CollisionRegionMax,
          deltaRMinTopSP=DeltaRMinTopSP,
          deltaRMaxTopSP=DeltaRMaxTopSP,
          deltaRMinBottomSP=DeltaRMinBottomSP,
          deltaRMaxBottomSP=DeltaRMaxBottomSP,
          deltaZMax=DeltaZMax,
        ),
        TruthEstimatedSeedingAlgorithmConfigArg(deltaR=(10.0 * u.mm, None)),
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared
        if truthSmearedSeeded
        else SeedingAlgorithm.TruthEstimated
        if truthEstimatedSeeded
        else SeedingAlgorithm.Orthogonal,
        geoSelectionConfigFile=geometrySelection,
        outputDirRoot=outputDir,
        rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
    )

#    addCKFTracks(
#        s,
#        trackingGeometry,
#        field,
#        outputDirRoot=outputDir,
#        outputDirCsv=outputDir / "csv" if outputCsv else None,
#    )

    return s


if "__main__" == __name__:
    options = getArgumentParser().parse_args()

    Inputdir = options.indir
    Outputdir = options.outdir

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    geo_dir = pathlib.Path("/Users/luiscoelho/lcoelho/acts/acts-itk")
    detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)

    field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))

    inputParticlePath = Path(Inputdir) / "pythia8_particles.root"
    if not inputParticlePath.exists():
        inputParticlePath = None
    else:
        print("--> Using particles from ", inputParticlePath)

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=geo_dir / "itk-hgtd/geoSelection-ITk.json",
        digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
        outputCsv=True,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=False,
        inputParticlePath=inputParticlePath,
        outputDir=Outputdir,
        NumEvents=options.nEvts,
        MaxSeedsPerSpM=options.sf_maxSeedsPerSpM,
        CotThetaMax=options.sf_cotThetaMax,
        SigmaScattering=options.sf_sigmaScattering,
        RadLengthPerSeed=options.sf_radLengthPerSeed,
        ImpactMax=options.sf_impactMax,
        MaxPtScattering=options.sf_maxPtScattering,
        DeltaPhiMax=options.sf_deltaPhiMax1,
        CollisionRegionMin=options.sf_collisionRegionMin,
        CollisionRegionMax=options.sf_collisionRegionMax,
        DeltaRMinTopSP=options.sf_deltaRMinTopSP,
        DeltaRMaxTopSP=options.sf_deltaRMaxTopSP,
        DeltaRMinBottomSP=options.sf_deltaRMinBottomSP,
        DeltaRMaxBottomSP=options.sf_deltaRMaxBottomSP,
        DeltaZMax=options.sf_deltaZMax,
    ).run()
