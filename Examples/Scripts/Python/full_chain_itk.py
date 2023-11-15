#!/usr/bin/env python3
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

from acts.examples import RootParticleReader

ttbar_pu200 = False
u = acts.UnitConstants
geo_dir = pathlib.Path("/afs/cern.ch/user/l/lfaldaul/work/actsODD/acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(events=100, numThreads=1)

s.addReader(
	RootParticleReader(
		level=acts.logging.INFO,
		filePath=str("/eos/home-l/lfaldaul/ACTS_v19Xv30/pythia8_particles.root"),
		particleCollection="particles_input",
		orderedEvents=False,
	)
)

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    preSelectParticles=ParticleSelectorConfig(
        rho=(0.0 * u.mm, 28.0 * u.mm),
        absZ=(0.0 * u.mm, 1000.0 * u.mm),
        eta=(-4.0, 4.0),
        pt=(150 * u.MeV, None),
        removeNeutral=True,
    ),
    outputDirRoot=outputDir,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)

addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None))
    if ttbar_pu200
    else TruthSeedRanges(),
    seedingAlgorithm=SeedingAlgorithm.Default,
    *acts.examples.itk.itkSeedingAlgConfig(
        acts.examples.itk.InputSpacePointsType.PixelSpacePoints
    ),
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
    outputDirRoot=outputDir,
)

s.run()
