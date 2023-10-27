#!/usr/bin/env python3
import pathlib, acts, acts.examples
from acts.examples.itk import buildITkGeometry

from acts.examples import RootParticleReader

u = acts.UnitConstants
geo_dir = pathlib.Path("/afs/cern.ch/user/l/lfaldaul/work/actsODD/acts-itk")
outputDir = pathlib.Path.cwd()

# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls
detector, trackingGeometry, decorators = buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

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
    TruthSeedRanges,
    SeedingAlgorithm,
    ParticleSmearingSigmas,
    addCKFTracks,
    CKFPerformanceConfig,
)

from acts.examples.itk import itkSeedingAlgConfig

s = acts.examples.Sequencer(events=100, numThreads=1)

s.addReader(
	RootParticleReader(
		level=acts.logging.INFO,
		filePath=str("/eos/home-l/lfaldaul/ACTS_v19Xv30/pythia8_particles.root"),
		particleCollection="particles_input",
		orderedEvents=False,
	)
)

s = addFatras(
    s,
    trackingGeometry,
    field,
    ParticleSelectorConfig(
        rho=(0.0 * u.mm, 28.0 * u.mm),
        absZ=(0.0 * u.mm, 1000.0 * u.mm),
        eta=(-4.0, 4.0),
        pt=(150 * u.MeV, None),
        removeNeutral=True,
    ),
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)
s = addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None)),
    # SeedingAlgorithm.TruthEstimated,
    # SeedingAlgorithm.TruthSmeared, ParticleSmearingSigmas(pRel=0.01), rnd=rnd,
    *itkSeedingAlgConfig("PixelSpacePoints"),
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
    outputDirRoot=outputDir,
)

s.run()
