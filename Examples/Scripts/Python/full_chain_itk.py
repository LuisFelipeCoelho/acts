#!/usr/bin/env python3
import pathlib, acts, acts.examples, acts.examples.itk
from acts.examples import RootParticleReader, CsvSpacePointReader
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
    CKFPerformanceConfig,
    TrackSelectorRanges,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
)

ttbar_pu200 = False
u = acts.UnitConstants
geo_dir = pathlib.Path("/Users/luiscoelho/lcoelho/acts/acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

events=200
s = acts.examples.Sequencer(events=events, numThreads=1, outputDir=str(outputDir))

#if not ttbar_pu200:
#    addParticleGun(
#        s,
#        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
#        EtaConfig(-4.0, 4.0, uniform=True),
#        ParticleConfig(200, acts.PdgParticle.eMuon, randomizeCharge=True),
#        rnd=rnd,
#    )
#else:
#    addPythia8(
#        s,
#        hardProcess=["Top:qqbar2ttbar=on"],
#        npileup=0,
#        vtxGen=acts.examples.GaussianVertexGenerator(
#            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
#            mean=acts.Vector4(0, 0, 0, 0),
#        ),
#        rnd=rnd,
#        outputDirRoot=outputDir,
#    )

s.addReader(
		RootParticleReader(
				level=acts.logging.INFO,
				filePath="/Users/luiscoelho/lcoelho/acts2/samples/ttbar_mu200_500evts/pythia8_particles.root",
				particleCollection="particles_input",
				orderedEvents=False,
		)
)

# https://github.com/acts-project/acts/pull/1714 https://github.com/acts-project/acts/pull/1757

addFatras(
    s,
    trackingGeometry,
    field,
    ParticleSelectorConfig(eta=(-4.0, 4.0), pt=(150 * u.MeV, None), removeNeutral=True)
    if ttbar_pu200
    else ParticleSelectorConfig(),
    outputDirRoot=outputDir,
    rnd=rnd,
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
	seedingAlgorithm=SeedingAlgorithm.Orthogonal,
	*acts.examples.itk.itkSeedingAlgConfig("PixelSpacePoints"),
	geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
	outputDirRoot=outputDir,
)

s.run()
