#!/usr/bin/env python3
import pathlib, acts, acts.examples, acts.examples.itk
from acts.examples import CsvSpacePointReader
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
    SeedFinderConfigArg,
		SeedFinderOptionsArg,
		SeedFilterConfigArg,
		SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
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
    addStandardSeeding
)

ttbar_pu200 = True
u = acts.UnitConstants
geo_dir = pathlib.Path("/Users/luiscoelho/lcoelho/acts/acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"
#acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

#detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(events=1, numThreads=1, outputDir=str(outputDir))

#s.addReader(
#		RootParticleReader(
#				level=acts.logging.INFO,
#				filePath="/Users/luiscoelho/lcoelho/acts2/samples/ttbar_mu200_500evts/pythia8_particles.root",
#				particleCollection="particles_input",
#				orderedEvents=False,
#		)
#)

pixel = True

if pixel:
	inputSpacePointsType = "PixelSpacePoints"
	inputCollection = "pixel"
	extendCollection = False
else:
	inputSpacePointsType = "StripSpacePoints"
	inputCollection = "strip"
	extendCollection = True

# Read input space points from input csv files
evReader = CsvSpacePointReader(
		level=acts.logging.INFO,
		inputStem="spacepoints",
		inputCollection=inputCollection,
		inputDir="/Users/luiscoelho/lcoelho/acts/data/CsvSpacePointsOutput_ttbar",
		outputSpacePoints=inputSpacePointsType,
		extendCollection=extendCollection,
)

s.addReader(evReader)

#addSeeding(
#    s,
#    trackingGeometry,
#    field,
#    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None))
#    if ttbar_pu200
#    else TruthSeedRanges(),
#    seedingAlgorithm=SeedingAlgorithm.Default,
#    *acts.examples.itk.itkSeedingAlgConfig(inputSpacePointsType),
#    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
#    outputDirRoot=outputDir,
#)

# run seeding
addStandardSeeding(
	s,
	evReader.config.outputSpacePoints,
	*acts.examples.itk.itkSeedingAlgConfig(inputSpacePointsType),
)

s.run()
