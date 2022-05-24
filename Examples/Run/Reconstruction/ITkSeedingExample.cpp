
// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointReader.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <limits>
#include <memory>

#include <boost/program_options.hpp>

#include "RecInput.hpp"

#include "Acts/Seeding/SeedConfirmationRange.hpp"

//#include "../../Tests/UnitTests/Core/Seeding/ATLASCuts.hpp"
#include "Acts/Seeding/SeedFilter.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);

  // Add specific options for this geometry
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  for (const auto& it : vm) {
    std::cout << it.first.c_str() << " ";
    auto& value = it.second.value();
    if (auto v = boost::any_cast<std::string>(&value))
      std::cout << *v;
    else if (auto a = boost::any_cast<float>(&value))
      std::cout << *a;
    else if (auto b = boost::any_cast<double>(&value))
      std::cout << *b;
    else if (auto c = boost::any_cast<int>(&value))
      std::cout << *c;
    else if (auto d = boost::any_cast<bool>(&value))
      std::cout << *d;
    else if (auto e = boost::any_cast<size_t>(&value))
      std::cout << *e;
    else
      std::cout << "wrong cast...";
    std::cout << std::endl;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // Setup the magnetic field
  Options::setupMagneticFieldServices(vm, sequencer);
  auto magneticField = Options::readMagneticField(vm);

  // ====== Pixel SP =======

  // read extend SP parameters
  bool extendCollection = false;  // pixel: false, strip: true

  // Read the space points and build the container
//  auto spReaderCfg =
//      setupSpacePointReading(vm, sequencer, "pixel", extendCollection);

	// Read truth hits from CSV files
	auto spReaderCfg = Options::readCsvSpacePointReaderConfig(vm);
	spReaderCfg.inputStem = "spacepoints";
	spReaderCfg.inputCollection = "pixel";
	// read extend SP parameters
	spReaderCfg.extendCollection = extendCollection;  // pixel: false, strip: true

	sequencer.addReader(
						std::make_shared<CsvSpacePointReader>(spReaderCfg, logLevel));

  //	spReaderCfg.outputSpacePoints = {"PixelSpacePoints",
  //"StripSpacePoints"};//, "OverlapSpacePoints"};
  spReaderCfg.outputSpacePoints = "PixelSpacePoints";

  // Seeding addAlgorithm
  SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {"PixelSpacePoints"};
  seedingCfg.outputSeeds = "PixelSeeds";
  seedingCfg.outputProtoTracks = "prototracks";

  //  seedingCfg.gridConfig.rMax = 1200._mm;
  seedingCfg.gridConfig.rMax = 320._mm;  // pixel: 320 mm, strip: 1000 mm
  //  seedingCfg.seedFinderConfig.rMax = seedingCfg.gridConfig.rMax;
  seedingCfg.seedFinderConfig.rMax = 320._mm;  // pixel: 320 mm, strip: 1200 mm

  seedingCfg.seedFilterConfig.deltaRMin = 20._mm;  // pixel: 20 mm, strip: 20._mm
  seedingCfg.seedFinderConfig.deltaRMin = seedingCfg.seedFilterConfig.deltaRMin;

  seedingCfg.seedFinderConfig.deltaRMinTopSP = 6._mm;  // pixel: 6 mm, strip: 20 mm
  seedingCfg.seedFinderConfig.deltaRMinBottomSP = 6._mm;  // pixel: 6 mm, strip: 20 mm

  //  seedingCfg.gridConfig.deltaRMax = 1000._mm;
  seedingCfg.gridConfig.deltaRMax = 280._mm;  // pixel: 280 mm, strip: 600 mm
  //  seedingCfg.seedFinderConfig.deltaRMax = seedingCfg.gridConfig.deltaRMax;
  seedingCfg.seedFinderConfig.deltaRMax = 280._mm;  // pixel: 280 mm, strip: 1000 mm

  // rMin = rMax - deltaRMax => pixel: 40 mm, strip: 600 mm

  seedingCfg.seedFinderConfig.deltaRMaxTopSP = 280._mm;  // pixel: 280 mm, strip:
  seedingCfg.seedFinderConfig.deltaRMaxBottomSP = 120._mm;  // pixel: 120 mm, strip: 300 mm

  seedingCfg.seedFinderConfig.collisionRegionMin = -200._mm;  // pixel: -200 mm, strip: -200 mm
  seedingCfg.seedFinderConfig.collisionRegionMax = 200._mm;  // pixel: 200 mm, strip: 200 mm

  seedingCfg.gridConfig.zMin = -3000._mm;
  seedingCfg.gridConfig.zMax = 3000._mm;
  seedingCfg.seedFinderConfig.zMin = seedingCfg.gridConfig.zMin;
  seedingCfg.seedFinderConfig.zMax = seedingCfg.gridConfig.zMax;

  seedingCfg.seedFilterConfig.maxSeedsPerSpM = 4;
  seedingCfg.seedFinderConfig.maxSeedsPerSpM =
      seedingCfg.seedFilterConfig.maxSeedsPerSpM;

  // cut on the compatibility between interaction point and SPs
  seedingCfg.seedFinderConfig.interactionPointCut =
      true;  // pixel: true, strip: false

  seedingCfg.gridConfig.cotThetaMax = 27.2899;  // pixel: 27.2899 (4.0 eta) ---> 27.2899 =
                // 1/np.tan(2*np.arctan(np.exp(-4)))
  seedingCfg.seedFinderConfig.cotThetaMax = seedingCfg.gridConfig.cotThetaMax;

  // number of standard deviations of Coulomb scattering angle that should be
  // considered
  seedingCfg.seedFinderConfig.sigmaScattering = 2;
  // radiation length in Highland equation
  seedingCfg.seedFinderConfig.radLengthPerSeed = 0.09804522341059585;

  // use arithmetic average in the calculation of the squared error on the
  // difference in tan(theta)
  seedingCfg.seedFinderConfig.arithmeticAverageCotTheta = false;  // pixel: false (uses geometric average), strip: true (uses
              // arithmetic average)

  seedingCfg.seedFinderConfig.deltaZMax = 10000000000000000._mm;

  seedingCfg.gridConfig.minPt = 900._MeV;
  seedingCfg.seedFinderConfig.minPt = seedingCfg.gridConfig.minPt;

  seedingCfg.gridConfig.bFieldInZ = 2_T;
  seedingCfg.seedFinderConfig.bFieldInZ = 1.997244192_T;

  seedingCfg.seedFinderConfig.beamPos = {0_mm, 0_mm};

  seedingCfg.gridConfig.impactMax = 2._mm;  // pixel: 2 mm, strip: 20 mm
  seedingCfg.seedFinderConfig.impactMax = seedingCfg.gridConfig.impactMax;

  seedingCfg.seedFinderConfig.maxPtScattering = std::numeric_limits<double>::infinity();

  // enable non equidistant binning in z, in case the binning is not defined the
  // edges are evaluated automatically using equidistant binning
  seedingCfg.gridConfig.zBinEdges = {-3000., -2500., -1400., -925.,
                                     -450.,  -250.,  250.,   450.,
                                     925.,   1400.,  2500.,  3000.};
  seedingCfg.seedFinderConfig.zBinEdges = seedingCfg.gridConfig.zBinEdges;
  // enable custom z looping when searching for SPs, must contain numbers from 1
  // to the total number of bin in zBinEdges
  seedingCfg.seedFinderConfig.zBinsCustomLooping = {1, 2, 3, 4, 11, 10,
                                                    9, 8, 6, 5, 7};
  // enable cotTheta sorting in SeedFinder
	seedingCfg.seedFinderConfig.skipPreviousTopSP = true;
	seedingCfg.seedFilterConfig.useDeltaRorTopRadius = true;
	seedingCfg.seedFilterConfig.maxSeedsPerSpMConf = 5;
	seedingCfg.seedFilterConfig.maxQualitySeedsPerSpMConf = 5;

  // LUT for building neighbors
  /* Guide for the following:
   * z == 6: central z region, |z|<250mm
   * [-3000, -2500., -1400., -925., -450., -250.,  250.,  450.,  925.,  1400.,
   * 2500.,  3000] 1       2       3      4      5      6      7      8      9
   * 10      11        z bin index
   * -------------------------------------------------------------------------------------------->
   * Z[mm] Z=-3000                                  IP,Z=0 Z=+3000
   */
  // allows to specify the number of neighbors desired for each bin
  // {-1,1} means one neighbor on the left and one on the right
  // vector containing the map of z bins for the top SP
  // for ITk pixel and strip: zBinNeighborsTop =
  // {{0,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,1},{0,1},{0,1},{0,1},{0,1},{0,0}};
  seedingCfg.zBinNeighborsTop = {{0, 0},  {-1, 0}, {-1, 0}, {-1, 0},
                                 {-1, 0}, {-1, 1}, {0, 1},  {0, 1},
                                 {0, 1},  {0, 1},  {0, 0}};
  // vector containing the map of z bins for the bottom SP
  // for ITk pixel: zBinNeighborsBottom =
  // {{0,1},{0,1},{0,1},{0,1},{0,1},{0,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0}};
  // for ITk strip: zBinNeighborsBottom =
  // {{0,1},{0,1},{0,1},{0,2},{0,1},{0,0},{-1,0},{-2,0},{-1,0},{-1,0},{-1,0}};
  seedingCfg.zBinNeighborsBottom = {{0, 1},  {0, 1},  {0, 1},  {0, 1},
                                    {0, 1},  {0, 0},  {-1, 0}, {-1, 0},
                                    {-1, 0}, {-1, 0}, {-1, 0}};
  // numPhiNeighbors for the Grid means "how many phiBin neighbors (plus the
  // current bin) should cover the full deflection of a minimum pT particle"
  // numPhiNeighbors for the BinFinder sets how many neighboring phi bins at
  // each side of the current bin are returned
  seedingCfg.gridConfig.phiBinDeflectionCoverage = 3;

  // radial range for middle SP -
  seedingCfg.seedFinderConfig.rRangeMiddleSP = {
      {40., 90.},  {40., 200.}, {46., 200.}, {46., 200.},
      {46., 250.}, {46., 250.}, {46., 250.}, {46., 200.},
      {46., 200.}, {40., 200.}, {40., 90.}};
  seedingCfg.seedFinderConfig.useVariableMiddleSPRange = true;
  seedingCfg.seedFinderConfig.deltaRMiddleMinSPRange = 10.;  // pixel: 10 mm, strip: 30 mm
  seedingCfg.seedFinderConfig.deltaRMiddleMaxSPRange = 10.;  // pixel: 10 mm, strip: 150 mm

  // ITk seed confirmation
  seedingCfg.seedFinderConfig.seedConfirmation = true;  // pixel: true, strip: true
  seedingCfg.seedFilterConfig.seedConfirmation = true;  // pixel: true, strip: false
  // contains parameters for central seed confirmation
  seedingCfg.seedFinderConfig.centralSeedConfirmationRange = Acts::SeedConfirmationRange(250., -250., 140., 1, 2);
  seedingCfg.seedFilterConfig.centralSeedConfirmationRange =
      seedingCfg.seedFinderConfig.centralSeedConfirmationRange;
  // contains parameters for forward seed confirmation
  seedingCfg.seedFinderConfig.forwardSeedConfirmationRange =
      Acts::SeedConfirmationRange(3000., -3000., 140., 1, 2);
  seedingCfg.seedFilterConfig.forwardSeedConfirmationRange =
      seedingCfg.seedFinderConfig.forwardSeedConfirmationRange;

  // parameters for the calculation of the weght of the seeds in the seed filter
  seedingCfg.seedFilterConfig.impactWeightFactor = 100.; // pixel: 100, strip: 1
  seedingCfg.seedFilterConfig.compatSeedWeight = 100.;
  // maximum number of seeds allowed after the filter
	seedingCfg.numPhiNeighbors = 1;
  seedingCfg.seedFilterConfig.compatSeedLimit = 3;  // pixel: 3, strip: 4

  // enable curvature sorting in SeedFilter
  seedingCfg.seedFilterConfig.curvatureSortingInFilter = true;

  // delete
//  seedingCfg.seedFinderConfig.inputCollectionTest = "pixel";
//  seedingCfg.seedFilterConfig.inputCollectionTest =
//      seedingCfg.seedFinderConfig.inputCollectionTest;

  sequencer.addAlgorithm(
      std::make_shared<SeedingAlgorithm>(seedingCfg, logLevel));

//  // ====== Strip SP =======
//
//  // read extend SP parameters
//  	bool extendCollection = true; // pixel: false, strip: true
//
////  	// Read the space points and build the container
////  	auto spReaderCfg = setupSpacePointReading(vm, sequencer, "strip",
////  extendCollection);
//
//	// Read truth hits from CSV files
//	auto spReaderCfg = Options::readCsvSpacePointReaderConfig(vm);
//	spReaderCfg.inputStem = "spacepoints";
//	spReaderCfg.inputCollection = "strip";
//	// read extend SP parameters
//	spReaderCfg.extendCollection = extendCollection;  // pixel: false, strip: true
//
//		sequencer.addReader(
//							std::make_shared<CsvSpacePointReader>(spReaderCfg, logLevel));
//
//  	//	spReaderCfg.outputSpacePoints = {"PixelSpacePoints",  "StripSpacePoints"};//, "OverlapSpacePoints"};
//	spReaderCfg.outputSpacePoints  = "StripSpacePoints";
//
//  	// Seeding addAlgorithm
//  	SeedingAlgorithm::Config seedingCfg;
//  	seedingCfg.inputSpacePoints = {"PixelSpacePoints"};
//  	seedingCfg.outputSeeds = "PixelSeeds";
//  	seedingCfg.outputProtoTracks = "prototracks";
//
//	seedingCfg.allowSeparateRMax = true;
//
//  	//  seedingCfg.gridConfig.rMax = 1200._mm;
//  	seedingCfg.gridConfig.rMax = 1000._mm; // pixel: 320 mm, strip: 1000 mm
//  	//  seedingCfg.seedFinderConfig.rMax = seedingCfg.gridConfig.rMax;
//  	seedingCfg.seedFinderConfig.rMax = 1200_mm; // pixel: 320 mm, strip: 1200 mm
//
//  	seedingCfg.seedFilterConfig.deltaRMin = 20._mm; // pixel: 20 mm, strip: 20._mm
//	seedingCfg.seedFinderConfig.deltaRMin =
//  seedingCfg.seedFilterConfig.deltaRMin;
//
//  	seedingCfg.seedFinderConfig.deltaRMinTopSP = 20._mm; // pixel: 6 mm, strip: 20 mm
//	seedingCfg.seedFinderConfig.deltaRMinBottomSP = 20._mm; // pixel: 6 mm, strip: 20 mm
//
//  	//  seedingCfg.gridConfig.deltaRMax = 1000._mm;
//  	seedingCfg.gridConfig.deltaRMax = 600._mm; // pixel: 280 mm, strip: 600  mm
//  	//  seedingCfg.seedFinderConfig.deltaRMax =  seedingCfg.gridConfig.deltaRMax;
//	seedingCfg.seedFinderConfig.deltaRMax =  600_mm; // pixel: 280 mm, strip: 1000 mm
//
//  	seedingCfg.seedFinderConfig.deltaRMaxTopSP = 300._mm; // pixel: 280 mm, strip: 300 mm
//	seedingCfg.seedFinderConfig.deltaRMaxBottomSP = 300._mm; //  pixel: 120 mm, strip: 300 mm
//
//  	seedingCfg.seedFinderConfig.collisionRegionMin = -200._mm; // pixel: 200  mm, strip: 200 mm
//	seedingCfg.seedFinderConfig.collisionRegionMax = 200._mm;  // pixel: 200 mm, strip: 200 mm
//
//  	seedingCfg.gridConfig.zMin = -3000._mm;
//  	seedingCfg.gridConfig.zMax = 3000._mm;
//  	seedingCfg.seedFinderConfig.zMin = seedingCfg.gridConfig.zMin;
//  	seedingCfg.seedFinderConfig.zMax = seedingCfg.gridConfig.zMax;
//
//  	seedingCfg.seedFilterConfig.maxSeedsPerSpM = 4;
//  	seedingCfg.seedFinderConfig.maxSeedsPerSpM =
//  	seedingCfg.seedFilterConfig.maxSeedsPerSpM;
//
//  	// cut on the compatibility between interaction point and SPs
//  	seedingCfg.seedFinderConfig.interactionPointCut = false; // pixel: true,  strip: false
//
//  	seedingCfg.gridConfig.cotThetaMax = 27.2899; // pixel: 27.2899 (4.0 eta)  ---> 27.2899 = 1/np.tan(2*np.arctan(np.exp(-4)))
//  	seedingCfg.seedFinderConfig.cotThetaMax =
//  seedingCfg.gridConfig.cotThetaMax;
//
//  	// number of standard deviations of Coulomb scattering angle that should  be considered
//	 seedingCfg.seedFinderConfig.sigmaScattering = 2;
//  	// radiation length in Highland equation
//  	seedingCfg.seedFinderConfig.radLengthPerSeed = 0.09804522341059585;
//
//  	// use arithmetic average in the calculation of the squared error on the difference in tan(theta)
//  	seedingCfg.seedFinderConfig.arithmeticAverageCotTheta = true; // pixel:  false (uses geometric average), strip: true (uses arithmetic average)
//
//	seedingCfg.seedFinderConfig.deltaZMax = 900._mm;
//
//  	seedingCfg.gridConfig.minPt = 900._MeV;
//  	seedingCfg.seedFinderConfig.minPt = seedingCfg.gridConfig.minPt;
//
//  	seedingCfg.gridConfig.bFieldInZ = 2_T;
//  	seedingCfg.seedFinderConfig.bFieldInZ = 1.997244311_T;
//
//  	seedingCfg.seedFinderConfig.beamPos = {0_mm, 0_mm};
//
//  	seedingCfg.gridConfig.impactMax = 20._mm; // pixel: 2 mm, strip: 20 mm
//  	seedingCfg.seedFinderConfig.impactMax = seedingCfg.gridConfig.impactMax;
//
//  	seedingCfg.seedFinderConfig.maxPtScattering =
//  std::numeric_limits<double>::infinity();
//
//  	// enable non equidistant binning in z, in case the binning is not defined the edges are evaluated automatically using equidistant binning
//  	seedingCfg.gridConfig.zBinEdges = {-3000., -2500., -1400., -925., -450.,
//  -250., 		250., 450., 925., 1400., 2500., 3000.};
//  	seedingCfg.seedFinderConfig.zBinEdges = seedingCfg.gridConfig.zBinEdges;
//  	// enable custom z looping when searching for SPs, must contain numbers from 1 to the total number of bin in zBinEdges
//  	seedingCfg.seedFinderConfig.zBinsCustomLooping = {6, 7, 5, 8, 4, 9, 3, 10, 2, 11, 1};
//  	// enable cotTheta sorting in SeedFinder
//
//	seedingCfg.seedFinderConfig.skipPreviousTopSP = false;
//
//  	// LUT for building neighbors
//  	/* Guide for the following:
//  	 * z == 6: central z region, |z|<250mm
//  	 * [-3000, -2500., -1400., -925., -450., -250.,  250.,  450.,  925.,
//  1400.,  2500.,  3000]
//  	 *       1       2       3      4      5      6      7      8      9 10
//  11        z bin index
//  	 *
//  -------------------------------------------------------------------------------------------->
//  Z[mm]
//  	 * Z=-3000                                  IP,Z=0 Z=+3000
//  	 */
//  	// allows to specify the number of neighbors desired for each bin
//  	// {-1,1} means one neighbor on the left and one on the right
//  	// vector containing the map of z bins for the top SP
//  	// for ITk pixel and strip: zBinNeighborsTop =  {{0,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,1},{0,1},{0,1},{0,1},{0,1},{0,0}};
//  	seedingCfg.zBinNeighborsTop =
//  {{0,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,1},{0,1},{0,1},{0,1},{0,1},{0,0}};
//  	// vector containing the map of z bins for the bottom SP
//  	// for ITk pixel: zBinNeighborsBottom =  {{0,1},{0,1},{0,1},{0,1},{0,1},{0,0},{-1,0},{-1,0},{-1,0},{-1,0},{-1,0}};
//  	// for ITk strip: zBinNeighborsBottom =  {{0,1},{0,1},{0,1},{0,2},{0,1},{0,0},{-1,0},{-2,0},{-1,0},{-1,0},{-1,0}};
//  	seedingCfg.zBinNeighborsBottom =
//  {{0,1},{0,1},{0,1},{0,2},{0,1},{0,0},{-1,0},{-2,0},{-1,0},{-1,0},{-1,0}};
//  	// numPhiNeighbors for the Grid means "how many phiBin neighbors (plus the current bin) should cover the full deflection of a minimum pT particle"
//  	// numPhiNeighbors for the BinFinder sets how many neighboring phi bins at each side of the current bin are returned
//  	seedingCfg.gridConfig.phiBinDeflectionCoverage = 3;
//	seedingCfg.numPhiNeighbors = 1;
//
//  	// radial range for middle SP -
//  	seedingCfg.seedFinderConfig.rRangeMiddleSP = {{40., 90.},{40.,
//  200.},{46., 200.},{46., 200.},{46., 250.},{46., 250.},{46., 250.},{46.,
//  200.},{46., 200.},{40., 200.},{40., 90.}};
//  	seedingCfg.seedFinderConfig.useVariableMiddleSPRange = true;
//  	seedingCfg.seedFinderConfig.deltaRMiddleMinSPRange = 30.; // pixel: 10  mm, strip: 30 mm
//	seedingCfg.seedFinderConfig.deltaRMiddleMaxSPRange = 150.;  // pixel: 10 mm, strip: 150 mm
//
//  	// ITk seed confirmation
//  	seedingCfg.seedFinderConfig.seedConfirmation = true; // pixel: true,  strip: true
//	seedingCfg.seedFilterConfig.seedConfirmation = false; // pixel:  true, strip: false
//  	// contains parameters for central seed confirmation
//  	seedingCfg.seedFinderConfig.centralSeedConfirmationRange =
//  Acts::SeedConfirmationRange(250., -250., 140., 1, 2);
//  	seedingCfg.seedFilterConfig.centralSeedConfirmationRange =
//  seedingCfg.seedFinderConfig.centralSeedConfirmationRange;
//  	// contains parameters for forward seed confirmation
//  	seedingCfg.seedFinderConfig.forwardSeedConfirmationRange =
//  Acts::SeedConfirmationRange(3000., -3000., 140., 1, 2);
//  	seedingCfg.seedFilterConfig.forwardSeedConfirmationRange =
//  seedingCfg.seedFinderConfig.forwardSeedConfirmationRange;
//
//	seedingCfg.seedFinderConfig.useDetailedDoubleMeasurementInfo = true;
//
//  	// parameters for the calculation of the weght of the seeds in the seed filter
//	seedingCfg.seedFilterConfig.impactWeightFactor = 1.; // pixel: 100, strip: 1
//	seedingCfg.seedFilterConfig.compatSeedWeight = 100.;
//  	// maximum number of seeds allowed after the filter
//  	seedingCfg.seedFilterConfig.compatSeedLimit = 4;
//	seedingCfg.seedFilterConfig.numSeedIncrement = 1;
//	seedingCfg.seedFilterConfig.seedWeightIncrement = 10100;
//
//		seedingCfg.seedFilterConfig.useDeltaRorTopRadius = false;
////		seedingCfg.seedFilterConfig.maxSeedsPerSpMConf = 5;
////		seedingCfg.seedFilterConfig.maxQualitySeedsPerSpMConf = 5;
//
//  	// enable curvature sorting in SeedFilter
//  	seedingCfg.seedFilterConfig.curvatureSortingInFilter = true;
//
////  	// delete
////  	seedingCfg.seedFinderConfig.inputCollectionTest = "strip";
////  	seedingCfg.seedFilterConfig.inputCollectionTest =
////  seedingCfg.seedFinderConfig.inputCollectionTest;
//
//  	sequencer.addAlgorithm(std::make_shared<SeedingAlgorithm>(seedingCfg,
//  logLevel));

  return sequencer.run();

  //   // Algorithm estimating track parameter from seed
  //   TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
  //   paramsEstimationCfg.inputSeeds = seedingCfg.outputSeeds;
  //   paramsEstimationCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  //   paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
  //   paramsEstimationCfg.outputTrackParametersSeedMap =
  //   "estimatedparams_seed_map"; paramsEstimationCfg.trackingGeometry =
  //   tGeometry; paramsEstimationCfg.magneticField = magneticField;
  //   sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
  //       paramsEstimationCfg, logLevel));
  //
  //   // Seeding performance Writers
  //   TrackFinderPerformanceWriter::Config tfPerfCfg;
  //   tfPerfCfg.inputProtoTracks = seedingCfg.outputProtoTracks;
  //   tfPerfCfg.inputParticles = inputParticles;
  //   tfPerfCfg.inputMeasurementParticlesMap =
  //       hitSmearingCfg.outputMeasurementParticlesMap;
  //   tfPerfCfg.outputDir = outputDir;
  //   tfPerfCfg.outputFilename = "performance_seeding_trees.root";
  //   sequencer.addWriter(
  //       std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));
  //
  //   SeedingPerformanceWriter::Config seedPerfCfg;
  //   seedPerfCfg.inputSeeds = seedingCfg.outputSeeds;
  //   seedPerfCfg.inputParticles = inputParticles;
  //   seedPerfCfg.inputMeasurementParticlesMap =
  //       hitSmearingCfg.outputMeasurementParticlesMap;
  //   seedPerfCfg.outputDir = outputDir;
  //   seedPerfCfg.outputFilename = "performance_seeding_hists.root";
  //   sequencer.addWriter(
  //       std::make_shared<SeedingPerformanceWriter>(seedPerfCfg, logLevel));
  //
  //   // The track parameters estimation writer
  //   RootTrackParameterWriter::Config trackParamsWriterCfg;
  //   trackParamsWriterCfg.inputSeeds = seedingCfg.outputSeeds;
  //   trackParamsWriterCfg.inputTrackParameters =
  //       paramsEstimationCfg.outputTrackParameters;
  //   trackParamsWriterCfg.inputTrackParametersSeedMap =
  //       paramsEstimationCfg.outputTrackParametersSeedMap;
  //   trackParamsWriterCfg.inputParticles = particleReader.outputParticles;
  //   trackParamsWriterCfg.inputSimHits = simHitReaderCfg.outputSimHits;
  //   trackParamsWriterCfg.inputMeasurementParticlesMap =
  //       hitSmearingCfg.outputMeasurementParticlesMap;
  //   trackParamsWriterCfg.inputMeasurementSimHitsMap =
  //       hitSmearingCfg.outputMeasurementSimHitsMap;
  //   trackParamsWriterCfg.outputDir = outputDir;
  //   trackParamsWriterCfg.outputFilename = "estimatedparams.root";
  //   trackParamsWriterCfg.outputTreename = "estimatedparams";
  //   sequencer.addWriter(std::make_shared<RootTrackParameterWriter>(
  //       trackParamsWriterCfg, logLevel));
  //
  //   return sequencer.run();
}
