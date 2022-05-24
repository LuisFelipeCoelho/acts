// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

ActsExamples::SeedingAlgorithm::SeedingAlgorithm(
    ActsExamples::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  for (const auto& i : m_cfg.inputSpacePoints) {
    if (i.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  if (m_cfg.gridConfig.rMax != m_cfg.seedFinderConfig.rMax) {
    ACTS_WARNING(
        "Inconsistent config rMax: using different values in gridConfig and "
        "seedFinderConfig. Ignore this warning if values are intentional")
  }

  if (m_cfg.seedFilterConfig.deltaRMin != m_cfg.seedFinderConfig.deltaRMin) {
    throw std::invalid_argument("Inconsistent config deltaRMin");
  }

  if (m_cfg.gridConfig.deltaRMax != m_cfg.seedFinderConfig.deltaRMax) {
    //    throw std::invalid_argument("Inconsistent config deltaRMax");
  }

  if (m_cfg.gridConfig.zMin != m_cfg.seedFinderConfig.zMin) {
    throw std::invalid_argument("Inconsistent config zMin");
  }

  if (m_cfg.gridConfig.zMax != m_cfg.seedFinderConfig.zMax) {
    throw std::invalid_argument("Inconsistent config zMax");
  }

  if (m_cfg.seedFilterConfig.maxSeedsPerSpM !=
      m_cfg.seedFinderConfig.maxSeedsPerSpM) {
    throw std::invalid_argument("Inconsistent config maxSeedsPerSpM");
  }

  if (m_cfg.gridConfig.cotThetaMax != m_cfg.seedFinderConfig.cotThetaMax) {
    throw std::invalid_argument("Inconsistent config cotThetaMax");
  }

  if (m_cfg.gridConfig.minPt != m_cfg.seedFinderConfig.minPt) {
    throw std::invalid_argument("Inconsistent config minPt");
  }

  if (m_cfg.gridConfig.bFieldInZ != m_cfg.seedFinderConfig.bFieldInZ) {
    //    throw std::invalid_argument("Inconsistent config bFieldInZ");
  }

  if (m_cfg.seedFinderConfig.cotThetaSorting == false and
      m_cfg.seedFinderConfig.enableCutsForSortedSP == true) {
    throw std::invalid_argument(
        "enableCutsForSortedSP cannot be true if cotThetaSorting is set to "
        "false");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 != m_cfg.zBinNeighborsTop.size() and
      m_cfg.zBinNeighborsTop.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsTop");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 !=
          m_cfg.zBinNeighborsBottom.size() and
      m_cfg.zBinNeighborsBottom.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsBottom");
  }

  //  if (m_cfg.gridConfig.zBinEdges.size()-1 !=
  //  m_cfg.seedFinderConfig.rRangeMiddleSP.size()) {
  //    throw std::invalid_argument("Inconsistent config rRangeMiddleSP");
  //  }

  if (m_cfg.seedFinderConfig.zBinsCustomLooping.size() != 0) {
    for (size_t i = 1; i != m_cfg.gridConfig.zBinEdges.size(); i++) {
      if (std::find(m_cfg.seedFinderConfig.zBinsCustomLooping.begin(),
                    m_cfg.seedFinderConfig.zBinsCustomLooping.end(),
                    i) == m_cfg.seedFinderConfig.zBinsCustomLooping.end()) {
        throw std::invalid_argument("Inconsistent config zBinsCustomLooping");
      }
    }
  }

  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(m_cfg.seedFilterConfig,
                                                        m_cfg.expCuts);
}

ActsExamples::ProcessCode ActsExamples::SeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // construct the combined input container of space point pointers from all
  // configured input sources.
  // pre-compute the total size required so we only need to allocate once
  size_t nSpacePoints = 0;
  for (const auto& isp : m_cfg.inputSpacePoints) {
    nSpacePoints += ctx.eventStore.get<SimSpacePointContainer>(isp).size();
  }

  //	float rMinMiddleSP;
  //	float rMaxMiddleSP;
  //	bool firstSP = true;

  //	std::vector<Acts::Vector3> rRangeSPVector;

  // extent used to store r range for middle spacepoint
  Acts::Extent rRangeSPExtent;

  std::vector<const SimSpacePoint*> spacePointPtrs;
  spacePointPtrs.reserve(nSpacePoints);
  for (const auto& isp : m_cfg.inputSpacePoints) {
    for (const auto& spacePoint :
         ctx.eventStore.get<SimSpacePointContainer>(isp)) {
      // since the event store owns the space points, their pointers should be
      // stable and we do not need to create local copies.
      spacePointPtrs.push_back(&spacePoint);

      //			std::cout << "|sorting| " << spacePoint.r() <<
      //std::endl; 			if (firstSP == true) { 				rMinMiddleSP = spacePoint.r();
      //				rMaxMiddleSP = spacePoint.r();
      //				firstSP = false;
      //			} else if (rMaxMiddleSP < spacePoint.r()) {
      //				rMaxMiddleSP = spacePoint.r();
      //			} else if (rMinMiddleSP > spacePoint.r()) {
      //				rMinMiddleSP = spacePoint.r();
      //			}
      //			rRangeSPVector.push_back({spacePoint.x(),
      //spacePoint.y(), spacePoint.z()});
      // store x,y,z values in extent
      rRangeSPExtent.check({spacePoint.x(), spacePoint.y(), spacePoint.z()});
    }
  }

  //	rRangeSPExtent.ranges[Acts::binR] = {rMinMiddleSP, rMaxMiddleSP};

  //	if (m_cfg.seedFinderConfig.useVariableMiddleSPRange == true ||
  //m_cfg.seedFinderConfig.rRangeMiddleSP.empty() == true) {
  //		m_cfg.seedFinderConfig.rMinMiddleSP = rMinMiddleSP;
  //		m_cfg.seedFinderConfig.rMaxMiddleSP = rMaxMiddleSP;
  //	}

  // construct the seeding tools
  // covariance tool, extracts covariances per spacepoint as required
  auto extractGlobalQuantities =
      [=](const SimSpacePoint& sp, float, float,
          float) -> std::pair<Acts::Vector3, Acts::Vector2> {
    Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
    Acts::Vector2 covariance{sp.varianceR(), sp.varianceZ()};
    return std::make_pair(position, covariance);
  };

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>(m_cfg.zBinNeighborsBottom,
                                     m_cfg.gridConfig.numPhiNeighbors));
  auto topBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>(m_cfg.zBinNeighborsTop,
                                     m_cfg.gridConfig.numPhiNeighbors));
  auto grid =
      Acts::SpacePointGridCreator::createGrid<SimSpacePoint>(m_cfg.gridConfig);
  std::cout << "phiBins, zBins = " << (grid->numLocalBins())[0] << ", "
            << (grid->numLocalBins())[1] << std::endl;
  std::cout << "phi min, z min = " << (grid->minPosition())[0] << ", "
            << (grid->minPosition())[1] << std::endl;
  std::cout << "phi max, z max = " << (grid->maxPosition())[0] << ", "
            << (grid->maxPosition())[1] << std::endl;
  auto spacePointsGrouping = Acts::BinnedSPGroup<SimSpacePoint>(
      spacePointPtrs.begin(), spacePointPtrs.end(), extractGlobalQuantities,
      bottomBinFinder, topBinFinder, std::move(grid), m_cfg.seedFinderConfig);
  auto finder = Acts::Seedfinder<SimSpacePoint>(m_cfg.seedFinderConfig);

  // run the seeding
  static thread_local SimSeedContainer seeds;
  seeds.clear();
  static thread_local decltype(finder)::State state;

  auto group = spacePointsGrouping.begin();
  auto groupEnd = spacePointsGrouping.end();

  std::string inputCollectionTest = m_cfg.seedFinderConfig.inputCollectionTest;

  int i = 0;
  int phibin = 1;
  int zbin = 1;
  bool pass = false;
  std::string s_mid = "|Groups Mid| (0) ==> |";

  for (; !(group == groupEnd); ++group) {
    auto bot = group.bottom();
    auto mid = group.middle();
    auto top = group.top();
    int n_bot = 0;
    int n_mid = 0;
    int n_top = 0;
    // auto mid_bin = grid->
    if (bot.begin() != bot.end() && mid.begin() != mid.end() &&
        top.begin() != top.end()) {
      std::cout << "------> Group " << i << " (phi,z)= " << phibin - 1 << ", "
                << zbin - 1 << std::endl;
      //          std::cout << "phi, z =" << oi << std::endl;
      for (auto sp : bot) {
        n_bot++;
        std::cout << std::setprecision(10) << "Bottom layer:  x= " << sp->x()
                  << " y= " << sp->y() << " z= " << sp->z()
                  << " r= " << sp->radius() << std::endl;
      }
      for (auto sp : mid) {
        n_mid++;
        std::cout << std::setprecision(10) << "Middle layer:  x= " << sp->x()
                  << " y= " << sp->y() << " z= " << sp->z()
                  << " r= " << sp->radius() << std::endl;
      }
      for (auto sp : top) {
        n_top++;
        std::cout << std::setprecision(10) << "Top layer:  x= " << sp->x()
                  << " y= " << sp->y() << " z= " << sp->z()
                  << " r= " << sp->radius() << std::endl;
      }
      i++;
      std::cout << "|Groups " << inputCollectionTest
                << "| n_group, n_mid, n_top, n_bot, phi_bin, z_bin: " << i
                << ", " << n_mid << ", " << n_top << ", " << n_bot << ", "
                << phibin - 1 << ", " << zbin - 1 << std::endl;
    }
    finder.createSeedsForGroup(state, std::back_inserter(seeds), group.bottom(),
                               group.middle(), group.top(), rRangeSPExtent);

    if (zbin % 11 == 0 && pass != true) {
      phibin++;
      zbin = 1;
      pass = true;
      s_mid += std::string("\n");
      s_mid += std::string("|Groups Mid| (") + std::to_string(phibin - 1) +
               std::string(") ==> |");
    } else if (zbin % 11 == 0 && pass == true) {
      zbin++;
      pass = false;
      s_mid += std::to_string(n_mid) + std::string("|");
    } else {
      zbin++;
      pass = false;
      s_mid += std::to_string(n_mid) + std::string("|");
    }
  }
  s_mid += std::string("|Groups Mid| ( z ) ==> |0|1|2|3|4|5|6|7|8|9|10|");
  //    std::cout<<s_mid<<std::endl;

  // extract proto tracks, i.e. groups of measurement indices, from tracks seeds
  size_t nSeeds = seeds.size();
  static thread_local ProtoTrackContainer protoTracks;
  protoTracks.clear();

// 	std::cout << "|Seeds Map| nSeeds_filter: " << nSeeds << " " << 0 <<  std::endl;

  protoTracks.reserve(nSeeds);
  for (const auto& seed : seeds) {
    ProtoTrack& protoTrack = protoTracks.emplace_back();
    protoTrack.reserve(seed.sp().size());
    for (auto spacePointPtr : seed.sp()) {
      protoTrack.push_back(spacePointPtr->measurementIndex());
    }

    // *******************

    auto lb = seed.sp()[0];
    auto lm = seed.sp()[1];
    auto lt = seed.sp()[2];

    float deltaXt = lt->x() - lm->x();
    float deltaYt = lt->y() - lm->y();
    float deltaZt = lt->z() - lm->z();
    float iDeltaR2t = 1. / (deltaXt * deltaXt + deltaYt * deltaYt);

    float deltaXb = lb->x() - lm->x();
    float deltaYb = lb->y() - lm->y();
    float deltaZb = lb->z() - lm->z();
    float iDeltaR2b = 1. / (deltaXb * deltaXb + deltaYb * deltaYb);

    //		float varianceZM = lm.varianceZ();
    //		float varianceRM = lm.varianceR();
    float cosPhiM = lm->x() / lm->r();
    float sinPhiM = lm->y() / lm->r();

    float Ut = (deltaXt * cosPhiM + deltaYt * sinPhiM) * iDeltaR2t;
    float Ub = (deltaXb * cosPhiM + deltaYb * sinPhiM) * iDeltaR2b;

    float Vt = (deltaYt * cosPhiM - deltaXt * sinPhiM) * iDeltaR2t;
    float Vb = (deltaYb * cosPhiM - deltaXb * sinPhiM) * iDeltaR2b;

    float dU = Ut - Ub;
    float A = (Vt - Vb) / dU;
    float S2 = 1. + A * A;
    float B = Vb - A * Ub;
    float B2 = B * B;

    float pT = 300. * 1.997244311 * std::sqrt(S2 / B2) / 2.;

    float cotThetaBt = deltaZt * std::sqrt(iDeltaR2t) * 1;
    float cotThetaBb = deltaZb * std::sqrt(iDeltaR2b) * (-1);

    float Im = std::abs((A - B * lm->r()) * lm->r());

    float eta = -std::log(
        std::tan(0.5 * std::atan(1. / std::sqrt(cotThetaBt * cotThetaBb))));

    std::cout << "|Seeds Map " << inputCollectionTest
              << "| pT_filter, eta_filter, dScore_filter, curvature_filter, "
                 "Im_filter: "
              << std::setprecision(10) << pT / 1000 << " " << eta << " "
              << 0 / 10 << " " << B / std::sqrt(S2) << " " << Im << std::endl;

    // *******************
  }

  std::cout << "|Seeds Map " << inputCollectionTest << "| seeds.size() "
            << seeds.size() << std::endl;

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePointPtrs.size() << " space points");

  ctx.eventStore.add(m_cfg.outputSeeds, SimSeedContainer{seeds});
  ctx.eventStore.add(m_cfg.outputProtoTracks, ProtoTrackContainer{protoTracks});
  return ActsExamples::ProcessCode::SUCCESS;
}
