// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

// Tools to make purity plots to show tracking purity.
// For the moment, the purity is taken as the fraction of successfully
// smoothed track over all tracks
class PurityPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("pT [GeV/c]", 40, 0, 100)}};
  };

  /// @brief Nested Cache struct
  struct PurityPlotCache {
		TEfficiency* trackPurity_vs_pT{nullptr};   ///< Tracking purity vs pT
		TEfficiency* trackPurity_vs_eta{nullptr};  ///< Tracking purity vs eta
		TEfficiency* trackPurity_vs_phi{nullptr};  ///< Tracking purity vs phi
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  PurityPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the purity plots
  ///
  /// @param purityPlotCache the cache for purity plots
  void book(PurityPlotCache& purityPlotCache) const;

  /// @brief fill purity plots
  ///
  /// @param purityPlotCache cache object for purity plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fill(PurityPlotCache& purityPlotCache,
            const ActsFatras::Particle& truthParticle, bool status) const;

  /// @brief write the purity plots to file
  ///
  /// @param purityPlotCache cache object for purity plots
  void write(const PurityPlotCache& purityPlotCache) const;

  /// @brief delete the purity plots
  ///
  /// @param purityPlotCache cache object for purity plots
  void clear(PurityPlotCache& purityPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
