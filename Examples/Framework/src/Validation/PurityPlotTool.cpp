// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/PurityPlotTool.hpp"

#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

ActsExamples::PurityPlotTool::PurityPlotTool(
    const ActsExamples::PurityPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("PurityPlotTool", lvl)) {}

void ActsExamples::PurityPlotTool::book(
    PurityPlotTool::PurityPlotCache& purityPlotCache) const {
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  ACTS_DEBUG("Initialize the histograms for purity plots");
  // purity vs pT
  purityPlotCache.trackPurity_vs_pT = PlotHelpers::bookEff(
      "trackPurity_vs_pT", "Tracking purity;Truth pT [GeV/c];Purity", bPt);
  // purity vs eta
  purityPlotCache.trackPurity_vs_eta = PlotHelpers::bookEff(
      "trackPurity_vs_eta", "Tracking purity;Truth #eta;Purity", bEta);
  // purity vs phi
  purityPlotCache.trackPurity_vs_phi = PlotHelpers::bookEff(
      "trackPurity_vs_phi", "Tracking purity;Truth #phi;Purity", bPhi);
}

void ActsExamples::PurityPlotTool::clear(
    PurityPlotCache& purityPlotCache) const {
  delete purityPlotCache.trackPurity_vs_pT;
  delete purityPlotCache.trackPurity_vs_eta;
  delete purityPlotCache.trackPurity_vs_phi;
}

void ActsExamples::PurityPlotTool::write(
    const PurityPlotTool::PurityPlotCache& purityPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  purityPlotCache.trackPurity_vs_pT->Write();
  purityPlotCache.trackPurity_vs_eta->Write();
  purityPlotCache.trackPurity_vs_phi->Write();
}

void ActsExamples::PurityPlotTool::fill(
    PurityPlotTool::PurityPlotCache& purityPlotCache,
    const ActsFatras::Particle& truthParticle, bool status) const {
  const auto t_phi = phi(truthParticle.unitDirection());
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillEff(purityPlotCache.trackPurity_vs_pT, t_pT, status);
  PlotHelpers::fillEff(purityPlotCache.trackPurity_vs_eta, t_eta, status);
  PlotHelpers::fillEff(purityPlotCache.trackPurity_vs_phi, t_phi, status);
}
