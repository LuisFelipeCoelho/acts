// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedFilter.hpp"

#include <cmath>
#include <numeric>
#include <type_traits>

namespace Acts {

template <typename external_spacepoint_t, typename platform_t>
Seedfinder<external_spacepoint_t, platform_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config)
    : m_config(config.toInternalUnits()) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);
  m_config.sigmapT2perRadius =
      m_config.pT2perRadius * std::pow(2 * m_config.sigmaScattering, 2);
}

template <typename external_spacepoint_t, typename platform_t>
template <template <typename...> typename container_t, typename sp_range_t>
void Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    State& state,
    std::back_insert_iterator<container_t<Seed<external_spacepoint_t>>> outIt,
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs,
    Extent rRangeSPExtent) const {
  //	std::sort(middleSPs.begin(), middleSPs.end(), [] (auto* a, auto* b) ->
  // bool { return (a->radius() < b->radius()); } );

  std::string inputCollectionTest = m_config.inputCollectionTest;

  int i_n = 0;
  for (auto spM : middleSPs) {
    // std::cout << std::endl;
    // std::cout << " rRangeSPExtent min | " << rRangeSPExtent.min(Acts::binR)
    // << std::endl; std::cout << " rRangeSPExtent max | " <<
    // rRangeSPExtent.max(Acts::binR)  << std::endl;

    float rM = spM->radius();
    float zM = spM->z();
    float varianceRM = spM->varianceR();
    float varianceZM = spM->varianceZ();

    std::cout << std::endl;
    std::cout << "|Seeds| rM: " << rM << " deltaRMin: "
              << std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
                     m_config.deltaRMiddleMinSPRange
              << std::endl;
    //		 std::cout << "|Seeds| rMinMiddleSP: " <<
    // m_config.rRanMiddleSP[0]
    //<< " rMaxMiddleSP: " << m_config.rRanMiddleSP[1] << std::endl;

    // ************ delete this --> just for printing values:
    size_t nSeeds = 0;
    size_t nSeeds_filter1 = 0;
    size_t nSeeds_filter2 = 0;
    size_t nSeeds_test1 = 0;
    size_t nSeeds_test2 = 0;
    size_t nSeeds_test3 = 0;

    // THESE ARE ONLY FOR DEBUGGING
    auto pVal = std::lower_bound(m_config.zBinEdges.begin(),
                                 m_config.zBinEdges.end(), zM);
    int zBin = std::distance(m_config.zBinEdges.begin(), pVal);
    zBin == 0 ? zBin : --zBin;

    // *************

    /// check if spM is outside our radial region of interest
    if (m_config.useVariableMiddleSPRange) {
      float rMinMiddleSP = std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
                           m_config.deltaRMiddleMinSPRange;
      float rMaxMiddleSP = std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2 -
                           m_config.deltaRMiddleMaxSPRange;
      std::cout << "rMaxMiddleSP, rMinMiddleSP: " << rMaxMiddleSP << " "
                << rMinMiddleSP << std::endl;
      if (rM < rMinMiddleSP || rM > rMaxMiddleSP) {
        continue;
      }
    } else if (not m_config.rRangeMiddleSP.empty()) {
      /// get zBin position of the middle SP
      auto zValMiddle = std::lower_bound(m_config.zBinEdges.begin(),
                                         m_config.zBinEdges.end(), zM);
      int zBinMiddle = std::distance(m_config.zBinEdges.begin(), zValMiddle);
      /// protects against zM at the limit of zBinEdges
      zBinMiddle == 0 ? zBinMiddle : --zBinMiddle;
      if (rM < m_config.rRangeMiddleSP[zBinMiddle][0] ||
          rM > m_config.rRangeMiddleSP[zBinMiddle][1]) {
        // std::cout << "|Seeds| !!! rM < rRangeMiddleSP || rM > rRangeMiddleSP
        // == TRUE !!!" << std::endl;
        continue;
      }
    }

    //		if (std::abs(zM) > 2700) {
    //			continue;
    //		}

    size_t nTopSeedConf = 0;
    if (m_config.seedConfirmation == true) {
      // check if middle SP is in the central or forward region
      SeedConfirmationRange seedConfRange =
          (zM > m_config.centralSeedConfirmationRange.zMaxSeedConf ||
           zM < m_config.centralSeedConfirmationRange.zMinSeedConf)
              ? m_config.forwardSeedConfirmationRange
              : m_config.centralSeedConfirmationRange;
      // set the minimum number of top SP depending on whether the middle SP is
      // in the central or forward region
      nTopSeedConf = rM > seedConfRange.rMaxSeedConf
                         ? seedConfRange.nTopForLargeR
                         : seedConfRange.nTopForSmallR;
    }

    state.compatTopSP.clear();

    std::cout << i_n << ") ==== TOP SP ====" << std::endl;

    for (auto topSP : topSPs) {
      float rT = topSP->radius();
      float deltaR = rT - rM;
      // std::cout << std::endl;
      std::cout << "|Seeds| --> rT: " << rT << " deltaR: " << deltaR
                << std::endl;
      std::cout << "|Seeds| deltaRMax: " << m_config.deltaRMaxTopSP
                << " deltaRMin: " << m_config.deltaRMinTopSP << std::endl;
      // if r-distance is too small, try next SP in bin
      if (deltaR < m_config.deltaRMinTopSP) {
        std::cout << "|Seeds| !!! deltaR < deltaRMin == TRUE !!!" << std::endl;
        continue;
      }
      // if r-distance is too big, try next SP in bin
      if (deltaR > m_config.deltaRMaxTopSP) {
        std::cout << "|Seeds| !!! deltaR > deltaRMax == TRUE !!!" << std::endl;
        continue;
      }
      // ratio Z/R (forward angle) of space point duplet
      float deltaZ = topSP->z() - zM;
			float iDeltaR = 1 / deltaR;
      float cotTheta = deltaZ * iDeltaR;
      //			topSP->setCotTheta(cotTheta);
      if (m_config.cotThetaMaxCut) {
        if (std::fabs(cotTheta) > m_config.cotThetaMax) {
          // std::cout << "|Seeds| !!! fabs(cotTheta) > cotThetaMax == TRUE !!!"
          // << std::endl;
          continue;
        }
      }
      // check if duplet origin on z axis within collision region
      float zOrigin = zM - rM * cotTheta;
      // std::cout << "|Seeds| zM - rM * cotTheta " << zM << " - " << rM << " *
      // " << cotTheta << std::endl; std::cout << "|Seeds| zOrigin: " << zOrigin
      // << " collisionRegionMin: " << m_config.collisionRegionMin << std::endl;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        // std::cout << "|Seeds| !!! zOrigin < collisionRegionMin || zOrigin >
        // collisionRegionMax == TRUE !!!" << std::endl;
        continue;
      }

      if (m_config.deltaZCut) {
        if (std::abs(deltaZ) > m_config.deltaZMax) {
          continue;
        }
      }
      // cut on the max curvature between top SP and interaction point
      // first transform the space point coordinates into a frame such that the
      // central space point SPm is in the origin of the frame and the x axis
      // points away from the interaction point in addition to a translation
      // transformation we also perform a rotation in order to keep the
      // curvature of the circle tangent to the x axis
      if (m_config.interactionPointCut) {
				float irM = 1 / rM;
        float xVal = (topSP->x() - spM->x()) * spM->x() * irM +
                     (topSP->y() - spM->y()) * spM->y() * irM;
        float yVal = (topSP->y() - spM->y()) * spM->x() * irM -
                     (topSP->x() - spM->x()) * spM->y() * irM;
        std::cout << std::abs(rM) << " * " << std::abs(yVal) << " > "
                  << -m_config.impactMax << " * " << xVal << std::endl;
        if (std::abs(rM * yVal) > m_config.impactMax * xVal) {
          // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the
          // circle into straight lines in the u/v plane the line equation can
          // be described in terms of aCoef and bCoef, where v = aCoef * u +
          // bCoef
          float uT = xVal / (xVal * xVal + yVal * yVal);
          float vT = yVal / (xVal * xVal + yVal * yVal);
          // in the rotated frame the interaction point is positioned at x = -rM
          // and y ~= impactParam
          float uIP = -1. / rM;
          float vIP = m_config.impactMax / (rM * rM);
          if (yVal > 0.)
            vIP = -vIP;
          // we can obtain aCoef as the slope dv/du of the linear function,
          // estimated using du and dv between the two SP bCoef is obtained by
          // inserting aCoef into the linear equation
          float aCoef = (vT - vIP) / (uT - uIP);
          float bCoef = vIP - aCoef * uIP;
          // the distance of the straight line from the origin (radius of the
          // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
          // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
          if ((bCoef * bCoef) >
              (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
            std::cout << bCoef << " * " << bCoef << " > 1+" << aCoef << " * "
                      << aCoef << " / " << m_config.minHelixDiameter2
                      << std::endl;
            std::cout << "|Seeds| !!!  impact parameter cut == TRUE !!!"
                      << std::endl;
            continue;
          }
        }
      }
      state.compatTopSP.push_back(topSP);
      std::cout << "|Seeds| # Fill SP" << std::endl;
    }
    if (state.compatTopSP.empty()) {
      // std::cout << "|Seeds| !!! top SP empty == TRUE !!!" << std::endl;
      continue;
    }
    // apply cut on the number of top SP if seedConfirmation is true
    if (m_config.seedConfirmation == true &&
        state.compatTopSP.size() < nTopSeedConf) {
      // std::cout << "|Seeds| !!! seedConfirmation top SP cut == TRUE !!!" <<
      // std::endl;
      continue;
    }

    state.compatBottomSP.clear();

    std::cout << i_n << ") ==== BOTTOM SP ====" << std::endl;
    i_n++;

    for (auto bottomSP : bottomSPs) {
      float rB = bottomSP->radius();
      float deltaR = rM - rB;
      // std::cout << std::endl;
      std::cout << "|Seeds| --> rB: " << rB << " deltaR: " << deltaR
                << std::endl;
      std::cout << "|Seeds| deltaRMax: " << m_config.deltaRMaxBottomSP
                << " deltaRMin: " << m_config.deltaRMinBottomSP << std::endl;
      // this condition is the opposite of the condition for top SP
      if (deltaR > m_config.deltaRMaxBottomSP) {
        std::cout << "|Seeds| !!! deltaR > deltaRMax == TRUE !!!" << std::endl;
        continue;
      }
      if (deltaR < m_config.deltaRMinBottomSP) {
        std::cout << "|Seeds| !!! deltaR < deltaRMin == TRUE !!!" << std::endl;
        continue;
      }
      // ratio Z/R (forward angle) of space point duplet
      float deltaZ = zM - bottomSP->z();
			float iDeltaR = 1 / deltaR;
			float cotTheta = deltaZ * iDeltaR;
      std::cout << "|Seeds| cotTheta: " << cotTheta
                << " cotThetaMax: " << m_config.cotThetaMax << std::endl;
      //			bottomSP->setCotTheta(cotTheta);
      if (m_config.cotThetaMaxCut) {
        if (std::fabs(cotTheta) > m_config.cotThetaMax) {
          //				//std::cout << "|Seeds| cotTheta: " <<
          // cotTheta << " cotThetaMax: " << m_config.cotThetaMax << std::endl;
          // std::cout << "|Seeds| !!! fabs(cotTheta) > cotThetaMax == TRUE !!!"
          // << std::endl;
          continue;
        }
      }
      // check if duplet origin on z axis within collision region
      std::cout << "|Seeds| zM - rM * cotTheta " << zM << " - " << rM << " * "
                << cotTheta << std::endl;
      float zOrigin = zM - rM * cotTheta;
      //
      //        //std::cout << "|Seeds| zOrigin: " << zOrigin << "
      //        collisionRegionMin: " << m_config.collisionRegionMin <<
      //        std::endl;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        // std::cout << "|Seeds| !!! zOrigin < collisionRegionMin || zOrigin >
        // collisionRegionMax == TRUE !!!" << std::endl;
        continue;
      }

      if (m_config.deltaZCut) {
        if (std::abs(deltaZ) > m_config.deltaZMax) {
          continue;
        }
      }
      // cut on the max curvature between bottom SP and interaction point
      // first transform the space point coordinates into a frame such that the
      // central space point SPm is in the origin of the frame and the x axis
      // points away from the interaction point in addition to a translation
      // transformation we also perform a rotation in order to keep the
      // curvature of the circle tangent to the x axis
      if (m_config.interactionPointCut) {
				float irM = 1 / rM;
        float xVal = (bottomSP->x() - spM->x()) * spM->x() * irM +
                     (bottomSP->y() - spM->y()) * spM->y() * irM;
        float yVal = (bottomSP->y() - spM->y()) * spM->x() * irM -
                     (bottomSP->x() - spM->x()) * spM->y() * irM;
        // std::cout << std::abs(rM) << " * " << std::abs(yVal) << " > " <<
        // -m_config.impactMax << " * " << xVal << std::endl;
        if (std::abs(rM * yVal) > -m_config.impactMax * xVal) {
          // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the
          // circle into straight lines in the u/v plane the line equation can
          // be described in terms of aCoef and bCoef, where v = aCoef * u +
          // bCoef
          float uB = xVal / (xVal * xVal + yVal * yVal);
          float vB = yVal / (xVal * xVal + yVal * yVal);
          // in the rotated frame the interaction point is positioned at x = -rM
          // and y ~= impactParam
          float uIP = -1. / rM;
          float vIP = m_config.impactMax / (rM * rM);
          if (yVal < 0.)
            vIP = -vIP;
          // we can obtain aCoef as the slope dv/du of the linear function,
          // estimated using du and dv between the two SP bCoef is obtained by
          // inserting aCoef into the linear equation
          float aCoef = (vB - vIP) / (uB - uIP);
          float bCoef = vIP - aCoef * uIP;
          // the distance of the straight line from the origin (radius of the
          // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
          // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
          std::cout << bCoef << "^2 > (1 - " << aCoef << "^2) /"
                    << m_config.minHelixDiameter2 << std::endl;
          std::cout << "|Seeds| !!! impact parameter cut == TRUE !!!"
                    << std::endl;
          if ((bCoef * bCoef) >
              (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
            continue;
          }
        }
      }
      state.compatBottomSP.push_back(bottomSP);
      std::cout << "|Seeds| # Fill SP" << std::endl;
    }
    // no bottom SP found -> try next spM
    if (state.compatBottomSP.empty()) {
      // std::cout << "|Seeds| !!! bottom SP empty == TRUE !!!" << std::endl;
      continue;
    }

    state.linCircleBottom.clear();
    state.linCircleTop.clear();

    transformCoordinates(state.compatBottomSP, *spM, true,
                         m_config.cotThetaSorting, state.linCircleBottom);
    transformCoordinates(state.compatTopSP, *spM, false,
                         m_config.cotThetaSorting, state.linCircleTop);

    state.topSpVec.clear();
    state.curvatures.clear();
    state.impactParameters.clear();
    state.seedsPerSpM.clear();
    //		state.topDeltaR.clear();

    size_t numBotSP = state.compatBottomSP.size();
    size_t numTopSP = state.compatTopSP.size();

    int m_nQualitySeeds = 0;
    int m_nSeeds = 0;

    std::cout << std::endl;
    std::cout << " ------- Filled SPs -------" << std::endl;

    if (spM->validDoubleMeasurementDetails() == false) {
      std::cout << "ERROR FALSE DETAiLS" << std::endl;
    }

    for (size_t b = 0; b < numBotSP; b++) {
      for (size_t t = 0; t < numTopSP; t++) {
        auto lb = state.linCircleBottom[b];
        auto lt = state.linCircleTop[t];
        std::cout << std::endl;
        std::cout << "--> rM: " << rM
                  << " rB: " << state.compatBottomSP[b]->radius()
                  << " rT: " << state.compatTopSP[t]->radius() << std::endl;
        std::cout << "xB: " << state.compatBottomSP[b]->x()
                  << " yB: " << state.compatBottomSP[b]->y()
                  << " zB: " << state.compatBottomSP[b]->z()
                  << " xT: " << state.compatTopSP[t]->x()
                  << " yT: " << state.compatTopSP[t]->y()
                  << " zT: " << state.compatTopSP[t]->z() << std::endl;
      }
    }

    std::cout << std::endl;
    std::cout << "|Seeds| --- Compare SP ---" << std::endl;

    size_t t0 = 0;

    for (size_t b = 0; b < numBotSP; b++) {
      auto lb = state.linCircleBottom[b];
      float Zob = lb.Zo;
//      float cotThetaB = lb.cotTheta;
			float cotThetaBB = lb.cotTheta;
      float Vb = lb.V;
      float Ub = lb.U;
      float ErB = lb.Er;
      float iDeltaRB = lb.iDeltaR;

      std::cout << " iSinTheta^-1 "
                << 1 / std::sqrt((1. + cotThetaBB * cotThetaBB)) << " iSinTheta2 "
                << (1. + cotThetaBB * cotThetaBB) << " cotThetaB " << cotThetaBB
                << std::endl;

      // 1+(cot^2(theta)) = 1/sin^2(theta)
      float iSinTheta2 = (1. + cotThetaBB * cotThetaBB);
      // calculate max scattering for min momentum at the seed's theta angle
      // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
      // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
      // scattering
      // but to avoid trig functions we approximate cot by scaling by
      // 1/sin^4(theta)
      // resolving with pT to p scaling --> only divide by sin^2(theta)
      // max approximation error for allowed scattering angles of 0.04 rad at
      // eta=infinity: ~8.5%
      float scatteringInRegion2 = m_config.maxScatteringAngle2 * iSinTheta2;
      // multiply the squared sigma onto the squared scattering
      scatteringInRegion2 *=
          m_config.sigmaScattering * m_config.sigmaScattering;

      // clear all vectors used in each inner for loop
      state.topSpVec.clear();
      state.curvatures.clear();
      state.impactParameters.clear();
			


      float iSinTheta = std::sqrt(iSinTheta2);
      float sinTheta = 1 / std::sqrt(iSinTheta2);  // problema com divisão??
      float cosTheta = cotThetaBB * sinTheta;
      float Ce = cotThetaBB * sinTheta;

      std::cout << "Ce = cotThetaB * iSinTheta^-1 " << Ce << " " << cotThetaBB
                << " " << 1 / iSinTheta << std::endl;

      float Sx = spM->x() * sinTheta / (spM->radius());
      float Sy = spM->y() * sinTheta / (spM->radius());

      //			state.topDeltaR.clear();
      for (size_t t = t0; t < numTopSP; t++) {
        auto lt = state.linCircleTop[t];

        float cotThetaT = lt.cotTheta;
//        float iDeltaRT = lt.iDeltaR;

        std::cout << std::endl;
        std::cout << "--> rM: " << rM
                  << " rB: " << state.compatBottomSP[b]->radius()
                  << " rT: " << state.compatTopSP[t]->radius() << std::endl;
        std::cout << "xB: " << state.compatBottomSP[b]->x()
                  << " yB: " << state.compatBottomSP[b]->y()
                  << " zB: " << state.compatBottomSP[b]->z()
                  << " xT: " << state.compatTopSP[t]->x()
                  << " yT: " << state.compatTopSP[t]->y()
                  << " zT: " << state.compatTopSP[t]->z() << std::endl;
        std::cout << "|Seeds Map| SP: rM, rB, rT, zM, zB, zT: " << rM << " "
                  << state.compatBottomSP[b]->radius() << " "
                  << state.compatTopSP[t]->radius() << " " << zM << " "
                  << state.compatBottomSP[b]->z() << " "
                  << state.compatTopSP[t]->z() << std::endl;

        nSeeds_test1 += 1;

        float rMCoord;
        float ub;
        float vb;

        float ut;
        float vt;

        if (spM->validDoubleMeasurementDetails() == true) {
          // protects against division by 0
          float dU = lt.U - Ub;
          if (dU == 0.) {
            continue;
          }
          // A and B are evaluated as a function of the circumference parameters
          // x_0 and y_0
          float A0 = (lt.V - Vb) / dU;

          // middle
          //          float iSinTheta = std::sqrt(iSinTheta2);
          //          float cosTheta = cotThetaB / iSinTheta;
          //					float rotationTermsUVtoXY[2] = {spM->x() /
          //(iSinTheta * spM->radius()), spM->y() / (iSinTheta *
          //spM->radius())};
          //					// position of Middle SP converted from UV to XY
          //assuming cotTheta evaluated from the Bottom and Middle SPs 					double
          //positionMiddle[3] = {rotationTermsUVtoXY[0] - rotationTermsUVtoXY[1]
          //* A0, 																			rotationTermsUVtoXY[0] * A0 + rotationTermsUVtoXY[1], 																			cosTheta
          //* std::sqrt(1 + A0 * A0)}; 					double dn[3] = {rotationTermsUVtoXY[0] -
          //rotationTermsUVtoXY[1] * A0, 						rotationTermsUVtoXY[0] * A0 +
          //rotationTermsUVtoXY[1], 						cosTheta * std::sqrt(1 + A0 * A0)};

          // remove this:
          //					float Ce = cotThetaB /
          //iSinTheta;
          //          float Sx = spM->x() / (iSinTheta * spM->radius());
          //          float Sy = spM->y() / (iSinTheta * spM->radius());
          float Cn = cosTheta * std::sqrt(1. + A0 * A0);
          double dn[3] = {Sx - Sy * A0, Sx * A0 + Sy, Cn};

          double rTestM[3];
          if (!coordinates(dn, rTestM, spM))
            continue;

          std::cout << " iSinTheta^-1 " << 1 / iSinTheta << " iSinTheta2 "
                    << iSinTheta2 << " cotThetaBB " << cotThetaBB << " cosTheta "
                    << cosTheta << " A0 " << A0 << std::endl;
          std::cout << lb.x << " " << lb.y << " " << lb.iDeltaR << std::endl;
          std::cout << lt.x << " " << lt.y << " " << lt.iDeltaR << std::endl;
          std::cout << Ce << " " << Sx << " " << Sy << " " << Cn << " " << A0
                    << " " << lt.U << " " << Ub << " " << lt.V << " " << Vb
                    << std::endl;
          std::cout << "rTest: " << rTestM[0] << " " << rTestM[1] << " " << rTestM[2]
                    << std::endl;

          // bottom
          float B0 = 2. * (Vb - A0 * Ub);
          float Cb = 1. - B0 * lb.y;
          float Sb = A0 + B0 * lb.x;
          double db[3] = {Sx * Cb - Sy * Sb, Sx * Sb + Sy * Cb, Cn};

          auto spB = state.compatBottomSP[b];
          double rTestB[3];
          if (!coordinates(db, rTestB, spB))
            continue;

          std::cout << B0 << " " << Cb << " " << Sb << std::endl;
          std::cout << "rTest: " << rTestB[0] << " " << rTestB[1] << " " << rTestB[2]
                    << std::endl;

          // top
          float Ct = 1. - B0 * lt.y;
          float St = A0 + B0 * lt.x;
          double dt[3] = {Sx * Ct - Sy * St, Sx * St + Sy * Ct, Cn};

          auto spT = state.compatTopSP[t];
          double rTestT[3];
          if (!coordinates(dt, rTestT, spT))
            continue;

          std::cout << Ct << " " << St << std::endl;
          std::cout << "rTest: " << rTestT[0] << " " << rTestT[1] << " " << rTestT[2]
                    << std::endl;

          float xB = rTestB[0] - rTestM[0];
          float yB = rTestB[1] - rTestM[1];
          float zB = rTestB[2] - rTestM[2];
          float xT = rTestT[0] - rTestM[0];
          float yT = rTestT[1] - rTestM[1];
          float zT = rTestT[2] - rTestM[2];

          float iDeltaRB2 = 1. / (xB * xB + yB * yB);
//          iDeltaRB = std::sqrt(iDeltaRB2);                 // ***** <======
          //					lb.iDeltaR = iDeltaRB;

          float iDeltaRT2 = 1. / (xT * xT + yT * yT);
//          iDeltaRT = std::sqrt(iDeltaRT2);
          //					lt.iDeltaR =
          // std::sqrt(iDeltaRT2);

          cotThetaBB = -zB * std::sqrt(iDeltaRB2);
          cotThetaT = zT * std::sqrt(iDeltaRT2);
          //					lb.cotTheta = cotThetaB;
          //					lt.cotTheta = zT * lt.iDeltaR;
          std::cout << xB << " " << yB << " " << zB << " " << xT << " " << yT
                    << " " << zT << " " << iDeltaRB2 << " " << iDeltaRT2 << " "
                    << cotThetaBB << " " << cotThetaT << std::endl;

          rMCoord = std::sqrt(rTestM[0] * rTestM[0] + rTestM[1] * rTestM[1]);
          float Ax = rTestM[0] / rMCoord;
          float Ay = rTestM[1] / rMCoord;

          ub = (xB * Ax + yB * Ay) * iDeltaRB2;
          vb = (yB * Ax - xB * Ay) * iDeltaRB2;
          ut = (xT * Ax + yT * Ay) * iDeltaRT2;
          vt = (yT * Ax - xT * Ay) * iDeltaRT2;

          std::cout << rMCoord << " " << Ax << " " << Ay << std::endl;

        } else {
//          cotThetaBB = cotThetaB;
        }

        float cotTheta2;
				float cotTheta;
        if (m_config.arithmeticAverageCotTheta) {
					cotTheta = (cotThetaBB + cotThetaT) / 2;
          cotTheta2 = std::pow((cotThetaBB + cotThetaT) / 2, 2);
        } else {
          cotTheta2 = cotThetaBB * cotThetaT;
        }
        // add errors of spB-spM and spM-spT pairs and add the correlation term
        // for errors on spM
        float error2 =
            lt.Er + ErB +
            2 * (cotTheta2 * varianceRM + varianceZM) * iDeltaRB * lt.iDeltaR;

        float deltaCotTheta = cotThetaBB - cotThetaT;
        float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
        float error;
        float dCotThetaMinusError2;
        if (m_config.cotThetaSorting) {
          // if the error is larger than the difference in theta, no need to
          // compare with scattering
          std::cout << "|Seeds| dT: " << deltaCotTheta2 << " - " << lt.Er + ErB
                    << " - "
                    << 2 * (cotThetaBB * cotThetaT * varianceRM + varianceZM) *
                           iDeltaRB * lt.iDeltaR
                    << std::endl;
          std::cout << "|Seeds| Ert + Erb: " << lt.Er << " + " << ErB
                    << std::endl;
          std::cout << "|Seeds| 3ºterm: " << 2 << " * (" << cotThetaBB << " * "
                    << cotThetaT << " * " << varianceRM << " + " << varianceZM
                    << ") * " << iDeltaRB << " * " << lt.iDeltaR << std::endl;
          std::cout << "|Seeds| ICSA: " << scatteringInRegion2 << " * "
                    << 134 * .05 * 9. * (1. / std::abs(m_config.minPt)) *
                           (1. / std::abs(m_config.minPt))
                    << std::endl;
          std::cout << "|Seeds| cotThetaB: " << cotThetaBB
                    << " cotThetaT: " << cotThetaT << " DT: " << deltaCotTheta
                    << std::endl;
          std::cout << "|Seeds| m_config.pTPerHelixRadius "
                    << m_config.pTPerHelixRadius << std::endl;
          if (deltaCotTheta2 - error2 > scatteringInRegion2) {
            // additional cut to skip top SPs when producing triplets
            if (m_config.enableCutsForSortedSP) {
              // break if cotThetaBB < cotThetaT because the SP are sorted by
              // cotTheta
              if (cotThetaBB - cotThetaT < 0) {
                std::cout << "** BREAK:  deltaCotTheta2 - error2 > "
                             "scatteringInRegion2"
                          << std::endl;
                break;
              }
              t0 = t + 1;
            }
            std::cout
                << "** Continue:  deltaCotTheta2 - error2 > scatteringInRegion2"
                << std::endl;
            continue;
          }
        } else {
          if (deltaCotTheta2 - error2 > 0) {
            deltaCotTheta = std::abs(deltaCotTheta);
            // if deltaTheta larger than the scattering for the lower pT cut,
            // skip
            error = std::sqrt(error2);
            dCotThetaMinusError2 =
                deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
            // avoid taking root of scatteringInRegion
            // if left side of ">" is positive, both sides of unequality can be
            // squared
            // (scattering is always positive)
            if (dCotThetaMinusError2 > scatteringInRegion2) {
              continue;
            }
          }
        }

        nSeeds_test2 += 1;

        // protects against division by 0
        float dU = lt.U - Ub;
        if (dU == 0.) {
          continue;
        }
        // A and B are evaluated as a function of the circumference parameters
        // x_0 and y_0
        float A = (lt.V - Vb) / dU;
        float S2 = 1. + A * A;
        float B = Vb - A * Ub;
        float B2 = B * B;

        if (spM->validDoubleMeasurementDetails() == true) {
          dU = ut - ub;
          if (dU == 0.) {
            continue;
          }
          A = (vt - vb) / dU;
          S2 = 1. + A * A;
          B = vb - A * ub;
          B2 = B * B;
        }

//        std::cout << "test IM: A " << A << " lt.V " << lt.V << " Vb " << Vb
//                  << " lt.U " << lt.U << " Ub " << Ub << " B " << B;
                std::cout << "test IM: A " << A << " lt.V " << vt << " Vb " <<
                vb
                          << " lt.U " << ut << " Ub " << ub << " B " << B
                          << std::endl;

        // sqrt(S2)/B = 2 * helixradius
        // calculated radius must not be smaller than minimum radius
        std::cout << "|Seeds| S2: " << S2 << " B2: " << B2
                  << " minHelixDiameter2: " << m_config.minHelixDiameter2 << " dT " << (deltaCotTheta2 - error2)
                  << std::endl;
        std::cout << "|Seeds| Vb - A * Ub: " << Vb << " - " << A << " * " << Ub
                  << std::endl;
        std::cout << "|Seeds| dT * S2 " << (deltaCotTheta2 - error2) * S2
                  << " CSA: "
                  << iSinTheta2 * 134 * .05 * 9 * 2 * 2 /
                         (m_config.pTPerHelixRadius * m_config.pTPerHelixRadius)
                  << " m_COFK " << 134 * .05 * 9 * 1000000 / (300 * 300)
                  << " iSinTheta2: " << iSinTheta2 << std::endl;
        // calculated radius must not be smaller than minimum radius
        if (S2 < B2 * m_config.minHelixDiameter2) {
          std::cout << "** Continue:  S2 < B2 * m_config.minHelixDiameter2"
                    << std::endl;
          continue;
        }

        nSeeds_test3 += 1;

        // refinement of the cut on the compatibility between the r-z slope of
        // the two seed segments using a scattering term scaled by the actual
        // measured pT
        float iHelixDiameter2 = B2 / S2;
        // calculate scattering for p(T) calculated from seed curvature
        float pT2scatterSigma = iHelixDiameter2 * m_config.sigmapT2perRadius;
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        float pT = m_config.pTPerHelixRadius * std::sqrt(S2 / B2) / 2.;
        if (pT > m_config.maxPtScattering) {
          float pTscatterSigma =
              (m_config.highland / m_config.maxPtScattering) *
              m_config.sigmaScattering;
          pT2scatterSigma = pTscatterSigma * pTscatterSigma;
        }
        // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
        // from rad to deltaCotTheta
        float p2scatterSigma = pT2scatterSigma * iSinTheta2;
        // if deltaTheta larger than allowed scattering for calculated pT, skip
        if (m_config.cotThetaSorting) {
          if (deltaCotTheta2 - error2 > p2scatterSigma) {
						std::cout << " deltaCotTheta2, error2, p2scatterSigma: " << deltaCotTheta2 << " " << error2 << " " << p2scatterSigma << std::endl;
            if (m_config.enableCutsForSortedSP) {
              if (cotThetaBB - cotThetaT < 0) {
                std::cout
                    << "** BREAK:  deltaCotTheta2 - error2 > p2scatterSigma"
                    << std::endl;
                break;
              }
              t0 = t;
            }
            std::cout
                << "** Continue:  deltaCotTheta2 - error2 > p2scatterSigma"
                << std::endl;
            continue;
          }
        } else {
          if ((deltaCotTheta2 - error2 > 0) &&
              (dCotThetaMinusError2 > p2scatterSigma)) {
            continue;
          }
        }

        // A and B allow calculation of impact params in U/V plane with linear
        // function
        // (in contrast to having to solve a quadratic function in x/y plane)
        float Im = std::abs((A - B * rM) * rM);

        if (spM->validDoubleMeasurementDetails() == true) {
          Im = std::abs((A - B * rMCoord) * rMCoord);
        }
//        std::cout << "|test ImpactPar|" << A << " " << B << " " << rM << " "
//                  << Im << std::endl;
        				std::cout << "|test ImpactPar|" << A << " " <<
        B << " " << rMCoord << " "
                          << Im << std::endl;

        if (Im <= m_config.impactMax) {
          nSeeds += 1;

          std::cout << std::endl;
          std::cout << "|Seeds Map| Seeds: rM, rB, rT, zM, zB, zT: " << rM
                    << " " << state.compatBottomSP[b]->radius() << " "
                    << state.compatTopSP[t]->radius() << " " << zM << " "
                    << state.compatBottomSP[b]->z() << " "
                    << state.compatTopSP[t]->z() << std::endl;

          state.topSpVec.push_back(state.compatTopSP[t]);
          // inverse diameter is signed depending if the curvature is
          // positive/negative in phi
          state.curvatures.push_back(B / std::sqrt(S2));
          state.impactParameters.push_back(Im);

          //					state.topDeltaR.push_back(lt.topDeltaR);
          //					std::cout << "TEEST " << lt.topDeltaR
          //<<
          //"
          //"
          //<< state.topDeltaR.size() << std::endl;

          if (B2 < 1e-8) {
            B2 = 1e-8;
            pT = m_config.pTPerHelixRadius * std::sqrt(S2 / B2) / 2.;
          }

          std::cout << "|test cotTheta2 geometric| " << cotThetaBB << " "
                    << cotThetaT << std::endl;

          std::cout << "|test cotTheta2 geometric| " << cotTheta2 << " "
                    << std::sqrt(cotTheta2) << " " << -std::sqrt(cotTheta2)
                    << " "
                    << -std::log(std::tan(
                           0.5 * std::atan(1. / -std::sqrt(cotTheta2))))
                    << std::endl;

          std::cout
              << "|test cotTheta2 arithmetic| "
              << std::pow((cotThetaBB + cotThetaT) / 2, 2) << " "
              << std::sqrt(std::pow((cotThetaBB + cotThetaT) / 2, 2)) << " "
              << -std::sqrt(std::pow((cotThetaBB + cotThetaT) / 2, 2)) << " "
              << -std::log(std::tan(
                     0.5 *
                     std::atan(1. / -std::sqrt(std::pow(
                                        (cotThetaBB + cotThetaT) / 2, 2)))))
              << std::endl;

//          cotTheta = std::sqrt(cotTheta2);
          
					if (std::abs(cotTheta) < 1e-8)
						cotTheta = 1e-8;

					if (state.compatTopSP[t]->z() < 0) {
							cotTheta = -cotTheta;
					}
					
          // evaluate eta and pT of the seed
          float theta = std::atan(1. / cotTheta);

          float eta = -std::log(std::tan(0.5 * theta));
          state.etaVec.push_back(eta);
          state.ptVec.push_back(pT);
          state.cotThetaVec.push_back(cotThetaT);

          // std::cout << "ACEPTED" << std::endl;
          std::cout << "|Seeds Map " << inputCollectionTest
                    << "| pT, eta, dScore, curvature, Im: "
                    << std::setprecision(10) << pT / 1000 << " " << eta << " "
                    << 0 << " " << B / std::sqrt(S2) << " " << Im << std::endl;
        }
      }

      std::cout << "|size1| " << (state.seedsPerSpM).size() << std::endl;

      if (!state.topSpVec.empty()) {
        m_config.seedFilter->filterSeeds_2SpFixed(
            *state.compatBottomSP[b], *spM, state.topSpVec, state.curvatures,
            state.impactParameters, state.cotThetaVec, Zob, m_nQualitySeeds,
            m_nSeeds, state.seedsPerSpM);

        //				m_config.seedFilter->filterSeeds_2SpFixedConfirmation(state.seedsPerSpM,
        //m_nQualitySeeds);
      }
    }
    std::cout << "|size2| " << (state.seedsPerSpM).size() << std::endl;

    nSeeds_filter1 = (state.seedsPerSpM).size();

    m_config.seedFilter->filterSeeds_1SpFixed(state.seedsPerSpM,
                                              m_nQualitySeeds, outIt);

    std::cout << "|size3| " << (state.seedsPerSpM).size() << std::endl;

    nSeeds_filter2 = (state.seedsPerSpM).size();

    // std::cout << std::endl;
    std::cout << "|Seeds Map " << inputCollectionTest
              << "| nSeeds, zBin, phiBin: " << nSeeds << " " << zBin << " "
              << std::abs(std::ceil(spM->phi() * 1 / (2 * 3.14159265359 / 138)))
              << std::endl;
    std::cout << "|Seeds Map " << inputCollectionTest
              << "| nSeeds_test: " << nSeeds_test1 << " " << nSeeds_test2 << " "
              << nSeeds_test3 << std::endl;
    //    std::cout << "|Seeds Map " << inputCollectionTest
    //              << "| nSeeds_filter: " << nSeeds_filter1 << " " <<
    //              nSeeds_filter2
    //              << std::endl;

    // std::cout << "|Seeds Map| nSeeds, zBin, phiBin: " <<
    // state.seedsPerSpM.size() << " " << zBin << " " <<
    // std::ceil(spM->phi()*1/(2*3.14159265359/138)) << std::endl;
  }
}

template <typename sp_range_t>
bool coordinates(const double* d, double* rTest, sp_range_t sp) {
  // check the compatibility of SPs coordinates in xyz assuming the
  // Bottom-Middle direction with the strip meassument details

  // add definition
  // d-> spacepointPosition
  // rTest-> outputCoordinates
  // return bool that says if SP is conpatible with being inside the detector
  // element

  std::cout << "|check strip coord| " << sp->topHalfStripLength() << std::endl;
  std::cout << "|check strip coord| " << sp->topStripDirection()[0] << " "
            << sp->topStripDirection()[1] << " " << sp->topStripDirection()[2]
            << std::endl;
  std::cout << "|check strip coord| " << sp->bottomStripDirection()[0] << " "
            << sp->bottomStripDirection()[1] << " "
            << sp->bottomStripDirection()[2] << std::endl;
  std::cout << "|check strip coord| " << sp->stripCenterDistance()[0] << " "
            << sp->stripCenterDistance()[1] << " "
            << sp->stripCenterDistance()[2] << std::endl;
  std::cout << "|check strip coord| " << sp->bottomStripCenterPosition()[0]
            << " " << sp->bottomStripCenterPosition()[1] << " "
            << sp->bottomStripCenterPosition()[2] << std::endl;
  std::cout << "|check strip coord| " << d[0] << " " << d[1] << " " << d[2]
            << std::endl;
//5.514911175
	// m_b1 = bottomHalfStripLength * bottomStripDirection
	// m_b0 = topHalfStripLength * topStripDirection
  // m_r0 = bottomStripCenterPosition
  // m_dr = stripCenterDistance
	
	// cross product between top strip vector and spacepointPosition
	double d1[3] = {
		(sp->topHalfStripLength() * sp->topStripDirection()[1]) * d[2] -
		(sp->topHalfStripLength() * sp->topStripDirection()[2]) * d[1],
		(sp->topHalfStripLength() * sp->topStripDirection()[2]) * d[0] -
		(sp->topHalfStripLength() * sp->topStripDirection()[0]) * d[2],
		(sp->topHalfStripLength() * sp->topStripDirection()[0]) * d[1] -
		(sp->topHalfStripLength() * sp->topStripDirection()[1]) * d[0]};

	std::cout << "|check strip coord| " << (sp->topHalfStripLength() * sp->topStripDirection()[0]) << " " << (sp->topHalfStripLength() * sp->topStripDirection()[1]) << " " << (sp->topHalfStripLength() * sp->topStripDirection()[2]) << std::endl;

	std::cout << "|check strip coord d1: | " << d[0] << " " << d[1] << " " << d[2] << std::endl;
	
	std::cout << "|check strip coord d1: | " << d1[0] << " " << d1[1] << " " << d1[2] << std::endl;
	
	// scalar product between bottom strip vector and d1
	double bd1 =
	(sp->bottomHalfStripLength() * sp->bottomStripDirection()[0]) * d1[0] +
	(sp->bottomHalfStripLength() * sp->bottomStripDirection()[1]) * d1[1] +
	(sp->bottomHalfStripLength() * sp->bottomStripDirection()[2]) * d1[2];
	// if vectors are perpendicular, spacepointPosition is not compatible with
	// strip directions
	//  if (bd1 == 0.)
	//    return false;

	// compatibility check using distance between strips to evaluate
	// spacepointPosition is inside the bottom detector element
	double s1 = (sp->stripCenterDistance()[0] * d1[0] +
							 sp->stripCenterDistance()[1] * d1[1] +
							 sp->stripCenterDistance()[2] * d1[2]);
	if (std::abs(s1) > std::abs(bd1) * 1.1)  // -> define tolerance parameter
		return false;
	
	std::cout << "|check strip coord s1 bd1: | " << s1 << " " << bd1 << std::endl;
	
	
	// cross product between bottom strip vector and spacepointPosition
  double d0[3] = {
      (sp->bottomHalfStripLength() * sp->bottomStripDirection()[1]) * d[2] -
          (sp->bottomHalfStripLength() * sp->bottomStripDirection()[2]) * d[1],
      (sp->bottomHalfStripLength() * sp->bottomStripDirection()[2]) * d[0] -
          (sp->bottomHalfStripLength() * sp->bottomStripDirection()[0]) * d[2],
      (sp->bottomHalfStripLength() * sp->bottomStripDirection()[0]) * d[1] -
          (sp->bottomHalfStripLength() * sp->bottomStripDirection()[1]) * d[0]};
	
	std::cout << "|check strip coord d0: | " << d0[0] << " " << d0[1] << " " << d0[2] << std::endl;

  // scalar product between top strip vector and d0
  double bd0 = (sp->topHalfStripLength() * sp->topStripDirection()[0]) * d0[0] +
               (sp->topHalfStripLength() * sp->topStripDirection()[1]) * d0[1] +
               (sp->topHalfStripLength() * sp->topStripDirection()[2]) * d0[2];
  // if vectors are perpendicular, spacepointPosition is not compatible with
  // strip directions
  //	if (bd0 == 0.)
  //    return false;

  // compatibility check using distance between strips to evaluate
  // spacepointPosition is inside the top detector element
  double s0 = (sp->stripCenterDistance()[0] * d0[0] +
                sp->stripCenterDistance()[1] * d0[1] +
                sp->stripCenterDistance()[2] * d0[2]);
	
//	double s0 = (sp->stripCenterDistance()[0] * ((sp->bottomHalfStripLength() * sp->bottomStripDirection()[1]) * d[2] -
//							 (sp->bottomHalfStripLength() * sp->bottomStripDirection()[2]) * d[1]) +
//							 sp->stripCenterDistance()[1] * ((sp->bottomHalfStripLength() * sp->bottomStripDirection()[2]) * d[0] -
//							 (sp->bottomHalfStripLength() * sp->bottomStripDirection()[0]) * d[2]) +
//							 sp->stripCenterDistance()[2] * ((sp->bottomHalfStripLength() * sp->bottomStripDirection()[0]) * d[1] -
//							 (sp->bottomHalfStripLength() * sp->bottomStripDirection()[1]) * d[0]));
	
  if (std::abs(s0) > std::abs(bd1) * 1.1)
    return false;
	
	std::cout << "|check strip coord s0 bd0: | " << s0 << " " << bd1 << std::endl;

  s0 = s0 / bd1;
  // if arive here spacepointPosition is compatible with strip directions and
  // detector elements

	std::cout << "|check strip coord s0/bd0 | " << s0 << std::endl;
	
  // spacepointPosition corected with respect to the top strip direction and the
  // distance between the strips
  rTest[0] = sp->bottomStripCenterPosition()[0] +
             (sp->topHalfStripLength() * sp->topStripDirection()[0]) * s0;
  rTest[1] = sp->bottomStripCenterPosition()[1] +
             (sp->topHalfStripLength() * sp->topStripDirection()[1]) * s0;
  rTest[2] = sp->bottomStripCenterPosition()[2] +
             (sp->topHalfStripLength() * sp->topStripDirection()[2]) * s0;
  return true;
}

}  // namespace Acts
