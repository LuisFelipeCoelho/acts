// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {
template <typename external_spacepoint_t>
inline LinCircle transformCoordinates(
    const InternalSpacePoint<external_spacepoint_t>& sp,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom) {
  auto extractFunction =
      [](const InternalSpacePoint<external_spacepoint_t>& obj)
      -> std::array<float, 6> {
    std::array<float, 6> output{obj.x(),      obj.y(),         obj.z(),
                                obj.radius(), obj.varianceR(), obj.varianceZ()};
    return output;
  };

  return transformCoordinates<InternalSpacePoint<external_spacepoint_t>>(
      sp, spM, bottom, std::move(extractFunction));
}

template <typename external_spacepoint_t, typename callable_t>
inline LinCircle transformCoordinates(const external_spacepoint_t& sp,
                                      const external_spacepoint_t& spM,
                                      bool bottom,
                                      callable_t&& extractFunction) {
  // The computation inside this function is exactly identical to that in the
  // vectorized version of this function, except that it operates on a single
  // spacepoint. Please see the other version of this function for more
  // detailed comments.

  auto [xM, yM, zM, rM, varianceRM, varianceZM] = extractFunction(spM);
  auto [xSP, ySP, zSP, rSP, varianceRSP, varianceZSP] = extractFunction(sp);

  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
  float deltaX = xSP - xM;
  float deltaY = ySP - yM;
  float deltaZ = zSP - zM;
  float x = deltaX * cosPhiM + deltaY * sinPhiM;
  float y = deltaY * cosPhiM - deltaX * sinPhiM;
  float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
  float iDeltaR = std::sqrt(iDeltaR2);
  int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
  float cot_theta = deltaZ * iDeltaR * bottomFactor;
  LinCircle l{};
  l.cotTheta = cot_theta;
  l.iDeltaR = iDeltaR;
  l.U = x * iDeltaR2;
  l.V = y * iDeltaR2;
  l.Er = ((varianceZM + varianceZSP) +
          (cot_theta * cot_theta) * (varianceRM + varianceRSP)) *
         iDeltaR2;
  return l;
}

inline LinCircle fillLineCircle(
    const std::array<float, 7>& lineCircleVariables) {
  auto [cotTheta, iDeltaR, Er, U, V, xNewFrame, yNewFrame] =
      lineCircleVariables;

  LinCircle l{};
  l.cotTheta = cotTheta;
  l.iDeltaR = iDeltaR;
  l.U = U;
  l.V = V;
  l.Er = Er;
  l.x = xNewFrame;
  l.y = yNewFrame;

  return l;
}

inline LinCircle fillLineCircleDetailed(
    const std::array<float, 7>& lineCircleVariables) {
  auto [cotTheta, iDeltaR, Er, U, V, xNewFrame, yNewFrame] =
      lineCircleVariables;

  LinCircle l{};
  l.cotTheta = cotTheta;
  l.iDeltaR = iDeltaR;
  l.U = U;
  l.V = V;
  l.Er = Er;
  l.x = xNewFrame;
  l.y = yNewFrame;

  return l;
}

template <typename external_spacepoint_t>
inline void transformCoordinates(
    Acts::SpacePointData& spacePointData,
    const std::vector<InternalSpacePoint<external_spacepoint_t>*>& vec,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    std::vector<LinCircle>& linCircleVec) {
  auto extractFunction =
      [](const InternalSpacePoint<external_spacepoint_t>& obj)
      -> std::array<float, 6> {
    std::array<float, 6> output{obj.x(),      obj.y(),         obj.z(),
                                obj.radius(), obj.varianceR(), obj.varianceZ()};
    return output;
  };

  transformCoordinates<InternalSpacePoint<external_spacepoint_t>>(
      spacePointData, vec, spM, bottom, linCircleVec,
      std::move(extractFunction));
}

template <typename external_spacepoint_t, typename callable_t>
inline void transformCoordinates(Acts::SpacePointData& spacePointData,
                                 const std::vector<external_spacepoint_t*>& vec,
                                 const external_spacepoint_t& spM, bool bottom,
                                 std::vector<LinCircle>& linCircleVec,
                                 callable_t&& extractFunction) {
  std::vector<std::size_t> indexes(vec.size());
  for (unsigned int i(0); i < indexes.size(); i++) {
    indexes[i] = i;
  }

  auto [xM, yM, zM, rM, varianceRM, varianceZM] = extractFunction(spM);

  // resize + operator[] is faster then reserve and push_back
  linCircleVec.resize(vec.size());

  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;

  for (std::size_t idx(0); idx < vec.size(); ++idx) {
    auto& sp = vec[idx];
    auto [xSP, ySP, zSP, rSP, varianceRSP, varianceZSP] = extractFunction(*sp);

    float deltaX = xSP - xM;
    float deltaY = ySP - yM;
    float deltaZ = zSP - zM;
    // calculate projection fraction of spM->sp vector pointing in same
    // direction as
    // vector origin->spM (x) and projection fraction of spM->sp vector pointing
    // orthogonal to origin->spM (y)
    float x = deltaX * cosPhiM + deltaY * sinPhiM;
    float y = deltaY * cosPhiM - deltaX * sinPhiM;
    // 1/(length of M -> SP)
    //
    float deltaR2 = (deltaX * deltaX + deltaY * deltaY);
    float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
    float iDeltaR = std::sqrt(iDeltaR2);
    //
    int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
    // cot_theta = (deltaZ/deltaR)
    float cot_theta = deltaZ * iDeltaR * bottomFactor;
    // VERY frequent (SP^3) access
    // LinCircle l{};
    // l.cotTheta = cot_theta;
    // l.iDeltaR = iDeltaR;
    // l.U = x * iDeltaR2;
    // l.V = y * iDeltaR2;
    // l.Er = ((varianceZM + sp->varianceZ()) +
    //        (cot_theta * cot_theta) * (varianceRM + sp->varianceR())) *
    //       iDeltaR2;

    // l.x = x;
    // l.y = y;

    // std::cout << "parm2 " << cot_theta << " " << iDeltaR << " " <<
    // ((varianceZM + sp->varianceZ()) +(cot_theta * cot_theta) * (varianceRM +
    // sp->varianceR())) *iDeltaR2 << " " << x * iDeltaR2 << " " << y * iDeltaR2
    // << " " << x << " " << y << std::endl;

    linCircleVec[idx] = fillLineCircle(
        {cot_theta, iDeltaR,
         ((varianceZM + sp->varianceZ()) +
          (cot_theta * cot_theta) * (varianceRM + sp->varianceR())) *
             iDeltaR2,
         x * iDeltaR2, y * iDeltaR2, x, y});
    // spacePointData.setDeltaR(sp->index(),
    //                         std::sqrt(deltaR2 + (deltaZ * deltaZ)));
  }
}

template <typename external_spacepoint_t>
inline bool xyzCoordinateCheck(
    Acts::SpacePointData& spacePointData,
    const Acts::SeedFinderConfig<external_spacepoint_t>& m_config,
    const Acts::InternalSpacePoint<external_spacepoint_t>& sp,
    const double* spacepointPosition, double* outputCoordinates) {
  // check the compatibility of SPs coordinates in xyz assuming the
  // Bottom-Middle direction with the strip measurement details
  bool hasValueStored = spacePointData.hasDynamicVariable();
  if (not hasValueStored) {
    return false;
  }

  std::size_t index = sp.index();

  // prepare variables
  const Acts::Vector3& topStripVector = spacePointData.getTopStripVector(index);
  const Acts::Vector3& bottomStripVector =
      spacePointData.getBottomStripVector(index);
  const Acts::Vector3& stripCenterDistance =
      spacePointData.getStripCenterDistance(index);

  const double& xTopStripVector = topStripVector[0];
  const double& yTopStripVector = topStripVector[1];
  const double& zTopStripVector = topStripVector[2];
  const double& xBottomStripVector = bottomStripVector[0];
  const double& yBottomStripVector = bottomStripVector[1];
  const double& zBottomStripVector = bottomStripVector[2];

  // cross product between top strip vector and spacepointPosition
  double d1[3] = {yTopStripVector * spacepointPosition[2] -
                      zTopStripVector * spacepointPosition[1],
                  zTopStripVector * spacepointPosition[0] -
                      xTopStripVector * spacepointPosition[2],
                  xTopStripVector * spacepointPosition[1] -
                      yTopStripVector * spacepointPosition[0]};

  // scalar product between bottom strip vector and d1
  double bd1 = xBottomStripVector * d1[0] + yBottomStripVector * d1[1] +
               zBottomStripVector * d1[2];

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  double s1 = (stripCenterDistance[0] * d1[0] + stripCenterDistance[1] * d1[1] +
               stripCenterDistance[2] * d1[2]);
  if (std::abs(s1) > std::abs(bd1) * m_config.toleranceParam) {
    return false;
  }

  // cross product between bottom strip vector and spacepointPosition
  double d0[3] = {yBottomStripVector * spacepointPosition[2] -
                      zBottomStripVector * spacepointPosition[1],
                  zBottomStripVector * spacepointPosition[0] -
                      xBottomStripVector * spacepointPosition[2],
                  xBottomStripVector * spacepointPosition[1] -
                      yBottomStripVector * spacepointPosition[0]};

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the top detector element
  double s0 = (stripCenterDistance[0] * d0[0] + stripCenterDistance[1] * d0[1] +
               stripCenterDistance[2] * d0[2]);
  if (std::abs(s0) > std::abs(bd1) * m_config.toleranceParam) {
    return false;
  }

  // if arive here spacepointPosition is compatible with strip directions and
  // detector elements

  const Acts::Vector3& topStripCenterPosition =
      spacePointData.getTopStripCenterPosition(index);

  // spacepointPosition corected with respect to the top strip position and
  // direction and the distance between the strips
  s0 = s0 / bd1;
  outputCoordinates[0] = topStripCenterPosition[0] + xTopStripVector * s0;
  outputCoordinates[1] = topStripCenterPosition[1] + yTopStripVector * s0;
  outputCoordinates[2] = topStripCenterPosition[2] + zTopStripVector * s0;
  return true;
}
}  // namespace Acts
