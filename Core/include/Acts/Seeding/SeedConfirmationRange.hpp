// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>

namespace Acts {
/// @brief contains parameters for  seed confirmation
struct SeedConfirmationRange {
  float zMinSeedConf =
      std::numeric_limits<float>::min() * Acts::UnitConstants::mm;
  float zMaxSeedConf =
      std::numeric_limits<float>::max() * Acts::UnitConstants::mm;
  float rMaxSeedConf =
      std::numeric_limits<float>::max() * Acts::UnitConstants::mm;
  size_t nTopForLargeR = 0;
  size_t nTopForSmallR = 0;
	
	SeedConfirmationRange(float zMin = std::numeric_limits<float>::min(),
												float zMax = std::numeric_limits<float>::max(),
												float rMax = std::numeric_limits<float>::max(),
												size_t nTopLargeR = 0, size_t nTopSmallR = 0)
	: zMinSeedConf(zMin),
	zMaxSeedConf(zMax),
	rMaxSeedConf(rMax),
	nTopForLargeR(nTopLargeR),
	nTopForSmallR(nTopSmallR) {}
};
}  // namespace Acts
