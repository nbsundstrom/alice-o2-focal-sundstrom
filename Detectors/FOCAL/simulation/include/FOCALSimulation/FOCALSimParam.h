// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef ALICEO2_FOCAL_SIMPARAM_H
#define ALICEO2_FOCAL_SIMPARAM_H

#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2::focal
{
struct FOCALSimParam : o2::conf::ConfigurableParamHelper<FOCALSimParam> {
  // Store one entry/exit segment per pixel-chip crossing in FOCPixelHit (input of the pixel digitizer)
  bool recordPixelHits = true;
  O2ParamDef(FOCALSimParam, "FOCALSimParam");
};
} // namespace o2::focal
#endif
