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
#ifndef ALICEO2_FOCAL_DIGITIZERPARAM_H
#define ALICEO2_FOCAL_DIGITIZERPARAM_H

#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2::focal
{

/// \struct FOCALDigitizerParam
/// \brief Configuration of the FoCal-E pixel digitizer (o2::itsmft::Digitizer with ALPIDE response)
///
/// Merges what the ITS/MFT digitizer reads from DPLAlpideParam and DPLDigitizerParam,
/// since these templates are restricted to ITS and MFT. Defaults are the ITS ones.
/// Keys are set with --configKeyValues "FOCALDigitizerParam.<member>=<value>".
struct FOCALDigitizerParam : public o2::conf::ConfigurableParamHelper<FOCALDigitizerParam> {
  // readout frame (as DPLAlpideParam)
  bool continuous = true;                                            ///< flag for continuous simulation
  int roFrameLengthInBC = o2::constants::lhc::LHCMaxBunches / 4;     ///< ROF length in BC for continuous mode
  float roFrameLengthTrig = 6000.;                                   ///< length of RO frame in ns for triggered mode
  float strobeDelay = o2::constants::lhc::LHCBunchSpacingNS * 4;     ///< strobe start (in ns) wrt ROF start
  float strobeLengthCont = -1.;                                      ///< if < 0, full ROF length - delay
  float strobeLengthTrig = 100.;                                     ///< length of the strobe in ns (sig. over threshold checked in this window only)
  int roFrameBiasInBC = 0;                                          ///< bias of the start of ROF wrt orbit start (in BC)
  // TODO: DETERMINE WHAT roFrameBiasInBC should be, temp set to 0

  // digitization (as DPLDigitizerParam)
  float noisePerPixel = 1e-8;                                        ///< ALPIDE noise per channel
  float strobeFlatTop = 7500.;                                       ///< strobe shape flat top
  float strobeMaxRiseTime = 1100.;                                   ///< strobe shape max rise time
  float strobeQRiseTime0 = 450.;                                     ///< q @ which strobe rise time is 0
  double timeOffset = 0.;                                            ///< time offset (in seconds!) to calculate ROFrame from hit time
  int chargeThreshold = 150;                                         ///< charge threshold in Nelectrons
  int minChargeToAccount = 15;                                       ///< minimum charge contribution to account (in Nelectrons)
  int nSimSteps = 7;                                                 ///< number of steps in response simulation
  float energyToNElectrons = 1. / 3.6e-9;                            ///< conversion of eloss (in GeV) to Nelectrons
  float Vbb = 0.0;                                                   ///< back bias absolute value (in Volt): 0 or 3 (for -3 V)

  O2ParamDef(FOCALDigitizerParam, "FOCALDigitizerParam");
};

} // namespace o2::focal

#endif
