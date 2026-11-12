// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PixelPatch.cxx
/// \brief Definition of the FOCAL PixelPatch class
/// \author nick.sundstrom@cern.ch
///
/// Geometry hierarchy built here (smallest -> largest):
///
///   ALPIDE chip  (one TGeoVolume, built once, reused everywhere)
///     -> String  (row of N chips; two flavours: type A / type B)
///       -> Patch (six strings on two carriers; one per tower and pixel layer)
///
/// Reference implementation with the same pattern:
///   Detectors/ITSMFT/ITS/simulation/src/V3Layer.cxx  (chip -> module -> stave -> layer)

#include "FOCALSimulation/PixelPatch.h"
#include "FOCALBase/PixelGeometryTGeo.h"

#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITSMFTSimulation/AlpideChip.h"

#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <fairlogger/Logger.h>

#include <cmath>

using namespace o2::focal;

namespace
{
using Segmentation = o2::itsmft::SegmentationAlpide;

// FoCal-E pixel patch layout. Units: cm. Source: <drawing number, revision, date>.
constexpr double ChipThickness = 100.e-4;                              // ALPIDE chip (November 2024 drawing, p. 6), as in the ITS outer barrel
constexpr double SensorThickness = Segmentation::SensorLayerThickness; // sensitive layer assumed by the ITSMFT digitizer
constexpr double ChipPitchX = Segmentation::SensorSizeCols + 0.;       // TODO(drawing): gap between the chips of a string
constexpr int ChipsPerString[2] = {15, 12};                            // type A, type B (beam-pipe row)
constexpr int StringsPerPatch = 6;
constexpr double CarrierThickness = 0.1;                            // two 1 mm aluminium carriers (November 2024 drawing, p. 6)
constexpr double CarrierZ[2] = {-0.1, 0.1};                         // TODO(drawing): carrier centres
constexpr double StringPitchY = Segmentation::ActiveMatrixSizeRows; // TODO(drawing): placeholder, active rows butted in projection
// chip centre on the outer face of a carrier
constexpr double StringZ = 0.5 * (CarrierZ[1] - CarrierZ[0] + CarrierThickness + ChipThickness);

/// Placement of a string in the patch frame
struct StringPlacement {
  double y;      ///< string centre (cm)
  double z;      ///< string centre (cm)
  bool backFace; ///< turned over about the string axis: collection side toward -z, row 0 toward -y
};

// TODO(drawing): placeholder that fits the envelope. Strings alternate between the outer faces of the two carriers.
constexpr std::array<StringPlacement, StringsPerPatch> StringPlacements = {{
  {-2.5 * StringPitchY, -StringZ, true},
  {-1.5 * StringPitchY, +StringZ, false},
  {-0.5 * StringPitchY, -StringZ, true},
  {+0.5 * StringPitchY, +StringZ, false},
  {+1.5 * StringPitchY, -StringZ, true},
  {+2.5 * StringPitchY, +StringZ, false},
}};
} // namespace

PixelPatch::PixelPatch(double sizeX, double sizeY, double sizeZ, double beamSideMargin, TGeoMedium* air, TGeoMedium* carrier)
  : mSize{sizeX, sizeY, sizeZ}, mBeamSideMargin(beamSideMargin), mAir(air)
{
  // ROOT's overlap check compares box corners, so it misses identical chips that are aligned in x and z.
  // Check directly that strings on the same face are at least one chip apart.
  for (int i = 0; i < StringsPerPatch; i++) {
    for (int j = i + 1; j < StringsPerPatch; j++) {
      if (StringPlacements[i].z == StringPlacements[j].z &&
          std::abs(StringPlacements[i].y - StringPlacements[j].y) < Segmentation::SensorSizeRows) {
        LOG(fatal) << "FOCAL pixel strings " << i << " and " << j << " overlap on the same carrier face";
      }
    }
  }

  mChip = createChip();
  mStrings = {createString(StringType::A), createString(StringType::B)};

  // TODO(drawing): carrier x/y size; flex and spacers. The 1.5 mm spacer regions include clearance and
  // components: do not model them as solid aluminium.
  mCarrier = gGeoManager->MakeBox("FOCALPixelCarrier", carrier, 0.5 * sizeX, 0.5 * sizeY, 0.5 * CarrierThickness);
  mCarrier->SetLineColor(kGray);
  mVolumeNames.emplace_back(mCarrier->GetName());
}

TGeoVolume* PixelPatch::createChip()
{
  // Chip with metal stack and sensitive layer, shared with ITS and MFT. The MetalStack volume inside is one
  // volume for all three detectors, so it is not listed in mVolumeNames.
  mVolumeNames.emplace_back(PixelGeometryTGeo::ChipVolumeName);
  mVolumeNames.emplace_back(PixelGeometryTGeo::SensorVolumeName);
  return o2::itsmft::AlpideChip::createChip(0.5 * ChipThickness, 0.5 * SensorThickness,
                                            PixelGeometryTGeo::ChipVolumeName, PixelGeometryTGeo::SensorVolumeName);
}

TGeoVolume* PixelPatch::createString(StringType type)
{
  const bool typeA = type == StringType::A;
  const int nChips = ChipsPerString[typeA ? 0 : 1];
  TGeoVolume* stringVolume = new TGeoVolumeAssembly(typeA ? "FOCALPixelStringA" : "FOCALPixelStringB");
  for (int iChip = 0; iChip < nChips; iChip++) {
    // chip rows -> string +y, chip thickness (collection side) -> string +z, chip columns -> string +x
    auto* rotation = new TGeoRotation();
    rotation->RotateX(90.);
    rotation->RotateZ(90.);
    stringVolume->AddNode(mChip, iChip, new TGeoCombiTrans((iChip - 0.5 * (nChips - 1)) * ChipPitchX, 0., 0., rotation));
  }
  // TODO(drawing): flex of the string
  return stringVolume;
}

TGeoVolume* PixelPatch::createPatch(Side side, StringType type)
{
  const bool left = side == Side::Left;
  const int typeIndex = type == StringType::A ? 0 : 1;
  const std::string name = std::string("FOCALPixelPatch") + (left ? "L" : "R") + (typeIndex == 0 ? "A" : "B");
  TGeoVolume* patch = gGeoManager->MakeBox(name.c_str(), mAir, 0.5 * mSize[0], 0.5 * mSize[1], 0.5 * mSize[2]);
  mVolumeNames.push_back(name);

  for (int iCarrier = 0; iCarrier < 2; iCarrier++) {
    patch->AddNode(mCarrier, iCarrier, new TGeoTranslation(0., 0., CarrierZ[iCarrier]));
  }

  // every string starts mBeamSideMargin from the beam-side edge: +x for a left patch, -x for a right patch
  const double stringLength = ChipsPerString[typeIndex] * ChipPitchX;
  const double x = (left ? 1. : -1.) * (0.5 * mSize[0] - mBeamSideMargin - 0.5 * stringLength);
  for (int iString = 0; iString < StringsPerPatch; iString++) {
    const auto& placement = StringPlacements[iString];
    auto* rotation = new TGeoRotation();
    if (placement.backFace) {
      rotation->RotateX(180.);
    }
    patch->AddNode(mStrings[typeIndex], iString, new TGeoCombiTrans(x, placement.y, placement.z, rotation));
  }
  LOG(debug) << name << ": " << StringsPerPatch << " strings of " << ChipsPerString[typeIndex] << " chips";
  return patch;
}
