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
#ifndef ALICEO2_FOCAL_PIXELGEOMETRYTGEO_H
#define ALICEO2_FOCAL_PIXELGEOMETRYTGEO_H

#include <map>
#include <memory>
#include <tuple>
#include <vector>
#include <TGeoMatrix.h>
#include "ITSMFTBase/GeometryTGeo.h"

namespace o2::focal
{

/// \class PixelGeometryTGeo
/// \brief ALPIDE chip view of the FoCal-E pixel layers for the ITSMFT digitizer
///
/// One chip per placement of the sensor volume built by o2::itsmft::AlpideChip, found below the ECAL as
/// tower / patch / string / chip / sensor. The chip frame is the sensor volume frame, as in ITS and MFT
/// (x: rows, y: thickness with +y on the collection side, z: columns), so o2::itsmft::Digitizer can be used as is.
/// Chip IDs are dense indices sorted by (layer, tower, string, chip); they are not readout addresses.
class PixelGeometryTGeo : public o2::itsmft::GeometryTGeo
{
 public:
  static constexpr const char* ChipVolumeName = "FOCALPixelChip";            ///< ALPIDE chip
  static constexpr const char* SensorVolumeName = "FOCALPixelSensor";        ///< sensitive layer of the chip
  static constexpr const char* ECALPath = "/cave_1/barrel_1/FOCAL_1/ECAL_1"; ///< common mother of the towers
  /// Depth of a sensor below the ECAL: tower / patch / string / chip / sensor
  static constexpr int SensorLevel = 5;

  PixelGeometryTGeo() : o2::itsmft::GeometryTGeo(o2::detectors::DetID::FOC) {}
  ~PixelGeometryTGeo() override = default;

  /// \brief Singleton built from gGeoManager on first access
  static PixelGeometryTGeo* Instance();

  /// \brief Collect the chips from gGeoManager and fill the requested matrix caches
  void Build(int loadTrans) override;
  void fillMatrixCache(int mask) override;

  /// \brief Pixel layer (0, 1, ...) of a chip
  int getLayer(int chipID) const override { return mChips[chipID].layer; }
  int getNumberOfPixelLayers() const { return mNumberOfLayers; }

  /// \brief Chip ID from the copy numbers of the tower, the patch (= FoCal layer), the string and the chip, -1 if unknown
  int getChipID(int towerCopy, int layerCopy, int stringCopy, int chipCopy) const;

 private:
  using ChipKey = std::tuple<int, int, int, int>; ///< (layer, tower, string, chip) copy numbers: sort order of the chip IDs

  struct Chip {
    ChipKey key;
    int layer = -1;            ///< dense pixel layer index
    TGeoHMatrix localToGlobal; ///< chip frame to global
  };

  std::vector<Chip> mChips;
  std::map<ChipKey, int> mChipIndex; ///< key -> chip ID
  int mNumberOfLayers = 0;

  static std::unique_ptr<PixelGeometryTGeo> sInstance;
};

} // namespace o2::focal

#endif
