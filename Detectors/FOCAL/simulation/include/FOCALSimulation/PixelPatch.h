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

#ifndef ALICEO2_FOCAL_PIXELPATCH_H_
#define ALICEO2_FOCAL_PIXELPATCH_H_

#include <array>
#include <string>
#include <vector>

class TGeoMedium;
class TGeoVolume;

namespace o2::focal
{
/// \class PixelPatch
/// \brief Builds the patch of one tower in a FoCal-E pixel layer from physical ALPIDE chips
///
/// ALPIDE chip (o2::itsmft::AlpideChip, built once)
///   -> string (chips butted along x; type A: 15 chips, type B: 12 chips for the beam-pipe row)
///     -> patch (six strings on two carriers in an air envelope; one per tower and pixel layer)
///
/// Patch frame: x along the strings, y transverse, z along the beam, origin at the envelope centre.
/// Same pattern as ITS (V3Layer: chip -> module -> stave) and MFT (Ladder).
class PixelPatch
{
 public:
  enum class StringType { A,
                          B };
  /// Tower column. The beam is at +x of a left patch and at -x of a right patch.
  enum class Side { Left,
                    Right };

  /// \param sizeX, sizeY, sizeZ full size of the patch envelope (cm)
  /// \param beamSideMargin distance from the beam-side edge of the envelope to the first chip (cm)
  /// \param air, carrier media of the envelope and of the carriers
  PixelPatch(double sizeX, double sizeY, double sizeZ, double beamSideMargin, TGeoMedium* air, TGeoMedium* carrier);

  /// \brief Envelope with the carriers and six strings of the given type, to be placed in a tower by the caller
  TGeoVolume* createPatch(Side side, StringType type);

  /// \brief Names of the volumes built so far, for the caller to register (see Detector::InitializeO2Detector)
  const std::vector<std::string>& getVolumeNames() const { return mVolumeNames; }

 private:
  TGeoVolume* createChip();
  TGeoVolume* createString(StringType type);

  std::array<double, 3> mSize; ///< full size of the patch envelope (cm)
  double mBeamSideMargin;      ///< distance from the beam-side edge of the envelope to the first chip (cm)
  TGeoMedium* mAir = nullptr;
  TGeoVolume* mChip = nullptr;                              ///< shared by all strings
  TGeoVolume* mCarrier = nullptr;                           ///< shared by all patches
  std::array<TGeoVolume*, 2> mStrings = {nullptr, nullptr}; ///< type A, type B; shared by all patches
  std::vector<std::string> mVolumeNames;                    ///< real (non-assembly) volumes built here
};
} // namespace o2::focal

#endif // ALICEO2_FOCAL_PIXELPATCH_H_
