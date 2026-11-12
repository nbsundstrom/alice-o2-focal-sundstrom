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

#include <algorithm>
#include <iterator>
#include <set>
#include <string_view>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <fairlogger/Logger.h>

#include "FOCALBase/PixelGeometryTGeo.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "MathUtils/Utils.h"

using namespace o2::focal;
using Segmentation = o2::itsmft::SegmentationAlpide;

std::unique_ptr<PixelGeometryTGeo> PixelGeometryTGeo::sInstance;

PixelGeometryTGeo* PixelGeometryTGeo::Instance()
{
  if (!sInstance) {
    sInstance = std::make_unique<PixelGeometryTGeo>();
    sInstance->Build(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  }
  return sInstance.get();
}

void PixelGeometryTGeo::Build(int loadTrans)
{
  if (isBuilt()) {
    LOG(warning) << "FOCAL pixel geometry is already built";
    return;
  }
  if (!gGeoManager) {
    LOG(fatal) << "Geometry is not loaded";
  }

  // the ECAL is the common mother of all towers
  gGeoManager->PushPath();
  if (!gGeoManager->cd(ECALPath)) {
    LOG(fatal) << "FOCAL ECAL volume " << ECALPath << " not found in the geometry";
  }
  const TGeoHMatrix ecalToGlobal(*gGeoManager->GetCurrentMatrix());
  TGeoVolume* ecal = gGeoManager->GetCurrentVolume();
  gGeoManager->PopPath();

  // As in ITS and MFT (GeometryTGeo::extractMatrixSensor), the chip frame is centred on the effective
  // sensitive layer, which shares the collection surface with the physical one
  const TGeoTranslation toEffectiveLayer(0., 0.5 * (Segmentation::SensorLayerThickness - Segmentation::SensorLayerThicknessEff), 0.);

  mChips.clear();
  TGeoIterator iter(ecal);
  while (auto node = iter()) {
    if (std::string_view(node->GetVolume()->GetName()) != SensorVolumeName) {
      continue;
    }
    const int level = iter.GetLevel();
    if (level != SensorLevel) {
      LOG(fatal) << "FOCAL pixel sensor found at depth " << level << " below the ECAL, expected tower/patch/string/chip/sensor";
    }
    Chip chip;
    chip.key = {iter.GetNode(level - 3)->GetNumber(), iter.GetNode(level - 4)->GetNumber(),
                iter.GetNode(level - 2)->GetNumber(), iter.GetNode(level - 1)->GetNumber()};
    chip.localToGlobal = ecalToGlobal;
    chip.localToGlobal.Multiply(iter.GetCurrentMatrix());
    chip.localToGlobal.Multiply(&toEffectiveLayer);
    mChips.push_back(chip);
  }
  if (mChips.empty()) {
    LOG(fatal) << "No FOCAL pixel chips (" << SensorVolumeName << ") found in the geometry";
  }
  std::sort(mChips.begin(), mChips.end(), [](const Chip& a, const Chip& b) { return a.key < b.key; });

  // dense layer index from the patch copy numbers: FoCal layers 4 and 9 become pixel layers 0 and 1
  std::set<int> layerCopies;
  for (const auto& chip : mChips) {
    layerCopies.insert(std::get<0>(chip.key));
  }
  mNumberOfLayers = layerCopies.size();

  mChipIndex.clear();
  for (int id = 0; id < static_cast<int>(mChips.size()); id++) {
    auto& chip = mChips[id];
    chip.layer = std::distance(layerCopies.begin(), layerCopies.find(std::get<0>(chip.key)));
    if (!mChipIndex.emplace(chip.key, id).second) {
      const auto& [layerCopy, towerCopy, stringCopy, chipCopy] = chip.key;
      LOG(fatal) << "Duplicate FOCAL pixel chip: layer " << layerCopy << ", tower " << towerCopy << ", string " << stringCopy << ", chip " << chipCopy;
    }
  }
  setSize(mChips.size());
  LOG(info) << "FOCAL pixel geometry: " << getSize() << " chips in " << mNumberOfLayers << " layers";
  fillMatrixCache(loadTrans);
}

void PixelGeometryTGeo::fillMatrixCache(int mask)
{
  if (mSize < 1) {
    LOG(warning) << "The method Build was not called yet";
    Build(mask);
    return;
  }
  // only local-to-global matrices are used by the digitizer
  if ((mask & o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G)) && !getCacheL2G().isFilled()) {
    auto& cacheL2G = getCacheL2G();
    cacheL2G.setSize(mSize);
    for (int chip = 0; chip < mSize; chip++) {
      cacheL2G.setMatrix(Mat3D(mChips[chip].localToGlobal), chip);
    }
  }
}

int PixelGeometryTGeo::getChipID(int towerCopy, int layerCopy, int stringCopy, int chipCopy) const
{
  const auto found = mChipIndex.find(ChipKey{layerCopy, towerCopy, stringCopy, chipCopy});
  return found == mChipIndex.end() ? -1 : found->second;
}
