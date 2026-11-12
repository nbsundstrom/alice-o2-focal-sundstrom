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

#ifndef STEER_DIGITIZERWORKFLOW_FOCALDIGITWRITER_H_
#define STEER_DIGITIZERWORKFLOW_FOCALDIGITWRITER_H_

#include <vector>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include "Framework/DataProcessorSpec.h"
#include "Framework/InputSpec.h"
#include "Framework/Logger.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "Headers/DataHeader.h"
#include "DataFormatsITSMFT/Digit.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "SimulationDataFormat/IOMCTruthContainerView.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2
{
namespace focal
{

template <typename T>
using BranchDefinition = framework::MakeRootTreeWriterSpec::BranchDefinition<T>;

/// Writer of the FoCal-E pixel digits (itsmft::Digit with ROF records and MC labels) to focaldigits.root
inline o2::framework::DataProcessorSpec getFOCALDigitWriterSpec(bool mctruth = true)
{
  using InputSpec = framework::InputSpec;
  using MakeRootTreeWriterSpec = framework::MakeRootTreeWriterSpec;
  auto detOrig = o2::header::gDataOriginFOC;

  auto logger = [](std::vector<o2::itsmft::Digit> const& inDigits) {
    LOG(info) << "RECEIVED DIGITS SIZE " << inDigits.size();
  };

  auto finishWriting = [](TFile* outputfile, TTree* outputtree) {
    const auto* brArr = outputtree->GetListOfBranches();
    int64_t nent = 0;
    for (const auto* brc : *brArr) {
      int64_t n = ((const TBranch*)brc)->GetEntries();
      if (nent && (nent != n)) {
        LOG(error) << "Branches have different number of entries";
      }
      nent = n;
    }
    outputtree->SetEntries(nent);
    outputfile->Write("", TObject::kOverwrite);
    outputfile->Close();
  };

  // Split label buffers that may exceed ROOT's entry size limit.
  auto fillLabels = [](TBranch& branch, std::vector<char> const& labelbuffer, framework::DataRef const& /*ref*/) {
    o2::dataformats::ConstMCTruthContainerView<o2::MCCompLabel> labels(labelbuffer);
    LOG(info) << "WRITING " << labels.getNElements() << " LABELS ";

    o2::dataformats::IOMCTruthContainerView outputcontainer;
    auto ptr = &outputcontainer;
    auto br = framework::RootTreeWriter::remapBranch(branch, &ptr);
    outputcontainer.adopt(labelbuffer);
    br->Fill();
    br->ResetAddress();
  };

  return MakeRootTreeWriterSpec("FOCALDigitWriter",
                                "focaldigits.root",
                                MakeRootTreeWriterSpec::TreeAttributes{"o2sim", "FOCAL Digits tree"},
                                MakeRootTreeWriterSpec::CustomClose(finishWriting),
                                BranchDefinition<std::vector<char>>{InputSpec{"digitsMCTR", detOrig, "DIGITSMCTR", 0},
                                                                    "FOCDigitMCTruth",
                                                                    (mctruth ? 1 : 0), fillLabels},
                                BranchDefinition<std::vector<itsmft::MC2ROFRecord>>{InputSpec{"digitsMC2ROF", detOrig, "DIGITSMC2ROF", 0},
                                                                                    "FOCDigitMC2ROF",
                                                                                    (mctruth ? 1 : 0)},
                                BranchDefinition<std::vector<itsmft::Digit>>{InputSpec{"digits", detOrig, "DIGITS", 0},
                                                                             "FOCDigit",
                                                                             logger},
                                BranchDefinition<std::vector<itsmft::ROFRecord>>{InputSpec{"digitsROF", detOrig, "DIGITSROF", 0},
                                                                                 "FOCDigitROF"})();
}

} // namespace focal
} // namespace o2

#endif /* STEER_DIGITIZERWORKFLOW_FOCALDIGITWRITER_H_ */
