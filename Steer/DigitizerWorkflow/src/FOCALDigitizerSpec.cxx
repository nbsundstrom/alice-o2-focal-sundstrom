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

#include "FOCALDigitizerSpec.h"
#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Task.h"
#include "Framework/CCDBParamSpec.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "DataFormatsITSMFT/Digit.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "DetectorsBase/BaseDPLDigitizer.h"
#include "DetectorsRaw/HBFUtils.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/SimTraits.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CommonConstants/LHCConstants.h"
#include "ITSMFTSimulation/Digitizer.h"
#include "ITSMFTSimulation/Hit.h"
#include "ITSMFTSimulation/AlpideSimResponse.h"
#include "FOCALBase/Hit.h"
#include "FOCALBase/PixelGeometryTGeo.h"
#include "FOCALSimulation/FOCALDigitizerParam.h"
#include <TChain.h>
#include <TStopwatch.h>
#include <TVector3.h>
#include <string>
#include <vector>

using namespace o2::framework;
using SubSpecificationType = o2::framework::DataAllocator::SubSpecificationType;

namespace o2::focal
{

using namespace o2::base;

class FOCALDPLDigitizerTask : BaseDPLDigitizer
{
 public:
  static constexpr o2::detectors::DetID ID{o2::detectors::DetID::FOC};
  static constexpr o2::header::DataOrigin Origin{o2::header::gDataOriginFOC};

  using BaseDPLDigitizer::init;

  FOCALDPLDigitizerTask(bool mctruth = true) : BaseDPLDigitizer(InitServices::FIELD | InitServices::GEOM), mWithMCTruth(mctruth) {}

  void initDigitizerTask(framework::InitContext& ic) override
  {
    mDisableQED = ic.options().get<bool>("disable-qed");
  }

  void run(framework::ProcessingContext& pc)
  {
    if (mFinished) {
      return;
    }

    const auto firstOrbitTF = pc.services().get<o2::framework::TimingInfo>().firstTForbit;
    const o2::InteractionRecord firstIR(0, firstOrbitTF);
    updateTimeDependentParams(pc);

    TStopwatch timer;
    timer.Start();
    LOG(info) << " CALLING FOCAL PIXEL DIGITIZATION ";

    auto context = pc.inputs().get<o2::steer::DigitizationContext*>("collisioncontext");
    context->initSimChains(ID, mSimChains);
    const bool withQED = context->isQEDProvided() && !mDisableQED;
    auto& timesview = context->getEventRecords(withQED);
    LOG(info) << "GOT " << timesview.size() << " COLLISSION TIMES";
    LOG(info) << "SIMCHAINS: " << mSimChains.size();

    if (timesview.empty()) {
      return;
    }

    // the chips of all pixel layers are digitized together (no staggering)
    mDigitizer.setDigits(&mDigits);
    mDigitizer.setROFRecords(&mROFRecords);
    mDigitizer.setMCLabels(&mLabels);
    mDigitizer.resetROFrameBounds();

    auto& digitsAccum = pc.outputs().make<std::vector<itsmft::Digit>>(Output{Origin, "DIGITS", 0});

    const auto& par = FOCALDigitizerParam::Instance();
    const int nROFsPerOrbit = o2::constants::lhc::LHCMaxBunches / par.roFrameLengthInBC;
    const int nROFsTF = nROFsPerOrbit * raw::HBFUtils::Instance().getNOrbitsPerTF();
    mROFRecordsAccum.reserve(nROFsTF);

    auto accumulate = [this, &digitsAccum]() {
      if (mDigits.empty()) {
        return;
      }
      auto ndigAcc = digitsAccum.size();
      digitsAccum.insert(digitsAccum.end(), mDigits.begin(), mDigits.end());

      auto nROFRecsOld = mROFRecordsAccum.size();

      for (int i = 0; i < mROFRecords.size(); i++) {
        auto& rof = mROFRecords[i];
        rof.setFirstEntry(ndigAcc + rof.getFirstEntry());
        rof.print();

        // Link each MC event to its first accumulated ROF record.
        for (int m2rid = mFixMC2ROF; m2rid < mMC2ROFRecordsAccum.size(); m2rid++) {
          auto& mc2rof = mMC2ROFRecordsAccum[m2rid];
          if (rof.getROFrame() == mc2rof.minROF) {
            mFixMC2ROF++;
            mc2rof.rofRecordID = nROFRecsOld + i;
            mc2rof.print();
          }
        }
      }

      mROFRecordsAccum.insert(mROFRecordsAccum.end(), mROFRecords.begin(), mROFRecords.end());
      if (mWithMCTruth) {
        mLabelsAccum.mergeAtBack(mLabels);
      }
      LOG(info) << "Added " << mDigits.size() << " digits";
      mLabels.clear();
      mDigits.clear();
      mROFRecords.clear();
    };

    const auto& eventParts = context->getEventParts(withQED);
    const int64_t bcShift = mDigitizer.getParams().getROFrameBiasInBC();
    for (int collID = 0; collID < timesview.size(); ++collID) {
      auto irt = timesview[collID];
      if (irt.toLong() < bcShift) { // due to the ROF misalignment the collision would go to negative ROF ID, discard
        continue;
      }
      irt -= bcShift; // account for the ROF start shift

      mDigitizer.setEventTime(irt);
      mDigitizer.resetEventROFrames(); // to estimate min/max ROF for this collID
      for (const auto& part : eventParts[collID]) {
        mFocalHits.clear();
        // FOCPixelHit: pixel-chip segments (entry/exit points) written by the FoCal simulation
        context->retrieveHits(mSimChains, o2::detectors::SimTraits::DETECTORBRANCHNAMES[ID][1].c_str(), part.sourceID, part.entryID, &mFocalHits);

        if (!mFocalHits.empty()) {
          LOG(debug) << "For collision " << collID << " eventID " << part.entryID << " found " << mFocalHits.size() << " hits ";
          convertHits();
          mDigitizer.process(&mHits, part.entryID, part.sourceID);
        }
      }
      mMC2ROFRecordsAccum.emplace_back(collID, -1, mDigitizer.getEventROFrameMin(), mDigitizer.getEventROFrameMax());
      accumulate();
    }
    mDigitizer.fillOutputContainer();
    accumulate();
    const uint64_t nDigits = digitsAccum.size();

    // Include empty ROFs so downstream consumers receive a continuous time structure starting at BC=0
    if (nROFsTF != mROFRecordsAccum.size()) {
      const size_t nROFsLayer = std::max((size_t)nROFsTF, mROFRecordsAccum.size());
      std::vector<o2::itsmft::ROFRecord> expDigitRofVec(nROFsLayer);
      for (int iROF{0}; iROF < nROFsLayer; ++iROF) {
        auto& rof = expDigitRofVec[iROF];
        int orb = iROF * par.roFrameLengthInBC / o2::constants::lhc::LHCMaxBunches + firstOrbitTF;
        int bc = iROF * par.roFrameLengthInBC % o2::constants::lhc::LHCMaxBunches;
        o2::InteractionRecord ir(bc, orb);
        rof.setBCData(ir);
        rof.setROFrame(iROF);
        rof.setNEntries(0);
        rof.setFirstEntry(-1);
      }
      for (const auto& rof : mROFRecordsAccum) {
        const auto& ir = rof.getBCData();
        const auto irToFirst = ir - firstIR;
        const int irROF = irToFirst.toLong() / par.roFrameLengthInBC;
        auto& expROF = expDigitRofVec[irROF];
        expROF.setFirstEntry(rof.getFirstEntry());
        expROF.setNEntries(rof.getNEntries());
        if (expROF.getBCData() != rof.getBCData()) {
          LOGP(fatal, "detected mismatch between expected {} and received {}", expROF.asString(), rof.asString());
        }
      }
      int prevFirst{0};
      for (auto& rof : expDigitRofVec) {
        if (rof.getFirstEntry() < 0) {
          rof.setFirstEntry(prevFirst);
        }
        prevFirst = rof.getFirstEntry();
      }
      // Drop ROFs that spilled into the next timeframe.
      expDigitRofVec.resize(nROFsTF);
      pc.outputs().snapshot(Output{Origin, "DIGITSROF", 0}, expDigitRofVec);
    } else {
      pc.outputs().snapshot(Output{Origin, "DIGITSROF", 0}, mROFRecordsAccum);
    }
    if (mWithMCTruth) {
      pc.outputs().snapshot(Output{Origin, "DIGITSMC2ROF", 0}, mMC2ROFRecordsAccum);
      auto& sharedlabels = pc.outputs().make<o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>>(Output{Origin, "DIGITSMCTR", 0});
      mLabelsAccum.flatten_to(sharedlabels);
      mLabels.clear_andfreememory();
      mLabelsAccum.clear_andfreememory();
    }

    LOG(info) << ID.getName() << ": Sending ROMode= " << mROMode << " to GRPUpdater";
    pc.outputs().snapshot(Output{Origin, "ROMode", 0}, mROMode);

    timer.Stop();
    LOG(info) << "Digitization took " << timer.CpuTime() << "s";
    LOG(info) << "Produced " << nDigits << " digits";

    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);

    mFinished = true;
  }

  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
  {
    // the response matching the requested Vbb is selected in updateTimeDependentParams
    if (matcher == ConcreteDataMatcher(Origin, "ALPIDERESPVbb0", 0)) {
      LOG(info) << ID.getName() << " loaded AlpideResponseData for Vbb=0V";
      mRespVbb0 = (o2::itsmft::AlpideSimResponse*)obj;
    }
    if (matcher == ConcreteDataMatcher(Origin, "ALPIDERESPVbbM3", 0)) {
      LOG(info) << ID.getName() << " loaded AlpideResponseData for Vbb=-3V";
      mRespVbbM3 = (o2::itsmft::AlpideSimResponse*)obj;
    }
  }

 protected:
  /// Convert the FoCal pixel hits of the current event to the ITSMFT hits expected by the digitizer
  void convertHits()
  {
    mHits.clear();
    mHits.reserve(mFocalHits.size());
    const int nChips = mGeometry->getNumberOfChips();
    // FoCal stores ns and eV; ITSMFT hits use seconds and GeV.
    constexpr double ns2sec = 1e-9;
    constexpr double eV2GeV = 1e-9;
    for (const auto& fhit : mFocalHits) {
      if (fhit.GetDetectorID() >= nChips) {
        LOG(warning) << "Skipping FoCal pixel hit with invalid chip ID " << fhit.GetDetectorID();
        continue;
      }
      const auto& pStart = fhit.getPosStart();
      const auto& pEnd = fhit.GetPos();
      // no momentum and track status information is available/needed by the digitizer
      mHits.emplace_back(fhit.GetTrackID(), fhit.GetDetectorID(),
                         TVector3(pStart.X(), pStart.Y(), pStart.Z()),
                         TVector3(pEnd.X(), pEnd.Y(), pEnd.Z()),
                         TVector3(0., 0., 0.), 0.,
                         fhit.GetTime() * ns2sec, fhit.GetEnergyLoss() * eV2GeV, 0, 0);
    }
  }

  void updateTimeDependentParams(ProcessingContext& pc)
  {
    pc.inputs().get<o2::itsmft::AlpideSimResponse*>("FOC_alpiderespvbb0");
    pc.inputs().get<o2::itsmft::AlpideSimResponse*>("FOC_alpiderespvbbm3");

    const auto& dopt = FOCALDigitizerParam::Instance();
    auto& digipar = mDigitizer.getParams();
    digipar.setContinuous(dopt.continuous);
    digipar.setROFrameBiasInBC(dopt.roFrameBiasInBC);
    digipar.setStrobeDelay(dopt.strobeDelay);
    if (dopt.continuous) {
      auto frameNS = dopt.roFrameLengthInBC * o2::constants::lhc::LHCBunchSpacingNS;
      digipar.setROFrameLengthInBC(dopt.roFrameLengthInBC);
      digipar.setROFrameLength(frameNS);
      digipar.setStrobeLength(dopt.strobeLengthCont > 0 ? dopt.strobeLengthCont : frameNS - dopt.strobeDelay);
    } else {
      digipar.setROFrameLength(dopt.roFrameLengthTrig);
      digipar.setStrobeLength(dopt.strobeLengthTrig);
    }
    digipar.getSignalShape().setParameters(dopt.strobeFlatTop, dopt.strobeMaxRiseTime, dopt.strobeQRiseTime0);
    // Round before DigiParams converts the minimum-charge fraction back to an integer.
    digipar.setChargeThreshold(dopt.chargeThreshold, dopt.chargeThreshold > 0 ? (dopt.minChargeToAccount + 0.5f) / dopt.chargeThreshold : 0.1f);
    digipar.setNoisePerPixel(dopt.noisePerPixel);
    digipar.setTimeOffset(dopt.timeOffset);
    digipar.setNSimSteps(dopt.nSimSteps);
    digipar.setEnergyToNElectrons(dopt.energyToNElectrons); // per GeV, see convertHits
    digipar.setVbb(dopt.Vbb);

    // The itsmft::Digitizer picks its response slot via the ITS/MFT Vbb parameters, which are irrelevant here:
    // provide the response matching the FoCal Vbb in both slots.
    const o2::itsmft::AlpideSimResponse* resp = nullptr;
    if (dopt.Vbb == 0.f) {
      resp = mRespVbb0;
      LOG(info) << ID.getName() << ": using ALPIDE response for Vbb=0V";
    } else if (dopt.Vbb == 3.f) {
      resp = mRespVbbM3;
      LOG(info) << ID.getName() << ": using ALPIDE response for Vbb=-3V";
    } else {
      LOG(fatal) << "Invalid FOCAL back-bias value " << dopt.Vbb << ", only 0 or 3 (-3V) are supported";
    }
    mDigitizer.setAlpideResponse(resp, 0);
    mDigitizer.setAlpideResponse(resp, 1);

    mROMode = digipar.isContinuous() ? o2::parameters::GRPObject::CONTINUOUS : o2::parameters::GRPObject::PRESENT;
    LOG(info) << ID.getName() << " simulated in "
              << ((mROMode == o2::parameters::GRPObject::CONTINUOUS) ? "CONTINUOUS" : "TRIGGERED")
              << " RO mode";

    mGeometry = o2::focal::PixelGeometryTGeo::Instance();
    LOG(info) << ID.getName() << ": " << mGeometry->getNumberOfPixelLayers() << " pixel layers, " << mGeometry->getNumberOfChips() << " chips";
    mDigitizer.setGeometry(mGeometry);
    mDigitizer.init();
  }

  bool mWithMCTruth = true;
  bool mFinished = false;
  bool mDisableQED = false;
  o2::itsmft::Digitizer mDigitizer;
  o2::focal::PixelGeometryTGeo* mGeometry = nullptr;
  std::vector<o2::itsmft::Digit> mDigits;
  std::vector<o2::itsmft::ROFRecord> mROFRecords;
  std::vector<o2::itsmft::ROFRecord> mROFRecordsAccum;
  std::vector<o2::focal::Hit> mFocalHits; // hits as read from the FOCPixelHit branch
  std::vector<o2::itsmft::Hit> mHits;     // same hits converted for the digitizer
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mLabels;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> mLabelsAccum;
  std::vector<o2::itsmft::MC2ROFRecord> mMC2ROFRecordsAccum;
  std::vector<TChain*> mSimChains;
  const o2::itsmft::AlpideSimResponse* mRespVbb0 = nullptr;
  const o2::itsmft::AlpideSimResponse* mRespVbbM3 = nullptr;

  int mFixMC2ROF = 0;                                                             // 1st entry in mc2rofRecordsAccum to be fixed for ROFRecordID
  o2::parameters::GRPObject::ROMode mROMode = o2::parameters::GRPObject::PRESENT; // readout mode
};

DataProcessorSpec getFOCALDigitizerSpec(int channel, bool mctruth)
{
  std::string detStr = o2::detectors::DetID::getName(FOCALDPLDigitizerTask::ID);
  auto detOrig = FOCALDPLDigitizerTask::Origin;
  std::vector<InputSpec> inputs;
  inputs.emplace_back("collisioncontext", "SIM", "COLLISIONCONTEXT", static_cast<SubSpecificationType>(channel), Lifetime::Timeframe);
  inputs.emplace_back("FOC_alpiderespvbb0", detOrig, "ALPIDERESPVbb0", 0, Lifetime::Condition, ccdbParamSpec("ITSMFT/Calib/ALPIDEResponseVbb0"));
  inputs.emplace_back("FOC_alpiderespvbbm3", detOrig, "ALPIDERESPVbbM3", 0, Lifetime::Condition, ccdbParamSpec("ITSMFT/Calib/ALPIDEResponseVbbM3"));

  std::vector<OutputSpec> outputs;
  outputs.emplace_back(detOrig, "DIGITS", 0, Lifetime::Timeframe);
  outputs.emplace_back(detOrig, "DIGITSROF", 0, Lifetime::Timeframe);
  if (mctruth) {
    outputs.emplace_back(detOrig, "DIGITSMC2ROF", 0, Lifetime::Timeframe);
    outputs.emplace_back(detOrig, "DIGITSMCTR", 0, Lifetime::Timeframe);
  }
  outputs.emplace_back(detOrig, "ROMode", 0, Lifetime::Timeframe);

  return DataProcessorSpec{.name = detStr + "Digitizer",
                           .inputs = inputs,
                           .outputs = outputs,
                           .algorithm = AlgorithmSpec{adaptFromTask<FOCALDPLDigitizerTask>(mctruth)},
                           .options = Options{{"disable-qed", o2::framework::VariantType::Bool, false, {"disable QED handling"}}}};
}

} // namespace o2::focal
