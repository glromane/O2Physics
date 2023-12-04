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
//
/// \brief Femto3D pair mixing task
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <TParameter.h>
#include <TH1F.h>

#include "Framework/ASoA.h"
#include "MathUtils/Utils.h"
#include "Framework/DataTypes.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"

#include "Framework/StaticFor.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include <vector>
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "PWGCF/Femto3D/Core/femto3dPairTask.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FemtoCorrelationsO2 {
  // using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe>; // aod::pidTPCPr
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> _min_P{"min_P", 0.0, "lower mometum limit"};
  Configurable<float> _max_P{"max_P", 100.0, "upper mometum limit"};
  Configurable<float> _eta{"eta", 100.0, "abs eta value limit"};
  Configurable<float> _dcaXY{"dcaXY", 10.0, "abs dcaXY value limit"};
  Configurable<float> _dcaZ{"dcaZ", 10.0, "abs dcaZ value limit"};
  Configurable<int16_t> _tpcNClsFound{"minTpcNClsFound", 0, "minimum allowed number of TPC clasters"};
  Configurable<float> _tpcChi2NCl{"tpcChi2NCl", 100.0, "upper limit for chi2 value of a fit over TPC clasters"};
  Configurable<float> _tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0, "lower limit of TPC CrossedRows/FindableCls value"};
  Configurable<int> _tpcNClsShared{"maxTpcNClsShared", 100, "maximum allowed number of TPC shared clasters"};
  Configurable<int> _itsNCls{"minItsNCls", 0, "minimum allowed number of ITS clasters"};
  Configurable<float> _itsChi2NCl{"itsChi2NCl", 100.0, "upper limit for chi2 value of a fit over ITS clasters"};
  Configurable<float> _vertexZ{"VertexZ", 10.0, "abs vertexZ value limit"};

  Configurable<int> _sign_1{"sign_1", 1, "sign of the first particle in a pair"};
  Configurable<int> _particlePDG_1{"particlePDG_1", 2212, "PDG code of the first particle in a pair to perform PID for (only proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma_1{"tpcNSigma_1", std::vector<float>{-3.0f, 3.0f}, "first particle PID: Nsigma range in TPC before the TOF is used"};
  Configurable<float> _PIDtrshld_1{"PIDtrshld_1", 10.0, "first particle PID: value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma_1{"tofNSigma_1", std::vector<float>{-3.0f, 3.0f}, "first particle PID: Nsigma range in TOF"};

  Configurable<int> _sign_2{"sign_2", 1, "sign of the second particle in a pair"};
  Configurable<int> _particlePDG_2{"particlePDG_2", 2212, "PDG code of the second particle in a pair to perform PID for (only proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma_2{"tpcNSigma_2", std::vector<float>{-3.0f, 3.0f}, "second particle PID: Nsigma range in TPC before the TOF is used"};
  Configurable<float> _PIDtrshld_2{"PIDtrshld_2", 10.0, "second particle PID: value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma_2{"tofNSigma_2", std::vector<float>{-3.0f, 3.0f}, "second particle PID: Nsigma range in TOF"};

  Configurable<int> _particlePDGtoReject{"particlePDGtoRejectFromSecond", 0, "applied only if the particles are non-identical and only to the second particle in the pair!!!"};
  Configurable<std::vector<float>> _rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-0.0f, 0.0f}, "TOF rejection Nsigma range for the particle specified with PDG to be rejected"};

  Configurable<float> _deta{"deta", 0.01, "minimum allowed defference in eta between two tracks in a pair"};
  Configurable<float> _dphi{"dphi", 0.01, "minimum allowed defference in phi_star between two tracks in a pair"};
  Configurable<float> _radiusTPC{"radiusTPC", 1.2, "TPC radius to calculate phi_star for"};

  ConfigurableAxis CFkStarBinning{"CFkStarBinning", {500, 0.005, 5.005}, "k* binning of the CF (Nbins, lowlimit, uplimit)"};

  Configurable<bool> _doME{"doME", false, "minimum allowed defference in eta between two tracks in a pair"};

  bool IsIdentical;

  o2::aod::singletrackselector::PIDcuts PIDsel_1;
  o2::aod::singletrackselector::PIDcuts PIDsel_2;

  using FilteredCollisions = aod::SingleCollSels;
  using FilteredTracks = aod::SingleTrackSels;

  typedef soa::Filtered<FilteredTracks>::iterator* trkType;
  std::unique_ptr<o2::aod::singletrackselector::FemtoPair<trkType>> Pair = std::make_unique<o2::aod::singletrackselector::FemtoPair<trkType>>();
  // o2::aod::singletrackselector::FemtoPair<trkType>* Pair = new o2::aod::singletrackselector::FemtoPair<trkType>();

  Filter pFilter = o2::aod::singletrackselector::p > _min_P&& o2::aod::singletrackselector::p < _max_P;
  Filter etaFilter = nabs(o2::aod::singletrackselector::eta) < _eta;
  Filter tpcTrkFilter = o2::aod::singletrackselector::tpcNClsFound >= _tpcNClsFound &&
                        o2::aod::singletrackselector::tpcChi2NCl < _tpcChi2NCl &&
                        o2::aod::singletrackselector::tpcCrossedRowsOverFindableCls > _tpcCrossedRowsOverFindableCls&&
                                                                                        o2::aod::singletrackselector::tpcNClsShared <= (uint8_t)_tpcNClsShared;
  Filter dcaFilter = nabs(o2::aod::singletrackselector::dcaXY) < _dcaXY && nabs(o2::aod::singletrackselector::dcaZ) < _dcaZ;
  Filter itsNClsFilter = o2::aod::singletrackselector::itsNCls >= (uint8_t)_itsNCls && o2::aod::singletrackselector::itsChi2NCl < _itsChi2NCl;

  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) < _vertexZ;

  Preslice<FilteredTracks> perCollId = o2::aod::singletrackselector::singleCollSelId;

  void init(o2::framework::InitContext&)
  {
    PIDsel_1.signedPDG = _particlePDG_1 * _sign_1;
    PIDsel_1.PIDtrshld = _PIDtrshld_1;
    PIDsel_1.TPClowLimit = _tpcNSigma_1.value[0] - o2::aod::singletrackselector::nsigma::binning::bin_width;
    PIDsel_1.TPCupLimit = _tpcNSigma_1.value[1] + o2::aod::singletrackselector::nsigma::binning::bin_width;
    PIDsel_1.TOFlowLimit = _tofNSigma_1.value[0] - o2::aod::singletrackselector::nsigma::binning::bin_width;
    PIDsel_1.TOFupLimit = _tofNSigma_1.value[1] + o2::aod::singletrackselector::nsigma::binning::bin_width;

    PIDsel_2.signedPDG = _particlePDG_2 * _sign_2;
    PIDsel_2.PIDtrshld = _PIDtrshld_2;
    PIDsel_2.TPClowLimit = _tpcNSigma_2.value[0] - o2::aod::singletrackselector::nsigma::binning::bin_width;
    PIDsel_2.TPCupLimit = _tpcNSigma_2.value[1] + o2::aod::singletrackselector::nsigma::binning::bin_width;
    PIDsel_2.TOFlowLimit = _tofNSigma_2.value[0] - o2::aod::singletrackselector::nsigma::binning::bin_width;
    PIDsel_2.TOFupLimit = _tofNSigma_2.value[1] - o2::aod::singletrackselector::nsigma::binning::bin_width;

    IsIdentical = (_sign_1 * _particlePDG_1 == _sign_2 * _particlePDG_2);

    Pair->SetIdentical(IsIdentical);
    Pair->SetPDG1(_particlePDG_1);
    Pair->SetPDG2(_particlePDG_2);

    const AxisSpec kStarAxis{CFkStarBinning, "k* (GeV/c)"};

    registry.add("SE", "SE", kTH1F, {kStarAxis});
    registry.add("ME", "ME", kTH1F, {kStarAxis});
    registry.add("p_first", "p", kTH1F, {{100, 0., 5., "p"}});
    registry.add("nsigmaTOF_first", Form("nsigmaTOF_%i", (int)_particlePDG_1), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    registry.add("nsigmaTPC_first", Form("nsigmaTPC_%i", (int)_particlePDG_1), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    if (!IsIdentical) {
      registry.add("p_second", Form("p_%i", (int)_particlePDG_2), kTH1F, {{100, 0., 5., "p"}});
      registry.add("nsigmaTOF_second", Form("nsigmaTOF_%i", (int)_particlePDG_2), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
      registry.add("nsigmaTPC_second", Form("nsigmaTPC_%i", (int)_particlePDG_2), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    }
  }

  SliceCache cacheTrk;
  ConfigurableAxis zBins{"binningVertex", {VARIABLE_WIDTH, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}, "vertex binning for mixing"};
  ConfigurableAxis MultBins{"binningMultiplicity", {VARIABLE_WIDTH, 0, 50, 100, 200, 300, 400, 500, 550, 600.1}, "multiplicity / centrality binning for mixing"};
  using BinningType = ColumnBinningPolicy<o2::aod::singletrackselector::PosZ, o2::aod::singletrackselector::Mult>;
  BinningType VertexMultBinning{{zBins, MultBins}, true}; // true is for 'ignore overflows' (true by default)

  //==================================================================================================================================

  template <typename ColType, typename TrackType>
  void fillSE(ColType const& collisions, TrackType& tracks)
  {
    if (!IsIdentical)
      LOGF(fatal, "Wrong 'fillSE' function has been called -- for non-identical use 2 sets of tracks instead of 1!!!");

    for (auto& collision : collisions) {
      auto SEtracks = tracks->sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision.index(), cacheTrk);

      Pair->SetMagField1(collision.magField());
      Pair->SetMagField2(collision.magField());

      for (auto& track : SEtracks) {
        registry.fill(HIST("p_first"), track.p());
        registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPr());
        registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPr());
      }

      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
      for (auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(SEtracks, SEtracks))) {
        Pair->SetPair(&trk1, &trk2);

        if (!Pair->IsClosePair(_deta, _dphi, _radiusTPC))
          registry.fill(HIST("SE"), Pair->GetKstar()); // close pair rejection and fillig the SE histo
        Pair->ResetPair();
      }
    }
  }

  template <typename ColType, typename TrackType>
  void fillSE(ColType const& collisions, TrackType& tracks1, TrackType& tracks2)
  {
    if (IsIdentical)
      LOGF(fatal, "Wrong 'fillSE' function has been called -- for identical use 1 set of tracks instead of 2!!!");

    for (auto& collision : collisions) {
      auto SEtracks1 = tracks1->sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision.index(), cacheTrk);
      auto SEtracks2 = tracks2->sliceByCached(o2::aod::singletrackselector::singleCollSelId, collision.index(), cacheTrk);

      Pair->SetMagField1(collision.magField());
      Pair->SetMagField2(collision.magField());

      for (auto& track : SEtracks1) {
        registry.fill(HIST("p_first"), track.p());
        registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPr());
        registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPr());
      }
      for (auto& track : SEtracks2) {
        registry.fill(HIST("p_second"), track.p());
        registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaPr());
        registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaPr());
      }

      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
      for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(SEtracks1, SEtracks2))) {
        Pair->SetPair(&trk1, &trk2);

        if (!Pair->IsClosePair(_deta, _dphi, _radiusTPC))
          registry.fill(HIST("SE"), Pair->GetKstar()); // close pair rejection and fillig the SE histo
        Pair->ResetPair();
      }
    }
  }

  template <typename ColType, typename TrackType>
  void fillME(ColType const& collisions, TrackType& tracks)
  {
    if (!IsIdentical)
      LOGF(fatal, "Wrong 'fillME' function has been called -- for non-identical use 2 sets of tracks instead of 1!!!");

    for (auto& [col1, col2] : selfPairCombinations(VertexMultBinning, 5, -1, collisions)) {

      auto tracksCol1 = tracks->sliceByCached(o2::aod::singletrackselector::singleCollSelId, col1.index(), cacheTrk);
      auto tracksCol2 = tracks->sliceByCached(o2::aod::singletrackselector::singleCollSelId, col2.index(), cacheTrk);

      Pair->SetMagField1(col1.magField());
      Pair->SetMagField2(col2.magField());

      for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(tracksCol1, tracksCol2))) {
        Pair->SetPair(&trk1, &trk2);

        if (!Pair->IsClosePair(_deta, _dphi, _radiusTPC))
          registry.fill(HIST("ME"), Pair->GetKstar()); // close pair rejection and fillig the SE histo
        Pair->ResetPair();
      }
    }
  }

  template <typename ColType, typename TrackType>
  void fillME(ColType const& collisions, TrackType& tracks1, TrackType& tracks2)
  {
    if (IsIdentical)
      LOGF(fatal, "Wrong 'fillME' function has been called -- for identical use 1 set of tracks instead of 2!!!");

    for (auto& [col1, col2] : selfPairCombinations(VertexMultBinning, 5, -1, collisions)) {

      auto tracksCol1 = tracks1->sliceByCached(o2::aod::singletrackselector::singleCollSelId, col1.index(), cacheTrk);
      auto tracksCol2 = tracks2->sliceByCached(o2::aod::singletrackselector::singleCollSelId, col2.index(), cacheTrk);

      Pair->SetMagField1(col1.magField());
      Pair->SetMagField2(col2.magField());

      for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(tracksCol1, tracksCol2))) {
        Pair->SetPair(&trk1, &trk2);

        if (!Pair->IsClosePair(_deta, _dphi, _radiusTPC))
          registry.fill(HIST("ME"), Pair->GetKstar()); // close pair rejection and fillig the SE histo
        Pair->ResetPair();
      }
    }
  }

  void process(soa::Filtered<FilteredCollisions> const& collisions, soa::Filtered<FilteredTracks> const& tracks)
  {
    if (_particlePDG_1 == 0 || _particlePDG_2 == 0)
      LOGF(fatal, "One of passed PDG is 0!!!");

    Partition<soa::Filtered<FilteredTracks>> groupedTracks1 = o2::aod::singletrackselector::PIDselection(PIDsel_1);
    groupedTracks1.bindTable(tracks);

    if (IsIdentical) {
      fillSE(collisions, groupedTracks1);
      if (_doME)
        fillME(collisions, groupedTracks1);
    } else {
      Partition<soa::Filtered<FilteredTracks>> groupedTracks2 = o2::aod::singletrackselector::PIDselection(PIDsel_2);
      groupedTracks2.bindTable(tracks);

      fillSE(collisions, groupedTracks1, groupedTracks2);
      if (_doME)
        fillME(collisions, groupedTracks1, groupedTracks2);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FemtoCorrelationsO2>(cfgc)};
}