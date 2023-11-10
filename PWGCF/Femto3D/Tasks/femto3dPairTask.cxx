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


struct FemtoCorrelations {
  //using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe>; // aod::pidTPCPr
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

  Configurable<int> _multbinwidth{"multbinwidth", 50, "width of multiplicity bins within which the mixing is done"};
  Configurable<int> _vertexbinwidth{"vertexbinwidth", 2, "width of vertexZ bins within which the mixing is done"};
  ConfigurableAxis CFkStarBinning{"CFkStarBinning", {500, 0.005, 5.005}, "k* binning of the CF (Nbins, lowlimit, uplimit)"};


  bool IsIdentical;

  std::pair<int, std::vector<float>> TPCcuts_1;
  std::pair<int, std::vector<float>> TOFcuts_1;

  std::pair<int, std::vector<float>> TPCcuts_2;
  std::pair<int, std::vector<float>> TOFcuts_2;

  std::vector<o2::aod::singletrackselector::FemtoParticle*> SEtracks1;
  std::vector<o2::aod::singletrackselector::FemtoParticle*> SEtracks2;

  std::map<int64_t, std::vector<o2::aod::singletrackselector::FemtoParticle*>> selectedtracks_1;
  std::map<int64_t, std::vector<o2::aod::singletrackselector::FemtoParticle*>> selectedtracks_2;
  std::map<std::pair<int, int>, std::vector<int>> mixbins;

  using FilteredCollisions = aod::SingleCollSels;
  using FilteredTracks = aod::SingleTrackSels;


  Filter pFilter = o2::aod::singletrackselector::p > _min_P && o2::aod::singletrackselector::p < _max_P;
  Filter etaFilter = nabs(o2::aod::singletrackselector::eta) < _eta;
  Filter tpcTrkFilter = o2::aod::singletrackselector::tpcNClsFound >= _tpcNClsFound &&
                        o2::aod::singletrackselector::tpcChi2NCl < _tpcChi2NCl &&
                        o2::aod::singletrackselector::tpcCrossedRowsOverFindableCls > _tpcCrossedRowsOverFindableCls &&
                        o2::aod::singletrackselector::tpcNClsShared <= (uint8_t)_tpcNClsShared;
  Filter dcaFilter = nabs(o2::aod::singletrackselector::dcaXY) < _dcaXY && nabs(o2::aod::singletrackselector::dcaZ) < _dcaZ;
  Filter itsNClsFilter = o2::aod::singletrackselector::itsNCls >= (uint8_t)_itsNCls && o2::aod::singletrackselector::itsChi2NCl < _itsChi2NCl;

  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) < _vertexZ;

  Preslice<FilteredTracks> perCollId = o2::aod::singletrackselector::singleCollSelId;


  void init(o2::framework::InitContext&)
  {
    IsIdentical = (_sign_1*_particlePDG_1 == _sign_2*_particlePDG_2);

    TPCcuts_1 = std::make_pair(_particlePDG_1, _tpcNSigma_1);
    TOFcuts_1 = std::make_pair(_particlePDG_1, _tofNSigma_1);
    TPCcuts_2 = std::make_pair(_particlePDG_2, _tpcNSigma_2);
    TOFcuts_2 = std::make_pair(_particlePDG_2, _tofNSigma_2);
    
    const AxisSpec kStarAxis{CFkStarBinning, "k* (GeV/c)"};

    registry.add("SE", "SE", kTH1F, {kStarAxis});
    registry.add("ME", "ME", kTH1F, {kStarAxis});
    registry.add("p_first", "p", kTH1F, {{100, 0., 5., "p"}});
    registry.add("nsigmaTOF_first", Form("nsigmaTOF_%i", (int)_particlePDG_1), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    registry.add("nsigmaTPC_first", Form("nsigmaTPC_%i", (int)_particlePDG_1), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    if(!IsIdentical){
      registry.add("p_second", Form("p_%i", (int)_particlePDG_2), kTH1F, {{100, 0., 5., "p"}});
      registry.add("nsigmaTOF_second", Form("nsigmaTOF_%i", (int)_particlePDG_2), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
      registry.add("nsigmaTPC_second", Form("nsigmaTPC_%i", (int)_particlePDG_2), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    }
  }


  void process(soa::Filtered<FilteredCollisions> const& collisions, soa::Filtered<FilteredTracks> const& tracks)
  {
    if(_particlePDG_1 == 0 || _particlePDG_2 == 0) LOGF(fatal, "One of passed PDG is 0!!!");

    o2::aod::singletrackselector::FemtoParticle* Particle;
    o2::aod::singletrackselector::FemtoPair* Pair = new o2::aod::singletrackselector::FemtoPair();

    for (auto& collision : collisions) {
      for (auto& track : tracks.sliceBy(perCollId, collision.index()) ) {

        if(track.sign() == _sign_1 && (track.p() < _PIDtrshld_1 ? o2::aod::singletrackselector::TPCselection(track, TPCcuts_1) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_1)) ){ // filling the map: eventID <-> selected particles1
          Particle = new o2::aod::singletrackselector::FemtoParticle( track.pt(), track.eta(), track.phi() );
          Particle->SetSignedPDG( _sign_1*_particlePDG_1 );
          Particle->SetMagField(collision.magField());
          SEtracks1.push_back(Particle);

          registry.fill(HIST("p_first"), track.p());
          if(_particlePDG_1 == 2212){
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPr());
          }
          if(_particlePDG_1 == 1000010020){
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaDe());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaDe());
          }
        }

        if(IsIdentical || track.sign() != _sign_2) continue;
        else if( !TOFselection(track, std::make_pair(_particlePDGtoReject, _rejectWithinNsigmaTOF)) && (track.p() < _PIDtrshld_2 ? o2::aod::singletrackselector::TPCselection(track, TPCcuts_2) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_2)) ){ // filling the map: eventID <-> selected particles2 if (see condition above ^)
          Particle = new o2::aod::singletrackselector::FemtoParticle( track.pt(), track.eta(), track.phi() );
          Particle->SetSignedPDG( _sign_2*_particlePDG_2 );
          Particle->SetMagField(collision.magField());
          SEtracks2.push_back(Particle);

          registry.fill(HIST("p_second"), track.p());
          if(_particlePDG_2 == 2212){
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaPr());
          }
          if(_particlePDG_2 == 1000010020){
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaDe());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaDe());
          }
        }
      }

      if(SEtracks1.size() || SEtracks2.size()) mixbins[std::pair<int, int>{round(collision.posZ()/_vertexbinwidth), floor(collision.mult()/_multbinwidth)}].push_back(collision.index());

      if(IsIdentical){ // start of same event identical

        for(int ii=0; ii<SEtracks1.size(); ii++){ // nested loop for all the combinations
          for(int iii=ii+1; iii<SEtracks1.size(); iii++){

            Pair->SetFirstParticle(SEtracks1[ii]);
            Pair->SetSecondParticle(SEtracks1[iii]);
            Pair->SetIdentical(IsIdentical);

            if(!Pair->IsClosePair(_deta, _dphi, _radiusTPC)) registry.fill(HIST("SE"), Pair->GetKstar());  // close pair rejection and fillig the SE histo
            Pair->Reset();
          }
        }
      } // end of same event identical
      else{ // start of same event non-identical

        for(auto& ii : SEtracks1 ){
          for(auto& iii : SEtracks2 ){

            Pair->SetFirstParticle(ii);
            Pair->SetSecondParticle(iii);
            Pair->SetIdentical(IsIdentical);

            if(!Pair->IsClosePair(_deta, _dphi, _radiusTPC)) registry.fill(HIST("SE"), Pair->GetKstar());  // close pair rejection and fillig the SE histo
            Pair->Reset();
          }
        }
      }// end of same event non-identical

      if(SEtracks1.size()) selectedtracks_1.insert( std::pair<int64_t, std::vector<o2::aod::singletrackselector::FemtoParticle*>>(collision.index(), SEtracks1) );
      SEtracks1.clear();

      if(!IsIdentical && SEtracks2.size()){
        selectedtracks_2.insert( std::pair<int64_t, std::vector<o2::aod::singletrackselector::FemtoParticle*>>(collision.index(), SEtracks2) );
        SEtracks2.clear();
      }
    }

    //====================================== mixing identical ======================================

    if(IsIdentical){ // start of same event identical

      for (auto i = mixbins.begin(); i != mixbins.end(); i++){ // start of mixed event identical

        for(int indx1=0; indx1<(i->second).size(); indx1++){ // nested loop for all the combinations of collisions in a chosen mult/vertex bin
          for(int indx2=indx1+1; indx2<(i->second).size(); indx2++){

            for(auto ii : selectedtracks_1[ (i->second)[indx1] ]){
              for(auto iii : selectedtracks_1[ (i->second)[indx2] ]){

                Pair->SetFirstParticle(ii);
                Pair->SetSecondParticle(iii);
                Pair->SetIdentical(IsIdentical);

                if(!Pair->IsClosePair(_deta, _dphi, _radiusTPC)) registry.fill(HIST("ME"), Pair->GetKstar());
                Pair->Reset();
              }
            }
          }
        }
      }// end of mixed event identical

    } //====================================== end of mixing identical ======================================

    else{ //====================================== mixing non-identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++){ // start of mixed event non-identical

        for(int indx1=0; indx1<(i->second).size(); indx1++){ // nested loop for all the combinations of collisions in a chosen mult/vertex bin
          for(int indx2=indx1+1; indx2<(i->second).size(); indx2++){

            for(auto ii : selectedtracks_1[ (i->second)[indx1] ]){
              for(auto iii : selectedtracks_2[ (i->second)[indx2] ]){

                Pair->SetFirstParticle(ii);
                Pair->SetSecondParticle(iii);
                Pair->SetIdentical(IsIdentical);

                if(!Pair->IsClosePair(_deta, _dphi, _radiusTPC)) registry.fill(HIST("ME"), Pair->GetKstar());
                Pair->Reset();
              }
            }
          }
        }
      } // end of mixed event non-idetical

    } //====================================== end of mixing non-identical ======================================


    // clearing up
    for(auto i = selectedtracks_1.begin(); i != selectedtracks_1.end(); i++){
      for(auto ii : i->second) delete ii;
      (i->second).clear();
    }
    selectedtracks_1.clear();

    if(!IsIdentical) for(auto i = selectedtracks_2.begin(); i != selectedtracks_2.end(); i++){
      for(auto ii : i->second) delete ii;
      (i->second).clear();
    }
    selectedtracks_2.clear();

    for(auto i = mixbins.begin(); i != mixbins.end(); i++){
      (i->second).clear();
    }
    mixbins.clear();

    delete Pair;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FemtoCorrelations>(cfgc)};
}