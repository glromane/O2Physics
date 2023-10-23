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
/// \brief utility functions for femto task
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 30 May 2023

#ifndef PWGCF_DATAMODEL_FEMTOTEST_H_
#define PWGCF_DATAMODEL_FEMTOTEST_H_

//#include "Framework/ASoA.h"
//#include "Framework/DataTypes.h"
//#include "Framework/AnalysisDataModel.h"
//#include "Common/DataModel/PIDResponse.h"
//#include "Framework/Logger.h"
//#include "Common/DataModel/Multiplicity.h"

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"

double particle_mass(int PDGcode){
  //if(PDGcode == 2212) return TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  if(PDGcode == 1000010020) return 1.87561294257;
  else return TDatabasePDG::Instance()->GetParticle(PDGcode)->Mass();
}

namespace o2::aod::singletrackselector
{
//====================================================================================

template <typename TrackType> class FemtoPair{
  public:
    FemtoPair(){}
    FemtoPair(TrackType const& first, TrackType const& second){_first = first; _second = second;}
    FemtoPair(TrackType const& first, TrackType const&  second, const bool& isidentical){_first = first; _second = second; _isidentical = isidentical;}

    FemtoPair(const FemtoPair& obj){ SetFirstParticle(obj.GetFirstParticle());   SetSecondParticle(obj.GetSecondParticle()); }
    FemtoPair(const FemtoPair* obj){ SetFirstParticle(obj->GetFirstParticle());   SetSecondParticle(obj->GetSecondParticle()); }
    ~FemtoPair(){}
    FemtoPair& operator=(const FemtoPair &obj){ if (this != &obj) {SetFirstParticle(obj.GetFirstParticle()); SetSecondParticle(obj.GetSecondParticle());}
                                                return *this;
                                              }

    void SetPair(TrackType const&  first, TrackType const&  second){_first = first; _second = second;}
    void SetFirstParticle(TrackType const&  first){_first = first;}
    void SetSecondParticle(TrackType const&  second){_second = second;}
    void SetIdentical(const bool& isidentical){_isidentical = isidentical;}
    void SetMagField1(const float& magfield1){_magfield1 = magfield1;}
    void SetMagField2(const float& magfield2){_magfield2 = magfield2;}
    void SetPDG1(const int& PDG1){_PDG1 = PDG1;}
    void SetPDG2(const int& PDG2){_PDG2 = PDG2;}
    void ResetPair();
    void ResetAll();

    TrackType* GetFirstParticle() const {return _first;}
    TrackType* GetSecondParticle() const {return _second;}
    bool IsIdentical(){return _isidentical;}

    bool IsClosePair(const float& deta = 0.01, const float& dphi = 0.01, const float& radius = 1.2);
    float GetEtaDiff() const {
      if(_first != NULL && _second != NULL) return _first->eta() - _second->eta();
      else return 1000;}
    float GetPhiStarDiff(const float& radius = 1.2) const {
      if(_first != NULL && _second != NULL) return _first->phiStar(_magfield1, radius) - _second->phiStar(_magfield2, radius);
      else return 1000;}
    float GetKstar() const;

  private:
    TrackType _first = NULL;
    TrackType _second = NULL;
    float _magfield1 = 0.0, _magfield2 = 0.0;
    int _PDG1 = 0, _PDG2 = 0;
    bool _isidentical = true;
};

template <typename TrackType> void FemtoPair<TrackType>::ResetPair(){
  _first = NULL;
  _second = NULL;
}

template <typename TrackType> void FemtoPair<TrackType>::ResetAll(){
  _first = NULL;
  _second = NULL;
  _magfield1 = 0.0;
  _magfield2 = 0.0;
  _PDG1 = 0;
  _PDG2 = 0;
  _isidentical = true;
}

template <typename TrackType> bool FemtoPair<TrackType>::IsClosePair(const float& deta, const float& dphi, const float& radius){
  if(_first == NULL || _second == NULL) return true;
  if(!(_magfield1*_magfield2)) return true;
  if(abs(GetEtaDiff()) < deta || abs(GetPhiStarDiff(radius)) < dphi) return true;

  return false;
}

template <typename TrackType> float FemtoPair<TrackType>::GetKstar() const {
  if(_first == NULL || _second == NULL) return -1000;
  if(!(_magfield1*_magfield2)) return -1000;
  if(!(_PDG1*_PDG2)) return -1000;

  TLorentzVector first4momentum;
  first4momentum.SetPtEtaPhiM(_first->pt(), _first->eta(), _first->phi(), particle_mass(_PDG1));
  TLorentzVector second4momentum;
  second4momentum.SetPtEtaPhiM(_second->pt(), _second->eta(), _second->phi(), particle_mass(_PDG2));

  if(_isidentical){
    TLorentzVector fourmomentadiff = first4momentum - second4momentum;
    return 0.5*abs(fourmomentadiff.Mag());
  }
  else{
    TLorentzVector fourmomentasum = first4momentum + second4momentum;

    first4momentum.Boost( (-1)*fourmomentasum.BoostVector() );
    second4momentum.Boost( (-1)*fourmomentasum.BoostVector() );

    TVector3 qinv = first4momentum.Vect() - second4momentum.Vect();
    return 0.5*abs(qinv.Mag());
  }
}
}// namespace singletrackselector

#endif // PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
