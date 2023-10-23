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

void comb(std::vector<std::vector<int>> &indxs, int N, int K)
{
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's

    // print integers and permute bitmask
    do {
        std::vector<int> temp;
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i]){
              temp.push_back(i);
            }
        }
        indxs.push_back(temp);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

//====================================================================================

class FemtoParticle {
  public:
    FemtoParticle(){}
    FemtoParticle(TLorentzVector* fourmomentum, const float& eta, const float& phi);
    //FemtoParticle(TLorentzVector* fourmomentum);
    FemtoParticle(const float& E, const float& px, const float& py, const float& pz, const float& eta, const float& phi);
    //FemtoParticle(const float& E, const float& px, const float& py, const float& pz);
    FemtoParticle(const FemtoParticle& obj);
    FemtoParticle(const FemtoParticle* obj);
    ~FemtoParticle();
    FemtoParticle& operator=(const FemtoParticle &obj);

    void Set4momentum(TLorentzVector* fourmomentum){_fourmomentum = fourmomentum;}
    void SetEta(const float& eta){_eta = eta;}
    void SetPhi(const float& phi){_phi = phi;}
    void SetSign(const float& sign){_sign = sign;}
    void SetMagField(const float& magField){_magField = magField;}

    TLorentzVector* Get4momentum() const {return _fourmomentum;}
    float GetEta() const {return _eta;}
    float GetPhi() const {return _phi;}
    float GetPhiStar(const float& radius = 1.2){return _phi + asin( -0.3*_magField*_sign*radius/(2.0*_fourmomentum->Pt()) );}
    float GetMagField() const {return _magField;}

  private:
    float _eta, _phi, _sign, _magField = 0.0;
    TLorentzVector* _fourmomentum;
};

FemtoParticle::FemtoParticle(TLorentzVector* fourmomentum, const float& eta, const float& phi){
  _fourmomentum = fourmomentum;
  _eta = eta;
  _phi = phi;
}

//FemtoParticle::FemtoParticle(TLorentzVector* fourmomentum){
//  _fourmomentum = fourmomentum;
//}

FemtoParticle::FemtoParticle(const float& E, const float& px, const float& py, const float& pz, const float& eta, const float& phi){
  _fourmomentum = new TLorentzVector(px, py, pz, E);
  _eta = eta;
  _phi = phi;
}

//FemtoParticle::FemtoParticle(const float& E, const float& px, const float& py, const float& pz){
//  _fourmomentum = new TLorentzVector(*px, *py, *pz, *E);
//}

FemtoParticle::FemtoParticle(const FemtoParticle& obj){
  Set4momentum(obj.Get4momentum());
  SetEta(obj.GetEta());
  SetPhi(obj.GetPhi());
}

FemtoParticle::FemtoParticle(const FemtoParticle* obj){
  Set4momentum(obj->Get4momentum());
  SetEta(obj->GetEta());
  SetPhi(obj->GetPhi());
}

FemtoParticle::~FemtoParticle()
{
}

FemtoParticle& FemtoParticle::operator=(const FemtoParticle &obj)
{
	if (this != &obj) {
		Set4momentum(obj.Get4momentum());
		SetEta(obj.GetEta());
		SetPhi(obj.GetPhi());
	}

	return *this;
}

//====================================================================================


class FemtoPair{
  public:
    FemtoPair(){};
    FemtoPair(FemtoParticle* first, FemtoParticle* second){_first = first; _second = second;}
    FemtoPair(FemtoParticle* first, FemtoParticle* second, const bool& isidentical){_first = first; _second = second; _isidentical = isidentical;}

    FemtoPair(const FemtoPair& obj){ SetFirstParticle(obj.GetFirstParticle());   SetSecondParticle(obj.GetSecondParticle()); }
    FemtoPair(const FemtoPair* obj){ SetFirstParticle(obj->GetFirstParticle());   SetSecondParticle(obj->GetSecondParticle()); }
    ~FemtoPair(){}
    FemtoPair& operator=(const FemtoPair &obj){ if (this != &obj) {SetFirstParticle(obj.GetFirstParticle()); SetSecondParticle(obj.GetSecondParticle());}
                                                return *this;
                                              }

    void SetFirstParticle(FemtoParticle* first){_first = first;}
    void SetSecondParticle(FemtoParticle* second){_second = second;}
    void SetIdentical(const bool& isidentical){_isidentical = isidentical;}

    FemtoParticle* GetFirstParticle() const {return _first;}
    FemtoParticle* GetSecondParticle() const {return _second;}
    bool IsIdentical(){return _isidentical;}

    bool IsClosePair(const float& deta = 0.01, const float& dphi = 0.01, const float& radius = 1.2);
    float GetEtaDiff() const {return _first->GetEta() - _second->GetEta();}
    float GetPhiStarDiff(const float& radius = 1.2) const {return _first->GetPhiStar(radius) - _second->GetPhiStar(radius);}
    float GetKstar() const;

  private:
    FemtoParticle* _first;
    FemtoParticle* _second;
    bool _isidentical = true;
};

bool FemtoPair::IsClosePair(const float& deta, const float& dphi, const float& radius){
  if(_first == NULL || _second == NULL) return true;
  if(!(_first->GetMagField()*_second->GetMagField())) return true;
  if(abs(GetEtaDiff()) < deta || abs(GetPhiStarDiff(radius)) < dphi) return true;

  return false;
}

float FemtoPair::GetKstar() const {
  if(_first == NULL || _second == NULL) return -1000;
  if(!(_first->GetMagField()*_second->GetMagField())) return -1000;

  TLorentzVector* first4momentum = new TLorentzVector( *(_first->Get4momentum()) );
  TLorentzVector* second4momentum = new TLorentzVector( *(_second->Get4momentum()) );

  if(_isidentical){
    TLorentzVector fourmomentadiff = *first4momentum - *second4momentum;
    return 0.5*abs(fourmomentadiff.Mag());
  }
  else{
    TLorentzVector fourmomentasum = *first4momentum + *second4momentum;

    first4momentum->Boost( (-1)*fourmomentasum.BoostVector() );
    second4momentum->Boost( (-1)*fourmomentasum.BoostVector() );

    TVector3 qinv = first4momentum->Vect() - second4momentum->Vect();
    return 0.5*abs(qinv.Mag());
  }
}

#endif // PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
