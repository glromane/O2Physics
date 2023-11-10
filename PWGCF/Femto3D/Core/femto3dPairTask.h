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

#include <memory>
#include "TLorentzVector.h"
#include "TVector3.h"

namespace o2::aod::singletrackselector
{

double particle_mass(int PDGcode)
{
  // if(PDGcode == 2212) return TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  if (PDGcode == 1000010020)
    return 1.87561294257;
  else
    return TDatabasePDG::Instance()->GetParticle(PDGcode)->Mass();
}

//====================================================================================

class FemtoParticle
{
 public:
  FemtoParticle() {}
  FemtoParticle(const float& pt, const float& eta, const float& phi);
  FemtoParticle(const float& pt, const float& eta, const float& phi, const int& signedPDG);
  FemtoParticle(const float& pt, const float& eta, const float& phi, const int& signedPDG, const float& magField);
  FemtoParticle(const FemtoParticle& obj);
  FemtoParticle(const FemtoParticle* obj);
  ~FemtoParticle();
  FemtoParticle& operator=(const FemtoParticle& obj);

  void Reset();

  void SetPt(const float& pt) { _pt = pt; }
  void SetEta(const float& eta) { _eta = eta; }
  void SetPhi(const float& phi) { _phi = phi; }
  void SetSignedPDG(const int& signedPDG) { _signedPDG = signedPDG; }
  void SetMagField(const float& magField) { _magField = magField; }

  float GetPt() const { return _pt; }
  float GetEta() const { return _eta; }
  float GetPhi() const { return _phi; }
  float GetSignedPDG() const { return _signedPDG; }
  float GetMagField() const { return _magField; }

  float GetPhiStar(const float& radius = 1.2) const;
  std::shared_ptr<TLorentzVector> Get4momentum() const;
  void Get4momentum(TLorentzVector& vec) const;

 private:
  float _pt = -1000.0, _eta = -1000.0, _phi = -1000.0, _magField = 0.0;
  int _signedPDG = 0;
};

FemtoParticle::FemtoParticle(const float& pt, const float& eta, const float& phi)
{
  _pt = pt;
  _eta = eta;
  _phi = phi;
}

FemtoParticle::FemtoParticle(const float& pt, const float& eta, const float& phi, const int& signedPDG)
{
  _pt = pt;
  _eta = eta;
  _phi = phi;
  _signedPDG = signedPDG;
}

FemtoParticle::FemtoParticle(const float& pt, const float& eta, const float& phi, const int& signedPDG, const float& magField)
{
  _pt = pt;
  _eta = eta;
  _phi = phi;
  _signedPDG = signedPDG;
  _magField = magField;
}

FemtoParticle::FemtoParticle(const FemtoParticle& obj)
{
  SetPt(obj.GetPt());
  SetEta(obj.GetEta());
  SetPhi(obj.GetPhi());
  SetSignedPDG(obj.GetSignedPDG());
  SetMagField(obj.GetMagField());
}

FemtoParticle::FemtoParticle(const FemtoParticle* obj)
{
  SetPt(obj->GetPt());
  SetEta(obj->GetEta());
  SetPhi(obj->GetPhi());
  SetSignedPDG(obj->GetSignedPDG());
  SetMagField(obj->GetMagField());
}

FemtoParticle::~FemtoParticle()
{
}

FemtoParticle& FemtoParticle::operator=(const FemtoParticle& obj)
{
  if (this != &obj) {
    SetPt(obj.GetPt());
    SetEta(obj.GetEta());
    SetPhi(obj.GetPhi());
    SetSignedPDG(obj.GetSignedPDG());
    SetMagField(obj.GetMagField());
  }

  return *this;
}

void FemtoParticle::Reset()
{
  _pt = -1000.0;
  _eta = -1000.0;
  _phi = -1000.0;
  _signedPDG = 0;
  _magField = 0.0;
}

float FemtoParticle::GetPhiStar(const float& radius) const
{
  if (_signedPDG && _magField)
    return _phi + asin(-0.3 * _magField * (_signedPDG / std::abs(_signedPDG)) * radius / (2.0 * _pt));
  else
    return _phi;
}

std::shared_ptr<TLorentzVector> FemtoParticle::Get4momentum() const
{
  std::shared_ptr<TLorentzVector> fourmomentum(new TLorentzVector(_pt * std::sin(_phi),
                                                                  _pt * std::cos(_phi),
                                                                  _pt * std::sinh(_eta),
                                                                  std::sqrt(_pt * std::cosh(_eta) * _pt * std::cosh(_eta) + particle_mass(std::abs(_signedPDG)) * particle_mass(std::abs(_signedPDG)))));
  return fourmomentum;
}

void FemtoParticle::Get4momentum(TLorentzVector& vec) const
{
  vec.SetPtEtaPhiM(_pt, _eta, _phi, particle_mass(std::abs(_signedPDG)));
}

//====================================================================================

class FemtoPair
{
 public:
  FemtoPair(){};
  FemtoPair(FemtoParticle* first, FemtoParticle* second)
  {
    _first = first;
    _second = second;
  }
  FemtoPair(FemtoParticle* first, FemtoParticle* second, const bool& isidentical)
  {
    _first = first;
    _second = second;
    _isidentical = isidentical;
  }

  FemtoPair(const FemtoPair& obj)
  {
    SetFirstParticle(obj.GetFirstParticle());
    SetSecondParticle(obj.GetSecondParticle());
    SetIdentical(obj.IsIdentical());
  }
  FemtoPair(const FemtoPair* obj)
  {
    SetFirstParticle(obj->GetFirstParticle());
    SetSecondParticle(obj->GetSecondParticle());
    SetIdentical(obj->IsIdentical());
  }
  ~FemtoPair() {}
  FemtoPair& operator=(const FemtoPair& obj)
  {
    if (this != &obj) {
      SetFirstParticle(obj.GetFirstParticle());
      SetSecondParticle(obj.GetSecondParticle());
      SetIdentical(obj.IsIdentical());
    }
    return *this;
  }

  void SetFirstParticle(FemtoParticle* first) { _first = first; }
  void SetSecondParticle(FemtoParticle* second) { _second = second; }
  void SetIdentical(const bool& isidentical) { _isidentical = isidentical; }
  void Reset();

  FemtoParticle* GetFirstParticle() const { return _first; }
  FemtoParticle* GetSecondParticle() const { return _second; }
  bool IsIdentical() const { return _isidentical; }

  bool IsClosePair(const float& deta = 0.01, const float& dphi = 0.01, const float& radius = 1.2);
  float GetEtaDiff() const { return _first->GetEta() - _second->GetEta(); }
  float GetPhiStarDiff(const float& radius = 1.2) const { return _first->GetPhiStar(radius) - _second->GetPhiStar(radius); }
  float GetKstar();

 private:
  FemtoParticle* _first = NULL;
  FemtoParticle* _second = NULL;
  bool _isidentical = true;
  TLorentzVector vec1, vec2;
};

void FemtoPair::Reset()
{
  _first = NULL;
  _second = NULL;
  _isidentical = true;
}

bool FemtoPair::IsClosePair(const float& deta, const float& dphi, const float& radius)
{
  if (_first == NULL || _second == NULL)
    return true;
  if (!(_first->GetMagField() * _second->GetMagField()))
    return true;
  if (abs(GetEtaDiff()) < deta || abs(GetPhiStarDiff(radius)) < dphi)
    return true;

  return false;
}

float FemtoPair::GetKstar()
{
  if (_first == NULL || _second == NULL)
    return -1000;

  _first->Get4momentum(vec1);
  _second->Get4momentum(vec2);

  if (_isidentical) {
    TLorentzVector fourmomentadiff = vec1 - vec2;
    return 0.5 * abs(fourmomentadiff.Mag());
  } else {
    TLorentzVector fourmomentasum = vec1 + vec2;

    vec1.Boost((-1) * fourmomentasum.BoostVector());
    vec1.Boost((-1) * fourmomentasum.BoostVector());

    TVector3 qinv = vec1.Vect() - vec2.Vect();
    return 0.5 * abs(qinv.Mag());
  }
}
} // namespace o2::aod::singletrackselector

#endif // PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
