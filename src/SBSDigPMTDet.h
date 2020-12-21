#ifndef SBSDIGPMTDET_H
#define SBSDIGPMTDET_H

#include <iostream>
#include <vector>
#include <map>
#include "gmn_tree.h"
#include "g4sbs_tree.h"
#include "SBSDigPMTSignal.h"

//________________________________
class SBSDigPMTDet {
 public:
  SBSDigPMTDet();
  SBSDigPMTDet(UShort_t uniqueid, UInt_t nchan);
  SBSDigPMTDet(UShort_t uniqueid, UInt_t nchan, double NpeChargeConv, double sigmapulse, double gatewidth);
  virtual ~SBSDigPMTDet();
  //void Digitize(gmn_tree* T, TRandom3* R);
  void Digitize(g4sbs_tree* T, TRandom3* R);
  void Clear(bool dosamples = false);
  void SetSamples(double sampsize);
  
  //private:
  UInt_t fNChan;
  UShort_t fUniqueID;
  //std::map<int, PMTSignal> PMTmap;
  std::vector<PMTSignal> PMTmap;
  SPEModel* fRefPulse;
  
  double fGain;
  double fPedestal;
  double fPedSigma;
  double fTrigOffset;
  double fThreshold;
  double fGateWidth;
  double fADCconv;
  double fADCbits;
  double fTDCconv;
  double fTDCbits;
};

#endif
