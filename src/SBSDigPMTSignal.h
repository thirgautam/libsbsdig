#ifndef SBSDIGPMTSIGNAL_H
#define SBSDIGPMTSIGNAL_H

#include <iostream>
#include <vector>
#include <map>
#include <TROOT.h>
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
//#include "gmn_tree.h"
#include "g4sbs_tree.h"

//
// classes for signal digitization
//
//________________________________
class SPEModel {
 public:
  SPEModel();
  SPEModel(UShort_t uniqueid, double sigma, double t0 = 0, double tmin = -50, double tmax = +50);
  virtual ~SPEModel();
  double Eval(double t){return fPulseHisto->Interpolate(t);};
  bool   PulseOverThr(double charge, double thr);
  bool   FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail);
  
 private:
  TH1D *fPulseHisto;//At least we'll have to setup and use one per detector
  double GetHistoX(double y, double x1, double x2);
  
  //ClassDef(SPEModel,1);//it actually doesn't like classdef... (not sure why...)
};

//_________________________________
class PMTSignal {
 public:
  PMTSignal();
  PMTSignal(double npechargeconv);
  void Fill(SPEModel *model, int npe, double thr, double evttime, int signal);
  void Fill(int npe, double thr, double evttime, double sigmatime, int signal);
  void Digitize(int chan, int detid, g4sbs_tree* T, //gmn_tree* T, 
		TRandom3* R, double ped, double ped_noise, double ADCconv, double ADCbits, double TDCconv, double TDCbits);
  void Clear(bool dosamples = false);
  ~PMTSignal(){Clear();};
  
  void AddSumEdep(double edep){
    fSumEdep+= edep;
  };
  void SetNpeChargeConv(double npechargeconv){fNpeChargeConv = npechargeconv;};
  void SetSamples(double tmin, double tmax, double sampsize);
    
  double SumEdep(){return fSumEdep;};
  UInt_t Npe(){return fNpe;};
  double Charge(){return fNpe*fNpeChargeConv;};
  UInt_t ADC(){return fADC;};

  double EventTime(){return fEventTime;};
  UInt_t LeadTimesSize(){return fLeadTimes.size();};
  double LeadTime(int i){return fLeadTimes.at(i);};
  UInt_t TrailTimesSize(){return fTrailTimes.size();};
  double TrailTime(int i){return fTrailTimes.at(i);};
  UInt_t TDCSize(){return fTDCs.size();};
  UInt_t TDC(int i){return fTDCs.at(i);};
  //SimEncoder::tdc_data TDCData() { return fTDCData; }

  UInt_t ADCSamples(int i){return fADCSamples[i];};
  
 private:
  //summing variables for dig...
  double fSumEdep;//Not forced to use it for everything
  UInt_t fNpe;
  double fNpeChargeConv;
  UInt_t fADC;// One unique ADC value ?

  double fEventTime;
  //TDCs: multiple values possible.
  std::vector<double> fLeadTimes;
  std::vector<double> fTrailTimes;
  std::vector<Int_t> fTDCs;
  //SimEncoder::tdc_data fTDCData;
  //TRndmManager* fRN;

  //let's try something for HCal
  TF1* f1;
  int fNADCSamps;
  int fNSamps;
  double fSampSize;
  double fADCSampSize;
  double fTmin;
  double* fSamples;
  double* fADCSamples;
};

#endif
