#ifndef SBSDIGGEMPLANE_H
#define SBSDIGGEMPLANE_H

#include <iostream>
#include <vector>
#include <map>
#include <TROOT.h>
//#include "TH1D.h"
//#include "TRandom3.h"

class SBSDigGEMPlane {
 public:
  SBSDigGEMPlane();
  SBSDigGEMPlane(short mod, int nstrips, int nsamples = 6, double thr = 100, double offset = 0, double roangle = 0);
  virtual ~SBSDigGEMPlane();
  void Clear();
  
  //void SetStripThreshold(double thr){Striphr = ;};

  Short_t Module(){return fModule;};
  Double_t dX(){return fdX;};
  Double_t Xoffset(){return fXoffset;};
  Double_t ROangle(){return fROangle;};
  Int_t GetNStrips(){return fNStrips;};
  Short_t GetADC(int strip, int samp){return fStripADC[strip*fNSamples+samp];};
  Int_t GetADCSum(int strip){return fStripADCsum[strip];};
  void SetADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripADCsum[strip]+= adc-fStripADC[strip*fNSamples+samp];
      fStripADC[strip*fNSamples+samp] = adc;
    }
  };
  void AddADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripADC[strip*fNSamples+samp]+=adc;
      fStripADCsum[strip]+=adc;
    }//else{
    //printf("strip = %d / %d, sample %d /%d \n", strip, fNStrips, samp, fNSamples);
    //}
  };
  
 private:
  // ADC sampled value of strip array of each axis
  Short_t fModule;
  Int_t fNStrips;
  Int_t fNSamples;
  Double_t fStripThr;//threshold for ADC sum
  Int_t* fStripADCsum;
  Short_t* fStripADC;
  double fdX;
  double fXoffset;
  double fROangle;
  //ClassDef(SBSDigGEMPlane, 1)
};
#endif
