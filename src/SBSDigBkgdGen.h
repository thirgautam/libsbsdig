#ifndef SBSDIGBKGDGEN_H
#define SBSDIGBKGDGEN_H

#include <iostream>
#include <vector>
#include <map>
#include "gmn_tree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"


//________________________________
class SBSDigBkgdGen {

 public:
  SBSDigBkgdGen();
  SBSDigBkgdGen(TFile* f_bkgd);
  ~SBSDigBkgdGen();
  void Initialize(TFile* f_bkgd);
  
  void GenerateBkgd(//double theta_sbs, double d_hcal, 
		    TRandom3* R, 
		    std::vector<SBSDigPMTDet*> pmtdets,
		    std::vector<int> detmap, 
		    std::vector<SBSDigGEMDet*> gemdets, 
		    std::vector<int> gemmap, 
		    double lumifrac);
  
 private:
  Double_t* NhitsBBGEMs;
  TH1D* h_EdephitBBGEMs;
  TH1D** h_xhitBBGEMs;
  TH1D** h_yhitBBGEMs;
  TH1D** h_dxhitBBGEMs;
  TH1D** h_dyhitBBGEMs;
  
  Double_t* NhitsHCal;
  TH1D* h_EdephitHCal;
  TH1D* h_zhitHCal;
  
  Double_t* NhitsBBPS;
  TH1D* h_EdephitBBPS;
  
  Double_t* NhitsBBSH;
  TH1D* h_EdephitBBSH;
  
  Double_t* NhitsBBHodo;
  TH1D* h_EdephitBBHodo;
  TH1D* h_xhitBBHodo;
  
  Double_t* P1hitGRINCH;
  Double_t* P2hitsGRINCH;
  TH1D* h_NpeGRINCH;
  
};

#endif // SBSDIGAUXI_H

