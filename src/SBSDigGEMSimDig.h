#ifndef SBSDIGGEMSIMDIG_H
#define SBSDIGGEMSIMDIG_H

#include "TRandom3.h"
#include "TVector3.h"
//#include "SBSDigGEMPlane.h"
//#include "SBSDigGEMDet.h"
//#include "TArrayS.h"
//#include "TArrayI.h"
//#include "TArrayD.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>
#include <vector>

class SBSDigGEMDet;
class SBSDigGEMPlane;
//class gmn_tree;
class g4sbs_tree;

class SBSDigGEMSimDig {
 public:
  //Constructor and destructor
  SBSDigGEMSimDig();
  SBSDigGEMSimDig(int nchambers, double* trigoffset, double zsup_thr, int napv = 0, double* commonmode_array = 0);
  virtual ~SBSDigGEMSimDig();
  void Print();
  
  Int_t Digitize (SBSDigGEMDet* gemdet, TRandom3* R);//, gmn_tree* T);
  //void CheckOut(SBSDigGEMDet* gemdet, TRandom3* R, gmn_tree* T);
  void CheckOut(SBSDigGEMDet* gemdet, const int uniqueid, TRandom3* R, g4sbs_tree* T);
  //void FillBBGEMTree(const SBSDigGEMPlane pl, gmn_tree* T, int j);
  void write_histos();
  
  struct IonPar_t {
    Double_t X;       // position of the point on the projection
    Double_t Y;
    Double_t Charge;  // Charge deposited by this ion
    Double_t SNorm;   // 3 x radius of ion diffusion area at readout
    Double_t R2;      // = SNorm^2 : radius of numerical integration area
    Double_t ggnorm;  // = Charge/R2/pi : charge per unit area
  };

  //private:
  void AvaModel(const int ic, //module number
		SBSDigGEMDet* gemdet,
		TRandom3* R,
		const TVector3& xi,
		const TVector3& xo,
		const Double_t t0);
  
  void IonModel (TRandom3* R,
		 const TVector3& xi,
		 const TVector3& xo,
		 const Double_t elost );
  
  std::vector<Double_t> fTriggerOffset; // trigger offset (ns), incl latency & readout offset
  //UInt_t fNChambers;  // # chambers
  //UInt_t* fNROPlanes;   // # planes in each chamber
  UInt_t   fRNIon;    // number of ions
  std::vector<IonPar_t> fRIon;
  Double_t fRSMax;
  Double_t fRTotalCharge;
  Double_t fRTime0;
  Double_t fTimeZero;
  
  std::vector<Double_t> fSumA;
  //std::vector<Short_t>  fDADC;

  //zero suppression and common mode
  Bool_t fDoZeroSup;
  Double_t fZeroSup;
  Bool_t fDoCommonMode;
  std::vector<Double_t> fCommonModeArray;
  
  Short_t ADCConvert(Double_t val, Double_t off, Double_t gain, Int_t bits);
  Double_t PulseShape(Double_t t, 
		      Double_t C,  // normalization factor
		      Double_t Tp); // shaping time 
  //ClassDef (SBSDigGEMSimDig, 0) 

  /*
  TH2D* h1_QvsX_ion;
  TH2D* h1_QvsY_ion;
  TH2D* h1_QnormvsX_ion;
  TH2D* h1_QnormvsY_ion;
  TH2D* h1_QareavsX_ion;
  TH2D* h1_QareavsY_ion;
  TH2D* h1_QintvsX_ion;
  TH2D* h1_QintvsY_ion;
  TH2D* h1_QvsX_ava;
  TH2D* h1_QvsY_ava;
  TH2D* h1_QintYvsX_ava;
  TH2D* h1_QintYvsY_ava;
  TH1D* h1_yGEM_preion;
  TH1D* h1_yGEM_preava;
  TH1D* h1_yGEM_inava;
  TH1D* h1_yGEM_inava_2;
  TH1D* h1_yGEM_inava_3;
  TH1D* h1_yGEM_inava_4;
  TH2D* h1_xGEMvsADC_inava_4;
  TH2D* h1_yGEMvsADC_inava_4;
  TH1D* h1_yGEM_incheckout;
  */
};

#endif


