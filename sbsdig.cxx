//includes: standard
#include <iostream>
#include <fstream>
#include <string>
#include <map>

//includes: root
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCut.h"
#include "TEventList.h"
#include "TMath.h"
#include "TRandom3.h"

//includes: specific
#include "G4SBSRunData.hh"
#include "g4sbs_types.h"
#include "gmn_tree.h"
#include "g4sbs_tree.h"
#include "SBSDigAuxi.h"
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"
#include "SBSDigGEMSimDig.h"
#include "SBSDigBkgdGen.h"

#ifdef __APPLE__
#include "unistd.h"
#endif

/*
//Defining here the parameters for the new detectors.
//TODO: write a list of parameters that are not "frozen" (e.g. gain, pedestal parameters, etc...) and switch them into databases...
#define NPlanes_bbgem 32 // modules...
#define NChan_bbps 52
#define NChan_bbsh 189
#define NChan_bbhodo 180 
#define NChan_grinch 510 
#define NChan_hcal 288

#define TriggerJitter 3.0 //ns
#define ADCbits 12 
#define gatewidth_PMT 100 //ns
#define gatewidth_GEM 400 //ns

#define FADC_sampsize 4.0 //ns

//DB???
#define sigmapulse_bbpsSH 3.0 //ns / 1.2 of 
#define sigmapulse_bbhodo 1.6 //ns
#define sigmapulse_grinch 3.75 //ns

#define gain_bbps 2.e6
#define ped_bbps 600.0 // ADC channel
#define pedsigma_bbps 3.0 // ADC channel
#define trigoffset_bbps 18.2 //ns
#define ADCconv_bbps 50 //fC/ch

#define gain_bbsh 7.5e5
#define ped_bbsh 500.0 // ADC channel
#define pedsigma_bbsh 4.5 // ADC channel
#define trigoffset_bbsh 18.5 //ns
#define ADCconv_bbsh 50 //fC/ch

#define gain_grinch 7.0e6
#define ped_grinch 0.0 
#define pedsigma_grinch 0.0
#define trigoffset_grinch 15.3 //ns
#define threshold_grinch 3.e-3 //V
#define ADCconv_grinch 100 //fC/ch
#define TDCconv_grinch 1.0 //ns/channel
#define TDCbits_grinch 16 //ns/channel

#define gain_bbhodo 1.0e5
#define ped_bbhodo 0.0 
#define pedsigma_bbhodo 0.0
#define trigoffset_bbhodo 18.6 //ns
#define threshold_bbhodo 3.e-3 //V
#define ADCconv_bbhodo 100 //fC/ch
#define TDCconv_bbhodo 0.1 //ns/channel
#define TDCbits_bbhodo 19 //ns/channel

#define gain_hcal 1.0e6
#define ped_hcal 0.0 
#define pedsigma_hcal 0.0
#define trigoffset_hcal 81.0 //ns
#define threshold_hcal 3.e-3 //V
#define ADCconv_hcal 1.0 //fC/ch //??
#define TDCconv_hcal 0.12 //ns/channel
#define TDCbits_hcal 16 //ns/channel
*/

using namespace std;
//____________________________________________________
int main(int argc, char** argv){
  
  // Step 0: read out arguments
  string db_file, inputsigfile, inputbkgdfile = "";//sources of files
  ULong64_t Nentries = -1;//number of events to process
  //UShort_t Nbkgd = 0;//number of background files to add to each event
  double LumiFrac = 0;
      
  if(argc<3){
    cout << "*** Not enough arguments! ***" << endl
	 << " Arguments: database (mandatory); " << endl
	 << "           list_of_sig_input_files (str, mandatory); " << endl
	 << "          nb_of_sig_evts_to_process (int, def=-1); " << endl
	 << "         list_of_bkgd_input_files (str, def=''); " << endl
	 << "        nb_of_bkgd_files_to_add_to_sig_evt (int, def=0); " << endl;
    return(-1);
  }
  
  db_file = argv[1];
  cout << " database file " << db_file << endl;
  inputsigfile = argv[2];
  cout << " Signal input files from: " << inputsigfile << endl;
  if(argc>3)Nentries = atoi(argv[3]);
  cout << " Number of (signal) events to process = " << Nentries << endl;
  if(argc>5){
    inputbkgdfile = argv[4];
    cout << " Background histgrams from: " << inputbkgdfile << endl;
    LumiFrac = max(0., atof(argv[5]));
    cout << " Fraction of background to superimpose to signal = " << LumiFrac << endl;
  }
  
  TFile* f_bkgd;
  SBSDigBkgdGen* BkgdGenerator;
  if(LumiFrac>0){
    f_bkgd = TFile::Open(inputbkgdfile.c_str());
    if(f_bkgd->IsZombie()){
      LumiFrac = 0;
    }else{
      BkgdGenerator = new SBSDigBkgdGen(f_bkgd);
    }
  }
  //f_bkgd->Close();
  
  // ------------------- // dev notes // ------------------- //
  // The loop on the input signal and background chains 
  // is going to happen here in the main I guess.
  //
  // First, we want to extend the input tree (for signal only!!!)
  // I guess in order to avoid adding extra layers of code, 
  // the tree extension might have to be coded in the custom tree class
  

  
  std::vector<SBSDigPMTDet*> PMTdetectors;
  std::vector<int> detmap;
  std::vector<SBSDigGEMDet*> GEMdetectors;
  std::vector<SBSDigGEMSimDig*> GEMsimDig;
  std::vector<int> gemdetmap;
  
  // Variable parameters. 
  // Can be configured with the database, but are provided with defaults.
  Int_t Rseed = 0;
  Double_t TriggerJitter = 3.0;
  
  std::vector<TString> detectors_list;
  
  Int_t NChan_bbps = 52;
  Double_t gatewidth_bbps = 100.;
  Double_t gain_bbps = 2.e6;
  Double_t ped_bbps = 600.;//
  Double_t pedsigma_bbps = 3.;//
  Double_t trigoffset_bbps = 18.2;//
  Double_t ADCconv_bbps = 50.;
  Int_t ADCbits_bbps = 12;
  Double_t sigmapulse_bbps = 3.0;
    
  Int_t NChan_bbsh = 189;
  Double_t gatewidth_bbsh = 100.;
  Double_t gain_bbsh = 7.5e5;
  Double_t ped_bbsh = 500.;
  Double_t pedsigma_bbsh = 4.5;
  Double_t trigoffset_bbsh = 18.5;
  Double_t ADCconv_bbsh = 50.;
  Int_t ADCbits_bbsh = 12;
  Double_t sigmapulse_bbsh = 3.0;
    
  Int_t NChan_grinch = 510;
  Double_t gatewidth_grinch = 100.;
  Double_t gain_grinch = 7.e6;
  Double_t ped_grinch = 0.;
  Double_t pedsigma_grinch = 0.;
  Double_t trigoffset_grinch = 15.3;
  Double_t threshold_grinch = 3.e-3;
  Double_t ADCconv_grinch = 100;
  Int_t ADCbits_grinch = 12;
  Double_t TDCconv_grinch = 1.;
  Int_t TDCbits_grinch = 16;
  Double_t sigmapulse_grinch = 3.75;
 
  Int_t NChan_bbhodo = 180;
  Double_t gatewidth_bbhodo = 100.;
  Double_t gain_bbhodo = 1.e5;
  Double_t ped_bbhodo = 0.;
  Double_t pedsigma_bbhodo = 0.;
  Double_t trigoffset_bbhodo = 18.6;
  Double_t threshold_bbhodo = 3.e3;
  Double_t ADCconv_bbhodo = 100.;
  Int_t ADCbits_bbhodo = 12;
  Double_t TDCconv_bbhodo = 0.1;
  Int_t TDCbits_bbhodo = 19;
  Double_t sigmapulse_bbhodo = 1.6;

  Int_t NChan_hcal = 288;
  Double_t gatewidth_hcal = 80;
  Double_t gain_hcal = 1.e6;
  Double_t ped_hcal = 0.;
  Double_t pedsigma_hcal = 0.;
  Double_t trigoffset_hcal = 81.;
  Double_t threshold_hcal = 3.e-3;
  Double_t ADCconv_hcal = 1.;
  Double_t TDCconv_hcal = 0.12;
  Int_t TDCbits_hcal = 16;
  Int_t FADC_ADCbits = 12;
  Double_t FADC_sampsize = 4.0;
 
  Int_t NPlanes_bbgem = 32;// number of planes/modules/readout
  Double_t gatewidth_bbgem = 400.;
  Double_t ZsupThr_bbgem = 240.;
  Int_t* nstrips_bbgem;
  Double_t* offset_bbgem;
  Double_t* RO_angle_bbgem;
  Double_t* triggeroffset_bbgem;
  Double_t* commonmode_array_bbgem;
  
  Int_t NPlanes_cepol_front = 32;// number of planes/modules/readout
  Double_t gatewidth_cepol_front = 400.;
  Double_t ZsupThr_cepol_front = 240.;
  Int_t* nstrips_cepol_front;
  Double_t* offset_cepol_front;
  Double_t* RO_angle_cepol_front;
  Double_t* triggeroffset_cepol_front;
  Double_t* commonmode_array_cepol_front;
  
  Int_t NPlanes_cepol_rear = 32;// number of planes/modules/readout
  Double_t gatewidth_cepol_rear = 400.;
  Double_t ZsupThr_cepol_rear = 240.;
  Int_t* nstrips_cepol_rear;
  Double_t* offset_cepol_rear;
  Double_t* RO_angle_cepol_rear;
  Double_t* triggeroffset_cepol_rear;
  Double_t* commonmode_array_cepol_rear;
  /*
  Int_t nstrips_bbgem[256];
  Double_t offset_bbgem[256];
  Double_t RO_angle_bbgem[256];
  Double_t triggeroffset_bbgem[256];
  Double_t commonmode_array_bbgem[65536];
  */
  UShort_t nAPV = 0;

  /*
  int nstrips_bbgem[32] = {1280, 1024, 1280, 1024, 1280, 1024, 
			   1280, 1024, 1280, 1024, 1280, 1024, 
			   1280, 1024, 1280, 1024, 1280, 1024, 
			   1280, 1024, 1280, 1024, 1280, 1024, 
			   1280, 1536, 1280, 1536, 
			   1280, 1536, 1280, 1536};
  double offset_bbgem[32] = {-0.512, 0., 0., 0., 0.512, 0., 
			     -0.512, 0., 0., 0., 0.512, 0., 
			     -0.512, 0., 0., 0., 0.512, 0., 
			     -0.512, 0., 0., 0., 0.512, 0., 
			     -0.768, 0., -0.256, 0., 
			     0.256, 0.,  0.768, 0.};
  double RO_angle_bbgem[32] = {0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 
			       0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 
			       0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 
			       0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 
			       0.0, 90.0, 0.0, 90.0, 
			       0.0, 90.0, 0.0, 90.0};
  for(int i = 0; i<NPlanes_bbgem; i++){
    RO_angle_bbgem[i]*= TMath::DegToRad();
    //cout << nstrips_bbgem[i] << " ";
  }//cout << endl;
  
  double triggeroffset_bbgem[16] = {121., 121., 121., 121.5, 121.5, 121.5,  
				    122., 122., 122., 122.5, 122.5, 122.5,  
				    126., 126., 126., 126.};
  
  double ZsupThr_bbgem = 240.;
  
  double commonmode_array_bbgem[1] = {1500.};nAPV = 1;
  */
  
  // ** How to add a new subsystem **
  // Add param for new detectors there...
 //polscint_bs
  Int_t NChan_polscint_bs = 48;
  Double_t gatewidth_polscint_bs = 30.;
  Double_t gain_polscint_bs = 3.e7;
  Double_t ped_polscint_bs = 300.;
  Double_t pedsigma_polscint_bs = 10.;
  Double_t trigoffset_polscint_bs = 37.6;
  Double_t threshold_polscint_bs = 3.e3;
  Double_t ADCconv_polscint_bs = 100.;
  Int_t ADCbits_polscint_bs = 12;
  Double_t TDCconv_polscint_bs = 0.1;
  Int_t TDCbits_polscint_bs = 19;
  Double_t sigmapulse_polscint_bs = 1.6;
 //polscint_fs
  Int_t NChan_polscint_fs = 48;
  Double_t gatewidth_polscint_fs = 30.;
  Double_t gain_polscint_fs = 3.e7;
  Double_t ped_polscint_fs = 300.;
  Double_t pedsigma_polscint_fs = 3.;
  Double_t trigoffset_polscint_fs = 37.6;
  Double_t threshold_polscint_fs = 3.e3;
  Double_t ADCconv_polscint_fs = 100.;
  Int_t ADCbits_polscint_fs = 12;
  Double_t TDCconv_polscint_fs = 0.1;
  Int_t TDCbits_polscint_fs = 19;
  Double_t sigmapulse_polscint_fs = 1.6;
  
 //activeana
  Int_t NChan_activeana = 32;
  Double_t gatewidth_activeana = 30.;
  Double_t gain_activeana = 3.e7;
  Double_t ped_activeana = 300.;
  Double_t pedsigma_activeana = 10.;
  Double_t trigoffset_activeana = 37.6;
  Double_t threshold_activeana = 3.e3;
  Double_t ADCconv_activeana = 100.;
  Int_t ADCbits_activeana = 19;
  Double_t TDCconv_activeana = 0.1;
  Int_t TDCbits_activeana = 19;
  Double_t sigmapulse_activeana = 1.6;

  //-----------------------------
  //  Read database
  //-----------------------------
  cout << "read database: " << db_file.c_str() << endl;
  ifstream in_db(db_file.c_str());
  if(!in_db.is_open()){
    cout << "database " << db_file.c_str() << " does not exist!!!" << endl;
    exit(-1);
  }
  
  TString currentline;
  while( currentline.ReadLine(in_db) && !currentline.BeginsWith("endconfig")){
    if( !currentline.BeginsWith("#") ){
      Int_t ntokens = 0;
      std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
      if( !tokens->IsEmpty() ) {
	ntokens = tokens->GetLast()+1;
      }
      //TObjArray *tokens = currentline.Tokenize(" ");//vg: def lost => versions prior to 6.06; should be fixed! ??? 
      //int ntokens = tokens->GetEntries();
      
      if( ntokens >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	
	if(skey=="Rseed"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Rseed = stemp.Atoi();
	}
	
	if(skey=="TriggerJitter"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TriggerJitter = stemp.Atof();
	}
	
	if(skey=="detectors_list"){
	  for(int k = 1; k<ntokens; k++){
	    TString sdet = ( (TObjString*) (*tokens)[k] )->GetString();
	    detectors_list.push_back(sdet);
	  }
	}
	
	//BBPS
	if(skey=="NChan_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbps = stemp.Atoi();
	}
	
	if(skey=="gatewidth_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbps = stemp.Atof();
	}
	
	if(skey=="gain_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_bbps = stemp.Atof();
	}
	
	if(skey=="ped_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_bbps = stemp.Atof();
	}
	
	if(skey=="pedsigma_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_bbps = stemp.Atof();
	}
	
	if(skey=="trigoffset_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbps = stemp.Atof();
	}
	
	if(skey=="ADCconv_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbps = stemp.Atof();
	}	

	if(skey=="ADCbits_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbps = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_bbps"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbps = stemp.Atof();
	}
	
	//BBSH
	if(skey=="NChan_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbsh = stemp.Atoi();
	}
	
	if(skey=="gatewidth_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbsh = stemp.Atof();
	}
	
	if(skey=="gain_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_bbsh = stemp.Atof();
	}
	
	if(skey=="ped_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_bbsh = stemp.Atof();
	}
	
	if(skey=="pedsigma_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_bbsh = stemp.Atof();
	}
	
	if(skey=="trigoffset_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbsh = stemp.Atof();
	}
	
	if(skey=="ADCconv_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbsh = stemp.Atof();
	}	

	if(skey=="ADCbits_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbsh = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_bbsh"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbsh = stemp.Atof();
	}

	//GRINCH
	if(skey=="NChan_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_grinch = stemp.Atoi();
	}
	
	if(skey=="gatewidth_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_grinch = stemp.Atof();
	}
	
	if(skey=="gain_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_grinch = stemp.Atof();
	}
	
	if(skey=="ped_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_grinch = stemp.Atof();
	}
	
	if(skey=="pedsigma_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_grinch = stemp.Atof();
	}
	
	if(skey=="trigoffset_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_grinch = stemp.Atof();
	}
	
	if(skey=="threshold_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_grinch = stemp.Atof();
	}
	
	if(skey=="ADCconv_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_grinch = stemp.Atof();
	}	

	if(skey=="ADCbits_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_grinch = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_grinch = stemp.Atof();
	}	
	
	if(skey=="TDCbits_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_grinch = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_grinch"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_grinch = stemp.Atof();
	}
	
	//BBHODO
	if(skey=="NChan_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_bbhodo = stemp.Atoi();
	}
	
	if(skey=="gatewidth_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbhodo = stemp.Atof();
	}
	
	if(skey=="gain_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_bbhodo = stemp.Atof();
	}
	
	if(skey=="ped_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_bbhodo = stemp.Atof();
	}
	
	if(skey=="pedsigma_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_bbhodo = stemp.Atof();
	}
	
	if(skey=="trigoffset_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_bbhodo = stemp.Atof();
	}
	
	if(skey=="threshold_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_bbhodo = stemp.Atof();
	}
	
	if(skey=="ADCconv_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_bbhodo = stemp.Atof();
	}	

	if(skey=="ADCbits_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_bbhodo = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_bbhodo = stemp.Atof();
	}	
	
	if(skey=="TDCbits_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_bbhodo = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_bbhodo"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_bbhodo = stemp.Atof();
	}
	
	//HCal
	if(skey=="NChan_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_hcal = stemp.Atoi();
	}
	
	if(skey=="gatewidth_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_hcal = stemp.Atof();
	}
	
	if(skey=="gain_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_hcal = stemp.Atof();
	}
	
	if(skey=="ped_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_hcal = stemp.Atof();
	}
	
	if(skey=="pedsigma_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_hcal = stemp.Atof();
	}
	
	if(skey=="trigoffset_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_hcal = stemp.Atof();
	}
	
	if(skey=="threshold_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_hcal = stemp.Atof();
	}
	
	if(skey=="ADCconv_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_hcal = stemp.Atof();
	}	
	
	if(skey=="TDCconv_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_hcal = stemp.Atof();
	}	
	
	if(skey=="TDCbits_hcal"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_hcal = stemp.Atoi();
	}	
	
	if(skey=="FADC_ADCbits"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  FADC_ADCbits = stemp.Atoi();
	}
	
	if(skey=="FADC_sampsize"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  FADC_sampsize = stemp.Atof();
	}
	
	// ** How to add a new subsystem **
	// Add reading of param from other detectors there...
	//GEn-RP Hodoscopes bs
	if(skey=="NChan_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_polscint_bs = stemp.Atoi();
	}
	
	if(skey=="gatewidth_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_polscint_bs = stemp.Atof();
	}
	
	if(skey=="gain_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_polscint_bs = stemp.Atof();
	}
	
	if(skey=="ped_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_polscint_bs = stemp.Atof();
	}
	
	if(skey=="pedsigma_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_polscint_bs = stemp.Atof();
	}
	
	if(skey=="trigoffset_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_polscint_bs = stemp.Atof();
	}
	
	if(skey=="threshold_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_polscint_bs = stemp.Atof();
	}
	
	if(skey=="ADCconv_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_polscint_bs = stemp.Atof();
	}	

	if(skey=="ADCbits_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_polscint_bs = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_polscint_bs = stemp.Atof();
	}	
	
	if(skey=="TDCbits_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_polscint_bs = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_polscint_bs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_polscint_bs = stemp.Atof();
	}
	
	//GEn-RP Hodoscopes fs
	if(skey=="NChan_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_polscint_fs = stemp.Atoi();
	}
	
	if(skey=="gatewidth_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_polscint_fs = stemp.Atof();
	}
	
	if(skey=="gain_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_polscint_fs = stemp.Atof();
	}
	
	if(skey=="ped_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_polscint_fs = stemp.Atof();
	}
	
	if(skey=="pedsigma_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_polscint_fs = stemp.Atof();
	}
	
	if(skey=="trigoffset_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_polscint_fs = stemp.Atof();
	}
	
	if(skey=="threshold_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_polscint_fs = stemp.Atof();
	}
	
	if(skey=="ADCconv_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_polscint_fs = stemp.Atof();
	}	

	if(skey=="ADCbits_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_polscint_fs = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_polscint_fs = stemp.Atof();
	}	
	
	if(skey=="TDCbits_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_polscint_fs = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_polscint_fs"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_polscint_fs = stemp.Atof();
	}

	//GEn-RP activeana
	if(skey=="NChan_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NChan_activeana = stemp.Atoi();
	}
	
	if(skey=="gatewidth_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_activeana = stemp.Atof();
	}
	
	if(skey=="gain_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gain_activeana = stemp.Atof();
	}
	
	if(skey=="ped_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ped_activeana = stemp.Atof();
	}
	
	if(skey=="pedsigma_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  pedsigma_activeana = stemp.Atof();
	}
	
	if(skey=="trigoffset_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  trigoffset_activeana = stemp.Atof();
	}
	
	if(skey=="threshold_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  threshold_activeana = stemp.Atof();
	}
	
	if(skey=="ADCconv_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCconv_activeana = stemp.Atof();
	}	

	if(skey=="ADCbits_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ADCbits_activeana = stemp.Atoi();
	}	
	
	if(skey=="TDCconv_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCconv_activeana = stemp.Atof();
	}	
	
	if(skey=="TDCbits_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  TDCbits_activeana = stemp.Atoi();
	}	
	
	if(skey=="sigmapulse_activeana"){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sigmapulse_activeana = stemp.Atof();
	}
	


	//GEMs
	if(skey=="NPlanes_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_bbgem = stemp.Atoi();
	  
	  nstrips_bbgem = new Int_t[NPlanes_bbgem];
	  offset_bbgem = new Double_t[NPlanes_bbgem];
	  RO_angle_bbgem = new Double_t[NPlanes_bbgem];
	  triggeroffset_bbgem = new Double_t[NPlanes_bbgem/2];
	}
	
	if(skey=="gatewidth_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_bbgem = stemp.Atof();
	}
	
	if(skey=="ZsupThr_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_bbgem = stemp.Atof();
	}
	
	if(skey=="nstrips_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_bbgem[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="offset_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_bbgem[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="RO_angle_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_bbgem[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="triggeroffset_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_bbgem/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_bbgem[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_bbgem = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_bbgem << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="commonmode_array_bbgem"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_bbgem = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_bbgem[k-1] = stemp.Atof();
	  }
	}

	//GEMs cepol_front...
	if(skey=="NPlanes_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_cepol_front = stemp.Atoi();
	  
	  nstrips_cepol_front = new Int_t[NPlanes_cepol_front];
	  offset_cepol_front = new Double_t[NPlanes_cepol_front];
	  RO_angle_cepol_front = new Double_t[NPlanes_cepol_front];
	  triggeroffset_cepol_front = new Double_t[NPlanes_cepol_front/2];
	}
	
	if(skey=="gatewidth_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_cepol_front = stemp.Atof();
	}
	
	if(skey=="ZsupThr_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_cepol_front = stemp.Atof();
	}
	
	if(skey=="nstrips_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_cepol_front[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="offset_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_cepol_front[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="RO_angle_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_cepol_front[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="triggeroffset_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_front/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_cepol_front[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_cepol_front = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_front << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="commonmode_array_cepol_front"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_cepol_front = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_cepol_front[k-1] = stemp.Atof();
	  }
	}

	//GEMs cepol_rear...
	if(skey=="NPlanes_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  NPlanes_cepol_rear = stemp.Atoi();
	  
	  nstrips_cepol_rear = new Int_t[NPlanes_cepol_rear];
	  offset_cepol_rear = new Double_t[NPlanes_cepol_rear];
	  RO_angle_cepol_rear = new Double_t[NPlanes_cepol_rear];
	  triggeroffset_cepol_rear = new Double_t[NPlanes_cepol_rear/2];
	}
	
	if(skey=="gatewidth_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  gatewidth_cepol_rear = stemp.Atof();
	}
	
	if(skey=="ZsupThr_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ZsupThr_cepol_rear = stemp.Atof();
	}
	
	if(skey=="nstrips_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      nstrips_cepol_rear[k-1] = stemp.Atoi();
	    }
	  }else{
	    cout << "number of entries for nstrips_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="offset_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      offset_cepol_rear[k-1] = stemp.Atof();
	    }
	  }else{
	    cout << "number of entries for offset_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="RO_angle_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      RO_angle_cepol_rear[k-1] = stemp.Atof()*TMath::DegToRad();
	    }
	  }else{
	    cout << "number of entries for RO_angle_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="triggeroffset_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  if(ntokens==NPlanes_cepol_rear/2+1){
	    for(int k = 1; k<ntokens; k++){
	      TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	      triggeroffset_cepol_rear[k-1] = stemp.Atof();
	      // if this is affected at k-1, program crashes... it should not and that's confusing...
	    }
	  }else{
	    cout << "number of entries for triggeroffset_cepol_rear = " << ntokens-1 
		 << " don't match Nplanes = " << NPlanes_cepol_rear << endl;
	    cout << "fix your db " << endl;
	    exit(-1);
	  }
	}
	
	if(skey=="commonmode_array_cepol_rear"){
	  cout << "reading " << skey.Data() << endl;
	  commonmode_array_cepol_rear = new Double_t[ntokens-1];
	  for(int k = 1; k<ntokens; k++){
	    nAPV++;
	    TString stemp = ( (TObjString*) (*tokens)[k] )->GetString();
	    commonmode_array_cepol_rear[k-1] = stemp.Atof();
	  }
	}

      }//end if( ntokens >= 2 )
      tokens->~TObjArray();// ineffective... :(
    }//end if( !currentline.BeginsWith("#"))
  }//end while
  
  //-----------------------------
  //  Declare detectors
  //-----------------------------
  cout << " declaring detectors " << endl;
  for(int k = 0; k<detectors_list.size(); k++){
    cout << "detector: " << detectors_list[k].Data() << "... " << endl;
    if(detectors_list[k] == "bbgem"){
      SBSDigGEMDet* bbgem = new SBSDigGEMDet(BBGEM_UNIQUE_DETID, NPlanes_bbgem, nstrips_bbgem, offset_bbgem, RO_angle_bbgem, 6, ZsupThr_bbgem);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_bbgem/2, triggeroffset_bbgem, ZsupThr_bbgem, nAPV, commonmode_array_bbgem);
      
      GEMdetectors.push_back(bbgem);
      gemdetmap.push_back(BBGEM_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    if(detectors_list[k] == "cepol_front"){
      SBSDigGEMDet* cepol_front = new SBSDigGEMDet(CEPOL_GEMFRONT_UNIQUE_DETID, NPlanes_cepol_front, nstrips_cepol_front, offset_cepol_front, RO_angle_cepol_front, 6, ZsupThr_cepol_front);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_cepol_front/2, triggeroffset_cepol_front, ZsupThr_cepol_front, nAPV, commonmode_array_cepol_front);
      
      GEMdetectors.push_back(cepol_front);
      gemdetmap.push_back(CEPOL_GEMFRONT_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "cepol_rear"){
      SBSDigGEMDet* cepol_rear = new SBSDigGEMDet(CEPOL_GEMREAR_UNIQUE_DETID, NPlanes_cepol_rear, nstrips_cepol_rear, offset_cepol_rear, RO_angle_cepol_rear, 6, ZsupThr_cepol_rear);
      SBSDigGEMSimDig* gemdig = new SBSDigGEMSimDig(NPlanes_cepol_rear/2, triggeroffset_cepol_rear, ZsupThr_cepol_rear, nAPV, commonmode_array_cepol_rear);
      
      GEMdetectors.push_back(cepol_rear);
      gemdetmap.push_back(CEPOL_GEMREAR_UNIQUE_DETID);
      GEMsimDig.push_back(gemdig);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbps"){
      SBSDigPMTDet* bbps = new SBSDigPMTDet(BBPS_UNIQUE_DETID, NChan_bbps, gain_bbps*qe, sigmapulse_bbps, gatewidth_bbps);

      bbps->fGain = gain_bbps;
      bbps->fPedestal = ped_bbps;
      bbps->fPedSigma = pedsigma_bbps;
      bbps->fTrigOffset = trigoffset_bbps;
      bbps->fGateWidth = gatewidth_bbps;
      bbps->fADCconv = ADCconv_bbps;
      bbps->fADCbits = ADCbits_bbps;
      
      PMTdetectors.push_back(bbps);
      detmap.push_back(BBPS_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbsh"){
      SBSDigPMTDet* bbsh = new SBSDigPMTDet(BBSH_UNIQUE_DETID, NChan_bbsh, gain_bbsh*qe, sigmapulse_bbsh, gatewidth_bbsh);
      
      bbsh->fGain = gain_bbsh;
      bbsh->fPedestal = ped_bbsh;
      bbsh->fPedSigma = pedsigma_bbsh;
      bbsh->fTrigOffset = trigoffset_bbsh;
      bbsh->fGateWidth = gatewidth_bbsh;
      bbsh->fADCconv = ADCconv_bbsh;
      bbsh->fADCbits = ADCbits_bbsh;    
      
      PMTdetectors.push_back(bbsh);
      detmap.push_back(BBSH_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "grinch"){
      SBSDigPMTDet* grinch = new SBSDigPMTDet(GRINCH_UNIQUE_DETID, NChan_grinch, gain_grinch*qe, sigmapulse_grinch, gatewidth_grinch);
  
      grinch->fGain = gain_grinch;
      grinch->fPedestal = ped_grinch;
      grinch->fPedSigma = pedsigma_grinch;
      grinch->fTrigOffset = trigoffset_grinch;
      grinch->fThreshold = threshold_grinch*spe_unit/ROimpedance;
      grinch->fGateWidth = gatewidth_grinch;
      grinch->fADCconv = ADCconv_grinch;
      grinch->fADCbits = ADCbits_grinch;
      grinch->fTDCconv = TDCconv_grinch;
      grinch->fTDCbits = TDCbits_grinch;
      
      PMTdetectors.push_back(grinch);
      detmap.push_back(GRINCH_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "bbhodo"){
      SBSDigPMTDet* bbhodo = new SBSDigPMTDet(HODO_UNIQUE_DETID, NChan_bbhodo, gain_bbhodo*qe, sigmapulse_bbhodo, gatewidth_bbhodo);
      
      bbhodo->fGain = gain_bbhodo;
      bbhodo->fPedestal = ped_bbhodo;
      bbhodo->fPedSigma = pedsigma_bbhodo;
      bbhodo->fTrigOffset = trigoffset_bbhodo;
      bbhodo->fThreshold = threshold_bbhodo*spe_unit/ROimpedance;
      bbhodo->fGateWidth = gatewidth_bbhodo;
      bbhodo->fADCconv = ADCconv_bbhodo;
      bbhodo->fADCbits = ADCbits_bbhodo;
      bbhodo->fTDCconv = TDCconv_bbhodo;
      bbhodo->fTDCbits = TDCbits_bbhodo; 
      
      PMTdetectors.push_back(bbhodo);
      detmap.push_back(HODO_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    if(detectors_list[k] == "hcal"){
      SBSDigPMTDet* hcal = new SBSDigPMTDet(HCAL_UNIQUE_DETID, NChan_hcal);
      
      hcal->fGain = gain_hcal;
      hcal->fPedestal = ped_hcal;
      hcal->fPedSigma = pedsigma_hcal;
      hcal->fTrigOffset = trigoffset_hcal;
      hcal->fThreshold = threshold_hcal*spe_unit/ROimpedance;
      hcal->fGateWidth = gatewidth_hcal;
      hcal->fADCconv = ADCconv_hcal;
      hcal->fADCbits = FADC_ADCbits;
      hcal->fTDCconv = TDCconv_hcal;
      hcal->fTDCbits = TDCbits_hcal; 
      hcal->SetSamples(FADC_sampsize);
      
      //ordered by increasing uinque id
      PMTdetectors.push_back(hcal);
      detmap.push_back(HCAL_UNIQUE_DETID);
      cout << " set up! " << endl;
    }
    
    // ** How to add a new subsystem **
    // Add the new detector here!
    if(detectors_list[k] == "prpolscint_bs"){
      SBSDigPMTDet* polscint_bs = new SBSDigPMTDet(PRPOLBS_SCINT_UNIQUE_DETID, NChan_polscint_bs, gain_polscint_bs*qe, sigmapulse_polscint_bs, gatewidth_polscint_bs);
      
      polscint_bs->fGain = gain_polscint_bs;
      polscint_bs->fPedestal = ped_polscint_bs;
      polscint_bs->fPedSigma = pedsigma_polscint_bs;
      polscint_bs->fTrigOffset = trigoffset_polscint_bs;
      polscint_bs->fThreshold = threshold_polscint_bs*spe_unit/ROimpedance;
      polscint_bs->fGateWidth = gatewidth_polscint_bs;
      polscint_bs->fADCconv = ADCconv_polscint_bs;
      polscint_bs->fADCbits = ADCbits_polscint_bs;
      polscint_bs->fTDCconv = TDCconv_polscint_bs;
      polscint_bs->fTDCbits = TDCbits_polscint_bs; 
      
      PMTdetectors.push_back(polscint_bs);
      detmap.push_back(PRPOLBS_SCINT_UNIQUE_DETID);
      cout << " set up! " << endl;
    } 
  
    if(detectors_list[k] == "prpolscint_fs"){
      SBSDigPMTDet* polscint_fs = new SBSDigPMTDet(PRPOLFS_SCINT_UNIQUE_DETID, NChan_polscint_fs, gain_polscint_fs*qe, sigmapulse_polscint_fs, gatewidth_polscint_fs);
      
      polscint_fs->fGain = gain_polscint_fs;
      polscint_fs->fPedestal = ped_polscint_fs;
      polscint_fs->fPedSigma = pedsigma_polscint_fs;
      polscint_fs->fTrigOffset = trigoffset_polscint_fs;
      polscint_fs->fThreshold = threshold_polscint_fs*spe_unit/ROimpedance;
      polscint_fs->fGateWidth = gatewidth_polscint_fs;
      polscint_fs->fADCconv = ADCconv_polscint_fs;
      polscint_fs->fADCbits = ADCbits_polscint_fs;
      polscint_fs->fTDCconv = TDCconv_polscint_fs;
      polscint_fs->fTDCbits = TDCbits_polscint_fs; 
      
      PMTdetectors.push_back(polscint_fs);
      detmap.push_back(PRPOLFS_SCINT_UNIQUE_DETID);
      cout << " set up! " << endl;
    } 

    if(detectors_list[k] == "activeana"){
      SBSDigPMTDet* activeana = new SBSDigPMTDet(ACTIVEANA_UNIQUE_DETID, NChan_activeana, gain_activeana*qe, sigmapulse_activeana, gatewidth_activeana);
      
      activeana->fGain = gain_activeana;
      activeana->fPedestal = ped_activeana;
      activeana->fPedSigma = pedsigma_activeana;
      activeana->fTrigOffset = trigoffset_activeana;
      activeana->fThreshold = threshold_activeana*spe_unit/ROimpedance;
      activeana->fGateWidth = gatewidth_activeana;
      activeana->fADCconv = ADCconv_activeana;
      activeana->fADCbits = ADCbits_activeana;
      activeana->fTDCconv = TDCconv_activeana;
      activeana->fTDCbits = TDCbits_activeana; 
      
      PMTdetectors.push_back(activeana);
      detmap.push_back(ACTIVEANA_UNIQUE_DETID);
      cout << " set up! " << endl;
    } 
  }
  /*  
  std::map<int, SBSDigPMTDet*> PMTdetectors;
  PMTdetectors[HCAL_UNIQUE_DETID] = hcal;
  PMTdetectors[HODO_UNIQUE_DETID] = bbhodo;
  PMTdetectors[BBPS_UNIQUE_DETID] = bbps;
  PMTdetectors[BBSH_UNIQUE_DETID] = bbsh;
  PMTdetectors[GRINCH_UNIQUE_DETID] = grinch;
  std::map<int, SBSDigGEMDet*> GEMdetectors;
  GEMdetectors[BBGEM_UNIQUE_DETID] = bbgem;
  */
  
  TRandom3* R = new TRandom3(Rseed);
  
  // Step 1: read input files build the input chains
  // build signal chain
  ifstream sig_inputfile(inputsigfile);
  TChain *C_s = new TChain("T");
  while( currentline.ReadLine(sig_inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C_s->Add(currentline.Data());
    }
  }
  TObjArray *fileElements_s=C_s->GetListOfFiles();
  TIter next_s(fileElements_s);
  TChainElement *chEl_s=0;
  
  /* need to change this... 
  // build background chain
  ifstream beam_inputfile(inputbkgdfile);
  TChain *C_b = new TChain("T");
  //TCut global_cut = "";
  //TEventList *eblist = new TEventList("eblist");
  if(Nbkgd!=0){
    while( currentline.ReadLine(beam_inputfile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	C_b->Add(currentline.Data());
	//global_cut += currentline.Data();
	//cout << currentline.Data() << endl;
      }
    }
    //C_b->Draw(">>eblist",global_cut);
  }
  //TObjArray *fileElements_b=C_b->GetListOfFiles();
  //TIter next_b(fileElements_b);
  //TChainElement *chEl_b=0;
  */
  
  //G4SBSRunData* run_data;
  
  double Theta_SBS, D_HCal;
  
  //gmn_tree *T_s_;//, *T_b;
  //g4sbs_tree *T_s;
  
  ULong64_t Nev_fs;//, Nev_fb;
  ULong64_t ev_s;//, ev_b;
  
  ULong64_t NEventsTotal = 0;
  //UShort_t nbkgd = 0;
  //int treenum = 0;
  //int oldtreenum = 0;
  
  int i_fs = 0;
  bool has_data;
  
  double timeZero;
  
  //T_b = new gmn_tree(C_b);
  //ev_b = 0;
  
  while (( chEl_s=(TChainElement*)next_s() )) {
    if(NEventsTotal>=Nentries){
      break;
    }
    TFile f_s(chEl_s->GetTitle(), "UPDATE");
    if(f_s.IsZombie())cout << "File " << chEl_s->GetTitle() << " cannot be found. Please check the path of your file." << endl; 
    //run_data = (G4SBSRunData*)f_s.Get("run_data");
    G4SBSRunData* run_data = (G4SBSRunData*)f_s.Get("run_data");
    Theta_SBS = run_data->fSBStheta;
    D_HCal = run_data->fHCALdist;
    //TFile fs_c(Form("digitized/simdigtest_%d.root", i_fs), "UPDATE");
    //f_s.Cp(Form("digitized/simdigtest_%d.root", i_fs));
    //if(fs_c.IsOpen())cout << "copy of file is open" << endl;
    //cout << fs_c->ReOpen("UPDATE") << endl;
    //C_s = (TChain*)fs_c.Get("T");
    C_s = (TChain*)f_s.Get("T");
    //T_s = new gmn_tree(C_s);
    //T_s = new g4sbs_tree(C_s, detectors_list);//vg: def lost
    g4sbs_tree *T_s = new g4sbs_tree(C_s, detectors_list);
    //g4sbs_tree T_s(C_s, detectors_list);
    
    // Expend tree here! (again, for signal only!!!)
    //T_s->AddDigBranches();
    
    Nev_fs = C_s->GetEntries();
    
    for(ev_s = 0; ev_s<Nev_fs; ev_s++, NEventsTotal++){
      if(NEventsTotal>=Nentries)break;
      if(NEventsTotal%100==0)
	cout << NEventsTotal << "/" << Nentries << endl;
      
      timeZero = R->Gaus(0.0, TriggerJitter);
      
      for(int k = 0; k<PMTdetectors.size(); k++){
	if(detmap[k]==HCAL_UNIQUE_DETID){
	  PMTdetectors[k]->Clear(true);
	}else{
	  PMTdetectors[k]->Clear();
	}
      }
      for(int k = 0; k<GEMdetectors.size(); k++){
	GEMdetectors[k]->Clear();
      }
      /*
      bbgem->Clear();
      //for(int i = 0; i<NPlanes_bbgem; i++){
      //cout << bbgem->GEMPlanes[i].GetNStrips() << " ";
      //}cout << endl;
      bbps->Clear();
      bbsh->Clear();
      grinch->Clear();
      bbhodo->Clear();
      hcal->Clear(true);
      */
      
      has_data = false;
      
      T_s->ClearDigBranches();
      T_s->GetEntry(ev_s);
      //T_s.ClearDigBranches();
      //T_s.GetEntry(ev_s);
      
      // unfold the thing then... but where???
      has_data = UnfoldData(T_s, Theta_SBS, D_HCal, R, PMTdetectors, detmap, GEMdetectors, gemdetmap, timeZero, 0);
      if(!has_data)continue;
      
      
      if(LumiFrac>0){
	BkgdGenerator->GenerateBkgd(R, PMTdetectors, detmap, GEMdetectors, gemdetmap, LumiFrac);
      }
      /*
      // loop here for background
      if(Nbkgd>0){
	nbkgd = 0;
	while( T_b->GetEntry(ev_b++) ){
	//while( T_b->GetEntry(eblist->GetEntry(ev_b++)) ){
	  treenum = C_b->GetTreeNumber();
	  if(treenum!=oldtreenum){
	    oldtreenum = treenum;
	    nbkgd++;
	    if(nbkgd>=Nbkgd)break;
	  }
	  timeZero = R->Uniform( -gatewidth_GEM-50., gatewidth_GEM/2.-50. );
	  
	  UnfoldData(T_b, Theta_SBS, D_HCal, R, PMTdetectors, detmap, GEMdetectors, gemdetmap, timeZero, 1);
	  //if(treenum)
	}
	
	while (( chEl_b=(TChainElement*)next_b() )) {
	  if(nbkgd>=Nbkgd)break;
	  cout << chEl_b->GetTitle() << endl;
	  TFile f_b(chEl_b->GetTitle());
	  C_b = (TChain*)f_b.Get("T");
	  T_b = new gmn_tree(C_b);
	  Nev_fb = C_b->GetEntries();
	  for(ev_b = 0; ev_b<Nev_fb; ev_b++){
	    T_b->GetEntry(ev_b);
	    UnfoldData(T_b, Theta_SBS, D_HCal, R);
	  }// end loop on background events
	  nbkgd++;
	}// end loop on background files
      }//end if Nbkgd>0
      */
      
      for(int k = 0; k<PMTdetectors.size(); k++){
	PMTdetectors[k]->Digitize(T_s,R);
	
	// bbps->Digitize(T_s,R);
	// bbsh->Digitize(T_s,R);
	// grinch->Digitize(T_s,R);
	// bbhodo->Digitize(T_s,R);
	// hcal->Digitize(T_s,R);
      }
      
      for(int k = 0; k<GEMdetectors.size(); k++){
	GEMsimDig[k]->Digitize(GEMdetectors[k], R);
	GEMsimDig[k]->CheckOut(GEMdetectors[k], gemdetmap[k], R, T_s);
      }
      //How come this function is so taxing in time??? 
      // Answer in the function... hope we've found a workaround
      //FillDigTree(T_s, PMTdetectors, GEMdetectors);

      T_s->FillDigBranches();
      //T_s.FillDigBranches();
      //T_s->fChain->Fill();
    }// end loop on signal events 
    /*
    for(int k = 0; k<GEMdetectors.size(); k++){
      GEMsimDig[k]->write_histos();
    }
    */
    T_s->fChain->Write("", TObject::kOverwrite);
    //T_s.fChain->Write("", TObject::kOverwrite);
    //fs_c.Write();
    //fs_c.Close();
    f_s.Write();
    f_s.Close();
    //T_s->~g4sbs_tree();
    //T_s.~g4sbs_tree();
    i_fs++;
  }// end loop on signal files
  
  
  exit(0);
}
