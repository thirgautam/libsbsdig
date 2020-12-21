
// Class handling digitization of GEMs

#ifndef __TGEMSBSSimDigitization__
#define __TGEMSBSSimDigitization__

#include "TRandom3.h"
#include "TVector3.h"
#include "TArrayS.h"
#include "TArrayI.h"
#include "TArrayD.h"


#include "THaAnalysisObject.h"

#include <vector>
#define MAXSAMP 6

class TFile;
class TTree;

class TGEMSBSGEMSimHitData;
class TGEMSBSGEMHit;
class TGEMSBSSpec;
//class TGEMSBSSimEvent;
class TGEMSBSDBManager;

// First an auxiliary class

// The whole strip plane; used to cumulate virtual strips charges
// and generate real strips

class TGEMSBSDigitizedPlane {
private:
  // ADC sampled value of strip array of each axis

  //TODO: make this a struct inside an STL vector or similar
  TArrayI fStripADC;  
  TArrayI fStripSimADC;  
  Short_t *fType;  // Type of track (primary, secondary) which left the hit for each strip
  Int_t   *fTotADC;  // number of ADC counts for each strip

  Float_t *fCharge;  // charge for each strip
  Float_t *fTime;   // time for each strip

  UShort_t fNSamples;   // number of ADC samples
  UShort_t fNStrips;   // number of strips in the plane
  Int_t    fThreshold;  // ADC threshold 

  UShort_t  fNOT;   // # strips over threshold
  Short_t*  fOverThr;  // # list of strips over threshold

  std::vector< std::vector<Short_t> > fStripClusters; // Clusters associated with each strip
  std::vector< std::vector<Double_t> > fStripWeightInCluster;
  std::vector< std::vector<Int_t> > fStripClusterADC[6];

  //used to simulate cross talk of APV25
  TRandom3 fRan;

public:
  // init and reset physics strips arrays
  TGEMSBSDigitizedPlane (UShort_t nstrip,
		      UShort_t nsample = 10,
		      Int_t    threshold = 0 );
  ~TGEMSBSDigitizedPlane();
  void Clear();

  // cumulate hits (strips signals)
  void Cumulate (const TGEMSBSGEMHit *vv, Short_t type, Short_t clusterID );
  
  //standard getters
  Short_t  GetType (Int_t n) const {return fType[n];}
  Int_t    GetTotADC (Int_t n) const {return fTotADC[n];}
  Float_t  GetTime (Int_t n) const {return fTime[n];}
  Float_t  GetCharge (Int_t n) const {return fCharge[n];}
  Int_t    GetADC (Int_t n, Int_t ks) const {return fStripADC[n*fNSamples+ks];}
  Int_t    GetSimADC (Int_t n, Int_t ks) const {return fStripSimADC[n*fNSamples+ks];}
  Int_t    GetSimADCSum (Int_t n) const {
    Int_t sum = 0;
    for(int i = 0; i<fNSamples; i++)sum+= GetSimADC (n, i);
    return sum;
  }
  void SetSimADC (Int_t n, Int_t ks, Int_t adc){fStripSimADC[n*fNSamples+ks] = adc;}
  UShort_t GetNSamples() const {return fNSamples;}
  UShort_t GetNStrips() const {return fNStrips;}

  UShort_t Threshold (Int_t thr);

  UShort_t GetNOverThr() const {return fNOT;}
  Short_t  GetIdxOverThr (Int_t n) const {return fOverThr[n];}

  const std::vector<Short_t>& GetStripClusters(UInt_t n) const { return fStripClusters[n]; }
  const std::vector<Double_t>& GetStripWeightInCluster(UInt_t n) const { return fStripWeightInCluster[n]; }
  const std::vector<Int_t>& GetStripClusterADC(UInt_t k, UInt_t n) const { return fStripClusterADC[k][n]; }
};



class TGEMSBSSimDigitization: public THaAnalysisObject
{
 public:
  //Constructor and destructor
  TGEMSBSSimDigitization( const TGEMSBSSpec& spect,
			  const char* name = "ratedig", TGEMSBSDBManager *manager = 0);
  virtual ~TGEMSBSSimDigitization();
  
  //full initialization of all parameters with database
  void Initialize(const TGEMSBSSpec& spect);
  Int_t ReadDatabase (const TDatime& date);
  
  //This is in those three functions that the job is done, more specifically in AddititveDigitize
  //Int_t Digitize (const TGEMSBSGEMSimHitData& gdata, const TGEMSBSSpec& spect); // digitize event
  void EventEnd() {};
  void EventStart();
  Int_t AdditiveDigitize (const TGEMSBSGEMSimHitData& gdata, const TGEMSBSSpec& spect); // add more digitized data, e.g. background
  void NoDigitize (const TGEMSBSGEMSimHitData& gdata, const TGEMSBSSpec& spect); // do not digitize event, just fill tree
  const TGEMSBSDigitizedPlane& GetDigitizedPlane (UInt_t ich, UInt_t ip) const {return *(fDP[ich][ip]);}; 
  void Print( Option_t* opt="" ) const;// print info
  void PrintCharges() const;
  void PrintSamples() const;
  
  Double_t GetGateWidth(){ return fGateWidth; }

  // Tree methods
  // To write a tree with digitization results:
  //   Call InitTree before main loop;
  //   Call SetTreeEvent in main loop (before or after Digitize)
  //   Call FillTree in main loop (after Digitize and SetTreeEvent)
  // Call WriteTree and CloseTree after main loop

  void InitTree (const TGEMSBSSpec& spect, const TString& ofile);
  //dpulication of the SetTreeEvent routine with G4SBS file input instead of EVIO file
  /*
  void SetTreeEvent (const TGEMSBSGEMSimHitData& tsgd,
		     const TGEMSBSGeant4File& f,
		     Int_t evnum = -1);
         */
  Short_t SetTreeHit (const UInt_t ih,
		      const TGEMSBSSpec& spect,
		      //TGEMSBSGEMHit* const *dh,
		      const TGEMSBSGEMSimHitData& tsgd,
		      Double_t t0 ); // called from Digitization
  void SetTreeStrips(); // called from Digitization
  //void FillTree ();
  //void WriteTree () const;
  //void CloseTree () const;

  // Access to results
  Short_t GetType (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetType (n);}
  Int_t   GetTotADC (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetTotADC (n);}
  Float_t GetTime (UInt_t ich, UInt_t ip, UInt_t n) const {return fDP[ich][ip]->GetTime (n);}
  Float_t GetCharge (UInt_t ich, UInt_t ip, UInt_t n) const {return fDP[ich][ip]->GetCharge (n);}
  Int_t   GetADC (UInt_t ich, UInt_t ip, Int_t n, Int_t ks) const {return fDP[ich][ip]->GetADC (n, ks);}
  Int_t   GetSimADC (UInt_t ich, UInt_t ip, Int_t n, Int_t ks) const {return fDP[ich][ip]->GetSimADC (n, ks);}
  Int_t   GetSimADCSum (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetSimADCSum (n);}
  void SetSimADC (UInt_t ich, UInt_t ip, Int_t n, Int_t ks, Int_t adc) {fDP[ich][ip]->SetSimADC (n, ks,adc);}
  UInt_t   GetNChambers() const {return fNChambers;};
  UInt_t   GetNPlanes (const UInt_t i) const {return fNPlanes[i];}
  UShort_t GetNSamples (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNSamples();}
  UShort_t GetNStrips (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNStrips();}
  UShort_t Threshold (UInt_t ich, UInt_t ip, Int_t thr) {return fDP[ich][ip]->Threshold (thr);}
  UShort_t GetNOverThr (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNOverThr();}
  Short_t  GetIdxOverThr (UInt_t ich, UInt_t ip, Int_t n) const
  { return fDP[ich][ip]->GetIdxOverThr(n); }

  const std::vector<Short_t>& GetStripClusters(UInt_t ich, UInt_t ip, UInt_t n) const
  { return fDP[ich][ip]->GetStripClusters(n); }
  const std::vector<Double_t>& GetStripWeightInCluster(UInt_t ich, UInt_t ip, UInt_t n) const
  { return fDP[ich][ip]->GetStripWeightInCluster(n); }
  const std::vector<Int_t>& GetStripClusterADC(UInt_t ich, UInt_t ip, UInt_t n, UInt_t k) const
  { return fDP[ich][ip]->GetStripClusterADC(k, n); }

  //TGEMSBSSimEvent* GetEvent() const { return fEvent; }
  
  Bool_t IsMapSector() const { return fDoMapSector; }
  void SetMapSector( Bool_t b = true ) { fDoMapSector = b; }
  void GetGEMData(TGEMSBSGEMSimHitData* gd);
  
  void SetTimeZero(double t0){fTimeZero = t0;}
  Double_t CommonMode(UInt_t i_mpd);
  Double_t ZeroSupThreshold(UInt_t i_mpd);
  
  // APV cross talk parameters
  static Int_t    fDoCrossTalk;  //whether we want to do cross talk simulation
  static Int_t    fNCStripApart; // # of strips the induced signal is away from the mean signal
  static Double_t fCrossFactor;  //reduction factor for the induced signal
  static Double_t fCrossSigma;   //uncertainty of the reduction factor
  
  //moved in "public" to allow it to compile with Root6/CentOS7
  struct IonPar_t {
    Double_t X;       // position of the point on the projection
    Double_t Y;
    Double_t Charge;  // Charge deposited by this ion
    Double_t SNorm;   // 3 x radius of ion diffusion area at readout
    Double_t R2;      // = SNorm^2 : radius of numerical integration area
    Double_t ggnorm;  // = Charge/R2/pi : charge per unit area
  };
  
 private:

  void IonModel (const TVector3& xi,
		 const TVector3& xo,
		 const Double_t elost );
  
  TGEMSBSGEMHit ** AvaModel (const Int_t ic,
			     const TGEMSBSSpec& spect,
			     const TVector3& xi,
			     const TVector3& xo,
			     const Double_t time_off);
  
  Double_t GetPedNoise(Double_t& phase, Double_t& amp, Int_t& isample);

  // Gas parameters
  Double_t fGasWion;               // eV
  Double_t fGasDiffusion;          // mm2/s
  Double_t fGasDriftVelocity;      // mm/s
  Double_t fAvalancheFiducialBand; // number of sigma defining the band around the avalanche in readout plane
  Int_t    fAvalancheChargeStatistics;  // 0 Furry, 1 Gaussian
  Double_t fGainMean;
  Double_t fGain0;
  UInt_t   fMaxNIon;               //maximum amount of ion pairs allowed in the digitization
  
  Double_t fSNormNsigma;           //fSNormNsigma is an arbitrary multiplicative factor for the avalance radius.
  Int_t    fAvaModel;              //0 for Heavyside, 1 for Gaussian, 2 for Cauchy-Lorentz
  Double_t fAvaGain;
  // Electronics parameters
  std::vector<Double_t> fTriggerOffset; // trigger offset (ns), incl latency & readout offset
  //Double_t fTriggerOffset;       // trigger offset (ns), incl latency & readout offset
  Double_t fTriggerJitter;       // trigger sigma jitter (ns)
  Double_t fAPVTimeJitter;       // time jitter associated with the APV internal clock
  Int_t    fEleSamplingPoints;
  Double_t fEleSamplingPeriod;   // ns
  Double_t fADCoffset;       // ADC offset
  Double_t fADCgain;         // ADC gain
  Int_t    fADCbits;         // ADC resolutions in bits
  Double_t fGateWidth;       // to be changed , ns - pulse shape width at ~1/10 max
  Int_t    fUseTrackerFrame;       // tracker frame is used in the original version, but not so in my version
                                   // Weizhi Xiong
  Double_t fEntranceRef;           // z position of the copper layer right before the first GEM gas layer,
                             // relative to the center of the GEM chamber
  Double_t fLateralUncertainty; // avalanche electrons can only pass through the holes of GEM foil
                                // which introduce additional uncertainty in the lateral direction

  //parameter for GEM pedestal noise
  Double_t fPulseNoiseSigma;  // additional sigma term of the pedestal noise
  Double_t fPulseNoisePeriod; // period of the pedestal noise, assuming sinusoidal function
  Double_t fPulseNoiseAmpConst;  // constant term of the pedestal noise amplitude
  Double_t fPulseNoiseAmpSigma;  // sigma term of the pedestal noise amplitude

  // Pulse shaping parameters
  Double_t fPulseShapeTau0;   // [ns] GEM model; = 50. in SiD model
  Double_t fPulseShapeTau1;   // [ns] GEM model only; if negative assume SiD model

  // Geometry
  Double_t fRoutZ;            // z-distance hit entrance to readout plane [mm]
  
  // Sector mapping // (???) are these relevant ? probably not (EFuchey 2016/12/05)
  Bool_t   fDoMapSector;
  Int_t    fSignalSector;

  //parameter for numerical integration
  UInt_t   fYIntegralStepsPerPitch;
  UInt_t   fXIntegralStepsPerPitch;
  
  TGEMSBSDigitizedPlane*** fDP; // 2D array of plane pointers indexed by chamber, plane #
  TGEMSBSGEMHit** fdh;// array of U & V GEM strips
  
  UInt_t fNChambers;  // # chambers
  UInt_t* fNPlanes;   // # planes in each chamber
  TRandom3 fTrnd;     // time randomizer
  UInt_t   fRNIon;    // number of ions
  std::vector<IonPar_t> fRIon;
  Double_t fRSMax;
  Double_t fRTotalCharge;
  Double_t fRTime0;
  Double_t fTimeZero;
  
  //std::vector<Double_t> fSumA;
  std::vector<Short_t>  fDADC;

  //zero suppression and common mode
  Bool_t fDoZeroSup;
  Double_t fZeroSup;
  Bool_t fDoCommonMode;
  std::vector<Double_t> fCommonModeArray;
  
  // Tree

  //TFile* fOFile;          // Output ROOT file
  //TTree* fOTree;          // Output tree
  //TGEMSBSSimEvent* fEvent;   // Output event structure, written to tree

  Bool_t fFilledStrips;   // True if no data changed since last SetTreeStrips

  void MakePrefix();
  void DeleteObjects();
  TGEMSBSDBManager *fManager;

 public:
  struct GEMCluster {
    Short_t   fID;        // Cluster number, cross-ref to GEMStrip
    // MC hit data
    Short_t   fSector;    // Sector number
    Short_t   fPlane;     // Plane number
    Short_t   fModule;    // Module number
    Short_t   fRealSector;// Real sector number (may be !=fSector if mapped)
    Short_t   fSource;    // MC data set source (0 = signal, >0 background)
    Int_t     fType;      // GEANT particle type (1 = primary)
    Int_t     fTRID;      // GEANT particle counter
    Int_t     fPID;       // PDG ID of particle generating the cluster
    TVector3  fP;         // Momentum of particle generating the cluster in the lab [GeV]
    TVector3  fPspec;     // Momentum of particle generating the cluster in the spec [GeV]
    TVector3  fXEntry;    // Track at chamber entrance in lab coords [m]
    TVector3  fMCpos;     // Approx. truth position of hit in lab [m]
    TVector3  fHitpos;    // fMCpos in Tracker frame [m]
    // Digitization results for this hit
    Float_t   fCharge;    // Charge of avalanche
    Float_t   fTime;      // Arrival time at electronics
    Int_t     fSize[2];   // Number of strips in cluster per axis
    Int_t     fStart[2];  // Number of first strip in cluster per axis
    Float_t   fXProj[2];  // fMCpos along projection axis [m]
    TVector3  fVertex;    // Vertex
  };

  std::vector<GEMCluster> fGEMClust;
  
  struct DigiGEMStrip {
    Short_t   fSector;    // Sector number
    Short_t   fPlane;     // Plane number
    Short_t   fModule;    // Module number
    Short_t   fProj;      // Readout coordinate ("x" = 0, "y" = 1)
    Short_t   fChan;      // Channel number
    Short_t   fSigType;   // Accumulated signal types (BIT(0) = signal)
    Float_t   fCharge;    // Total charge in strip
    Float_t   fTime1;     // Time of first sample
                          //   relative to event start in target (TBC)
    UShort_t  fNsamp;     // Number of ADC samples
    Int_t     fADC[MAXSAMP]; // [fNsamp] ADC samples
    //Int_t     fCommonMode[MC_MAMSAMP];// Real Common mode added to strip----going to work on next, needs to digitize all strips// seems no need to add strip to strip offset since it is just adding a constant in digitization and subtracting same known constant in analysis. "Strip specific pedestal rms and common mode are important"
    TArrayS   fClusters;  // Clusters ID contributing to signal on this strip
    TArrayI   fClusterRatio[MAXSAMP];
    TArrayD   fStripWeightInCluster;
  };


  ClassDef (TGEMSBSSimDigitization, 0)
};

#endif

