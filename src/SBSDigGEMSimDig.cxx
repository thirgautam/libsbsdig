#include "SBSDigGEMSimDig.h"
#include "SBSDigGEMDet.h"
//#include "gmn_tree.h"
#include "g4sbs_tree.h"

//gas parameters
#define fGasWion 26             // eV
#define fGasDiffusion 1.e5       // mm2/s
#define fGasDriftVelocity 5.5e7   // mm/s
#define fAvalancheFiducialBand 10. // number of sigma defining the band around the avalanche in readout plane
#define fAvalancheChargeStatistics 1 // 0 Furry, 1 Gaussian
#define fGainMean 8.e3
#define fGain0 20.
#define fMaxNIon 1.e4               //maximum amount of ion pairs allowed in the digitization
  
#define fSNormNsigma 18.          //fSNormNsigma is an arbitrary multiplicative fact  
#define fAvaGain 20.
#define fLateralUncertainty 0.

//electronics parameters
#define fAPVTimeJitter 25.    // time jitter associated with the APV internal clock
  
#define fEleSamplingPoints 6
#define fEleSamplingPeriod 25. // ns
#define fADCoffset 0.         // ADC offset
#define fADCgain 1.          // ADC gain
#define fADCbits 12         // ADC resolutions in bits
#define fGateWidth 400.    // to be changed , ns - pulse shape width at ~1/10 max

//parameter for GEM pedestal noise
#define fPulseNoiseSigma 20.   // additional sigma term of the pedestal noise
#define fPulseNoisePeriod 200. // period of the pedestal noise, assuming sinusoidal function
#define fPulseNoiseAmpConst 0  // constant term of the pedestal noise amplitude
#define fPulseNoiseAmpSigma 0  // sigma term of the pedestal noise amplitude

//paramters for cross-talk
#define fNCStripApart 32  // # of strips the induced signal is away from the mean signal
#define fCrossFactor 0.1  //reduction factor for the induced signal
#define fCrossSigma 0.03  //uncertainty of the reduction factor

// Pulse shaping parameters
#define fPulseShapeTau 56.   // [ns] GEM model 0 = 50. in SiD model
//#define fPulseShapeTau1 0.   // [ns] GEM model only; if negative assume SiD model

#define fEntranceRef -1.5  // z position of the copper layer right before the first GEM gas layer,             // relative to the center of the GEM chamber
                         // which introduce additional uncertainty in the lateral direction
#define fRoutZ 9.185   // z-distance hit entrance to readout plane [mm]

//numerical integration parameters
#define fYIntegralStepsPerPitch 4
#define fXIntegralStepsPerPitch 4

#define fNROPlanes 2
#define fStripPitch 4.e-4

using namespace std;

//stupid but: let's do the defautl constructor:
SBSDigGEMSimDig::SBSDigGEMSimDig()// :
  //fGasWion(), fGasDiffusion(), fGasDriftVelocity(), fAvalancheFiducialBand(), fAvalancheChargeStatistics(), fGainMean(), fGain0(), fMaxNIon(), fSNormNsigma(), fAvaGain(), fAPVTimeJitter() 
{
}

SBSDigGEMSimDig::SBSDigGEMSimDig(int nchambers, double* trigoffset, double zsup_thr, int napv, double* commonmode_array) : fZeroSup(zsup_thr) 
{
  for(int i = 0; i<nchambers; i++){
    fTriggerOffset.push_back(trigoffset[i]);
    cout << i << "/" << nchambers << ": " << fTriggerOffset[i] << endl;
  }
  if(fZeroSup>0)fDoZeroSup = true;
  if(napv){
    fDoCommonMode = true;
    for(int i = 0; i<napv; i++){
      fCommonModeArray.push_back(commonmode_array[i]);
      cout << i << "/" << napv << ": " << fCommonModeArray[i] << endl;
    }
  }
  fRIon.resize((int)fMaxNIon);
  
  /*
  h1_QvsX_ion = new TH2D("h1_QvsX_ion", "", 250, -0.25, 0.25, 200, 0, 2.e4);
  h1_QvsY_ion = new TH2D("h1_QvsY_ion", "", 200, -0.2, 0.2, 200, 0, 2.e4);
  h1_QnormvsX_ion = new TH2D("h1_QnormvsX_ion", "", 250, -0.25, 0.25, 200, 0, 2.e4);
  h1_QnormvsY_ion = new TH2D("h1_QnormvsY_ion", "", 200, -0.2, 0.2, 200, 0, 2.e4);
  h1_QareavsX_ion = new TH2D("h1_QareavsX_ion", "", 750, -0.75, 0.75, 100, 0, 0.01);
  h1_QareavsY_ion = new TH2D("h1_QareavsY_ion", "", 200, -0.2, 0.2, 100, 0, 0.01);
  h1_QintvsX_ion = new TH2D("h1_QintvsX_ion", "", 250, -0.25, 0.25, 1000, 0, 1.e6);
  h1_QintvsY_ion = new TH2D("h1_QintvsY_ion", "", 200, -0.2, 0.2, 1000, 0, 1.e6);
  h1_QvsX_ava = new TH2D("h1_QvsX_ava", "", 750, -0.75, 0.75, 1000, 0, 1.e5);
  h1_QvsY_ava = new TH2D("h1_QvsY_ava", "", 200, -0.2, 0.2, 1000, 0, 1.e5);
  h1_QintYvsX_ava = new TH2D("h1_QintYvsX_ava", "", 750, -0.75, 0.75, 1000, 0, 1.e5);
  h1_QintYvsY_ava = new TH2D("h1_QintYvsY_ava", "", 200, -0.2, 0.2, 1000, 0, 1.e5);
  h1_yGEM_preion = new TH1D("h1_yGEM_preion", "", 200, -0.2, 0.2);
  h1_yGEM_preava = new TH1D("h1_yGEM_preava", "", 200, -0.2, 0.2);
  h1_yGEM_inava = new TH1D("h1_yGEM_inava", "", 200, -0.2, 0.2);
  h1_yGEM_inava_2 = new TH1D("h1_yGEM_inava_2", "", 200, -0.2, 0.2);
  h1_yGEM_inava_3 = new TH1D("h1_yGEM_inava_3", "", 200, -0.2, 0.2);
  h1_yGEM_inava_4 = new TH1D("h1_yGEM_inava_4", "", 200, -0.2, 0.2);
  h1_xGEMvsADC_inava_4 = new TH2D("h1_xGEMvsADC_inava_4", "", 750, -0.75, 0.75, 1000, -19, 20000-19);
  h1_yGEMvsADC_inava_4 = new TH2D("h1_yGEMvsADC_inava_4", "", 200, -0.2, 0.2, 1000, -19, 20000-19);
  h1_yGEM_incheckout = new TH1D("h1_yGEM_incheckout", "", 200, -0.2, 0.2);
  */
}


SBSDigGEMSimDig::~SBSDigGEMSimDig()
{
}


//.......................................................
// ionization Model
//
void
SBSDigGEMSimDig::IonModel(TRandom3* R,
			  const TVector3& xi,
			  const TVector3& xo,
			  const Double_t elost ) // eV
{
#define DBG_ION 0

  TVector3 vseg = xo-xi; // mm
  
  // ---- extract primary ions from Poisson
  fRNIon = R->Poisson(elost/fGasWion);

  if (fRNIon <=0)
    return;

#if DBG_ION > 0
  cout << "E lost = " << elost << ", " << fRNIon << " ions" << endl;
#endif
  if (fRNIon > fMaxNIon) {
#if DBG_ION > 0
    cout << __FUNCTION__ << ": WARNING: too many primary ions " << fRNIon << " limit to "
	 << fMaxNIon << endl;
#endif
    fRNIon = fMaxNIon;
  }

  fRSMax = 0.;
  fRTotalCharge = 0;
  fRTime0 = 999999.; // minimum time of drift

  for (UInt_t i=0;i<fRNIon;i++) { // first loop used to evaluate quantities
    IonPar_t ip;

    Double_t lion = R->Uniform(0.,1.); // position of the hit along the track segment (fraction)

    //In principle, the lateral uncertainty should have been put in the Ava model, but not here
    //But since we are not simulating the details of the avalanche, I think it is ok (Weizhi)
    ip.X = vseg.X()*lion+xi.X() + R->Gaus(0., fLateralUncertainty);
    ip.Y = vseg.Y()*lion+xi.Y() + R->Gaus(0., fLateralUncertainty);

    // Note the definition of fRoutZ is the distance from xi.Z() to xrout.Z():
    //        xi               xo   xrout
    // |<-LD->|<-----vseg----->|    |
    // |<-------fRoutZ---------|--->|
    // |      |<-lion*vseg->   |    |
    // |      |             <--LL-->|
    
    Double_t LD = TMath::Abs(xi.Z() - fEntranceRef);//usually should be 0,
                                            //unless particle is produced inside the gas layer

    Double_t LL = TMath::Abs(fRoutZ - LD - vseg.Z()*lion);
    Double_t ttime = LL/fGasDriftVelocity; // traveling time from the drift gap to the readout
    
    //cout << " rout Z  (mm?) " << fRoutZ << ", LD (mm?) " << LD << " vseg Z (mm?) " << vseg.Z()  << endl;
    //cout << " travelling length (mm?) " << LL << ", travelling time:  " <<  ttime << endl;
    
    fRTime0 = TMath::Min(ttime, fRTime0); // minimum traveling time [s]

    ip.SNorm = TMath::Sqrt(2.*fGasDiffusion*ttime); // spatial sigma on readout [mm]
    // cout<<"ip.SNorm: "<<ip.SNorm<<endl;
    //  cout<<"vseg: "<<vseg.X()<<" deltaX: "<<vseg.X()*lion<<endl;
    if( fAvalancheChargeStatistics == 1 ) {
      Double_t gnorm = fGainMean/TMath::Sqrt(fGain0); // overall gain TBC
      ip.Charge = R->Gaus(fGainMean, gnorm); // Gaussian distribution of the charge
    }
    else {
      ip.Charge = R->Exp(fGainMean); // Furry distribution
    }

    if( ip.Charge > 0 )
      fRTotalCharge += ip.Charge;
    else
      ip.Charge = 0;

    fRSMax = TMath::Max(ip.SNorm, fRSMax);

    // Derived quantities needed by the numerical integration in AvaModel
    ip.SNorm *= fSNormNsigma;
    ip.R2 = ip.SNorm * ip.SNorm;
    ip.ggnorm = ip.Charge * TMath::InvPi() / ip.R2; // normalized charge

#if DBG_ION > 1
    printf("z coords %f %f %f %f lion %f LL %lf\n",
	   xi.Z(), xo.Z(), vseg.Z(), lion, LL);
    printf("ttime = %e\n", ttime);
#endif
#if DBG_ION > 0
    cout << " x, y = " << ip.X << ", " << ip.Y << " snorm = "
	 << ip.SNorm/fSNormNsigma << " charge " << ip.Charge << endl;
    cout << "fRTime0 = " << fRTime0 << endl;
    cout << "fRion size " << fRIon.size() << " " << i << endl;
#endif
    /*
    h1_QvsX_ion->Fill(ip.X, ip.Charge);
    h1_QvsY_ion->Fill(ip.Y, ip.Charge);
    h1_QvsX_ion->Fill(ip.X, ip.ggnorm);
    h1_QvsY_ion->Fill(ip.Y, ip.ggnorm);
    */
    fRIon[i] = ip;
  }
  return;
}

Short_t 
SBSDigGEMSimDig::ADCConvert(Double_t val, Double_t off, Double_t gain, Int_t bits)
{
  // Convert analog value 'val' to integer ADC reading with 'bits' resolution
  assert( bits >= 0 && bits <= 12 );

  if( val < 0. )
    val = 0.;
  Double_t vvv = (val - off)/gain;
  //std::cout<<val<<" : "<<vvv<<std::endl;
  //printf("offset = %1.3f, gain = %1.3f, input value = %1.3f , output value = %1.3f  \n", off, gain, val, vvv);
  Double_t saturation = static_cast<Double_t>( (1<<bits)-1 );
  if( vvv > saturation )
    vvv = saturation;

  Short_t dval =
    static_cast<Short_t>( TMath::Floor( (vvv>saturation) ? saturation : vvv ));

  //  cerr << val << " dval = " << dval << endl;
  if( dval < 0 ) dval = 0;
  return dval;

}

// Pulse Shape SiD model
// APV25 time function from  M. Friedl et al NIMA 572 (2007) pg 385-387 (APV25 for silicon!)
//

Double_t 
SBSDigGEMSimDig::PulseShape(Double_t t, 
			    Double_t C,  // normalization factor
			    Double_t Tp) // shaping time 
{

  Double_t v;
  Double_t x;
  x = t/Tp;
  v = C/Tp * x * TMath::Exp(-x);
  
  return ( v>0. ) ? v : 0.;

}



//.......................................................
// avalanche model
//

//TGEMSBSGEMHit **
void SBSDigGEMSimDig::AvaModel(const int ic,
			       SBSDigGEMDet* gemdet, 
			       TRandom3* R,
			       const TVector3& xi,
			       const TVector3& xo,
			       const Double_t t0)
{
#define DBG_AVA 0
#if DBG_AVA > 0
  cout << "Chamber " << ic << "----------------------------------" << endl;
  cout << "In  " << xi.X() << " " << xi.Y() << " " << xi.Z() << endl;
  cout << "Out " << xo.X() << " " << xo.Y() << " " << xo.Z() << endl;
#endif

  // xi, xo are in chamber frame, in mm

  Double_t nsigma = fAvalancheFiducialBand; // coverage factor

#if DBG_AVA > 0
  cout << "fRSMax, nsigma " << fRSMax << " " << nsigma << endl;
#endif

  Double_t x0,y0,x1,y1; // lower and upper corners of avalanche diffusion area
  
  if (xi.X()<xo.X()) {
    x0 = xi.X()-nsigma*fRSMax;
    x1 = xo.X()+nsigma*fRSMax;
  } else {
    x1 = xi.X()+nsigma*fRSMax;
    x0 = xo.X()-nsigma*fRSMax;
  }

  if (xi.Y()< xo.Y()) {
    y0 = xi.Y()-nsigma*fRSMax;
    y1 = xo.Y()+nsigma*fRSMax;
  } else {
    y1 = xi.Y()+nsigma*fRSMax;
    y0 = xo.Y()-nsigma*fRSMax;
  }
  
  // Check if any part of the avalanche region is in the active area of the sector.
  // Here, "active area" means the chamber's *bounding box*, which is
  // larger than the wedge's active area (section of a ring)
  
  //const TGEMSBSGEMChamber& chamber = spect.GetChamber(ic);
  Double_t glx = (-gemdet->GEMPlanes[ic*2].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(0).GetStripLowerEdge(0)+chamber.GetPlane(0).GetSPitch()/2.0) * 1000.0;
  Double_t gly = (-gemdet->GEMPlanes[ic*2+1].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;
  Double_t gux = (gemdet->GEMPlanes[ic*2].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(0).GetStripUpperEdge(chamber.GetPlane(0).GetNStrips()-1) -chamber.GetPlane(0).GetSPitch()/2.0) * 1000.0;
  Double_t guy = (gemdet->GEMPlanes[ic*2+1].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(1).GetStripUpperEdge(chamber.GetPlane(1).GetNStrips()-1) -chamber.GetPlane(1).GetSPitch()/2.0) * 1000.0;
  
  if (x1<glx || x0>gux ||
      y1<gly || y0>guy) { // out of the sector's bounding box
    cerr << __FILE__ << " " << __FUNCTION__ << ": out of sector, "
	 << "chamber " << ic << " sector " << ic/30 << " plane " << ic%30 << endl
	 << "Following relations should hold:" << endl
	 << "(x1 " << x1 << ">glx " << glx << ") (x0 " << x0 << "<gux " << gux << ")" << endl
	 << "(y1 " << y1 << ">gly " << gly << ") (y0 " << y0 << "<guy " << guy << ")" << endl;
    //return 0;
  }
  
  //bool bb_clipped = (x0<glx||y0<gly||x1>gux||y1>guy);
  if(x0<glx) x0=glx;
  if(y0<gly) y0=gly;
  if(x1>gux) x1=gux;
  if(y1>guy) y1=guy;

  // Loop over chamber planes
  double roangle_mod, dx_mod, xoffset_mod;
  int GEMstrips;
  
  double xt_factor;
  int isLeft;
  //TGEMSBSGEMHit **virs;
  //virs = new TGEMSBSGEMHit *[fNROPlanes[ic]];
  for (UInt_t ipl = 0; ipl < fNROPlanes; ++ipl){
#if DBG_AVA > 0
    cout << "coordinate " << ipl << " =========================" << endl;
#endif
     
    xt_factor = R->Gaus(fCrossFactor, fCrossSigma);
    isLeft = R->Uniform(1.) < 0.5 ? -1 : 1;
    
    //cout << xt_factor << endl;
    
    // Compute strips affected by the avalanche
    //const SBSDigGEMPlane& pl = gemdet->GEMPlanes[ic*2+ipl];
    // strip angle = - strip angle for plane -> strip rotation purposes.
    roangle_mod = -gemdet->GEMPlanes[ic*2+ipl].ROangle();
    dx_mod = gemdet->GEMPlanes[ic*2+ipl].dX();
    xoffset_mod = gemdet->GEMPlanes[ic*2+ipl].Xoffset();
    GEMstrips = gemdet->GEMPlanes[ic*2+ipl].GetNStrips();
    
    // Positions in strip frame
    Double_t xs0 = x0*cos(roangle_mod) - y0*sin(roangle_mod);
    Double_t ys0 = x0*sin(roangle_mod) + y0*cos(roangle_mod);
    Double_t xs1 = x1*cos(roangle_mod) - y1*sin(roangle_mod); 
    Double_t ys1 = x1*sin(roangle_mod) + y1*cos(roangle_mod);
#if DBG_AVA > 0
    cout << "glx gly gux guy " << glx << " " << gly << " " << gux << " " << guy << endl;
    cout //<< ic << " " << ipl << " " << roangle_mod 
      << " xs0 ys0 xs1 ys1 " << xs0 << " " << ys0 << " " << xs1 << " " << ys1 << endl;
#endif
     //if(ipl==1 && ic<12)h1_yGEM_inava->Fill(xs0*1.e-3);

    Int_t iL = max(0, Int_t((xs0*1.e-3+dx_mod/2.)/fStripPitch) );
    iL = min(iL, GEMstrips);
    //pl.GetStrip (xs0 * 1e-3, ys0 * 1e-3);
    Int_t iU = min(Int_t((xs1*1.e-3+dx_mod/2.)/fStripPitch), GEMstrips);
    iU = max(0, iU);
    //pl.GetStrip (xs1 * 1e-3, ys1 * 1e-3);
     
    // Check for (part of) the avalanche area being outside of the strip region
    //if( (iL <= 0 && iU <= 0) || (iL>=pl.GetNStrips() && iU>=pl.GetNStrips()) ) {
      // All of the avalanche outside -> nothing to do
      // TODO: what if this happens for only one strip coordinate (ipl)?
// #if DBG_AVA > 0
//       cerr << __FILE__ << " " << __FUNCTION__ << ": out of active area, "
// 	   << "chamber " << ic << " sector " << ic%30 << " plane " << ic/30 << endl
// 	   << "iL_raw " << pl.GetStripUnchecked(xs0*1e-3) << " "
// 	   << "iU_raw " << pl.GetStripUnchecked(xs1*1e-3) << endl
// 	   << endl << endl;
// #endif
    if(iL==iU){//nothing to do
      return;
    }

    if(iU<iL)swap(iU, iL);


    //
    // Bounds of rectangular avalanche region, in strip frame
    //

    // Limits in x are low edge of first strip to high edge of last
    Double_t xl = (iL*fStripPitch-dx_mod/2.)*1.e3;//pl.GetStripLowerEdge (iL) * 1000.0;
    Double_t xr = (iU*fStripPitch-dx_mod/2.)*1.e3;//pl.GetStripUpperEdge (iU) * 1000.0;

#if DBG_AVA > 0
    cout << "iL gsle " << iL << " " << xl << endl;
    cout << "iU gsue " << iU << " " << xr << endl;
#endif
    
    //if(ipl==1 && ic<12)h1_yGEM_inava_2->Fill((xl+xr)*5.e-4);

    // Limits in y are y limits of track plus some reasonable margin
    // We do this in units of strip pitch for convenience (even though
    // this is the direction orthogonal to the pitch direction)

    // Use y-integration step size of 1/10 of strip pitch (in mm)
    Double_t yq = fStripPitch * 1000.0 / fYIntegralStepsPerPitch;
    Double_t yb = ys0, yt = ys1;
    if (yb > yt)
      swap( yb, yt );
    yb = yq * TMath::Floor (yb / yq);
    yt = yq * TMath::Ceil  (yt / yq);

    // We should also allow x to have variable bin size based on the db
    // the new avalanche model (Cauchy-Lorentz) has a very sharp full width
    // half maximum, so if the bin size is too large, it can introduce
    // fairly large error on the charge deposition. Setting fXIntegralStepsPerPitch
    // to 1 will go back to the original version -- Weizhi Xiong

    Int_t nstrips = iU - iL + 1;
    Int_t nx = (iU - iL + 1) * fXIntegralStepsPerPitch;
    Int_t ny = TMath::Nint( (yt - yb)/yq );
#if DBG_AVA > 0
    cout << "xr xl yt yb nx ny "
	 << xr << " " << xl << " " << yt << " " << yb
	 << " " << nx << " " << ny << endl;
#endif
    assert( nx > 0 && ny > 0 );

    // define function, gaussian and sum of gaussian

    Double_t xbw = (xr - xl) / nx;
    Double_t ybw = (yt - yb) / ny;
#if DBG_AVA > 0
    cout << "xbw ybw " << xbw << " " << ybw << endl;
#endif
    
    Int_t sumASize = nx * ny;
#if DBG_AVA > 0
    cout<<nx<<" : "<<ny<< ", nx*ny " << sumASize <<" Nstrips: "<<nstrips<<endl;
#endif
    fSumA.resize(sumASize);
    memset (&fSumA[0], 0, fSumA.size() * sizeof (Double_t));
    //Double_t fSumA[sumASize];
    //memset (fSumA, 0, sumASize * sizeof (Double_t));
    //for(Int_t i=0; i<sumASize; i++){
    //fSumA[i] = 0;
    //}
    
    //    fSumA.resize(nx*ny);
    //memset (&fSumA[0], 0, fSumA.size() * sizeof (Double_t));
#if DBG_AVA > 0
    cout << fRNIon << " " << fRIon.size() << endl;
#endif
    
    for (UInt_t i = 0; i < fRNIon; i++){
      Double_t frxs = fRIon[i].X*cos(roangle_mod) - fRIon[i].Y*sin(roangle_mod);
      Double_t frys = fRIon[i].X*sin(roangle_mod) + fRIon[i].Y*cos(roangle_mod);
      //pl.PlaneToStrip (frxs, frys);
      //frxs *= 1e3; frys *= 1e3;
      //  cout<<"IonStrip: "<<pl.GetStrip(frxs*1e-3,frys*1e-3)<<endl;
      // bin containing center and # bins each side to process
      Int_t ix = (frxs-xl) / xbw;
      Int_t iy = (frys-yb) / ybw;
      Int_t dx = fRIon[i].SNorm / xbw  + 1;
      Int_t dy = fRIon[i].SNorm / ybw  + 1;
#if DBG_AVA > 1
       cout << "ix dx iy dy " << ix << " " << dx << " " << iy << " " << dy << endl;
#endif

       //if(ipl==1 && ic<12)h1_yGEM_inava_3->Fill(frxs*1.e-3);

    
      //
      // NL change:
      //
      // ggnorm is the avalance charge for the i^th ion, and R2 is the square of the radius of the diffusion 
      // circle, mutiplied by the kSNormNsigma factor: (ip.SNorm * ip.SNorm)*kSNormNsigma*kSNormNsigma. All 
      // strips falling within this circle are considered in charge summing. 
      //
      // The charge contribution to a given strip by the i^th ion is evaluated by a Lorentzian (or Gaussian)
      // distribution; the sigma for this distribution is eff_sigma, which is the actual avalance sigma. 
      //
      Double_t ggnorm = fRIon[i].ggnorm;
      Double_t r2 = fRIon[i].R2;
      Double_t eff_sigma_square = r2/(fSNormNsigma*fSNormNsigma);
      Double_t eff_sigma = TMath::Sqrt(eff_sigma_square);
      Double_t current_ion_amplitude = fAvaGain*ggnorm*(1./(TMath::Pi()*eff_sigma))*(eff_sigma*eff_sigma);
      
      //if(ipl==0 && ic<12)h1_QnormvsX_ion->Fill(frxs*1.e-3, current_ion_amplitude);
      //if(ipl==1 && ic<12)h1_QnormvsY_ion->Fill(frxs*1.e-3, current_ion_amplitude);
      
      // xc and yc are center of current bin
      double sumA = 0;
      // Loop over bins
      Int_t min_binNb_x = max(ix-dx,0);
      Int_t min_binNb_y = max(iy-dy,0);
      Int_t max_binNb_x = min(ix+dx+1,nx);
      Int_t max_binNb_y = min(iy+dy+1,ny);
      Int_t jx = min_binNb_x;
      Double_t xc = xl + (jx+0.5) * xbw;
      for (; jx < max_binNb_x; ++jx, xc += xbw){
	Double_t xd2 = frxs-xc; xd2 *= xd2;
	//if( xd2 > r2 ){
	//  if( (xc - frxs)>0 )
	//    break;
	//  else
	//    continue;
	//}
	Int_t jx_base = jx * ny;
	Int_t jy = min_binNb_y;
	Double_t yc = yb + (jy+0.5) * ybw;
	
	for (; jy < max_binNb_y; ++jy, yc += ybw){
	  Double_t yd2 = frys-yc; yd2 *= yd2;
	  //if( yd2 > r2 ){
	  //  if( (yc - frys)>0 )
	  //    break;
	  //  else
	  //    continue;
	  // }
	  if( xd2 + yd2 <= r2 ) {
	    //cout << current_ion_amplitude << " " << xd2+yd2 << " " << eff_sigma_square << endl;
	    fSumA[jx_base+jy] += current_ion_amplitude / ((xd2+yd2)+eff_sigma_square);
	    sumA+= current_ion_amplitude / ((xd2+yd2)+eff_sigma_square);
	  }
	}//cout<<endl;
      }//cout<<"##########################################################################"<<endl<<endl;getchar();
      /*
      if(ipl==0 && ic<12){
	h1_QintvsX_ion->Fill(frxs*1.e-3, sumA);
      }
      if(ipl==1 && ic<12){
	h1_QintvsY_ion->Fill(frxs*1.e-3, sumA);
      }
      */
    }
    
#if DBG_AVA > 0
    cout << "t0 = " << t0 << " plane " << ipl 
	 << endl;
#endif

    //virs[ipl] = new TGEMSBSGEMHit(nx,fEleSamplingPoints);
    //virs[ipl]->SetTime(t0);
    //virs[ipl]->SetHitCharge(fRTotalCharge);
    
    //Int_t ai=0;
    Double_t area = xbw * ybw;

//when we integrate in order to get the signal pulse, we want all charge
    //deposition on the area of a single strip -- Weizhi
    
    //cout << "number of strips: " << nstrips << ", number of samples " << fEleSamplingPoints << " area: " << area << endl;
    
    // if(nstrips>0){cout<<"nstrips: "<<nstrips<<" Nion: "<<fRNIon<<endl;}
    
    for (Int_t j = 0; j < nstrips; j++){
      //  cout<<"strip: "<<iL+j<<":    ";
      //Int_t posflag = 0;
      Double_t us = 0.;
      for (UInt_t k=0; k<fXIntegralStepsPerPitch; k++){
	Double_t integralY_tmp = 0;
	int kx = (j * fXIntegralStepsPerPitch + k) * ny;
	for( Int_t jy = ny; jy != 0; --jy )
	  integralY_tmp += fSumA[kx++];
	/*
	if(ipl==0 && ic<12)h1_QintYvsX_ava->Fill((iL+j)*fStripPitch-dx_mod/2.+xoffset_mod, integralY_tmp);
	if(ipl==1 && ic<12)h1_QintYvsY_ava->Fill((iL+j)*fStripPitch-dx_mod/2.+xoffset_mod, integralY_tmp);
	*/
	/*
	if(integralY_tmp==0){
	  cout << " ipl " << ipl << " ny " << ny << " k " << k << " kx_0 " << (j * fXIntegralStepsPerPitch + k) * ny << " -> kx " << kx << " => " << kx -(j * fXIntegralStepsPerPitch + k) * ny  << endl;
	}
	*/
	us += integralY_tmp * area;
	//	us += IntegralY( fSumA, j * fXIntegralStepsPerPitch + k, nx, ny ) * area;
	//if(us>0)cout << "k " << k << ", us " << us << endl;
      }
      /*
      if(ipl==0 && ic<12){
	h1_QvsX_ava->Fill((iL+j)*fStripPitch-dx_mod/2.+xoffset_mod, us);
	h1_QareavsX_ion->Fill((iL+j)*fStripPitch-dx_mod/2.+xoffset_mod, area);
      }
      if(ipl==1 && ic<12){
	h1_QvsY_ava->Fill((iL+j)*fStripPitch-dx_mod/2.+xoffset_mod, us);
	h1_QareavsY_ion->Fill((iL+j)*fStripPitch-dx_mod/2.+xoffset_mod, area);
      }
      */
#if DBG_AVA > 0
      cout << "strip " << iL+j << " us " << us << endl;
#endif
      //if(ipl==1 && ic<12)h1_yGEM_inava_4->Fill( (iL+j)*fStripPitch-dx_mod/2. );
      
      // cout <<setw(6)<< (Int_t)(us/100);
      // cout<<iL+j<<" : "<<us<<endl;
      //generate the random pedestal phase and amplitude
      // Double_t phase = fTrnd.Uniform(0., fPulseNoisePeriod);
      // Double_t amp = fPulseNoiseAmpConst + fTrnd.Gaus(0., fPulseNoiseAmpSigma);
      
      for (Int_t b = 0; b < fEleSamplingPoints; b++){
	Double_t pulse = PulseShape (fEleSamplingPeriod * b - t0,
				     us,
				     fPulseShapeTau);
	//fPulseShapeTau0, fPulseShapeTau1 );
	
	Short_t dadc = ADCConvert( pulse,
				   0,// fADCoffset,
				   fADCgain,
				   fADCbits );
	/*
#if DBG_AVA > 0
	if(pulse>0)
	  cout << "strip number " << iL+j << ", sampling number " << b << ", t0 = " << t0 << endl
	       << "pulse = " << pulse << ", (val - off)/gain = " 
	       << (pulse-fADCoffset)/fADCgain << ", dadc = " << dadc << endl;
#endif
	*/
	//fDADC[b] = dadc;
	//cout << ic*2+ipl << " " << iL+j << " " << gemdet->GEMPlanes[ic*2+ipl].GetADC(iL+j, b) << " " << gemdet->GEMPlanes[ic*2+ipl].GetADCSum(iL+j) << " ==> ";
	gemdet->GEMPlanes[ic*2+ipl].AddADC(iL+j, b, dadc);
	//if(gemdet->GEMPlanes[ic*2+ipl].GetADC(iL+j, b)<0 || gemdet->GEMPlanes[ic*2+ipl].GetADC(iL+j, b)>4096)cout << " hou  " << iL+j << " " << b << " " << gemdet->GEMPlanes[ic*2+ipl].GetADC(iL+j, b) << endl;
	//posflag += dadc;
	//if(dadc>0)cout << t0 << " " << pulse << " " << dadc << endl;
	//cross talk here ?
	if(xt_factor>0){
	  if(iL+j+isLeft*fNCStripApart>=0 && iL+j+isLeft*fNCStripApart<GEMstrips){
	    //cout << "induced: " << iL+j+isLeft*fNCStripApart << " " << gemdet->GEMPlanes[ic*2+ipl].GetADC(iL+j+isLeft*fNCStripApart, b) << " " << gemdet->GEMPlanes[ic*2+ipl].GetADCSum(iL+j+isLeft*fNCStripApart) << " ==> ";
	    gemdet->GEMPlanes[ic*2+ipl].AddADC(iL+j+isLeft*fNCStripApart, b, TMath::Nint(dadc*xt_factor));
	    //cout << gemdet->GEMPlanes[ic*2+ipl].GetADC(iL+j+isLeft*fNCStripApart, b) << " " << gemdet->GEMPlanes[ic*2+ipl].GetADCSum(iL+j+isLeft*fNCStripApart) << endl;
	  }
	}
	//cout << ic*2+ipl << " " << iL+j << " " << b << " " << gemdet->GEMPlanes[ic*2+ipl].GetADC(iL+j, b) << "; " << iL+j+isLeft*fNCStripApart << " " << gemdet->GEMPlanes[ic*2+ipl].GetADC(iL+j+isLeft*fNCStripApart, b)<< endl;
      }//cout <<"  "<<t0<< endl;

      /*
      if(ipl==0 && ic<12)h1_xGEMvsADC_inava_4->Fill( (iL+j)*fStripPitch-dx_mod/2.+xoffset_mod , gemdet->GEMPlanes[ic*2+ipl].GetADCSum(iL+j) );
      if(ipl==1 && ic<12)h1_yGEMvsADC_inava_4->Fill( (iL+j)*fStripPitch-dx_mod/2.+xoffset_mod , gemdet->GEMPlanes[ic*2+ipl].GetADCSum(iL+j) );
      */
    }//end loop on strips
        
    // if(ic*2+ipl==30){
    //   for(int j = 570; j<580; j++){
    // 	cout << j << "   ";
    // 	for(int b = 0; b<6; b++){
    // 	  cout << " " << gemdet->GEMPlanes[ic*2+ipl].GetADC(j, b);
    // 	}cout << endl;
    //   }
    // }
    
  }//end loop on planes
  
  //return virs;
}


Int_t
SBSDigGEMSimDig::Digitize (SBSDigGEMDet* gemdet,
			   TRandom3* R)//, 
//gmn_tree* T)
{
  // Digitize event. Add results to any existing digitized data.

  //UInt_t nh = gdata.GetNHit();
  bool is_background = false;
  Float_t event_time=0,time_zero=0;
  Double_t trigger_jitter = R->Uniform(-fAPVTimeJitter/2, fAPVTimeJitter/2);
  
  // For signal data, determine the sector of the primary track
  
  for(size_t ih = 0; ih<gemdet->fGEMhits.size(); ih++){
    is_background = (gemdet->fGEMhits[ih].source==0);
    UInt_t igem = gemdet->fGEMhits[ih].module;
    //UInt_t igem = iplane/2;
    
    //if(igem>=16)cout << igem << endl;
    //cout<<igem<<":"<<imodule<<":"<<iplane<<endl;
    //Short_t itype = (gdata.GetParticleType(ih)==1) ? 1 : 2; // primary = 1, secondaries = 2
    // if(gdata.GetParticleType(ih)!=1){cout<<"x"<<endl;getchar();}

    // These vectors are in the spec frame, we need them in the chamber frame
    TVector3 vv1(gemdet->fGEMhits[ih].xin*1.e3,
		 gemdet->fGEMhits[ih].yin*1.e3, 
		 gemdet->fGEMhits[ih].zin*1.e3);
    TVector3 vv2(gemdet->fGEMhits[ih].xout*1.e3,
		 gemdet->fGEMhits[ih].yout*1.e3, 
		 gemdet->fGEMhits[ih].zout*1.e3);
    
    if(abs(vv1.X()-vv2.X())>50 || abs(vv1.Y()-vv2.Y())>50){//in mm
      //cout<<abs(vv1.X()-vv2.X())<<endl;
      //getchar();
      continue;
    }
    //if(igem<12)h1_yGEM_preion->Fill(vv1.Y()*1.e-3);
    IonModel (R, vv1, vv2, gemdet->fGEMhits[ih].edep );
    
    // Get Signal Start Time 'time_zero'
    if( is_background ) {
      // For background data, uniformly randomize event time between
      // -fGateWidth to +75 ns (assuming 3 useful 25 ns samples).
      // Not using HitTime from simulation file but randomize HitTime to cycle use background files
      //event_time = m(-fGateWidth, 6*fEleSamplingPeriod);
      event_time = fTimeZero;//fTrnd.Uniform(-fGateWidth/2.-fEleSamplingPeriod, fGateWidth-fEleSamplingPeriod);
      //event_time = fTrnd.Uniform((-fGateWidth+2*fEleSamplingPeriod), 8*fEleSamplingPeriod);
    } else {
      // Signal events occur at t = 0, 
      event_time = fTimeZero+gemdet->fGEMhits[ih].t;
    }
    //  cout<<event_time<<"  "<<ih<<endl;
    // Adding drift time and trigger_jitter
    time_zero = event_time - fTriggerOffset[igem] + fRTime0*1e9 - trigger_jitter;
    
    //cout << time_zero << " " << fTimeZero << " " << gemdet->fGEMhits[ih].t 
    //<< " " << trigger_jitter << " " << fRTime0*1e9 << endl;
    
#if DBG_AVA > 0
    if(time_zero>200.0)
      cout << "time_zero " << time_zero 
	   << "; evt time " << event_time 
	   << "; hit time " << gemdet->fGEMhits[ih].t
	   << "; drift time " << fRTime0*1e9
	   << endl;
#endif
    if (fRNIon > 0) {
      //cout << "AvaModel..." << endl;
      //if(igem<12)h1_yGEM_preava->Fill(vv1.Y()*1.e-3);
      AvaModel (igem, gemdet, R, vv1, vv2, time_zero);
      //cout << "Done!" << endl;
      //cout << " hou " << gemdet->GEMPlanes[4].GetADCSum(400) << endl;
      //CheckOut(gemdet, R, T);
    }
    
  }//end loop on hits
  //fFilledStrips = false;
  
  //Cumulate(gemdet, R);
    
  return 0;
}

//-------------------------------------------------------
// Helper functions for integration in AvaModel
inline static
Double_t IntegralY( Double_t* a, Int_t ix, Int_t nx, Int_t ny )
{
  register double sum = 0.;
  register int kx = ix*ny;
  for( Int_t jy = ny; jy != 0; --jy )
    sum += a[kx++];

  return sum;
}


void SBSDigGEMSimDig::CheckOut(SBSDigGEMDet* gemdet,
			       const int uniqueid, 
			       TRandom3* R, 
			       g4sbs_tree* T)
			       //gmn_tree* T)
{
  //int test = gemdet->GEMPlanes[4].GetADCSum(400);
  //cout << " hou hou " << test << endl;
  // for(int j = 570; j<580; j++){
  //   cout << j << "   ";
  //   for(int b = 0; b<6; b++){
  //     cout << " " << gemdet->GEMPlanes[30].GetADC(j, b);
  //   }cout << endl;
  // }
  
  short strip;
  double commonmode = 0;
  if(fDoCommonMode){commonmode = fCommonModeArray[0];}
  //cout << "commonmode " << commonmode << endl;
  int apv_ctr;
  for(size_t i = 0; i<gemdet->GEMPlanes.size(); i++){
    
    for(int j = 0; j<gemdet->GEMPlanes[i].GetNStrips(); j++){
      //if(gemdet->GEMPlanes[4].GetADCSum(400)!=test){
      //cout << gemdet->GEMPlanes[4].GetADCSum(400) << "!=" << test << ": " << i << " " << j << endl;
      //test = gemdet->GEMPlanes[4].GetADCSum(400);
      //}
      //cout << fDoCommonMode << endl;
      if(fDoCommonMode){
	//cout << " hou " << endl;
	if(j%128==0 && apv_ctr<fCommonModeArray.size()){
	  commonmode = fCommonModeArray[apv_ctr++];
	  //cout << commonmode << endl;
	}
      }
      //if(i==4 && j==400){
      //cout << gemdet->GEMPlanes[i].GetADCSum(j) << endl;
      //}
      if(gemdet->GEMPlanes[i].GetADCSum(j)>0){
	//if(i%2==1 && i<24)h1_yGEM_incheckout->Fill(j*fStripPitch-gemdet->GEMPlanes[i].dX()/2.);
	for(int k = 0; k<6; k++){
	  //cout << i << " " << j << " " << k << " " << gemdet->GEMPlanes[i].GetADC(j, k) << " " << gemdet->GEMPlanes[i].GetADCSum(j) << " = = > ";
	  //int ped = TMath::Nint(R->Gaus(commonmode, fPulseNoiseSigma));
	  //if(gemdet->GEMPlanes[i].GetADC(j, k)<0 || gemdet->GEMPlanes[i].GetADC(j, k)>4096)cout << i << " " << j << " " << k << " " << gemdet->GEMPlanes[i].GetADC(j, k) << " " << ped << " " << commonmode << " " << fPulseNoiseSigma << endl;
	  gemdet->GEMPlanes[i].AddADC(j, k, TMath::Nint(R->Gaus(commonmode, fPulseNoiseSigma)));
	  //cout << gemdet->GEMPlanes[i].GetADC(j, k) << " " << gemdet->GEMPlanes[i].GetADCSum(j)<< endl;
	  //handle saturation
	  if(gemdet->GEMPlanes[i].GetADC(j, k)>pow(2, fADCbits) ){
	    //cout << gemdet->GEMPlanes[i].GetADC(j, k) << " => ";
	    gemdet->GEMPlanes[i].SetADC(j, k, TMath::Nint(pow(2, fADCbits)) );
	    //cout << gemdet->GEMPlanes[i].GetADC(j, k) << endl;
	  }
	}
	//#ifdef 
	if( (fDoZeroSup && gemdet->GEMPlanes[i].GetADCSum(j)-commonmode*6>fZeroSup) || !fDoZeroSup) {
	  //if(i<4)cout << i << " " << gemdet->GEMPlanes[i].GetNStrips() << " " << commonmode << endl;
	  if(uniqueid==BBGEM_UNIQUE_DETID){
	     T->Earm_BBGEM_Dig.nstrips++;
	     T->Earm_BBGEM_Dig.module->push_back(i);
	     T->Earm_BBGEM_Dig.strip->push_back(j);
	     T->Earm_BBGEM_Dig.adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	     T->Earm_BBGEM_Dig.adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	     T->Earm_BBGEM_Dig.adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	     T->Earm_BBGEM_Dig.adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	     T->Earm_BBGEM_Dig.adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	     T->Earm_BBGEM_Dig.adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	  }
	  //if(gemdet->GEMPlanes[i].GetADCSum(j)-commonmode*6>fZeroSup){
	  //FillBBGEMTree(gemdet->GEMPlanes[i], T, j);
	  //#ifdef QUENOUILLE	 
	  //for(int k = 0; k<6; k++){
	  //if(gemdet->GEMPlanes[i].GetADC(j, k)>4096 || gemdet->GEMPlanes[i].GetADC(j, k)<0)cout << i << " " << j << " " << k << " " << gemdet->GEMPlanes[i].GetADC(j, k) << endl;
	  //}
	  
	  /*
	  if(gemdet->GEMPlanes[i].Module()<3){
	    strip = j+gemdet->GEMPlanes[i].GetNStrips()*gemdet->GEMPlanes[i].Module();
	    if(gemdet->GEMPlanes[i].ROangle()==0){
	      T->Earm_BBGEM_1x_dighit_nstrips++;
	      T->Earm_BBGEM_1x_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_1x_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_1x_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_1x_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_1x_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_1x_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_1x_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	      //cout << gemdet->GEMPlanes[i].GetADC(j, 2) << " " << T->Earm_BBGEM_1x_dighit_nstrips << " " << T->Earm_BBGEM_1x_dighit_adc_2->size()-1 << " " << T->Earm_BBGEM_1x_dighit_adc_2->at(T->Earm_BBGEM_1x_dighit_adc_2->size()-1) << endl;
	    }else{
	      T->Earm_BBGEM_1y_dighit_nstrips++;
	      T->Earm_BBGEM_1y_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_1y_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_1y_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_1y_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_1y_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_1y_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_1y_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }
	  }else if(gemdet->GEMPlanes[i].Module()<6){
	    strip = j+gemdet->GEMPlanes[i].GetNStrips()*(gemdet->GEMPlanes[i].Module()-3);
	    if(gemdet->GEMPlanes[i].ROangle()==0){
	      T->Earm_BBGEM_2x_dighit_nstrips++;
	      T->Earm_BBGEM_2x_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_2x_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_2x_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_2x_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_2x_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_2x_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_2x_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }else{
	      T->Earm_BBGEM_2y_dighit_nstrips++;
	      T->Earm_BBGEM_2y_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_2y_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_2y_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_2y_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_2y_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_2y_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_2y_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }
	  }else if(gemdet->GEMPlanes[i].Module()<9){
	    strip = j+gemdet->GEMPlanes[i].GetNStrips()*(gemdet->GEMPlanes[i].Module()-6);
	    if(gemdet->GEMPlanes[i].ROangle()==0){
	      T->Earm_BBGEM_3x_dighit_nstrips++;
	      T->Earm_BBGEM_3x_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_3x_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_3x_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_3x_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_3x_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_3x_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_3x_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }else{
	      T->Earm_BBGEM_3y_dighit_nstrips++;
	      T->Earm_BBGEM_3y_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_3y_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_3y_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_3y_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_3y_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_3y_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_3y_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }
	  }else if(gemdet->GEMPlanes[i].Module()<12){
	    strip = j+gemdet->GEMPlanes[i].GetNStrips()*(gemdet->GEMPlanes[i].Module()-9);
	    if(gemdet->GEMPlanes[i].ROangle()==0){
	      T->Earm_BBGEM_4x_dighit_nstrips++;
	      T->Earm_BBGEM_4x_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_4x_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_4x_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_4x_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_4x_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_4x_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_4x_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }else{
	      T->Earm_BBGEM_4y_dighit_nstrips++;
	      T->Earm_BBGEM_4y_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_4y_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_4y_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_4y_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_4y_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_4y_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_4y_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }
	  }else{
	    strip = j+gemdet->GEMPlanes[i].GetNStrips()*(gemdet->GEMPlanes[i].Module()-12);
	    if(gemdet->GEMPlanes[i].ROangle()==0){
	      T->Earm_BBGEM_5x_dighit_nstrips++;
	      T->Earm_BBGEM_5x_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_5x_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_5x_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_5x_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_5x_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_5x_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_5x_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }else{
	      T->Earm_BBGEM_5y_dighit_nstrips++;
	      T->Earm_BBGEM_5y_dighit_strip->push_back(strip);
	      T->Earm_BBGEM_5y_dighit_adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->Earm_BBGEM_5y_dighit_adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->Earm_BBGEM_5y_dighit_adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->Earm_BBGEM_5y_dighit_adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->Earm_BBGEM_5y_dighit_adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->Earm_BBGEM_5y_dighit_adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	    }
	  }
	  */
	}//end if(...)
	//#endif
	//}else{
	//FillBBGEMTree(gemdet->GEMPlanes[i], T, j); 
	//}
      }
    }
  }  
}

/*
void SBSDigGEMSimDig::FillBBGEMTree(const SBSDigGEMPlane pl, gmn_tree* T, int j)
{
  short strip;
  if(pl.Module()<3){
    strip = j+pl.GetNStrips()*pl.Module();
    if(pl.ROangle()==0){
      T->Earm_BBGEM_1x_dighit_nstrips++;
      T->Earm_BBGEM_1x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_1x_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_1x_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_1x_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_1x_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_1x_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_1x_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }else{
      T->Earm_BBGEM_1y_dighit_nstrips++;
      T->Earm_BBGEM_1y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_1y_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_1y_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_1y_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_1y_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_1y_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_1y_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }
  }else if(pl.Module()<6){
    strip = j+pl.GetNStrips()*(pl.Module()-3);
    if(pl.ROangle()==0){
      T->Earm_BBGEM_2x_dighit_nstrips++;
      T->Earm_BBGEM_2x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_2x_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_2x_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_2x_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_2x_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_2x_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_2x_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }else{
      T->Earm_BBGEM_2y_dighit_nstrips++;
      T->Earm_BBGEM_2y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_2y_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_2y_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_2y_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_2y_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_2y_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_2y_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }
  }else if(pl.Module()<9){
    strip = j+pl.GetNStrips()*(pl.Module()-6);
    if(pl.ROangle()==0){
      T->Earm_BBGEM_3x_dighit_nstrips++;
      T->Earm_BBGEM_3x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_3x_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_3x_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_3x_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_3x_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_3x_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_3x_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }else{
      T->Earm_BBGEM_3y_dighit_nstrips++;
      T->Earm_BBGEM_3y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_3y_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_3y_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_3y_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_3y_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_3y_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_3y_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }
  }else if(pl.Module()<12){
    strip = j+pl.GetNStrips()*(pl.Module()-9);
    if(pl.ROangle()==0){
      T->Earm_BBGEM_4x_dighit_nstrips++;
      T->Earm_BBGEM_4x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_4x_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_4x_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_4x_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_4x_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_4x_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_4x_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }else{
      T->Earm_BBGEM_4y_dighit_nstrips++;
      T->Earm_BBGEM_4y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_4y_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_4y_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_4y_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_4y_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_4y_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_4y_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }
  }else{
    strip = j+pl.GetNStrips()*(pl.Module()-12);
    if(pl.ROangle()==0){
      T->Earm_BBGEM_5x_dighit_nstrips++;
      T->Earm_BBGEM_5x_dighit_strip->push_back(strip);
      T->Earm_BBGEM_5x_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_5x_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_5x_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_5x_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_5x_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_5x_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }else{
      T->Earm_BBGEM_5y_dighit_nstrips++;
      T->Earm_BBGEM_5y_dighit_strip->push_back(strip);
      T->Earm_BBGEM_5y_dighit_adc_0->push_back(pl.GetADC(j, 0));
      T->Earm_BBGEM_5y_dighit_adc_1->push_back(pl.GetADC(j, 1));
      T->Earm_BBGEM_5y_dighit_adc_2->push_back(pl.GetADC(j, 2));
      T->Earm_BBGEM_5y_dighit_adc_3->push_back(pl.GetADC(j, 3));
      T->Earm_BBGEM_5y_dighit_adc_4->push_back(pl.GetADC(j, 4));
      T->Earm_BBGEM_5y_dighit_adc_5->push_back(pl.GetADC(j, 5));
    }
  }
  
}
*/

//___________________________________________________________________________________
void SBSDigGEMSimDig::Print()
{
  cout << "GEM digitization:" << endl;
  cout << "  Gas parameters:" << endl;
  cout << "    Gas ion width: " << fGasWion << endl;
  cout << "    Gas diffusion: " << fGasDiffusion << endl;
  cout << "    Gas drift velocity: " << fGasDriftVelocity << endl;
  cout << "    Avalanche fiducial band: " << fAvalancheFiducialBand << endl;
  cout << "    Avalanche charge statistics: " << fAvalancheChargeStatistics << endl;
  cout << "    Gain mean: " << fGainMean << endl;
  cout << "    Gain 0: " << fGain0 << endl;

  cout << "  Electronics parameters:" << endl;
  //cout << "    Trigger offsets: "; //<< fTriggerOffset 
  //for(int i = 0; i<fManager->GetNChamber(); i++)cout << fTriggerOffset[i] << " ";
  //cout << endl;
  cout << "    APV time jitter: " << fAPVTimeJitter << endl;
  cout << "    Sampling Period: " << fEleSamplingPeriod << endl;
  cout << "    Sampling Points: " << fEleSamplingPoints   << endl;
  cout << "    Pulse Noise width: " << fPulseNoiseSigma << endl;
  cout << "    ADC offset: " << fADCoffset << endl;
  cout << "    ADC gain: " << fADCgain << endl;
  cout << "    ADC bits: " << fADCbits << endl;
  cout << "    Gate width: " << fGateWidth << endl;

  cout << "  Pulse shaping parameters:" << endl;
  cout << "    Pulse shape tau: " << fPulseShapeTau << endl;
  //cout << "    Pulse shape tau1: " << fPulseShapeTau1 << endl;
}

void SBSDigGEMSimDig::write_histos()
{
  /*
  h1_QvsX_ion->Write();
  h1_QvsY_ion->Write();
  h1_QnormvsX_ion->Write();
  h1_QnormvsY_ion->Write();
  h1_QareavsX_ion->Write();
  h1_QareavsY_ion->Write();
  h1_QintvsX_ion->Write();
  h1_QintvsY_ion->Write();
  h1_QvsX_ava->Write();
  h1_QvsY_ava->Write();
  h1_QintYvsX_ava->Write();
  h1_QintYvsY_ava->Write();
  h1_yGEM_preion->Write();
  h1_yGEM_preava->Write();
  h1_yGEM_inava->Write();
  h1_yGEM_inava_2->Write();
  h1_yGEM_inava_3->Write();
  h1_yGEM_inava_4->Write();
  h1_xGEMvsADC_inava_4->Write();
  h1_yGEMvsADC_inava_4->Write();
  h1_yGEM_incheckout->Write();
  */
}

/*
Double_t
TGEMSBSSimDigitization::CommonMode(UInt_t i_mpd)
{
  if(fCommonModeArray.size() && fDoCommonMode){
    i_mpd = (i_mpd<fCommonModeArray.size() ? i_mpd: 0);
    return fCommonModeArray[i_mpd];
  }else{
    return 0;
  }
}

Double_t
TGEMSBSSimDigitization::ZeroSupThreshold(UInt_t i_mpd)
{
  //i_mpd
  if(fDoZeroSup){
    if(fDoCommonMode){
      return fZeroSup+CommonMode(i_mpd)*fEleSamplingPoints;
    }else{
      return fZeroSup;
    }
  }else{
    return -1000;
  }
}
*/
