#include "SBSDigPMTDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigPMTDet::SBSDigPMTDet()
{
}

SBSDigPMTDet::SBSDigPMTDet(UShort_t uniqueid, UInt_t nchan):
  fUniqueID(uniqueid), fNChan(nchan)
{
  //for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal();
  for(int i = 0; i<fNChan; i++)PMTmap.push_back(PMTSignal());
}

SBSDigPMTDet::SBSDigPMTDet(UShort_t uniqueid, UInt_t nchan, double NpeChargeConv, double sigmapulse, double gatewidth):
  fUniqueID(uniqueid), fNChan(nchan)
{
  //for(int i = 0; i<fNChan; i++)PMTmap[i] = PMTSignal(NpeChargeConv);
  for(int i = 0; i<fNChan; i++)PMTmap.push_back(PMTSignal(NpeChargeConv));
  fRefPulse = new SPEModel(fUniqueID, sigmapulse, 0, -gatewidth/2., gatewidth/2.);
}

SBSDigPMTDet::~SBSDigPMTDet()
{
  
}

//void SBSDigPMTDet::Digitize(gmn_tree* T, TRandom3* R)
void SBSDigPMTDet::Digitize(g4sbs_tree* T, TRandom3* R)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Digitize(i, fUniqueID, T, R, fPedestal, fPedSigma, fADCconv, fADCbits, fTDCconv, fTDCbits);
}
  
void SBSDigPMTDet::SetSamples(double sampsize)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].SetSamples(-fGateWidth/2, fGateWidth/2, sampsize);
}

void SBSDigPMTDet::Clear(bool dosamples)
{
  for(int i = 0; i<fNChan; i++)PMTmap[i].Clear(dosamples);
}
