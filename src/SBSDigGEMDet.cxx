#include "SBSDigGEMDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

SBSDigGEMDet::SBSDigGEMDet()
{
}

SBSDigGEMDet::SBSDigGEMDet(UShort_t uniqueid, UInt_t nplanes, int* nstrips, double* offset, double* roangle, int nsamp, double zsup_thr):
  fUniqueID(uniqueid), fNPlanes(nplanes)
{
  //for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i] = SBSDigGEMPlane(nstrips[i], nsamp, zsup_thr);
  for(uint i = 0; i<fNPlanes; i++){
    cout << i << " " << nstrips[i] << " " << offset[i] << " " << roangle[i] << endl; 
    GEMPlanes.push_back(SBSDigGEMPlane(i/2, nstrips[i], nsamp, zsup_thr, offset[i], roangle[i]));
  }
}

SBSDigGEMDet::~SBSDigGEMDet()
{
  
}

void SBSDigGEMDet::Clear()
{
  for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i].Clear();
  fGEMhits.clear();
}
