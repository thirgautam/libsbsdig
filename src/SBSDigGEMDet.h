#ifndef SBSDIGGEMDET_H
#define SBSDIGGEMDET_H

#include <iostream>
#include <vector>
#include <map>
#include "SBSDigGEMPlane.h"

//________________________________
class SBSDigGEMDet {
 public:
  SBSDigGEMDet();
  SBSDigGEMDet(UShort_t uinqueid, UInt_t nplanes, int* nstrips, double* offset, double* roangle, int nsamp, double zsup_thr);
  virtual ~SBSDigGEMDet();
  void Clear();

  struct gemhit{
    int source;
    int module;
    double edep;
    //double tmin;
    //double tmax;
    double t;
    double xin;
    double yin;
    double zin;
    double xout;
    double yout;
    double zout;
  };

  std::vector<gemhit> fGEMhits;
    
  //private:
  UShort_t fUniqueID;
  UInt_t fNPlanes;
  //std::map<int, SBSDigGEMPlane> GEMPlanes;
  std::vector<SBSDigGEMPlane> GEMPlanes;
};

#endif
