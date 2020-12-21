#ifndef SBSDIGAUXI_H
#define SBSDIGAUXI_H

#include <iostream>
#include <vector>
#include <map>
#include "gmn_tree.h"
#include "g4sbs_tree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "SBSDigPMTDet.h"
#include "SBSDigGEMDet.h"

bool UnfoldData(g4sbs_tree* T, double theta_sbs, double d_hcal, TRandom3* R, 
		std::vector<SBSDigPMTDet*> pmtdets, 
		std::vector<int> detmap, 
		std::vector<SBSDigGEMDet*> gemdets, 
		std::vector<int> gemmap, 
		//std::map<int, SBSDigPMTDet*> pmtdets, 
		//std::map<int, SBSDigGEMDet*> gemdets, 
		double tzero,
		int signal);

//bool FillDigTree(gmn_tree* T, std::map<int, SBSDigPMTDet*> pmtdets, std::map<int, SBSDigGEMDet*> gemdets);

#endif // SBSDIGAUXI_H

