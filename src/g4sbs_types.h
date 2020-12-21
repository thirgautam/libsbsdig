#ifndef __G4SBS_TYPES_H
#define __G4SBS_TYPES_H

#include <vector>
#include <string>

////////////////////////////////////////////////////////
//  Data for extracting things from GEMC
//
//  we'll hardcode them here, but it would be nice to
//  maybe get them into a database
//  I guess this could also be done through a mysql
//  interface, but I think that makes it more complicated
//  and breakable

#define qe 1.602e-19
#define spe_unit 1.0e-9 //to convert ns to s...
#define m_e 511.e-6
#define n_lg 1.68
#define ROimpedance 50.0 //Ohm

// List of detector unique IDs: 
// by (proposed) convention: DetUniqueID = DetType*10+DetID
// DetType of type det_type defined in g4sbs_types: kHCal(0), kECal(1), kCher(2), kScint(3), kGEM(4);
// TODO: also put those in the detector DB, and have DB manager ensure no detector share and indentical unique ID
#define HCAL_UNIQUE_DETID 0
#define BBPS_UNIQUE_DETID 10
#define BBSH_UNIQUE_DETID 11
#define ECAL_UNIQUE_DETID 12
#define GRINCH_UNIQUE_DETID 20
#define RICH_UNIQUE_DETID 21
#define HODO_UNIQUE_DETID 30
#define CDET_UNIQUE_DETID 31
#define ACTIVEANA_UNIQUE_DETID 32
#define PRPOLBS_SCINT_UNIQUE_DETID 33
#define PRPOLFS_SCINT_UNIQUE_DETID 34
#define BBGEM_UNIQUE_DETID 40
#define SBSGEM_UNIQUE_DETID 41
#define FT_UNIQUE_DETID 42
#define FPP1_UNIQUE_DETID 43
#define FPP2_UNIQUE_DETID 44
#define CEPOL_GEMFRONT_UNIQUE_DETID 45
#define CEPOL_GEMREAR_UNIQUE_DETID 46
#define PRPOLBS_GEM_UNIQUE_DETID 47
#define PRPOLFS_GEM_UNIQUE_DETID 48

/*
enum exp_type{
  kNeutronExp, kGEp, kGEnRP, 
  kSIDIS, kA1n, kTDIS, kDVCS
};

enum det_type{
  kHCal, kECal,
  kCher, kScint,
  kGEM
};
*/

const std::string kProj_str[2] = {"x", "y"};
const double bbgem_z[5] = {0.85, 1.0, 1.15, 1.30, 2.33223544};
const double cepol_front_z[4] = {0.85, 1.0, 1.15, 1.30};
const double cepol_rear_z[4] = {0.85, 1.0, 1.15, 1.30};

#endif//__GEMC_TYPES_H
