#ifndef TGEMSBSDBMANAGER_H
#define TGEMSBSDBMANAGER_H

#include "THaAnalysisObject.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "Rtypes.h"
#include "VarDef.h"
#include "TMath.h"
#include "types.h"

class TGEMSBSSpec;

class TGEMSBSDBManager : public THaAnalysisObject {
public:
    ~TGEMSBSDBManager();
    TGEMSBSDBManager(const char *spec = "bb", const char *det = "gem");
    //TODO: Remove this static instance!
    //static TGEMSBSDBManager* GetInstance() {
    //    if (fManager == NULL) fManager = new TGEMSBSDBManager();
    //    return fManager;
    //}
    
    void LoadGeneralInfo(const std::string& fileName);
    void LoadGeoInfo(const std::string& fileName);
    std::vector<int>& GetChanMap() { return fChanMap; }
   
    /* int       DoMapSector() const          { return fDoMapSector;         } */
    /* int       DoSelfDefineSector() const   { return fDoSelfDefinedSector; } */
    /* int       GetSectorMapped() const      { return fMappedSector;        } */
    int       GetNChamber() const          { return fNChamber;            }
    int       GetNSector() const           { return fNSector;             }
    int       GetNGEMPlane() const         { return fNGEMPlane;           }
    int       GetNModule(int plane) const  { if(plane<0) return 0; else return fNModule[plane]; }
    int       GetPlaneID(int igem)         { return fmIgemtoPlane[igem];  }
    int       GetModuleID(int igem)        { return fmIgemtoModule[igem]; }
    int       GetGEMID(int ip,int im)      { return fmPMtoIgem[ip][im];   }
    /* // see comment l. 93-95. */
    /* int      GetNChamber2() const          { return fNChamber2;            } */
    /* int      GetNSector2() const           { return fNSector2;             } */
    int       GetNReadOut() const          { return fNReadOut;            }
    /* int       GetGEMDriftID() const        { return fGEMDriftID;          } */
    /* int       GetGEMCopperFrontID() const  { return fGEMCopperFrontID;    } */
    /* int       GetGEMCopperBackID() const   { return fGEMCopperBackID;     } */
    /* int       GetGEMStripID() const        { return fGEMStripID;          } */
    /* int       GetNSigParticle() const      { return fNSigParticle;        } */
    /* int       GetFAECID() const            { return fFAECID;              } */
    /* int       GetLAECID() const            { return fLAECID;              } */
    int       GetChanPerSlot() const       { return fChanPerSlot;         }
    int       GetModulesPerReadOut() const { return fModulesPerReadOut;   }
    int       GetModulesPerChamber() const { return fModulesPerChamber;   }
    int       GetChambersPerCrate() const  { return fChambersPerCrate;    }
    
    int       GetSigPID(unsigned int i) const;
    int       GetSigTID(unsigned int i) const;
    
    int       Getg4sbsDetectorType() const { return fg4sbsDetectorType;   }
    double    Getg4sbsZSpecOffset() const  { return fg4sbsZSpecOffset;    }
    double    GetZ0() const                { return fgZ0;                 }
    
    //double    GetCaloThreshold() const     { return fCaloThr;             }
    //double    GetCaloZ() const             { return fgCaloZ;              }
    //double    GetCaloRes() const           { return fgCaloRes;            }
    //int       DoCalo() const               { return fgDoCalo;             }

    void     SetZ0( Double_t z0 ) { fgZ0 = z0; }
    // Support for calorimeter emulation. Static functions to allow script access
    //void     EmulateCalorimeter( Bool_t f = true ) { fgDoCalo = f; }
    //void     SetCaloZ( Double_t z )     { fgCaloZ   = z; }
    //void     SetCaloRes( Double_t res ) { fgCaloRes = res; }
    
    double    GetDMag(int i, int j);
    double    GetThetaV(int i, int j);
    double    GetD0(int i, int j);
    double    GetXOffset(int i, int j);
    double    GetDepth(int i, int j);
    double    GetDX(int i, int j);
    double    GetDY(int i, int j);
    double    GetStripAngle(int i, int j, int k);
    double    GetPitch(int i, int j, int k);
    
    int GetModuleIDFromPos(int iplane, double x, double y = 0);
    double GetPosFromModuleStrip(int iproj, int iplane, int isector, int istrip);
    UInt_t GetGlobalStripPlane(uint lstrip, int plane, int module, int proj);
    void SetDBFileName(const std::string &filename) { fDBFileName = filename; }
    const std::string& GetDBFileName() { return fDBFileName; }

    const std::string& GetPrefix() { return fPrefix; }
    const std::string& GetSpecName() { return fSpecName; }
    const std::string& GetDetName() { return fDetName; }

    TGEMSBSSpec& GetSpec() { return *fSpec; }
    int GetNChan() { return fNChan; }

    void InitializeGEMs();
     
    void GetPMfromGlobalPlaneNum(uint gplanenum, int& plane, int& module);
    
protected:
    //int    LoadDB(std::ifstream& inp, DBRequest* request, const std::string& prefix);
    //std::string FindKey( std::ifstream& inp, const std::string& key ) const;
    bool   CheckIndex(int i, int j=0, int k=0) const;
    
    //static TGEMSBSDBManager* fManager;

    //variable for data base information
    /* int fDoMapSector; */
    /* int fMappedSector; */
    /* int fDoSelfDefinedSector; */
    
    int    fNChamber;
    std::vector<std::string> fChambers;
    int    fNSector;
    int    fNGEMPlane;
    std::vector<Int_t> fNModule;
    std::vector<std::vector<std::string> > fModules;
    std::map<Int_t,std::map<Int_t, Int_t> > fmPMtoIgem;
    std::map<Int_t, Int_t> fmIgemtoPlane;
    std::map<Int_t, Int_t> fmIgemtoModule;

    int    fNReadOut;
    int    fNSigParticle;
    int    fChanPerSlot;
    int    fModulesPerReadOut;
    int    fModulesPerChamber;
    int    fChambersPerCrate;
    
    // Parameters for TGEMSBSGeant4File
    int fg4sbsDetectorType;// flag to determine which type of GEM should be read.
    //Choices are: 1 - BB GEMs
    //Choices are: 2 - SIDIS SBS GEMs
    //Choices are: 3 - GEP SBS GEMs: FT
    //Choices are: 3 - GEP SBS GEMs: FPP
    double fg4sbsZSpecOffset;
    double fgZ0;        // z position of first chamber plane
    //Offset between the local z value recorded in g4sbs and the actual distance 
    //of the GEM from the midplane pivot

    // Parameters for simple calorimeter analysis
    // Parameters for TGEMSBSSimDecoder
    // Calorimeter emulation
    // double fCaloThr;    
    // double fgCaloZ;     // z position of emulated calorimeter
    // double fgCaloRes;   // Resolution (sigma) of emulated calorimeter (m)
    // int    fgDoCalo;    // Enable calorimeter emulation
    
    int    fErrID;
    double fErrVal;
    
    std::vector<int>    fSigPID;
    std::vector<int>    fSigTID;
    
    /* vector<double> fChamberZ; */
    //std::map< int, std::vector<GeoInfo> > fGeoInfo;
    std::map< int, std::vector<GeoInfo> > fPMGeoInfo; //plane module format geo info
    std::string fSpecName;
    std::string fDetName;
    std::string fDBFileName;
    std::string fPrefix;

    TGEMSBSSpec *fSpec;
    int fNChan;
    std::vector<int> fChanMap;

    ClassDef(TGEMSBSDBManager,1);
};

#endif
