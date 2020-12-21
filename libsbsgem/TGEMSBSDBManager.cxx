#include "TGEMSBSDBManager.h"
#include "TGEMSBSSpec.h"
#include "TGEMSBSGEMChamber.h"
#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector2.h"
#include "TRandom3.h"
#include "TSystem.h"

ClassImp(TGEMSBSDBManager);

using namespace std;

//TGEMSBSDBManager * TGEMSBSDBManager::fManager = NULL;

TGEMSBSDBManager::TGEMSBSDBManager(const char *spec, const char* det)
  : //fDoMapSector(0), fMappedSector(0), fDoSelfDefinedSector(0),
    fNChamber(0), fNSector(0), fNGEMPlane(0), fNReadOut(0), fNSigParticle(0),
    fChanPerSlot(2048), fModulesPerReadOut(1), fModulesPerChamber(1), fChambersPerCrate(1),
    fg4sbsDetectorType(0), fg4sbsZSpecOffset(0), fgZ0(0),
    //fCaloThr(0), fgCaloZ(0), fgCaloRes(0), fgDoCalo(0), 
    fErrID(-999), fErrVal(-999.), fSpecName(spec), fDetName(det)
{
  fPrefix = fSpecName+"."+fDetName;
}
//______________________________________________________________
TGEMSBSDBManager::~TGEMSBSDBManager()
{
}

//______________________________________________________________
void TGEMSBSDBManager::InitializeGEMs()
{
  if(fDebug>=2)cout << "TGEMSBSDBManager::InitializeGEMs()" << endl;
  // Make spectrometer
  fSpec = new TGEMSBSSpec(fPrefix.c_str(),
      Form("Temporary GEM spectrometer for %s",fDetName.c_str()));
  // And make all the GEM chambers for this spectrometer
  TGEMSBSGEMChamber *dGEM;
  fNChan = 0;
  for(Int_t plane = 0; plane < GetNGEMPlane(); plane++) {
    for(Int_t mod = 0; mod < GetNModule(plane); mod++) {
      dGEM = new TGEMSBSGEMChamber(Form("%s.%s",
					/*GetPrefix().c_str(),*/
					fChambers[plane].c_str(),
					fModules[plane][mod].c_str()),
				   Form("%s: SimGEMChamber on %s.%s",
					fDetName.c_str(),fChambers[plane].c_str(),
					fModules[plane][mod].c_str())
				   );
      dGEM->SetGeometry(GetD0(plane, mod),
			GetXOffset(plane, mod),
			GetDepth(plane, mod),
			GetDX(plane, mod),
			GetDY(plane, mod),
			GetDMag(plane, mod),
			GetThetaV(plane, mod));
      
      //dGEM->Print();
      string coord[2] = {"x", "y"};
      for(int k = 0; k<fNReadOut; k++){
	if(fDebug>=2)cout << fChambers[plane].c_str() << " " << fModules[plane][mod].c_str() << " " << coord[k].c_str() << endl;
	dGEM->InitPlane(k, 
			Form("%s.%s.%s",
			     /*GetPrefix().c_str(),*/
			     fChambers[plane].c_str(),
			     fModules[plane][mod].c_str(),
			     coord[k].c_str()),
			Form("%s: SimGEMChamber on %s.%s.%s",
			     fDetName.c_str(),fChambers[plane].c_str(),
			     fModules[plane][mod].c_str(), 
			     coord[k].c_str()),
			GetStripAngle(plane, mod, k), GetPitch(plane, mod, k));
	
      }
      if(fDebug>=2)dGEM->Print();
      dGEM->SetApparatus(fSpec);
      if(dGEM->Init()) { // true == error
        if(fDebug>=1)std::cerr << "ERROR!: TGEMSBSDBManager::InitializeGEMs() error error initializing GEM: " << Form("%s.%s.%s", fPrefix.c_str(), fChambers[plane].c_str(), fModules[plane][mod].c_str()) << std::endl;
      } else {
        // Get total number of strips (channels)
        fNChan += dGEM->GetNStripTotal();
        fSpec->AddGEM(dGEM);
      }
    }
  }
   if(fDebug>=1)std::cout << "Initialized " << fSpecName << "." << fDetName << " with " << fNChan << " GEM strips." << std::endl;
}


//______________________________________________________________
/*
static bool OpenInput( const string& filename, ifstream& ifs )
{
  // Open input stream 'ifs' for file 'filename'.
  // Look first in current directory, then in $DB_DIR, then in
  // $LIBSBSGEM/db

  ifs.close();

  struct FileLoc {
    FileLoc() : env(0) {}
    const char* env;
    string subdir;
  };
  const int N = 4;
  FileLoc fileloc[N];
  fileloc[0].env = "";
  fileloc[1].env = gSystem->Getenv("DB_DIR");
  fileloc[2].env = gSystem->Getenv("LIBSBSGEM");
  fileloc[2].subdir = "db";
  fileloc[3].env = gSystem->Getenv("SBS_DIGI_DB");

  for( int i=0; i<N; i++ ) {
    FileLoc& f = fileloc[i];
    if( !f.env )
      continue;
    string path(f.env);
    if( !path.empty() )
      path += "/";
    if( !f.subdir.empty() )
      path += f.subdir + "/";
    path += filename;

    ifs.clear();
    ifs.open(path.c_str());
    if( ifs.good() )
      return true;
  }
  return false;
}
*/

//______________________________________________________________
void TGEMSBSDBManager::LoadGeneralInfo(const string& fileName)
{  
  //The "sector-plane" concept is not suitable for SBS GEMTrackers. 
  //Instead, "Plane-Module" is introduced. "Plane" means tracking plane and 
  //"Module" means a independent GEM module which is a sub division of the "Plane"
  
    FILE *input = OpenFile( fileName.c_str(), GetInitDate() );
    //if ( !OpenInput(fileName,input, GetInitDate()) ){
    if(!input) {
        cout<<"cannot find general information file "<<fileName
            <<". Exiting the program"<<endl;
        exit(0);
    }
    //const string prefix = "generalinfo.";
    const string prefix = fPrefix + ".info.";

    //std::vector<Int_t>* NModule = new vector<Int_t>;
    DBRequest request[] = {
        // {"do_map_sector",       &fDoMapSector         , kInt,    0, 1},
        // {"self_define_sector",  &fDoSelfDefinedSector , kInt,    0, 1},
        // {"sector_mapped",       &fMappedSector        , kInt,    0, 1},
	{"nchamber",            &fNChamber            , kInt,    0, 1},
        {"nsector",             &fNSector             , kInt,    0, 1},
	//{"nplane",              &fNGEMPlane           , kInt,    0, 1},
	//{"nmodule",             NModule               , kIntV        },
        {"nreadout",            &fNReadOut            , kInt,    0, 1},
        {"nsignal",             &fNSigParticle        , kInt,    0, 1},
        {"chan_per_slot",       &fChanPerSlot         , kInt,    0, 1},
        {"modules_per_readout", &fModulesPerReadOut   , kInt,    0, 1},
	{"g4sbs_detectortype",  &fg4sbsDetectorType   , kInt,    0, 1},
	{"g4sbs_z_specoffset",  &fg4sbsZSpecOffset    , kDouble, 0, 1},
	{"z0",                  &fgZ0                 , kDouble, 0, 1},
	// {"calo_thr",            &fCaloThr             , kDouble, 0, 1},
	// {"calo_z",              &fgCaloZ              , kDouble, 0, 1},
	// {"calo_res",            &fgCaloRes            , kDouble, 0, 1},
	// {"docalo",              &fgDoCalo             , kInt,    0, 1},
        { 0 }
    };
    int pid, tid;
    DBRequest signalRequest[] = {
        {"pid",                 &pid,                   kInt, 0, 1},
        {"tid",                 &tid,                   kInt, 0, 1},
        { 0 }
    };

    int err = LoadDB( input, GetInitDate(), request,  prefix.c_str());
    if( err ) {cout<<"Load DB error"<<endl;exit(2);} 
    
    /*
    if(fNGEMPlane!=(int)NModule->size()) {
      cout<<"Check consistency of number of GEM Planes"<<endl;
      exit(2);
    }
    int nGEMtot=0;
    for(int i=0;i<fNGEMPlane;i++)
      {
	int nmodule = NModule->at(i);
	fNModule.push_back(nmodule);
	for(int j=0;j<nmodule;j++)
	  {
	    fmPMtoIgem[i][j]=nGEMtot;
	    fmIgemtoPlane[nGEMtot]=i;
	    fmIgemtoModule[nGEMtot]=j;
	    nGEMtot++;
	  }
      }
    for (int i=0; i<GetNGEMPlane(); i++){
      vector<GeoInfo> thisInfo;
      thisInfo.clear();
      fPMGeoInfo[i] = thisInfo;
    }

    //delete NModule;

    */
   
    
    for (int i=0; i<fNSigParticle; i++){
        ostringstream signal_prefix(prefix, ios_base::ate);
        signal_prefix<<"signal"<<i+1<<".";
        
        err = LoadDB(input, GetInitDate(), signalRequest, signal_prefix.str().c_str());
        
        fSigPID.push_back(pid);
        fSigTID.push_back(tid);
	
	if( err ) exit(2); 
    }
        

    //fModulesPerChamber = fModulesPerReadOut * fNReadOut;
    
    // fChambersPerCrate = 
    // (TGEMSBSSimDecoder::GetMAXSLOT()/fModulesPerChamber/fNChamber) * fNChamber;
    //input.close();
    fclose(input);
}

void TGEMSBSDBManager::LoadGeoInfo(const string& fileName)
{
  //const string& fileName = "db_"+prefix+".dat";
  //const string prefix = "geo."+fSpecName+"."+fDetName+".";
  const string prefix = fSpecName+"."+fDetName+".";
    
  //ifstream input;
  FILE *input = OpenFile(fileName.c_str(), GetInitDate());;
  //if( !OpenInput(fileName,input) ) {
  if(!input) {
    cout<<"cannot find geometry file "<<fileName
	<<". Exiting the program"<<endl;
    exit(0);
  } else {
    cerr << "Opened file " << fileName << " yay!!!" << std::endl;
  }

  // First, get the number of chambers

  std::string chambers = "";
  DBRequest chambers_request[] = {
    // Chambers are called "Tracker" planes in this simulation, apparently
    {"chambers", &chambers, kString}, ///< REQUIRED!
    { 0 }
  };
  int err = LoadDB(input, GetInitDate(), chambers_request, prefix.c_str());
  if(err) {
    std::cerr << "Got this error: " << err << std::endl;
    exit(2);
    return;
  }
  fChambers.clear();
  fChambers = THaAnalysisObject::vsplit(chambers);

  GeoInfo thisGeo;
  
  DBRequest request[] = {
    {"dmag",        &thisGeo.dmag,         kDouble, 0, 1},
    {"d0",          &thisGeo.d0,           kDouble, 0, 1},
    {"xoffset",     &thisGeo.xoffset,      kDouble, 0, 1},
    {"dx",          &thisGeo.dx,           kDouble, 0, 1},
    {"dy",          &thisGeo.dy,           kDouble, 0, 1},
    {"thetaV",      &thisGeo.thetaV,       kDouble, 0, 1},
    {"depth",       &thisGeo.depth,        kDouble, 0, 1},
    { 0 }
  };
  
  std::vector<Int_t> x_chanmap;
  std::vector<Int_t> y_chanmap;
  fChanMap.clear();
  DBRequest plane_request[] = {
    { "x.stripangle",     &thisGeo.stripangle_u,   kDouble, 0, 1},
    { "x.pitch",          &thisGeo.pitch_u,        kDouble, 0, 1},
    { "y.stripangle",     &thisGeo.stripangle_v,   kDouble, 0, 1},
    { "y.pitch",          &thisGeo.pitch_v,        kDouble, 0, 1},
    { "x.chanmap",        &x_chanmap,      kIntV, 0, true},
    { "y.chanmap",        &y_chanmap,      kIntV, 0, true},
    { 0 }
  };
  std::string modules;
  DBRequest modules_request[] = {
    {"modules", &modules, kString, 0, false }, ///< REQUIRED!
    { 0 }
  };
  fNGEMPlane = fChambers.size();
  fModules.resize(fNGEMPlane);
  int nGEMtot=0;
  
  fChanMap.clear();
  for (int i=0; i<fNGEMPlane; i++){
    vector<GeoInfo> thisInfoVec;
    fPMGeoInfo[i] = thisInfoVec;
    ostringstream plane_prefix(prefix, ios_base::ate);
    plane_prefix<<fChambers[i]<<".";
    err = LoadDB(input, GetInitDate(), modules_request, plane_prefix.str().c_str());
    if( err) exit(2);
    fModules[i] = THaAnalysisObject::vsplit(modules);
    fNModule.push_back(fModules[i].size());
    
    for (int j=0; j<fNModule[i]; j++){
      fmPMtoIgem[i][j]=nGEMtot;
      fmIgemtoPlane[nGEMtot]=i;
      fmIgemtoModule[nGEMtot]=j;
      nGEMtot++;

      ostringstream module_prefix(plane_prefix.str(), ios_base::ate);
      //      int idx = j;
      module_prefix<<fModules[i][j]<<".";
      
      int err = LoadDB(input, GetInitDate(), request, module_prefix.str().c_str());
      if( err ) exit(2);
      thisGeo.thetaV*= TMath::DegToRad();
      x_chanmap.clear();
      y_chanmap.clear();
      err = LoadDB(input, GetInitDate(), plane_request, module_prefix.str().c_str());
      thisGeo.stripangle_u*= TMath::DegToRad();
      thisGeo.stripangle_v*= TMath::DegToRad();
      
      if (err) exit(2);

      // Now process the channel map
      for(std::vector<Int_t>::iterator it = x_chanmap.begin(); it !=
          x_chanmap.end(); it++) {
        fChanMap.push_back(*it);
      }
      for(std::vector<Int_t>::iterator it = y_chanmap.begin(); it !=
          y_chanmap.end(); it++) {
        fChanMap.push_back(*it);
      }
      
      if(fDebug>=2)cout << "plane " << i << " module " << j << endl 
			<< " thisGeo.dmag: " << thisGeo.dmag 
			<< " thisGeo.d0 " << thisGeo.d0 << endl
			<< " thisGeo.xoffset " << thisGeo.xoffset
			<< " thisGeo.dx " << thisGeo.dx
			<< " thisGeo.dy " << thisGeo.dy << endl
			<< " thisGeo.thetaV " << thisGeo.thetaV
			<< " thisGeo.z " << thisGeo.z
			<< " thisGeo.depth " << thisGeo.depth << endl
			<< " thisGeo.stripangle_u " << thisGeo.stripangle_u
			<< " thisGeo.stripangle_v " << thisGeo.stripangle_v
			<< " thisGeo.pitch_u " << thisGeo.pitch_u
			<< " thisGeo.pitch_v " << thisGeo.pitch_v << endl;
      
      fPMGeoInfo[i].push_back(thisGeo);
    }
  }
  // And now that we are done, process the channel map
  //input.close();
  fclose(input);
}



/*
//______________________________________________________________
string TGEMSBSDBManager::FindKey( ifstream& inp, const string& key ) const
{
  static const string empty("");
  string line;
  string::size_type keylen = key.size();
  inp.seekg(0); // could probably be more efficient, but it's fast enough
  while( getline(inp,line) ) {
    if( line.size() <= keylen )
      continue;
    if( line.compare(0,keylen,key) == 0 ) {
      if( keylen < line.size() ) {
	string::size_type pos = line.find_first_not_of(" \t=", keylen);
	if( pos != string::npos )
	  return line.substr(pos);
      }
      break;
    }
  }
  return empty;
}
*/
//_________________________________________________________________________
bool TGEMSBSDBManager::CheckIndex(int i, int j, int k) const//(plane, module, readoutAxis)
{
    if (i >= fNChamber || i < 0){
        cout<<"invalid chamber ID requested: "<<i<<endl;
        return false;
    }
    else if(j>=fNModule[i]|| j<0){
      cout<<"invalid module id requested: "<<j<<endl;
      return false;
    }
    else if (k >= fNReadOut || k < 0){
        cout<<"invalid readout id requested: "<<k<<endl;
    }
    return true;
}
/*
//_________________________________________________________________
int TGEMSBSDBManager::LoadDB( ifstream& inp, DBRequest* request, const string& prefix )
{
  DBRequest* item = request;
  while( item->name ) {
    ostringstream sn(prefix, ios_base::ate);
    sn << item->name;
    const string& key = sn.str();
    string val = FindKey(inp,key);
    Int_t tempval;
    if( !val.empty() ) {
      istringstream sv(val);
      switch(item->type){
        case kString:
          *((std::string*)item->var) = val;
          break;
        case kDouble:
          sv >> *((double*)item->var);
          break;
        case kInt:
          sv >> *((Int_t*)item->var);
          break;
        case kIntV:
	  while(1){
	    if(!sv.good()) break;
	    sv >> tempval;
            ((std::vector<Int_t>*)item->var)->push_back(tempval);
	  }
  	  break;
  
        default:
          return 1;
        break;
      }
      if( !sv ) {
	cerr << "Error converting key/value = " << key << "/" << val << endl;
	return 1;
      }
    } else {
      cerr << "key \"" << key << "\" not found" << endl;
      return 2;
    }
    ++item;
  }
  return 0;
}
*/
//_____________________________________________________________________
int TGEMSBSDBManager::GetSigPID(unsigned int i) const
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigPID[i];
}
//______________________________________________________________________
int TGEMSBSDBManager::GetSigTID(unsigned int i) const
{
    if ( i >= fSigPID.size() ){ 
        cout<<"only "<<fSigPID.size()<<" signal particle registered"<<endl;
        return fErrID;
    }
    return fSigTID[i];
}

//______________________________________________________________________
double TGEMSBSDBManager::GetDMag(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).dmag;
}
//______________________________________________________________________
double TGEMSBSDBManager::GetThetaV(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).thetaV;
}
//______________________________________________________________________
double TGEMSBSDBManager::GetD0(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).d0;
}
//______________________________________________________________________
double TGEMSBSDBManager::GetXOffset(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).xoffset;
}
//______________________________________________________________________
double TGEMSBSDBManager::GetDepth(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).depth;
}
//______________________________________________________________________
double TGEMSBSDBManager::GetDX(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).dx;
}
//______________________________________________________________________
double TGEMSBSDBManager::GetDY(int i, int j)
{
  if (!CheckIndex(i, j)) return fErrVal;
  return fPMGeoInfo[i].at(j).dy;
}
//_________________________________________________________________________
double TGEMSBSDBManager::GetStripAngle(int i, int j, int k)
{
  if (!CheckIndex(i, j, k)) return fErrVal;
  if (k == 0) return fPMGeoInfo[i].at(j).stripangle_u;
  else return fPMGeoInfo[i].at(j).stripangle_v;
}
//_________________________________________________________________________
double TGEMSBSDBManager::GetPitch(int i, int j, int k)
{
   if (!CheckIndex(i, j, k)) return fErrVal;
   if (k == 0) return fPMGeoInfo[i].at(j).pitch_u;
   else return fPMGeoInfo[i].at(j).pitch_v;
}



int TGEMSBSDBManager::GetModuleIDFromPos(int iplane, double x, double /*y*/)
{
  if (!CheckIndex(iplane)) return fErrVal;
  
  int module = -1;
  for(size_t k = 0; k<fPMGeoInfo[iplane].size(); k++){
    if(fPMGeoInfo[iplane].at(k).xoffset-fPMGeoInfo[iplane].at(k).dx/2.0<=x && 
       x<=fPMGeoInfo[iplane].at(k).xoffset+fPMGeoInfo[iplane].at(k).dx/2.0)
      {
        module = k;
      }
  }

  return module;
}

//__________________________________________________________________________

double TGEMSBSDBManager::GetPosFromModuleStrip(int iproj, int iplane,
					    int imodule, int istrip)
{
  if (!CheckIndex(iplane, imodule)) return fErrVal;

  double pos = fErrVal;
  if(iproj==0){
    pos = fPMGeoInfo[iplane].at(imodule).pitch_u*istrip
         -fPMGeoInfo[iplane].at(imodule).dx/2.0
         +fPMGeoInfo[iplane].at(imodule).xoffset;
    }
  
  if(iproj==1){
    pos = fPMGeoInfo[iplane].at(imodule).pitch_v*istrip
      -fPMGeoInfo[iplane].at(imodule).dy/2.0;
  }
  
  //cout << " " << pos << endl;
  return pos;
}



void TGEMSBSDBManager::GetPMfromGlobalPlaneNum(uint gplanenum, 
					       int& plane, 
					       int& module)
{
  if(gplanenum>fmIgemtoPlane.size()){
    cout << "TGEMSBSDBManager::GetPMfromGlobalPlaneNum(uint, int&, int&): "
	 << " global plane num = " << gplanenum
	 << " should be < total number of planes " << fmIgemtoPlane.size() << endl;
    return;
  }
  plane = fmIgemtoPlane[gplanenum];
  module = fmIgemtoModule[gplanenum];
}


UInt_t TGEMSBSDBManager::GetGlobalStripPlane(uint lstrip, int plane, int module, int proj)
{
  if (!CheckIndex(plane, module)) return fErrVal;
  uint nstrips = 0;
  for(int i = 0; i<module; i++){
    if(proj==0){
      nstrips+= round(fPMGeoInfo[plane].at(module).dx/fPMGeoInfo[plane].at(module).pitch_u);
    }
    if(proj==1){
      nstrips+= round(fPMGeoInfo[plane].at(module).dy/fPMGeoInfo[plane].at(module).pitch_v);
    }
  }
  return nstrips+lstrip;
}

