#include "TGEMSBSSpec.h"
#include "TGEMSBSGEMChamber.h"
#include <iostream>

using namespace std;

TGEMSBSSpec::TGEMSBSSpec(const char* name, const char* desc )
    :THaSpectrometer(name,desc) {

  // We don't need run db (not yet at least)
  fProperties &= ~kNeedsRunDB;
  return;
}

TGEMSBSSpec::~TGEMSBSSpec()
{
  // Destructor: delete all plane objects

  for( vector<TGEMSBSGEMChamber*>::iterator it = fChambers.begin();
       it != fChambers.end(); ++it ) {
    delete *it;
  }
}

Int_t 
TGEMSBSSpec::AddGEM (TGEMSBSGEMChamber* pdet)
{
  // Add a detector to the internal lists of spectrometer detectors.
  // The detector object must be allocated and deleted by the caller.
  // Duplicate detector names are not allowed.

  fChambers.push_back(pdet);
  return 0;
}

Int_t TGEMSBSSpec::CoarseTrack(){
    // Needs work
    return 0;
}

Int_t TGEMSBSSpec::CoarseReconstruct(){
    return 0;
}

Int_t TGEMSBSSpec::Track(){
    return 0;
}

Int_t TGEMSBSSpec::Reconstruct(){
    return 0;
}

Int_t TGEMSBSSpec::FindVertices(TClonesArray &){
    return 0;
}

void
TGEMSBSSpec::Print(Option_t*) const
{
  cout << "Hello, I'm a spectrometer named " << GetName() << endl;
	
  for( vector<TGEMSBSGEMChamber*>::const_iterator it = fChambers.begin();
       it != fChambers.end(); ++it ) {
    (*it)->Print();
  }
}
