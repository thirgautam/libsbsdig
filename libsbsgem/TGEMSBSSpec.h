#ifndef __TGEMSBSSPEC_H
#define __TGEMSBSSPEC_H

#include "THaSpectrometer.h"

#include "types.h"
#include <vector>

class TGEMSBSGEMChamber;

// Class TGEMSBSSpec is more or less a "container" to store the information of all GEM chambers.
// Ultimately, it will also contain information on reconstrcuted tracks, 
// but that needs more work.
// This class inherits from class THaSpectrometer, 
// which grants it all the functions from its class
// (see http://hallaweb.jlab.org/podd/doc/html_v16/ClassIndex.html for more info).

class TGEMSBSSpec : public THaSpectrometer {
    public:
        //Constructor and destructor
	TGEMSBSSpec( const char *name, const char *desc );
        virtual ~TGEMSBSSpec();

	Int_t AddGEM (TGEMSBSGEMChamber* pdet);

	// Useless: the actual job is done by TreeSearch.
	// However, those methods seem to have to be declared, 
	// perhaps for reasons of inheritence from THaSpectrometer
	Int_t CoarseTrack();
	Int_t CoarseReconstruct();
	Int_t Track();
	Int_t Reconstruct();
	
	Int_t TrackCalc() { return 0; }
	
	Int_t FindVertices(TClonesArray &);
	/* void MakePrefix(){ return; } */
	
	//Access to GEM chambers info 
	UInt_t GetNChambers() const { return fChambers.size(); }
	TGEMSBSGEMChamber &GetChamber(Int_t i) const { return *(fChambers.at(i)); }
	
	//Print spectrometer info, with each individual GEM chamber
	void Print(Option_t* opt="") const;

    private:
        std::vector<TGEMSBSGEMChamber*>  fChambers;

    public:
	ClassDef(TGEMSBSSpec,0)
};

#endif//__TGEMSBSSPEC_H
