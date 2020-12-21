#include "SBSDigPMTSignal.h"
#include "TMath.h"
#include "TFormula.h"
#include "g4sbs_types.h"

using namespace std;

//
// Class SPEModel
//
SPEModel::SPEModel()
{
  fPulseHisto = new TH1D("fPulseHisto", "", 1000, -50, 50);
}

SPEModel::SPEModel(UShort_t uniqueid, double sigma, 
		   double t0, double tmin, double tmax)
{
  //TF1 fFunc("fFunc", "landaun", tmin, tmax);//garbage (sorry)
  // power law x exp decay...
  TF1 fFunc("fFunc", 
	    "TMath::Max(0., [0]*((x-[1]+[2]*0.4)/([2]*[2]*0.16))*TMath::Exp(-(x-[1]+[2]*0.4)/([2]*0.4)) )", 
	    tmin, tmax); 
  
  fFunc.SetParameters(1., t0, sigma);
  const int NbinsTotal = int(tmax-tmin)*20;// 20 bins/ns should do... since we will extrapolate after...
  //let's try 20...
  fPulseHisto = new TH1D(Form("fPulseHisto_%d", uniqueid), "", NbinsTotal, tmin, tmax);
  double t_i;//, t_j;
  double ps_i;//, g_j;
  for(int i = 1; i<=NbinsTotal; i++){
    t_i = fPulseHisto->GetBinCenter(i+1);
    ps_i = fFunc.Eval(t_i);
    fPulseHisto->Fill(t_i, ps_i);
    //cout << fPulseHisto->GetBinContent(i) << " ";
  }
  //cout << endl;
}

SPEModel::~SPEModel()
{
  fPulseHisto->Delete();
}

bool SPEModel::PulseOverThr(double charge, double thr)
{
  if(fPulseHisto->GetMaximum()<thr/charge){
    return false;
  }else{
    return true;
  }
};
 
bool SPEModel::FindLeadTrailTime(double charge, double thr, double &t_lead, double &t_trail)
{
  //cout << charge << " " << thr << endl;  
  if(!PulseOverThr(charge, thr)){
    t_lead = 1.0e38;
    t_trail = 1.0e38;
    return false;
  }else{
    double xmax = fPulseHisto->GetBinCenter(fPulseHisto->GetMaximumBin());
    //cout << fPulseHisto->GetBinContent(1)<< " " << thr/charge << " " << fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX())<< endl; 
    if(fPulseHisto->GetBinContent(1)<thr/charge){
      t_lead = GetHistoX(thr/charge, fPulseHisto->GetBinLowEdge(1), xmax);
    }else{
      t_lead = 1.0e38;
      return false;
    }
    if(fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX())<thr/charge){
      t_trail = GetHistoX(thr/charge, xmax, fPulseHisto->GetBinLowEdge(fPulseHisto->GetNbinsX()+1));
    }else{
      t_trail = 1.0e38;
      return false;
    }
    //why does this function seem to oblitarate fMCHit containers when t_trail is calculated to be 1.e38... GetHistoX?
    if(t_trail>1.e30) cout << "thr/charge" << thr/charge << " >? " <<  fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX()) << " or " << fPulseHisto->GetBinContent(fPulseHisto->GetNbinsX()+1) << endl;
    return true;
  }
}

double SPEModel::GetHistoX(double y, double x1, double x2)
{
  double splineslope;
  for(int k = fPulseHisto->FindBin(x1); k<=fPulseHisto->FindBin(x2); k++){
    if(  ( (fPulseHisto->GetBinContent(k+1)-y)*(fPulseHisto->GetBinContent(k)-y) )<0 ){
      // threshold crossed if diff(TH1::GetBinContent-y) changes sign
      splineslope = (fPulseHisto->GetBinContent(k+1)-fPulseHisto->GetBinContent(k))/(fPulseHisto->GetBinCenter(k+1)-fPulseHisto->GetBinCenter(k));
      
      return fPulseHisto->GetBinCenter(k)+(y-fPulseHisto->GetBinContent(k))/splineslope;
    }
  }
  return 1.0e38;
}

//
// Class PMTSignal
//
PMTSignal::PMTSignal()
  : fSumEdep(0), fNpe(0), fNpeChargeConv(1.0), fADC(0), fEventTime(0), fNADCSamps(0), fNSamps(0)
{ 
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
  
  f1 = 0;
  //R = TRndmManager::GetInstance();
}

PMTSignal::PMTSignal(double npechargeconv)
  : fSumEdep(0), fNpe(0), fNpeChargeConv(npechargeconv), fADC(0), fEventTime(0), fNADCSamps(0), fNSamps(0)
{ 
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
  
  f1 = 0;
}

void PMTSignal::Fill(SPEModel *model, int npe, double thr, double evttime, int signal)
{
  if(signal==0)fEventTime = evttime;
  fNpe+= npe;
  //if(model->PulseOverThr(fCharge, thr))fNpe//fADC = model->GetCharge()*model->GetADCconversion();
  
  //determine lead and trail times
  double t_lead, t_trail;
  // find the lead and trail time for *this* pulse, not the total pulse
  bool goodtime = model->FindLeadTrailTime(npe*fNpeChargeConv, thr, t_lead, t_trail);
  
  //if(goodtime)cout << evttime << " " << t_lead << " " << t_trail << endl;
  
  t_lead+=evttime;
  t_trail+=evttime;
  
  if(goodtime){
    //Filter here the lead and trail times
    if(fLeadTimes.size()>0){
      // Check if the lead and trail times are inside an existing lead time- trail time pair
      // we assume here that fLeadTimes and fTrailTimes are same size 
      // *if we do things correctly, that should be the case*
      // we shall keep lead-trail times pair in timing order
      // we neglect pulse overlaps ftm.
      if(fLeadTimes.size()!=fTrailTimes.size()){
	cout << " B - Warning: size of lead times container: " << fLeadTimes.size() 
	     << " != size of trail times container: " << fTrailTimes.size() << endl;
      }
      
      for(size_t i = 0; i<fLeadTimes.size(); i++){
	// possibility of the current pair straddling with others.... :/
	// treat those separately to simplify...
	// tL < tT_i-1 < tL_i < tT
	if(i>0){
	  if(t_lead < fTrailTimes.at(i-1) && fLeadTimes.at(i) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i-1 => tL *replaces* tL_i-1
	    if(t_lead < fLeadTimes.at(i-1)){
	      fLeadTimes.erase(fLeadTimes.begin()+i-1);
	      fLeadTimes.insert(fLeadTimes.begin()+i-1, t_lead);
	    }
	    // tT_i < tT => tT *replaces* tT_i
	    if(fTrailTimes.at(i) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i);
	      fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i);
	    fTrailTimes.erase(fTrailTimes.begin()+i-1);
	    break;
	  }
	}//end if(i>0)
	// tL < tT_i < tL_i+1 < tT
	if(i<fLeadTimes.size()-1){
	  if(t_lead < fTrailTimes.at(i) && fLeadTimes.at(i+1) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i => tL *replaces* tL_i
	    if(t_lead < fLeadTimes.at(i)){
	      fLeadTimes.erase(fLeadTimes.begin()+i);
	      fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	    }
	    // tT_i+1 < tT => tT *replaces* tT_i+1
	    if(fTrailTimes.at(i+1) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i+1);
	      fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i+1);
	    fTrailTimes.erase(fTrailTimes.begin()+i);
	    break;
	  }
	}//end if(i<fLeadTimes.size()-1)
	
	// if not, 6 cases to consider:
	// tL < tT < tL_i < tT_i => both inserted *before* existing pair 
	if(t_trail < fLeadTimes.at(i)){
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail); 
	  break;
	}
	// tL < tL_i < tT < tT_i => tL *replaces* tL_i
	if(t_lead < fLeadTimes.at(i) && fLeadTimes.at(i) < t_trail && t_trail < fTrailTimes.at(i)){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  break;
	}
	// tL_i < tL < tT < tT_i => tL *replaces* tL_i AND tT *replaces* tT_i
	if(t_lead < fLeadTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_trail < fTrailTimes.at(i)){
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_lead < fTrailTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	}
	if(fTrailTimes.at(i) < t_lead){
	  fLeadTimes.insert(fLeadTimes.begin()+i+1, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	  break;
	}
	if(fLeadTimes.size()!=fTrailTimes.size()){
	  cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	       << " != size of trail times container: " << fTrailTimes.size() << endl;
	}
	if(i>=fLeadTimes.size())cout << "Warning: i = " << i << " >= size of containers = " << fLeadTimes.size() << endl;
      }
    }else{
      //of course, if initial size was 0, just psuh it back
      //hopefully it will be the case most of the time
      fLeadTimes.push_back(t_lead);
      fTrailTimes.push_back(t_trail);
    }
    if(fLeadTimes.size()!=fTrailTimes.size()){
      cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	   << " != size of trail times container: " << fTrailTimes.size() << endl;
    }
  }//end if(t_lead && t_trail<30)
}

void PMTSignal::Fill(int npe, double thr, double evttime, double sigmatime, int signal)
{
  if(signal==0)fEventTime = evttime;
  fNpe+= npe;
  //if(model->PulseOverThr(fCharge, thr))fNpe//fADC = model->GetCharge()*model->GetADCconversion();
  
  if(evttime>=fTmin+fNSamps*fSampSize)return;
  
  f1->SetParameters(npe*fNpeChargeConv, evttime, sigmatime);
  //determine lead and trail times
  double t_lead, t_trail;
  bool goodtime = false;//model->FindLeadTrailTime(npe*fNpeChargeConv, thr, t_lead, t_trail);
  if(fNSamps){
    fSamples[0]+= f1->Eval(fTmin+(0.5)*fSampSize);//*fSampSize;
    //Evaluate this function might be a bit of a time drain!
    for(int i = 1; i<fNSamps; i++){
      fSamples[i]+= f1->Eval(fTmin+(i+0.5)*fSampSize);//*fSampSize;
      //if(i>0){
      if(fSamples[i-1]<=thr && thr<fSamples[i]){
	t_lead = fTmin+(i-0.5)*fSampSize+fSampSize*(thr-fSamples[i-1])/(fSamples[i]-fSamples[i-1]);
	goodtime = true;
      }
      if(fSamples[i-1]>=thr && thr>fSamples[i]){
	t_trail = fTmin+(i-0.5)*fSampSize+fSampSize*(thr-fSamples[i-1])/(fSamples[i]-fSamples[i-1]);
	goodtime = true;
      }
    }
  }
    
  // find the lead and trail time for *this* pulse, not the total pulse
  //t_lead+=evttime;
  //t_trail+=evttime;
  //if(goodtime)cout << evttime << " " << t_lead << " " << t_trail << endl;
  
  if(goodtime){
    //Filter here the lead and trail times
    if(fLeadTimes.size()>0){
      // Check if the lead and trail times are inside an existing lead time- trail time pair
      // we assume here that fLeadTimes and fTrailTimes are same size 
      // *if we do things correctly, that should be the case*
      // we shall keep lead-trail times pair in timing order
      // we neglect pulse overlaps ftm.
      if(fLeadTimes.size()!=fTrailTimes.size()){
	cout << " B - Warning: size of lead times container: " << fLeadTimes.size() 
	     << " != size of trail times container: " << fTrailTimes.size() << endl;
      }
      
      for(size_t i = 0; i<fLeadTimes.size(); i++){
	// possibility of the current pair straddling with others.... :/
	// treat those separately to simplify...
	// tL < tT_i-1 < tL_i < tT
	if(i>0){
	  if(t_lead < fTrailTimes.at(i-1) && fLeadTimes.at(i) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i-1 => tL *replaces* tL_i-1
	    if(t_lead < fLeadTimes.at(i-1)){
	      fLeadTimes.erase(fLeadTimes.begin()+i-1);
	      fLeadTimes.insert(fLeadTimes.begin()+i-1, t_lead);
	    }
	    // tT_i < tT => tT *replaces* tT_i
	    if(fTrailTimes.at(i) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i);
	      fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i);
	    fTrailTimes.erase(fTrailTimes.begin()+i-1);
	    break;
	  }
	}//end if(i>0)
	// tL < tT_i < tL_i+1 < tT
	if(i<fLeadTimes.size()-1){
	  if(t_lead < fTrailTimes.at(i) && fLeadTimes.at(i+1) < t_trail){
	    //do necessary substitutions first, then erase
	    // tL < tL_i => tL *replaces* tL_i
	    if(t_lead < fLeadTimes.at(i)){
	      fLeadTimes.erase(fLeadTimes.begin()+i);
	      fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	    }
	    // tT_i+1 < tT => tT *replaces* tT_i+1
	    if(fTrailTimes.at(i+1) < t_trail){
	      fTrailTimes.erase(fTrailTimes.begin()+i+1);
	      fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	    }
	    fLeadTimes.erase(fLeadTimes.begin()+i+1);
	    fTrailTimes.erase(fTrailTimes.begin()+i);
	    break;
	  }
	}//end if(i<fLeadTimes.size()-1)
	
	// if not, 6 cases to consider:
	// tL < tT < tL_i < tT_i => both inserted *before* existing pair 
	if(t_trail < fLeadTimes.at(i)){
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail); 
	  break;
	}
	// tL < tL_i < tT < tT_i => tL *replaces* tL_i
	if(t_lead < fLeadTimes.at(i) && fLeadTimes.at(i) < t_trail && t_trail < fTrailTimes.at(i)){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  break;
	}
	// tL_i < tL < tT < tT_i => tL *replaces* tL_i AND tT *replaces* tT_i
	if(t_lead < fLeadTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fLeadTimes.erase(fLeadTimes.begin()+i);
	  fLeadTimes.insert(fLeadTimes.begin()+i, t_lead);
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_trail < fTrailTimes.at(i)){
	  break;
	}
	if(fLeadTimes.at(i) < t_lead && t_lead < fTrailTimes.at(i) && fTrailTimes.at(i) < t_trail){
	  fTrailTimes.erase(fTrailTimes.begin()+i);
	  fTrailTimes.insert(fTrailTimes.begin()+i, t_trail);
	}
	if(fTrailTimes.at(i) < t_lead){
	  fLeadTimes.insert(fLeadTimes.begin()+i+1, t_lead);
	  fTrailTimes.insert(fTrailTimes.begin()+i+1, t_trail);
	  break;
	}
	if(fLeadTimes.size()!=fTrailTimes.size()){
	  cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	       << " != size of trail times container: " << fTrailTimes.size() << endl;
	}
	if(i>=fLeadTimes.size())cout << "Warning: i = " << i << " >= size of containers = " << fLeadTimes.size() << endl;
      }
    }else{
      //of course, if initial size was 0, just psuh it back
      //hopefully it will be the case most of the time
      fLeadTimes.push_back(t_lead);
      fTrailTimes.push_back(t_trail);
    }
    if(fLeadTimes.size()!=fTrailTimes.size()){
      cout << " A - Warning: size of lead times container: " << fLeadTimes.size() 
	   << " != size of trail times container: " << fTrailTimes.size() << endl;
    }
  }//end if(t_lead && t_trail<30)
}



void PMTSignal::Digitize(int chan, int detid, g4sbs_tree* T, //gmn_tree* T, 
			 TRandom3* R, double ped, double ped_noise, double ADCconv, double ADCbits, double TDCconv, double TDCbits)
{
  if(fNpe<=0){
    //fADC = R->Gaus(ped, ped_noise);
    return;
  }

  fADC = TMath::Nint(Charge()*1.0e15/ADCconv+R->Gaus(ped, ped_noise));

  if( fADC>UInt_t(TMath::Nint( TMath::Power(2, ADCbits) )) ){
    fADC = TMath::Nint( TMath::Power(2, ADCbits) );
  }
  
  Int_t tdc_value;
  if(fLeadTimes.size()){
    for(size_t i = 0; i<fLeadTimes.size(); i++){
      //cout << "detid " << detid << " fLeadTimes.at(" << i << ") " << fLeadTimes.at(i) << " fTrailTimes.at(" << i << ") " << fTrailTimes.at(i) << endl;
      // trim "all" bits that are above the number of TDC bits - a couple to speed it up
      // (since TDC have a revolving clock, as far as I understand)
      // let's use an arbitrary reference time offset of 1us before the trigger
      tdc_value = TMath::Nint((fLeadTimes.at(i))/TDCconv)+1000;
      for(int j = 30; j>=TDCbits; j--){
	tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (j) );
      }
      tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (31) );
      //fTDCs.insert(fTDCs.begin()+0, TMath::Nint(fLeadTimes.at(0)*diginfo.TDCConversion()));//bug!!!!
      fTDCs.push_back(tdc_value);//they're already sorted in order, presumably
      // also mark the traling time with setting bin 31 to 1 // need to reconvert then
      tdc_value = TMath::Nint((fTrailTimes.at(i))/TDCconv)+1000;
      for(int j = 30; j>=TDCbits; j--){
	tdc_value ^= ( -0 ^ tdc_value) & ( 1 << (j) );
      }
      tdc_value ^= ( -1 ^ tdc_value) & ( 1 << (31) );
      fTDCs.push_back(tdc_value);
    }
  }
  //cout << "detid " << detid << " TDC size " << fTDCs.size() << endl;
  
  if(fNSamps){
    fADC = 0;
    Int_t Nconv = fNSamps/fNADCSamps;
    for(int i = 0; i<fNADCSamps; i++){
      for(int j = 0; j<Nconv; j++)fADCSamples[i]+=fSamples[i*Nconv+j]*fSampSize;//renormalize the sample for the integration;
      fADCSamples[i]*=ADCconv;
      fADCSamples[i]+=R->Gaus(ped, ped_noise);
      if(fADCSamples[i]>pow(2, ADCbits))
      fADC+=fADCSamples[i];
    }
  }
  //fSumEdep*=1.0e9;// store in eV.
  
  //Fill in directly (hoping it takes less time...)
  //switch(detid){
  //case(BBPS_UNIQUE_DETID):
  //}
  if(detid==BBPS_UNIQUE_DETID){
    // T->Earm_BBPS_dighit_nchan++;
    // T->Earm_BBPS_dighit_chan->push_back(chan);
    // T->Earm_BBPS_dighit_adc->push_back(fADC);
    T->Earm_BBPS_Dig.nchan++;
    T->Earm_BBPS_Dig.chan->push_back(chan);
    T->Earm_BBPS_Dig.adc->push_back(fADC);
  }
  
  if(detid==BBSH_UNIQUE_DETID){
    // T->Earm_BBSH_dighit_nchan++;
    // T->Earm_BBSH_dighit_chan->push_back(chan);
    // T->Earm_BBSH_dighit_adc->push_back(fADC);
    T->Earm_BBSH_Dig.nchan++;
    T->Earm_BBSH_Dig.chan->push_back(chan);
    T->Earm_BBSH_Dig.adc->push_back(fADC);
  }
  
  if(detid==HODO_UNIQUE_DETID){
    // T->Earm_BBHodo_dighit_nchan++;
    // T->Earm_BBHodo_dighit_chan->push_back(chan);
    // T->Earm_BBHodo_dighit_adc->push_back(fADC);
    T->Earm_BBHodo_Dig.nchan++;
    T->Earm_BBHodo_Dig.chan->push_back(chan);
    T->Earm_BBHodo_Dig.adc->push_back(fADC);
    if(fTDCs.size()==2){
      for(int j = 0;j<fTDCs.size(); j++){
	if(fTDCs[j] & ( 1 << (31) )){
	  fTDCs[j] ^= ( -0 ^ fTDCs[j] ) & ( 1 << (31) );
	  //T->Earm_BBHodo_dighit_tdc_t->push_back(fTDCs[j]-1000);
	  T->Earm_BBHodo_Dig.tdc_t->push_back(fTDCs[j]-1000);
	}else{
	  //T->Earm_BBHodo_dighit_tdc_l->push_back(fTDCs[j]-1000);
	  T->Earm_BBHodo_Dig.tdc_l->push_back(fTDCs[j]-1000);
	}
	/*
	if(j>3 && j%2==0){
	  // T->Earm_BBHodo_dighit_nchan++;
	  // T->Earm_BBHodo_dighit_chan->push_back(chan);
	  // T->Earm_BBHodo_dighit_adc->push_back(-1000000);
	  T->Earm_BBHodo_Dig.nchan++;
	  T->Earm_BBHodo_Dig.chan->push_back(chan);
	  T->Earm_BBHodo_Dig.adc->push_back(-1000000);
	}
	*/
      }
    }else{
      // T->Earm_BBHodo_dighit_tdc_l->push_back(-1000000);
      // T->Earm_BBHodo_dighit_tdc_t->push_back(-1000000);
      T->Earm_BBHodo_Dig.tdc_l->push_back(-1000000);
      T->Earm_BBHodo_Dig.tdc_t->push_back(-1000000);
    }
  }
  
  if(detid==GRINCH_UNIQUE_DETID){
    // T->Earm_GRINCH_dighit_nchan++;
    // T->Earm_GRINCH_dighit_chan->push_back(chan);
    // T->Earm_GRINCH_dighit_adc->push_back(fADC);
    T->Earm_GRINCH_Dig.nchan++;
    T->Earm_GRINCH_Dig.chan->push_back(chan);
    T->Earm_GRINCH_Dig.adc->push_back(fADC);
    if(fTDCs.size()==2){
      for(int j = 0;j<fTDCs.size(); j++){
	if(fTDCs[j] & ( 1 << (31) )){
	  fTDCs[j] ^= ( -0 ^ fTDCs[j] ) & ( 1 << (31) );
	  //T->Earm_GRINCH_dighit_tdc_t->push_back(fTDCs[j]-1000);
	  T->Earm_GRINCH_Dig.tdc_t->push_back(fTDCs[j]-1000);
	}else{
	  //T->Earm_GRINCH_Dig.tdc_l->push_back(fTDCs[j]-1000);
	  T->Earm_GRINCH_Dig.tdc_l->push_back(fTDCs[j]-1000);
	}
	/*
	if(j>3 && j%2==0){
	  // T->Earm_GRINCH_dighit_nchan++;
	  // T->Earm_GRINCH_dighit_chan->push_back(chan);
	  // T->Earm_GRINCH_dighit_adc->push_back(-1000000);
	  T->Earm_GRINCH_Dig.nchan++;
	  T->Earm_GRINCH_Dig.chan->push_back(chan);
	  T->Earm_GRINCH_Dig.adc->push_back(-1000000);
	}
	*/
      }
    }else{
      // T->Earm_GRINCH_dighit_tdc_l->push_back(-1000000);
      // T->Earm_GRINCH_dighit_tdc_t->push_back(-1000000);
      T->Earm_GRINCH_Dig.tdc_l->push_back(-1000000);
      T->Earm_GRINCH_Dig.tdc_t->push_back(-1000000);
    }
  }
  
  if(detid==HCAL_UNIQUE_DETID){
    //T->Harm_HCal_dighit_nchan++;
    //T->Harm_HCal_dighit_chan->push_back(chan);
    T->Harm_HCal_Dig.nchan++;
    T->Harm_HCal_Dig.chan->push_back(chan);
    // T->Harm_HCal_dighit_adc_0->push_back(fADCSamples[0]);
    // T->Harm_HCal_dighit_adc_1->push_back(fADCSamples[1]);
    // T->Harm_HCal_dighit_adc_2->push_back(fADCSamples[2]);
    // T->Harm_HCal_dighit_adc_3->push_back(fADCSamples[3]);
    // T->Harm_HCal_dighit_adc_4->push_back(fADCSamples[4]);
    // T->Harm_HCal_dighit_adc_5->push_back(fADCSamples[5]);
    // T->Harm_HCal_dighit_adc_6->push_back(fADCSamples[6]);
    // T->Harm_HCal_dighit_adc_7->push_back(fADCSamples[7]);
    // T->Harm_HCal_dighit_adc_8->push_back(fADCSamples[8]);
    // T->Harm_HCal_dighit_adc_9->push_back(fADCSamples[9]);
    // T->Harm_HCal_dighit_adc_10->push_back(fADCSamples[10]);
    // T->Harm_HCal_dighit_adc_11->push_back(fADCSamples[11]);
    // T->Harm_HCal_dighit_adc_12->push_back(fADCSamples[12]);
    // T->Harm_HCal_dighit_adc_13->push_back(fADCSamples[13]);
    // T->Harm_HCal_dighit_adc_14->push_back(fADCSamples[14]);
    // T->Harm_HCal_dighit_adc_15->push_back(fADCSamples[15]);
    // T->Harm_HCal_dighit_adc_16->push_back(fADCSamples[16]);
    // T->Harm_HCal_dighit_adc_17->push_back(fADCSamples[17]);
    // T->Harm_HCal_dighit_adc_18->push_back(fADCSamples[18]);
    // T->Harm_HCal_dighit_adc_19->push_back(fADCSamples[19]);
    T->Harm_HCal_Dig.adc_0->push_back(fADCSamples[0]);
    T->Harm_HCal_Dig.adc_1->push_back(fADCSamples[1]);
    T->Harm_HCal_Dig.adc_2->push_back(fADCSamples[2]);
    T->Harm_HCal_Dig.adc_3->push_back(fADCSamples[3]);
    T->Harm_HCal_Dig.adc_4->push_back(fADCSamples[4]);
    T->Harm_HCal_Dig.adc_5->push_back(fADCSamples[5]);
    T->Harm_HCal_Dig.adc_6->push_back(fADCSamples[6]);
    T->Harm_HCal_Dig.adc_7->push_back(fADCSamples[7]);
    T->Harm_HCal_Dig.adc_8->push_back(fADCSamples[8]);
    T->Harm_HCal_Dig.adc_9->push_back(fADCSamples[9]);
    T->Harm_HCal_Dig.adc_10->push_back(fADCSamples[10]);
    T->Harm_HCal_Dig.adc_11->push_back(fADCSamples[11]);
    T->Harm_HCal_Dig.adc_12->push_back(fADCSamples[12]);
    T->Harm_HCal_Dig.adc_13->push_back(fADCSamples[13]);
    T->Harm_HCal_Dig.adc_14->push_back(fADCSamples[14]);
    T->Harm_HCal_Dig.adc_15->push_back(fADCSamples[15]);
    T->Harm_HCal_Dig.adc_16->push_back(fADCSamples[16]);
    T->Harm_HCal_Dig.adc_17->push_back(fADCSamples[17]);
    T->Harm_HCal_Dig.adc_18->push_back(fADCSamples[18]);
    T->Harm_HCal_Dig.adc_19->push_back(fADCSamples[19]);
    if(fTDCs.size()){
      //T->Harm_HCal_dighit_tdc->push_back(fTDCs[0]);
      T->Harm_HCal_Dig.tdc->push_back(fTDCs[0]-1000);
    }else{
      //T->Harm_HCal_dighit_tdc->push_back(-1000000);
      T->Harm_HCal_Dig.tdc->push_back(-1000000);
    }
  }
  
  // ** How to add a new subsystem **
  // fill the new detector output here
  //genrp Hoda Harm_PRPolScintBeamSide
  if(detid==PRPOLBS_SCINT_UNIQUE_DETID){
    T->Harm_PRPolScintBeamSide_Dig.nchan++;
    T->Harm_PRPolScintBeamSide_Dig.chan->push_back(chan);
    T->Harm_PRPolScintBeamSide_Dig.adc->push_back(fADC);

    if(fTDCs.size()==2){
      for(int j = 0;j<fTDCs.size(); j++){
	if(fTDCs[j] & ( 1 << (31) )){
	  fTDCs[j] ^= ( -0 ^ fTDCs[j] ) & ( 1 << (31) );
	  T->Harm_PRPolScintBeamSide_Dig.tdc_t->push_back(fTDCs[j]-1000);
	}else{
	  T->Harm_PRPolScintBeamSide_Dig.tdc_l->push_back(fTDCs[j]-1000);
	}
      }
    }else{
      T->Harm_PRPolScintBeamSide_Dig.tdc_l->push_back(-1000000);
      T->Harm_PRPolScintBeamSide_Dig.tdc_t->push_back(-1000000);
    }

  }
  //genrp Hoda Harm_PRPolScintFarSide
  if(detid==PRPOLFS_SCINT_UNIQUE_DETID){
    T->Harm_PRPolScintFarSide_Dig.nchan++;
    T->Harm_PRPolScintFarSide_Dig.chan->push_back(chan);
    T->Harm_PRPolScintFarSide_Dig.adc->push_back(fADC);
    if(fTDCs.size()==2){
      for(int j = 0;j<fTDCs.size(); j++){
	if(fTDCs[j] & ( 1 << (31) )){
	  fTDCs[j] ^= ( -0 ^ fTDCs[j] ) & ( 1 << (31) );
	  T->Harm_PRPolScintFarSide_Dig.tdc_t->push_back(fTDCs[j]-1000);
	}else{
	  T->Harm_PRPolScintFarSide_Dig.tdc_l->push_back(fTDCs[j]-1000);
	}
      }
    }else{
      T->Harm_PRPolScintFarSide_Dig.tdc_l->push_back(-1000000);
      T->Harm_PRPolScintFarSide_Dig.tdc_t->push_back(-1000000);
    }
  }
  //genrp ActiveAna
  if(detid==ACTIVEANA_UNIQUE_DETID){
    T->Harm_ActAn_Dig.nchan++;
    T->Harm_ActAn_Dig.chan->push_back(chan);
    T->Harm_ActAn_Dig.adc->push_back(fADC);
    if(fTDCs.size()==2){
      for(int j = 0;j<fTDCs.size(); j++){
	if(fTDCs[j] & ( 1 << (31) )){
	  fTDCs[j] ^= ( -0 ^ fTDCs[j] ) & ( 1 << (31) );
	  T->Harm_ActAn_Dig.tdc_t->push_back(fTDCs[j]-1000);
	}else{
	  T->Harm_ActAn_Dig.tdc_l->push_back(fTDCs[j]-1000);
	}
      }
    }else{
      T->Harm_ActAn_Dig.tdc_l->push_back(-1000000);
      T->Harm_ActAn_Dig.tdc_t->push_back(-1000000);
    }
  }
  
}//

void PMTSignal::SetSamples(double tmin, double tmax, double sampsize)
{
  fTmin = tmin;
  fADCSampSize = sampsize;
  //fSampSize = sampsize/10;//the bin size is too large for the tdc size 
  fSampSize = sampsize/64;
  fNADCSamps = round((tmax-tmin)/fADCSampSize);
  fNSamps = round((tmax-tmin)/fSampSize);
  fADCSamples = new double[fNADCSamps];
  fSamples = new double[fNSamps];
  
  //cout << fNADCSamps << " " << fNSamps << " " << fADCSampsSize << " " << fSampSize << endl;
  
  memset(fSamples, 0, fNSamps*sizeof(double));
  memset(fADCSamples, 0, fNADCSamps*sizeof(double));
  //f1 = new TF1("f1", "landaun", tmin, tmax);
  f1 = new TF1("fFunc", 
	       "TMath::Max(0., [0]*((x-[1]+[2]*0.4)/([2]*[2]*0.16))*TMath::Exp(-(x-[1]+[2]*0.4)/([2]*0.4)) )", 
	       tmin, tmax); 
}

void PMTSignal::Clear(bool dosamples)
{
  //cout << " PMTSignal::Clear() " << endl;
  
  fSumEdep = 0;
  fNpe = 0;
  fADC = 0;
  
  fEventTime = 0;
  fLeadTimes.clear();
  fTrailTimes.clear();
  fTDCs.clear();
  
  if(dosamples){
    memset(fSamples, 0, fNSamps*sizeof(double));
    memset(fADCSamples, 0, fNADCSamps*sizeof(double));
  }
}



//ClassImp(SPEModel)

