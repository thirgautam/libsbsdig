#include "g4sbs_types.h"
#include "SBSDigAuxi.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

bool UnfoldData(g4sbs_tree* T, double theta_sbs, double d_hcal, TRandom3* R, 
		std::vector<SBSDigPMTDet*> pmtdets, 
		std::vector<int> detmap,
		std::vector<SBSDigGEMDet*> gemdets, 
		std::vector<int> gemmap,
		//std::map<int, SBSDigPMTDet*> pmtdets, 
		//std::map<int, SBSDigGEMDet*> gemdets, 
		double tzero,
		int signal)
{
  bool has_data = false;
  
  int Npe;
  double t;
  
  double x_ref = -d_hcal*sin(theta_sbs);
  double z_ref = d_hcal*cos(theta_sbs);
  
  double z_hit, Npe_Edep_ratio, sigma_tgen;

  //double z_det, pz, E, 
  double beta, sin2thetaC;
  
  int chan;
  int mod;
  
  int idet = 0;
  if(!detmap.empty()){
    //ordering by increasing unique det ID
    while(idet<(int)detmap.size()){
      if(detmap[idet]!=HCAL_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(detmap[idet]!=HCAL_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;  
    if(idet>=0){// && T->Harm_HCalScint.sumedep) {
      for(size_t k = 0; k < T->Harm_HCalScint.sumedep->size(); k++) {
	chan = T->Harm_HCalScint.cell->at(k);
	//T->Harm_HCalScint_hit_sumedep->at(k);
      
	z_hit = -(T->Harm_HCalScint.xhitg->at(k)-x_ref)*sin(theta_sbs)+(T->Harm_HCalScint.zhitg->at(k)-z_ref)*cos(theta_sbs);
      
	// Evaluation of number of photoelectrons from energy deposit documented at:
	// https://sbs.jlab.org/DocDB/0000/000043/002/Harm_HCal_Digi_EdepOnly_2.pdf
	// TODO: put that stuff in DB...
	Npe_Edep_ratio = 5.242+11.39*z_hit+10.41*pow(z_hit, 2);
	Npe = R->Poisson(Npe_Edep_ratio*T->Harm_HCalScint.sumedep->at(k)*1.0e3);
	t = tzero+R->Gaus(T->Harm_HCalScint.tavg->at(k)+10.11, 1.912)-pmtdets[idet]->fTrigOffset;
      
	sigma_tgen = 0.4244+11380/pow(Npe+153.4, 2);
	//Generate here,...
	//cout << " HCal : t = " << t << ", t_zero = " << tzero << ", t_avg = " << T->Harm_HCalScint.tavg->at(k) << ", -t_offset = " << -pmtdets[idet]->fTrigOffset << endl;
	pmtdets[idet]->PMTmap[chan].Fill(Npe, pmtdets[idet]->fThreshold, t, sigma_tgen, signal);
      }
      has_data = true;
    }
  
    while(idet<(int)detmap.size()){
      if(idet<0)idet++;
      if(detmap[idet]!=BBPS_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(detmap[idet]!=BBPS_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;  
    // Process BB PS data
    if(idet>=0){// && T->Earm_BBPSTF1.nhits){
      for(int i = 0; i<T->Earm_BBPSTF1.nhits; i++){
	// Evaluation of number of photoelectrons and time from energy deposit documented at:
	if(T->Earm_BBPSTF1.sumedep->at(i)<1.e-4)continue;
	//check probability to generate p.e. yield
	bool genpeyield = true;
	if(T->Earm_BBPSTF1.sumedep->at(i)<1.e-2)
	  genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*T->Earm_BBPSTF1.sumedep->at(i));
	//if we're go, let's generate
	if(genpeyield){
	  beta = sqrt( pow(m_e+T->Earm_BBPSTF1.sumedep->at(i), 2)-m_e*m_e )/(m_e + T->Earm_BBPSTF1.sumedep->at(i));
	  sin2thetaC = TMath::Max(1.-1./pow(n_lg*beta, 2), 0.);
	  //1500. Used to be 454.: just wrong
	  Npe = R->Poisson(300.0*T->Earm_BBPSTF1.sumedep->at(i)*sin2thetaC/(1.-1./(n_lg*n_lg)) );
	  t = tzero+T->Earm_BBPSTF1.tavg->at(i)+R->Gaus(3.2-5.805*T->Earm_BBPSTF1.zhit->at(i)-17.77*pow(T->Earm_BBPSTF1.zhit->at(i), 2), 0.5)-pmtdets[idet]->fTrigOffset;
	  chan = T->Earm_BBPSTF1.cell->at(i);
	  //T->Earm_BBPSTF1_hit_sumedep->at(i);
	
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, signal);
	}
      }
      has_data = true;
    }
  
    while(idet<(int)detmap.size()){
      if(idet<0)idet++;
      if(detmap[idet]!=BBSH_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(detmap[idet]!=BBSH_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;
    if(idet>=0){// && T->Earm_BBSHTF1.nhits){
      for(int i = 0; i<T->Earm_BBSHTF1.nhits; i++){
	// Evaluation of number of photoelectrons and time from energy deposit documented at:
	// 
	if(T->Earm_BBSHTF1.sumedep->at(i)<1.e-4)continue;
	//check probability to generate p.e. yield
	bool genpeyield = true;
	if(T->Earm_BBSHTF1.sumedep->at(i)<1.e-2)genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*T->Earm_BBSHTF1.sumedep->at(i));
	//if we're go, let's generate
	if(genpeyield){
	  beta = sqrt( pow(m_e+T->Earm_BBSHTF1.sumedep->at(i), 2)-m_e*m_e )/(m_e + T->Earm_BBSHTF1.sumedep->at(i));
	  sin2thetaC = TMath::Max(1.-1./pow(n_lg*beta, 2), 0.);
	  //1800. Used to be 932.: just wrong
	  Npe = R->Poisson(360.0*T->Earm_BBSHTF1.sumedep->at(i)*sin2thetaC/(1.-1./(n_lg*n_lg)) );
	  t = tzero+T->Earm_BBSHTF1.tavg->at(i)+R->Gaus(2.216-8.601*T->Earm_BBSHTF1.zhit->at(i)-7.469*pow(T->Earm_BBSHTF1.zhit->at(i), 2), 0.8)-pmtdets[idet]->fTrigOffset;
	  chan = T->Earm_BBSHTF1.cell->at(i);
	  //T->Earm_BBSHTF1_hit_sumedep->at(i);
		
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, signal);
	}
      }
      has_data = true;
    }

    while(idet<(int)detmap.size()){
      if(idet<0)idet++;
      if(detmap[idet]!=GRINCH_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(detmap[idet]!=GRINCH_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;
    // Process GRINCH data
    if(idet>=0){// && T->Earm_GRINCH.nhits){
      for(int i = 0; i<T->Earm_GRINCH.nhits; i++){
	chan = int(T->Earm_GRINCH.PMT->at(i)/5)-1;
	t = tzero+T->Earm_GRINCH.Time_avg->at(i)+pmtdets[idet]->fTrigOffset;
	Npe = T->Earm_GRINCH.NumPhotoelectrons->at(i);
      
	//if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
      }
      has_data = true;
    }
  
    while(idet<(int)detmap.size()){
      if(idet<0)idet++;
      if(detmap[idet]!=HODO_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(detmap[idet]!=HODO_UNIQUE_DETID && idet<(int)detmap.size())idet++;
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;
    // Process hodoscope data
    if(idet>=0){// && T->Earm_BBHodoScint.nhits){
      for(int i = 0; i<T->Earm_BBHodoScint.nhits; i++){
	for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	  // Evaluation of number of photoelectrons and time from energy deposit documented at:
	  // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	  Npe = R->Poisson(1.0e7*T->Earm_BBHodoScint.sumedep->at(i)*0.113187*exp(-(0.3+pow(-1, j)*T->Earm_BBHodoScint.xhit->at(i))/1.03533)* 0.24);
	  t = tzero+T->Earm_BBHodoScint.tavg->at(i)+(0.55+pow(-1, j)*T->Earm_BBHodoScint.xhit->at(i))/0.15-pmtdets[idet]->fTrigOffset;
	  chan = T->Earm_BBHodoScint.cell->at(i)*2+j;
	  //T->Earm_BBHodoScint_hit_sumedep->at(i);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
	}
      }
      has_data = true;
    }

    // ** How to add a new subsystem **
    // Unfold here the data for the new detector
    //genrp detectors beam side
    while(idet<(int)detmap.size()){
      if(idet<0)idet++;
      if(detmap[idet]!=PRPOLBS_SCINT_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(detmap[idet]!=PRPOLBS_SCINT_UNIQUE_DETID && idet<(int)detmap.size()){
    if(idet>=detmap.size())idet = -1;
    //cout << " " << idet;
    // Process hodoscope data
    if(idet>=0){// && T->Harm_PRPolScintBeamSide.nhits){
      for(int i = 0; i<T->Harm_PRPolScintBeamSide.nhits; i++){
	for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	  Npe = R->Poisson(1.0e7*T->Harm_PRPolScintBeamSide.sumedep->at(i)*0.113187*exp(-(0.3+pow(-1, j)*T->Harm_PRPolScintBeamSide.xhit->at(i))/1.03533)* 0.24);
	  t = tzero+T->Harm_PRPolScintBeamSide.tavg->at(i)+(0.55+pow(-1, j)*T->Harm_PRPolScintBeamSide.xhit->at(i))/0.15-pmtdets[idet]->fTrigOffset;
	  chan = T->Harm_PRPolScintBeamSide.cell->at(i)*2+j;
	  pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
	}
      }
      has_data = true;
    }
    //genrp detectors far side scint
    while(idet<(int)detmap.size()){
      if(idet<0)idet++;
      if(detmap[idet]!=PRPOLFS_SCINT_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    if(idet>=detmap.size())idet = -1;
    if(idet>=0){// && T->Harm_PRPolScintFarSide.nhits){
      for(int i = 0; i<T->Harm_PRPolScintFarSide.nhits; i++){
	for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	  Npe = R->Poisson(1.0e7*T->Harm_PRPolScintFarSide.sumedep->at(i)*0.113187*exp(-(0.3+pow(-1, j)*T->Harm_PRPolScintFarSide.xhit->at(i))/1.03533)* 0.24);
	  t = tzero+T->Harm_PRPolScintFarSide.tavg->at(i)+(0.55+pow(-1, j)*T->Harm_PRPolScintFarSide.xhit->at(i))/0.15-pmtdets[idet]->fTrigOffset;
	  chan = T->Harm_PRPolScintFarSide.cell->at(i)*2+j;
	  pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
	}
      }
      has_data = true;
    }

    //genrp detectors Activeana
    while(idet<(int)detmap.size()){
      if(idet<0)idet++;
      if(detmap[idet]!=ACTIVEANA_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    if(idet>=detmap.size())idet = -1;
    if(idet>=0){// && T->Harm_ActAnScint.nhits){
      for(int i = 0; i<T->Harm_ActAnScint.nhits; i++){
	//for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
     //Npe = R->Poisson(Npe_edep_unit*T->Harm_ActAnScint.sumedep->at(i)); //Find Npe_edep_unit: Ave. amount of light produced per energy deposit 
// The number of photoelectrons for each PMT is the product of the raw number of photoelectrons produced
//(which depends on the energy deposit sumedep)times the light attenuation
//(which depends on the distance between the light production and the PMT)
// => Npe = (Npe_edep_unit*sumedep)*exp(-(|x_hit-x_PMT|)/Lambda)
     
       Npe = R->Poisson(1.0e7*T->Harm_ActAnScint.sumedep->at(i)*0.113187*exp(-(0.3+pow(-1,0)*T->Harm_ActAnScint.xhit->at(i))/1.03533)* 0.24);
	  t = tzero+T->Harm_ActAnScint.tavg->at(i)+(0.55+pow(-1,0)*T->Harm_ActAnScint.xhit->at(i))/0.15-pmtdets[idet]->fTrigOffset;
	  chan = T->Harm_ActAnScint.cell->at(i)*2;
	  pmtdets[idet]->PMTmap[chan].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, signal);
	//}
      }
      has_data = true;
    }
    
}//end if(!detmap.empty())
  
  //GEMs
  if(!gemmap.empty()){
    idet = 0;
    //genrp detectors
    while(idet<(int)gemmap.size()){
      if(gemmap[idet]!=BBGEM_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(gemmap[idet]!=BBGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
    if(idet>=gemmap.size())idet = -1;
    //cout << " gem " << idet << endl;
    // Now process the GEM data
    if(idet>=0){// && T->Earm_BBGEM.nhits){
      for(int k = 0; k<T->Earm_BBGEM.nhits; k++){
	if(T->Earm_BBGEM.edep->at(k)>0){
	  SBSDigGEMDet::gemhit hit; 
	  hit.source = signal;
	  if(T->Earm_BBGEM.plane->at(k)==5){
	    if(fabs(T->Earm_BBGEM.xin->at(k))>=1.024)continue;
	    mod = 12 + floor((T->Earm_BBGEM.xin->at(k)+1.024)/0.512);
	  }else{
	    if(fabs(T->Earm_BBGEM.xin->at(k))>=0.768)continue;
	    mod = (T->Earm_BBGEM.plane->at(k)-1)*3 + floor((T->Earm_BBGEM.xin->at(k)+0.768)/0.512);
	  }
	  //if(mod<2)cout << mod << " " << T->Earm_BBGEM.plane->at(k) << " " << T->Earm_BBGEM.xin->at(k) << endl;
	  hit.module = mod; 
	  hit.edep = T->Earm_BBGEM.edep->at(k)*1.0e9;//eV! not MeV!!!!
	  //hit.tmin = T->Earm_BBGEM_hit_tmin->at(k);
	  //hit.tmax = T->Earm_BBGEM_hit_tmax->at(k);
	  hit.t = tzero+T->Earm_BBGEM.t->at(k);
	  //cout << mod*2 << " " << gemdets[idet]->GEMPlanes[mod*2].Xoffset() << endl;
	  hit.xin = T->Earm_BBGEM.xin->at(k)-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	  hit.yin = T->Earm_BBGEM.yin->at(k);
	  hit.zin = T->Earm_BBGEM.zin->at(k)-bbgem_z[T->Earm_BBGEM.plane->at(k)-1]+0.8031825;
	  hit.xout = T->Earm_BBGEM.xout->at(k)-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	  hit.yout = T->Earm_BBGEM.yout->at(k);
	  hit.zout = T->Earm_BBGEM.zout->at(k)-bbgem_z[T->Earm_BBGEM.plane->at(k)-1]+0.8031825;
	  //cout << mod << " " << hit.xin << " " << hit.xout << endl;
	  gemdets[idet]->fGEMhits.push_back(hit);
	}//end if(sumedep>0)
	
      }
      has_data = true;  
    }
 }//end if(!gemmap.empty())...  
    
//CEPolFront GEMs
  if(!gemmap.empty()){
    idet = 0;
    //genrp detectors
    while(idet<(int)gemmap.size()){
      if(gemmap[idet]!=CEPOL_GEMFRONT_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(gemmap[idet]!=BBGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
    if(idet>=gemmap.size())idet = -1;
    //cout << " gem " << idet << endl;
    // Now process the GEM data
    if(idet>=0){// && T->Earm_BBGEM.nhits){
      for(int k = 0; k<T->Harm_CEPolFront.nhits; k++){
	if(T->Harm_CEPolFront.edep->at(k)>0){
	  SBSDigGEMDet::gemhit hit; 
	  hit.source = signal;
/*	  if(T->Harm_CEPolFront.plane->at(k)==5){
	    if(fabs(T->Harm_CEPolFront.xin->at(k))>=1.024)continue;
	    mod = 12 + floor((T->Harm_CEPolFront.xin->at(k)+1.024)/0.512);
	  }else{
	    if(fabs(T->Harm_CEPolFront.xin->at(k))>=0.768)continue;
	    mod = (T->Harm_CEPolFront->at(k)-1)*3 + floor((T->Harm_CEPolFront.xin->at(k)+0.768)/0.512);
	  }*/
	  //if(mod<2)cout << mod << " " << T->Earm_BBGEM.plane->at(k) << " " << T->Earm_BBGEM.xin->at(k) << endl;
	  hit.module = mod; 
	  hit.edep = T->Harm_CEPolFront.edep->at(k)*1.0e9;//eV! not MeV!!!!
	  //hit.tmin = T->Earm_BBGEM_hit_tmin->at(k);
	  //hit.tmax = T->Earm_BBGEM_hit_tmax->at(k);
	  hit.t = tzero+T->Harm_CEPolFront.t->at(k);
	  //cout << mod*2 << " " << gemdets[idet]->GEMPlanes[mod*2].Xoffset() << endl;
	  hit.xin = T->Harm_CEPolFront.xin->at(k);//-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	  hit.yin = T->Harm_CEPolFront.yin->at(k);
	  hit.zin = T->Harm_CEPolFront.zin->at(k)-cepol_front_z[T->Harm_CEPolFront.plane->at(k)-1]+0.8031825;
	  hit.xout = T->Harm_CEPolFront.xout->at(k);//-gemdets[idet]->Harm_CEPolFront[mod*2].Xoffset();
	  hit.yout = T->Harm_CEPolFront.yout->at(k);
	  hit.zout = T->Harm_CEPolFront.zout->at(k)-cepol_front_z[T->Harm_CEPolFront.plane->at(k)-1]+0.8031825;
	  //cout << mod << " " << hit.xin << " " << hit.xout << endl;
	  gemdets[idet]->fGEMhits.push_back(hit);
	}//end if(sumedep>0)
	
      }
      has_data = true;  
    }
 }//end if(!gemmap.empty())...

//CEPolRear GEMs
  if(!gemmap.empty()){
    idet = 0;
    //genrp detectors
    while(idet<(int)gemmap.size()){
      if(gemmap[idet]!=CEPOL_GEMREAR_UNIQUE_DETID){
	idet++;
      }else{
	break;
      }
    }
    //while(gemmap[idet]!=BBGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
    if(idet>=gemmap.size())idet = -1;
    //cout << " gem " << idet << endl;
    // Now process the GEM data
    if(idet>=0){// && T->Earm_BBGEM.nhits){
      for(int k = 0; k<T->Harm_CEPolRear.nhits; k++){
	if(T->Harm_CEPolRear.edep->at(k)>0){
	  SBSDigGEMDet::gemhit hit; 
	  hit.source = signal;
	  /*if(T->Harm_CEPolRear.plane->at(k)==5){
	    if(fabs(T->Harm_CEPolRear.xin->at(k))>=1.024)continue;
	    mod = 12 + floor((T->Harm_CEPolRear.xin->at(k)+1.024)/0.512);
	  }else{
	    if(fabs(T->Harm_CEPolRear.xin->at(k))>=0.768)continue;
	    mod = (T->Harm_CEPolRear->at(k)-1)*3 + floor((T->Harm_CEPolRear.xin->at(k)+0.768)/0.512);
	  }*/
	  //if(mod<2)cout << mod << " " << T->Earm_BBGEM.plane->at(k) << " " << T->Earm_BBGEM.xin->at(k) << endl;
	  hit.module = mod; 
	  hit.edep = T->Harm_CEPolRear.edep->at(k)*1.0e9;//eV! not MeV!!!!
	  //hit.tmin = T->Earm_BBGEM_hit_tmin->at(k);
	  //hit.tmax = T->Earm_BBGEM_hit_tmax->at(k);
	  hit.t = tzero+T->Harm_CEPolRear.t->at(k);
	  //cout << mod*2 << " " << gemdets[idet]->GEMPlanes[mod*2].Xoffset() << endl;
	  hit.xin = T->Harm_CEPolRear.xin->at(k);//-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	  hit.yin = T->Harm_CEPolRear.yin->at(k);
	  hit.zin = T->Harm_CEPolRear.zin->at(k)-cepol_rear_z[T->Harm_CEPolRear.plane->at(k)-1]+0.8031825;
	  hit.xout = T->Harm_CEPolRear.xout->at(k);//-gemdets[idet]->Harm_CEPolRear[mod*2].Xoffset();
	  hit.yout = T->Harm_CEPolRear.yout->at(k);
	  hit.zout = T->Harm_CEPolRear.zout->at(k)-cepol_rear_z[T->Harm_CEPolRear.plane->at(k)-1]+0.8031825;
	  //cout << mod << " " << hit.xin << " " << hit.xout << endl;
	  gemdets[idet]->fGEMhits.push_back(hit);
	}//end if(sumedep>0)
	
      }
      has_data = true;  
    }
 
  }//end if(!gemmap.empty())...
  return has_data;
}

