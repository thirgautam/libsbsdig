#include "g4sbs_types.h"
#include "SBSDigBkgdGen.h"

#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

SBSDigBkgdGen::SBSDigBkgdGen()
{
  NhitsBBGEMs = new Double_t[5];
  h_xhitBBGEMs = new TH1D*[5];
  h_yhitBBGEMs = new TH1D*[5];
  h_dxhitBBGEMs = new TH1D*[5];
  h_dyhitBBGEMs = new TH1D*[5];
  NhitsHCal = new Double_t[288];
  NhitsBBPS = new Double_t[52];
  NhitsBBSH = new Double_t[189];
  NhitsBBHodo = new Double_t[90];
  P1hitGRINCH = new Double_t[510];
  P2hitsGRINCH = new Double_t[510];
}

SBSDigBkgdGen::SBSDigBkgdGen(TFile* f_bkgd)
{
  NhitsBBGEMs = new Double_t[5];
  h_xhitBBGEMs = new TH1D*[5];
  h_yhitBBGEMs = new TH1D*[5];
  h_dxhitBBGEMs = new TH1D*[5];
  h_dyhitBBGEMs = new TH1D*[5];
  NhitsHCal = new Double_t[288];
  NhitsBBPS = new Double_t[52];
  NhitsBBSH = new Double_t[189];
  NhitsBBHodo = new Double_t[90];
  P1hitGRINCH = new Double_t[510];
  P2hitsGRINCH = new Double_t[510];
  
  Initialize(f_bkgd);
}

SBSDigBkgdGen::~SBSDigBkgdGen()
{
}

void SBSDigBkgdGen::Initialize(TFile* f_bkgd)
{
  double mu, sigma;
  
  // GEMs
  TH1D* h1_BBGEM_nhits_[5];
  TF1* f1_gemnhits_[5];
  TH2D* h1_BBGEM_yVsx_[5];
  TH2D* h1_BBGEM_dyVsdx_[5];
  TH1D* h1_BBGEM_Edep_[5];
  
  cout << "GEMs" << endl;
  
  for(int m = 0; m<5; m++){
    //Nhits
    h1_BBGEM_nhits_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_nhits_%d",m));
    f1_gemnhits_[m] = new TF1(Form("f1_gemnhits_%d", m), "gaus", 100., 400.);
    h1_BBGEM_nhits_[m]->Fit(f1_gemnhits_[m], "QRN");
    mu = f1_gemnhits_[m]->GetParameter(1);
    sigma = f1_gemnhits_[m]->GetParameter(2);
    if(mu>=0){
      f1_gemnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
      h1_BBGEM_nhits_[m]->Fit(f1_gemnhits_[m], "QRN");
    }
    NhitsBBGEMs[m] = max(1.0, f1_gemnhits_[m]->GetParameter(1));
    
    h1_BBGEM_yVsx_[m] = (TH2D*)f_bkgd->Get(Form("h1_BBGEM_yVsx_%d",m));
    h_xhitBBGEMs[m] = h1_BBGEM_yVsx_[m]->ProjectionX(Form("h1_xhitBBGEM_%d",m));
    h_yhitBBGEMs[m] = h1_BBGEM_yVsx_[m]->ProjectionY(Form("h1_yhitBBGEM_%d",m));
    
    h1_BBGEM_dyVsdx_[m] = (TH2D*)f_bkgd->Get(Form("h1_BBGEM_dyVsdx_%d",m));
    h_dxhitBBGEMs[m] = h1_BBGEM_dyVsdx_[m]->ProjectionX(Form("h1_dxhitBBGEM_%d",m));
    h_dyhitBBGEMs[m] = h1_BBGEM_dyVsdx_[m]->ProjectionY(Form("h1_dyhitBBGEM_%d",m));
    
    h1_BBGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_Edep_%d",m));
    //h1_BBGEM_Edep_[m] = (TH1D*)f_bkgd->Get(Form("h1_BBGEM_Edep_log_%d",m));
    
    if(m==0){
      h_EdephitBBGEMs = h1_BBGEM_Edep_[m];
    }else{
      h_EdephitBBGEMs->Add(h1_BBGEM_Edep_[m]);
    }
    
    //cout << m << " " << NhitsBBGEMs[m] << endl;
  }
  
  //HCal
  cout << "HCal" << endl;
  TH2D *h1_HCal_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_HCal_nhitsVsChan");
  TH1D* h1_HCal_nhits_[288];
  TF1* f1_hcalnhits_[288];
  TH2D *h1_HCal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_EdepHitVsChan");
  //TH2D *h1_HCal_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_EdepHitVsChan_log");
  TH2D *h1_HCal_zHitVsChan = (TH2D*)f_bkgd->Get("h1_HCal_zHitVsChan");

  for(int m = 0; m<288; m++){
    h1_HCal_nhits_[m] = h1_HCal_nhitsVsChan->ProjectionY(Form("h1_HCal_nhits_%d", m), m+1, m+1);
    f1_hcalnhits_[m] = new TF1(Form("f1_hcalnhits_%d", m), "gaus", 0, 50);
    h1_HCal_nhits_[m]->Fit(f1_hcalnhits_[m], "QR0");
    mu = f1_hcalnhits_[m]->GetParameter(1);
    sigma = f1_hcalnhits_[m]->GetParameter(2);
    if(mu>=0){
      f1_hcalnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
      h1_HCal_nhits_[m]->Fit(f1_hcalnhits_[m], "QR0");
    }
    NhitsHCal[m] = max(1.0, f1_hcalnhits_[m]->GetParameter(1));
    //cout << m << " " << NhitsHCal[m] << endl;
  }
  h_EdephitHCal = h1_HCal_EdepHitVsChan->ProjectionY("h_EdephitHCal");
  h_zhitHCal = h1_HCal_zHitVsChan->ProjectionY("h_zhitHCal");

  //PS
  cout << "PS" << endl;
  TH2D *h1_BBPS_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_nhitsVsChan");
  TH1D* h1_BBPS_nhits_[52];
  TF1* f1_bbpsnhits_[52];
  TH2D *h1_BBPS_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_EdepHitVsChan");
  //TH2D *h1_BBPS_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBPS_EdepHitVsChan_log");
  
  for(int m = 0; m<52; m++){
    h1_BBPS_nhits_[m] = h1_BBPS_nhitsVsChan->ProjectionY(Form("h1_BBPS_nhits_%d", m), m+1, m+1);
    f1_bbpsnhits_[m] = new TF1(Form("f1_bbpsnhits_%d", m), "gaus", 0, 150);
    h1_BBPS_nhits_[m]->Fit(f1_bbpsnhits_[m], "QR0");
    mu = f1_bbpsnhits_[m]->GetParameter(1);
    sigma = f1_bbpsnhits_[m]->GetParameter(2);
    if(mu>=0){
      f1_bbpsnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
      h1_BBPS_nhits_[m]->Fit(f1_bbpsnhits_[m], "QR0");
    }
    
    NhitsBBPS[m] = max(1.0, f1_bbpsnhits_[m]->GetParameter(1));
    //cout << m << " " << NhitsBBPS[m] << endl;
  }
  h_EdephitBBPS = h1_BBPS_EdepHitVsChan->ProjectionY("h_EdephitBBPS");
  
  //SH
  cout << "SH" << endl;
  TH2D *h1_BBSH_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_nhitsVsChan");
  TH1D* h1_BBSH_nhits_[189];
  TF1* f1_bbshnhits_[189];
  TH2D *h1_BBSH_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_EdepHitVsChan");
  //TH2D *h1_BBSH_EdepHitVsChan = (TH2D*)f_bkgd->Get("h1_BBSH_EdepHitVsChan_log");
  
  for(int m = 0; m<189; m++){
    h1_BBSH_nhits_[m] = h1_BBSH_nhitsVsChan->ProjectionY(Form("h1_BBSH_nhits_%d", m), m+1, m+1);
    f1_bbshnhits_[m] = new TF1(Form("f1_bbshnhits_%d", m), "gaus", 0, 150);
    h1_BBSH_nhits_[m]->Fit(f1_bbshnhits_[m], "QR0");
    mu = f1_bbshnhits_[m]->GetParameter(1);
    sigma = f1_bbshnhits_[m]->GetParameter(2);
    if(mu>=0){
      f1_bbshnhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
      h1_BBSH_nhits_[m]->Fit(f1_bbshnhits_[m], "QR0");
    }
    NhitsBBSH[m] = max(1.0, f1_bbshnhits_[m]->GetParameter(1));
    //cout << m << " " << NhitsBBSH[m] << endl;
  }
  h_EdephitBBSH = h1_BBSH_EdepHitVsChan->ProjectionY("h_EdephitBBSH");
  
  //BB Hodo
  cout << "Hodo" << endl;
  TH2D *h1_BBHodo_nhitsVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_nhitsVsSlat");
  TH1D* h1_BBHodo_nhits_[90];
  TF1* f1_bbhodonhits_[90];
  
  TH2D *h1_BBHodo_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_EdepHitVsSlat");
  //TH2D *h1_BBHodo_EdepHitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_EdepHitVsSlat_log");
  TH2D *h1_BBHodo_xhitVsSlat = (TH2D*)f_bkgd->Get("h1_BBHodo_xhitVsSlat");
  
  for(int m = 0; m<90; m++){
    h1_BBHodo_nhits_[m] = h1_BBHodo_nhitsVsSlat->ProjectionY(Form("h1_BBHodo_nhits_%d", m), m+1, m+1);
    f1_bbhodonhits_[m] = new TF1(Form("f1_bbhodonhits_%d", m), "gaus", 0, 100);
    h1_BBHodo_nhits_[m]->Fit(f1_bbhodonhits_[m], "QR0");
    mu = f1_bbhodonhits_[m]->GetParameter(1);
    sigma = f1_bbhodonhits_[m]->GetParameter(2);
    if(mu>=0){
      f1_bbhodonhits_[m]->SetRange(mu-2*sigma, mu+2*sigma);
      h1_BBHodo_nhits_[m]->Fit(f1_bbhodonhits_[m], "QR0");
    }
    NhitsBBHodo[m] = max(1.0, f1_bbhodonhits_[m]->GetParameter(1));
    //cout << m << " " << NhitsBBHodo[m] << endl;
  }
  
  h_EdephitBBHodo = h1_BBHodo_EdepHitVsSlat->ProjectionY("h_EdephitBBHodo");
  h_xhitBBHodo = h1_BBHodo_xhitVsSlat->ProjectionY("h_xhitBBHodo");
  
  //GRINCH
  cout << "GRINCH" << endl;
  TH2D *h1_GRINCH_nhitsVsChan = (TH2D*)f_bkgd->Get("h1_GRINCH_nhitsVsChan");
  // TH1D* h1_GRINCH_Chan_1hit = h1_GRINCH_nhitsVsChan->ProjectionX("h1_GRINCH_Chan_1hit", 2, 2);
  // TH1D* h1_GRINCH_Chan_2hits = h1_GRINCH_nhitsVsChan->ProjectionX("h1_GRINCH_Chan_2hit", 3, -1);
  TH2D *h1_GRINCH_NpeVsChan = (TH2D*)f_bkgd->Get("h1_GRINCH_NpeVsChan");
  
  for(int m = 1; m<=510; m++){
    P1hitGRINCH[m] = h1_GRINCH_nhitsVsChan->Integral(m, m, 2, 2)/h1_GRINCH_nhitsVsChan->Integral(m, m, 0, -1);
    P2hitsGRINCH[m] = h1_GRINCH_nhitsVsChan->Integral(m, m, 3, -1)/h1_GRINCH_nhitsVsChan->Integral(m, m, 0, -1);
  }
  h_NpeGRINCH = h1_GRINCH_NpeVsChan->ProjectionY("h1_GRINCH_Npe");  
}



void SBSDigBkgdGen::GenerateBkgd(//double theta_sbs, double d_hcal, 
				 TRandom3* R, 
				 std::vector<SBSDigPMTDet*> pmtdets,
				 std::vector<int> detmap, 
				 std::vector<SBSDigGEMDet*> gemdets, 
				 std::vector<int> gemmap, 
				 double lumifrac)
{
  int nhits;
  int mod;
  double edep;
  int Npe;
  double t, p;
  double z_hit, Npe_Edep_ratio, sigma_tgen;
  double beta, sin2thetaC;
  double x_hit, y_hit;
  int idet = 0;
  //ordering by increasing unique det ID
  while(detmap[idet]!=HCAL_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "hcal" << endl;
    for(int m = 0; m<288; m++){
      nhits = R->Poisson(NhitsHCal[m]*lumifrac);
      //cout << m << " " << NhitsHCal[m]*lumifrac << " " << nhits << endl;
      for(int i = 0; i<nhits; i++){
	edep = h_EdephitHCal->GetRandom();//R);
	z_hit = h_zhitHCal->GetRandom();//R); //R->Uniform(-0.91, 0.91);//for the time being
	//cout << " " << i << " " << edep << " " << z_hit << endl;
	Npe_Edep_ratio = 5.242+11.39*z_hit+10.41*pow(z_hit, 2);
	Npe = R->Poisson(Npe_Edep_ratio*edep*1.0e3);
	sigma_tgen = 0.4244+11380/pow(Npe+153.4, 2);
	
	t = R->Uniform(-40.,40.);//by hand for the time being
	
	pmtdets[idet]->PMTmap[m].Fill(Npe, pmtdets[idet]->fThreshold, t, sigma_tgen, 1);
      }
    }
  }
  
  while(detmap[idet]!=BBPS_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "ps" << endl;
    for(int m = 0; m<52; m++){
      nhits = R->Poisson(NhitsBBPS[m]*lumifrac);
      //cout << m << " " << NhitsBBPS[m]*lumifrac << " " << nhits << endl;
      for(int i = 0; i<nhits; i++){
	edep = h_EdephitBBPS->GetRandom();//R);
	
	if(edep<1.e-4)continue;
	//check probability to generate p.e. yield
	bool genpeyield = true;
	if(edep<1.e-2)
	  genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*edep);
	//if we're go, let's generate
	if(genpeyield){
	  beta = sqrt( pow(m_e+edep, 2)-m_e*m_e )/(m_e + edep);
	  sin2thetaC = TMath::Max(1.-1./pow(n_lg*beta, 2), 0.);
	  //1500. Used to be 454.: just wrong
	  Npe = R->Poisson(300.0*edep*sin2thetaC/(1.-1./(n_lg*n_lg)) );
	  t = R->Uniform(-50.,50.);
	  
	  //cout << " " << i << " " << edep << " " << Npe << endl;
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, 1);
	}
      }
    }
  }
  
  while(detmap[idet]!=BBSH_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "sh" << endl;
    for(int m = 0; m<189; m++){
      nhits = R->Poisson(NhitsBBSH[m]*lumifrac);
      for(int i = 0; i<nhits; i++){
	edep = h_EdephitBBSH->GetRandom();//R);
	
	if(edep<1.e-4)continue;
	//check probability to generate p.e. yield
	bool genpeyield = true;
	if(edep<1.e-2)
	  genpeyield = R->Uniform(0, 1)<=1-exp(0.29-950.*edep);
	//if we're go, let's generate
	if(genpeyield){
	  beta = sqrt( pow(m_e+edep, 2)-m_e*m_e )/(m_e + edep);
	  sin2thetaC = TMath::Max(1.-1./pow(n_lg*beta, 2), 0.);
	  //1500. Used to be 454.: just wrong
	  Npe = R->Poisson(360.0*edep*sin2thetaC/(1.-1./(n_lg*n_lg)) );
	  t = R->Uniform(-50.,50.);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, 0, t, 1);
	}
      }
    }
  }
  
  while(detmap[idet]!=GRINCH_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  if(idet>=0){
    //cout << "grinch" << endl;
    for(int m = 0; m<510; m++){
      p = R->Uniform(0, 1);
      if(p<P2hitsGRINCH[m]*lumifrac){
	nhits = 2;
      }else if(p<P1hitGRINCH[m]*lumifrac)nhits = 1;
      
      for(int i = 0; i<nhits; i++){
	t = R->Uniform(-50.,50.);
	Npe = h_NpeGRINCH->GetRandom();
	
	pmtdets[idet]->PMTmap[m].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
      }
    }
  }
  
  while(detmap[idet]!=HODO_UNIQUE_DETID && idet<(int)detmap.size())idet++;
  if(idet>=detmap.size())idet = -1;
  //cout << " " << idet;
  // Process hodoscope data
  if(idet>=0){
    //cout << "hodo" << endl;
    for(int m = 0; m<90; m++){
      nhits = R->Poisson(NhitsBBHodo[m]*lumifrac);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitBBHodo->GetRandom()*1.e6;
	x_hit =  h_xhitBBHodo->GetRandom();
	
	p = R->Uniform(-50.,50.);
	for(int j = 0; j<2; j++){//j = 0: close beam PMT, j = 1: far beam PMT
	  // Evaluation of number of photoelectrons and time from energy deposit documented at:
	  // https://hallaweb.jlab.org/dvcslog/SBS/170711_172759/BB_hodoscope_restudy_update_20170711.pdf
	  Npe = R->Poisson(1.0e7*edep*0.113187*exp(-(0.3+pow(-1, j)*x_hit)/1.03533)* 0.24);
	  t = p+(0.55+pow(-1, j)*x_hit)/0.15;
	  //T->Earm_BBHodoScint_hit_sumedep->at(i);
	  //if(chan>pmtdets[idet]->fNChan)cout << chan << endl;
	  pmtdets[idet]->PMTmap[m*2+j].Fill(pmtdets[idet]->fRefPulse, Npe, pmtdets[idet]->fThreshold, t, 1);
	}
      }
    }
  }
  
  idet = 0;
  while(gemmap[idet]!=BBGEM_UNIQUE_DETID && idet<(int)gemmap.size())idet++;
  if(idet>=gemmap.size())idet = -1;
  
  if(idet>=0){
    //    cout << "bbgems" << endl;
    for(int m = 0; m<5; m++){
      nhits = R->Poisson(NhitsBBHodo[m]*lumifrac);
      for(int i = 0; i<nhits; i++){
	edep =  h_EdephitBBGEMs->GetRandom();
	x_hit =  h_xhitBBGEMs[m]->GetRandom();
	y_hit =  h_yhitBBGEMs[m]->GetRandom();
	
	//x_hit =  h_dxhitBBGEMs[m]->GetRandom();
	//y_hit =  h_dyhitBBGEMs[m]->GetRandom();
	SBSDigGEMDet::gemhit hit; 
	hit.source = 1;
	if(m==4){
	  if(fabs(x_hit)>=1.024)continue;
	  mod = 12 + floor((x_hit+1.024)/0.512);
	}else{
	  if(fabs(x_hit)>=0.768)continue;
	  mod = m*3 + floor((x_hit+0.768)/0.512);
	}
	hit.module = mod; 
	hit.edep = edep;
	
	hit.xin = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset();
	hit.yin = y_hit;
	hit.zin = -0.0015;
	hit.xout = x_hit-gemdets[idet]->GEMPlanes[mod*2].Xoffset()+h_dxhitBBGEMs[m]->GetRandom();
	hit.yout = y_hit+h_dyhitBBGEMs[m]->GetRandom();
	hit.zout = 0.0015;
	
	gemdets[idet]->fGEMhits.push_back(hit);
      }
    }
  }
  
}
