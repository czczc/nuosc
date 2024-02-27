// Barger, Marfatia and Whisnant PRD 65 073023 (2002)
Double_t numu2nue_Barger(Double_t E, Double_t delta=3.1415926/2., Double_t mh=1, Double_t type=1, Double_t L=1300)
{
  // E GeV 
  // delta: CP phase
  // mh: normal 1; inverted -1 
  // type: numu 1; antinumu -1
  // L km 
    
  Double_t s13 =sin(asin(sqrt(0.09))/2.);
  Double_t c13 = sqrt(1.-s13*s13);
  // s23
  // s23^2 = 0.5 +0.07 - 0.06 
  Double_t s23 = sin(45./180.*3.1415926);
  Double_t c23 = sqrt(1.-s23*s23);
  //s12
  // s12 = 0.304 +0.022-0.016
  Double_t s12 = sqrt(0.304);
  Double_t c12 = sqrt(1.-s12*s12);
  
  // Delta_m_31 = 2.4 +0.12 - 0.11  * e-3 eV^2
  Double_t Dm31 = 2.41e-3+7.59e-5;
  // Delta_m_12 = 7.65 + 0.23 - 0.2 * e-5 eV^2
  Double_t Dm21 = 7.59e-5;
  
  Dm31 = mh * Dm31;
  
  Double_t Dm12 = Dm21 * (-1.);
  Double_t Dm13 = Dm31 * (-1.);
  
  Double_t Dm32 = Dm31 - Dm21;
  Double_t Dm23 = Dm32 * (-1.);

  Double_t prob = 0;

  // rho  = 3.0 g/cm^3
  Double_t a = 1.54e-4 * E * 2.7  * 0.5;
  //Double_t a = 1.54e-4 * E * 5.5  * 0.5;


  Double_t Ahat =  fabs(a / Dm31);
  
  Double_t Delta = fabs(1.267*Dm31*L/E);
  Double_t alpha = fabs(Dm21/Dm31);
  Double_t x = s23*2.*s13*c13;
  Double_t y =alpha * c23*2.*s12*c12;
  Double_t f = sin((1-Ahat)*Delta)/(1-Ahat);
  Double_t fbar = sin((1+Ahat)*Delta)/(1+Ahat);
  Double_t g =  sin(Ahat * Delta)/Ahat;
  if (type == 1 && Dm31 >0){
    prob = x*x*f*f + 2*x*y*f*g*(cos(delta)*cos(Delta)-sin(delta)*sin(Delta)) + y*y*g*g;
  }else if (type==-1 && Dm31 > 0){
    prob = x*x*fbar*fbar + 2*x*y*fbar*g*(cos(delta)*cos(Delta)+sin(delta)*sin(Delta)) + y*y*g*g;
  }else if (type==1 && Dm31 < 0){
    prob = x*x*fbar*fbar - 2*x*y*fbar*g*(cos(delta)*cos(Delta)+sin(delta)*sin(Delta)) + y*y*g*g;
  }else if (type==-1 && Dm31 < 0){
    prob = x*x*f*f - 2*x*y*f*g*(cos(delta)*cos(Delta) - sin(delta)*sin(Delta)) + y*y*g*g;
  }
  
  
  return prob;
  
}


Double_t numu2nue(Double_t E, Double_t delta=3.1415926/2., Double_t mh=1, Double_t type=1, Double_t L=1300)
{
  // E GeV 
  // delta: CP phase
  // mh: normal 1; inverted -1 
  // type: numu 1; antinumu -1
  // L km 
  
  Double_t a1=2.7; // density g/cm^3
    
  Double_t s13 =sin(asin(sqrt(0.09))/2.);
  Double_t c13 = sqrt(1.-s13*s13);
  // s23
  // s23^2 = 0.5 +0.07 - 0.06 
  Double_t s23 = sin(45./180.*3.1415926);
  Double_t c23 = sqrt(1.-s23*s23);
  //s12
  // s12 = 0.304 +0.022-0.016
  Double_t s12 = sqrt(0.304);
  Double_t c12 = sqrt(1.-s12*s12);
  
  // Delta_m_31 = 2.4 +0.12 - 0.11  * e-3 eV^2
  Double_t Dm31 = 2.41e-3+7.59e-5;
  // Delta_m_12 = 7.65 + 0.23 - 0.2 * e-5 eV^2
  Double_t Dm21 = 7.59e-5;
  Dm31 = mh * Dm31;
  
  Double_t Dm12 = Dm21 * (-1.);
  Double_t Dm13 = Dm31 * (-1.);
  
  Double_t Dm32 = Dm31 - Dm21;
  Double_t Dm23 = Dm32 * (-1.);

  Double_t prob = 0;

  // rho  = 3.6 g/cm^3
  Double_t a = 1.54e-4 * E * a1 * 0.5;
  
  Double_t x=Dm21 + Dm31 + a;
  Double_t y=Dm21*Dm31+a*(Dm21*(1-s12*s12*c13*c13)+Dm31*(1-s13*s13));
  Double_t z=cos(1./3.*acos((2*pow(x,3)-9*x*y+27*a*Dm21*Dm31*c12*c12*c13*c13)/(2.*pow(x*x-3*y,3./2.))));

  Double_t Dmt1 = 1./3.*x-1./3.*sqrt(x*x-3*y)*(z+sqrt(3.*(1-z*z)));
  Double_t Dmt2 = 1./3.*x-1./3.*sqrt(x*x-3*y)*(z-sqrt(3.*(1-z*z)));
  Double_t Dmt3 = 1./3.*x+2./3.*z*sqrt(x*x-3*y);
  
  Double_t Dm1 = 0.;
  Double_t Dm2 = Dm21;
  Double_t Dm3 = Dm31;

  Double_t Dh21 = Dm2 - Dmt1;
  Double_t Dh31 = Dm3 - Dmt1;
  Double_t Dh32 = Dm3 - Dmt2;
  
  Double_t Dh12 = Dm1 - Dmt2;
  Double_t Dh13 = Dm1 - Dmt3;
  Double_t Dh23 = Dm2 - Dmt3;
  
  Double_t Dh11 = Dm1 - Dmt1;
  Double_t Dh22 = Dm2 - Dmt2;
  Double_t Dh33 = Dm3 - Dmt3;


  Double_t Dt21 = Dmt2 - Dmt1;
  Double_t Dt31 = Dmt3 - Dmt1;
  Double_t Dt32 = Dmt3 - Dmt2;

  Double_t Dt12 = Dmt1 - Dmt2;
  Double_t Dt13 = Dmt1 - Dmt3;
  Double_t Dt23 = Dmt2 - Dmt3;


  TComplex Ve1(c12*c13,0);
  TComplex Ve2(s12*c13,0);
  TComplex Ve3(s13*cos(delta),-s13*sin(delta));

  TComplex Vmu1(-s12*c23-c12*s23*s13*cos(delta),-c12*s23*s13*sin(delta));
  TComplex Vmu2(c12*c23-s12*s23*s13*cos(delta),-s12*s23*s13*sin(delta));
  TComplex Vmu3(s23*c13,0);


  TComplex Ve1s = TComplex::Conjugate(Ve1);
  TComplex Ve2s = TComplex::Conjugate(Ve2);
  TComplex Ve3s = TComplex::Conjugate(Ve3);

  TComplex Vmu1s = TComplex::Conjugate(Vmu1);
  TComplex Vmu2s = TComplex::Conjugate(Vmu2);
  TComplex Vmu3s = TComplex::Conjugate(Vmu3);
  
  
  TComplex Ve1Vmu1s = Dh21*Dm31/Dt21/Dt31*Ve1*Vmu1s + Dh11*Dm32/Dt21/Dt31*Ve2*Vmu2s;
  TComplex Ve2Vmu2s = Dh32*Dm21/Dt32/Dt21*Ve2*Vmu2s + Dh22*Dm31/Dt32/Dt21*Ve3*Vmu3s;
  TComplex Ve3Vmu3s = Dh13*Dm23/Dt13/Dt23*Ve3*Vmu3s + Dh33*Dm21/Dt13/Dt23*Ve1*Vmu1s;

  TComplex Vmu1Ve1s = TComplex::Conjugate(Ve1Vmu1s);
  TComplex Vmu2Ve2s = TComplex::Conjugate(Ve2Vmu2s);
  TComplex Vmu3Ve3s = TComplex::Conjugate(Ve3Vmu3s);

  TComplex v12 = Vmu1Ve1s*Ve2Vmu2s;
  TComplex v13 = Vmu1Ve1s*Ve3Vmu3s;
  TComplex v23 = Vmu2Ve2s*Ve3Vmu3s;

  prob -= 4*v12.Re()*pow(sin(1.267*Dt21*L/E),2);
  prob -= 4*v13.Re()*pow(sin(1.267*Dt31*L/E),2);
  prob -= 4*v23.Re()*pow(sin(1.267*Dt32*L/E),2);

  Double_t J=0;

  TComplex vv = Vmu1Ve1s*Ve2Vmu2s;
  J = -vv.Im();

  prob -= 8*J*sin(1.267*Dt21*L/E)*sin(1.267*Dt31*L/E)*sin(1.267*Dt32*L/E);
  
  return prob;
}

void SetStyle()
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleOffset(1.0, "x");
    gStyle->SetTitleOffset(1.0, "y");
    gStyle->SetTitleFont(42, "P");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);
    gStyle->SetLabelSize(0.05, "x");
    gStyle->SetLabelSize(0.05, "y");
    gStyle->SetHistLineWidth(2);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetMarkerSize(0.3);
}

void plot_numu2nue() 
{
  gROOT->ForceStyle();
  SetStyle();

  int N = 400;
  TH1F *h1 = new TH1F("h1", "", N, 0.3, 6.3);
  TH1F *h2 = new TH1F("h2", "", N, 0.3, 6.3);
  TH1F *h3 = new TH1F("h3", "", N, 0.3, 6.3);
  TH1F *hb = new TH1F("hb", "", N, 0.3, 6.3);
  
  const double pi=3.1415926;
  double e, prob;
  double delta1 = 0;
  double delta2 = pi/2;
  double delta3 = -pi/2;

  for (int i=1; i<=N; i++){
    e = h1->GetBinCenter(i);
    prob = numu2nue(e, delta1);
    h1->SetBinContent(i, prob);
    prob = numu2nue(e, delta2);
    h2->SetBinContent(i, prob);
    prob = numu2nue(e, delta3);
    h3->SetBinContent(i, prob);

    prob = numu2nue_Barger(e, delta2);
    hb->SetBinContent(i, prob);
  }
  h1->GetXaxis()->SetTitle("E (GeV)");
  h1->GetYaxis()->SetTitle("Probability");
  h1->Draw();
  h2->SetLineColor(kRed);
  h2->Draw("same");
  h3->SetLineColor(kBlue);
  h3->Draw("same");

  // hb->SetLineColor(kGreen);
  // hb->Draw("same");

  TLegend *leg = new TLegend(0.65,0.60, 0.87,0.87);
  leg->AddEntry(h1, " #delta=0", "l");
  leg->AddEntry(h2,   " #delta=#pi/2", "l");
  leg->AddEntry(h3,   " #delta=-#pi/2", "l");
  // leg->AddEntry(hb,   " Barger", "l");
  leg->SetFillColor(kWhite);
  leg->Draw();
  gPad->SetGridx();
  gPad->SetGridy();


  TCanvas *c2 = new TCanvas();
  TFile *f = new TFile("flux/flux_dune_neutrino_FD.root");
  TH1D *numu_cceventrate = (TH1D*)f->Get("numu_cceventrate");

  double binW = 0.125; //GeV
  N = 50;
  TH1D *hOsc1 = (TH1D*)numu_cceventrate->Clone("hOsc1");
  TH1D *hOsc2 = (TH1D*)numu_cceventrate->Clone("hOsc2");
  TH1D *hOsc3 = (TH1D*)numu_cceventrate->Clone("hOsc3");
  // scale to 40 kt.year (from /kt/POT; assume 1.2 MW, 120 GeV proton)
  // 1.2 MW . 86400 sec . 365 d / 120 GeV / 1.6e-19 = 2e21
  const double scale = 2e21*40; 
  for (int i=1; i<=N; i++){
    e = numu_cceventrate->GetBinCenter(i);
    prob = numu2nue(e, delta1);
    hOsc1->SetBinContent(i, prob*numu_cceventrate->GetBinContent(i)*scale);
    prob = numu2nue(e, delta2);
    hOsc2->SetBinContent(i, prob*numu_cceventrate->GetBinContent(i)*scale);
    prob = numu2nue(e, delta3);
    hOsc3->SetBinContent(i, prob*numu_cceventrate->GetBinContent(i)*scale);
  }
  hOsc3->GetXaxis()->SetRangeUser(0, 6);
  hOsc3->GetXaxis()->SetTitle("E (GeV)");
  hOsc3->GetYaxis()->SetTitle("#nu_{e}CC Events / (40 kTon #upoint y #upoint 1.2 MW)");
  hOsc3->SetTitle("");
  hOsc3->SetLineColor(kBlue);
  hOsc3->Draw("hist");

  hOsc1->Draw("hist same");
  hOsc2->SetLineColor(kRed);
  hOsc2->Draw("hist same");

  TLegend *leg2 = new TLegend(0.65,0.60, 0.87,0.87);
  leg2->AddEntry(h1, " #delta=0", "l");
  leg2->AddEntry(h2,   " #delta=#pi/2", "l");
  leg2->AddEntry(h3,   " #delta=-#pi/2", "l");
  leg2->SetFillColor(kWhite);
  leg2->Draw();
  gPad->SetGridx();
  gPad->SetGridy();

}