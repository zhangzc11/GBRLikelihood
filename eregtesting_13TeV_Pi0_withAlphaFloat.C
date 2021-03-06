#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHybridBDTAutoPdf.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "TRandom.h"
#include "TGraph.h"
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "TFile.h"
#include "TNtuple.h"
#include "HybridGBRForest.h"
#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
//#include "HZZ2L2QRooPdfs.h"
#include "RooDoubleCBFast.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TAxis.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
#include "RooLinearVar.h"
#include <TStyle.h>


using namespace RooFit;

Double_t FWHM(TH1 * hist)
{
int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
Double_t fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
return fwhm;
}
//effsigma function from Chris
Double_t effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
//   if(total < 100.) {
//     cout << "effsigma: Too few entries " << total << endl;
//     return 0.;
//   }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}

void eregtesting_13TeV_Pi0(bool dobarrel=true, bool doele=false,int gammaID=0) {
  
  //output dir
  TString EEorEB = "EE";
  if(dobarrel)
	{
	EEorEB = "EB";
	}
  TString gammaDir = "bothGammas";
  if(gammaID==1)
  {
   gammaDir = "gamma1";
  }
  else if(gammaID==2)
  {
   gammaDir = "gamma2";
  }
  TString dirname = TString::Format("ereg_test_plots/%s_%s",gammaDir.Data(),EEorEB.Data());
  
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);    
  
  //read workspace from training
  TString fname;
  if (doele && dobarrel) 
    fname = "wereg_ele_eb.root";
  else if (doele && !dobarrel) 
    fname = "wereg_ele_ee.root";
  else if (!doele && dobarrel) 
    fname = "wereg_ph_eb.root";
  else if (!doele && !dobarrel) 
    fname = "wereg_ph_ee.root";
  
  TString infile = TString::Format("../../ereg_ws/%s/%s",gammaDir.Data(),fname.Data());
  
  TFile *fws = TFile::Open(infile); 
  RooWorkspace *ws = (RooWorkspace*)fws->Get("wereg");
  
  //read variables from workspace
  RooGBRTargetFlex *meantgt = static_cast<RooGBRTargetFlex*>(ws->arg("sigmeant"));  
  RooRealVar *tgtvar = ws->var("tgtvar");
  
  
  RooArgList vars;
  vars.add(meantgt->FuncVars());
  vars.add(*tgtvar);
   
  //read testing dataset from TTree
  RooRealVar weightvar("weightvar","",1.);

  TTree *dtree;
  
  if (doele) {
    //TFile *fdin = TFile::Open("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
    TFile *fdin = TFile::Open("/data/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");

    TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterSingleInvert");
    dtree = (TTree*)ddir->Get("hPhotonTreeSingle");       
  }
  else {
    if(dobarrel)
    {
    TFile *fdin = TFile::Open("/afs/cern.ch/work/z/zhicaiz/public/ECALpro_MC_TreeForRegression/Gun_Pi0_Pt1To15_FlatPU0to50RAW_withHLT_80X_mcRun2_GEN-SIM-RAW_ALL_EcalNtp_ALL_EB_combine_test.root");//("root://eoscms.cern.ch///eos/cms/store/cmst3/user/bendavid/idTreesAug1/hgg-2013Final8TeV_ID_s12-h124gg-gf-v7n_noskim.root");
   // TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterPreselNoSmear");
	if(gammaID==0)
	{
	dtree = (TTree*)fdin->Get("Tree_Optim_gamma");
	}
	else if(gammaID==1)
	{
	dtree = (TTree*)fdin->Get("Tree_Optim_gamma1");
	}
	else if(gammaID==2)
	{
	dtree = (TTree*)fdin->Get("Tree_Optim_gamma2");
	}
    }      
   else
    {
  TFile *fdin = TFile::Open("/afs/cern.ch/work/z/zhicaiz/public/ECALpro_MC_TreeForRegression/Gun_Pi0_Pt1To15_FlatPU0to50RAW_withHLT_80X_mcRun2_GEN-SIM-RAW_ALL_EcalNtp_ALL_EE_combine_test.root");//("root://eoscms.cern.ch///eos/cms/store/cmst3/user/bendavid/idTreesAug1/hgg-2013Final8TeV_ID_s12-h124gg-gf-v7n_noskim.root");
   // TDirectory *ddir = (TDirectory*)fdin->FindObjectAny("PhotonTreeWriterPreselNoSmear");
   	if(gammaID==0)
	{
	dtree = (TTree*)fdin->Get("Tree_Optim_gamma");
	}
	else if(gammaID==1)
	{
	dtree = (TTree*)fdin->Get("Tree_Optim_gamma1");
	}
	else if(gammaID==2)
	{
	dtree = (TTree*)fdin->Get("Tree_Optim_gamma2");
	}
    } 
  }
  
  //selection cuts for testing
  //TCut selcut = "(STr2_enG1_true/cosh(STr2_Eta_1)>1.0) && (STr2_S4S9_1>0.75)";
  TCut selcut = "(STr2_enG_nocor/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.75) && (STr2_isMerging < 2) && (STr2_DeltaR < 0.03)";
  //TCut selcut = "(STr2_enG_nocor/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.75) && (STr2_isMerging < 2) && (STr2_DeltaR < 0.03) && (abs(STr2_iEtaiX)<60)";
  //TCut selcut = "(STr2_enG_nocor/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.75) && (STr2_isMerging < 2) && (STr2_DeltaR < 0.03) && (abs(STr2_iEtaiX)>60)";
  //TCut selcut = "(STr2_enG_nocor/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.9) && (STr2_S2S9>0.85)&& (STr2_isMerging < 2) && (STr2_DeltaR < 0.03) && (abs(STr2_iEtaiX)<60)";
  //TCut selcut = "(STr2_enG_nocor/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.9) && (STr2_S2S9>0.85)&& (STr2_isMerging < 2) && (STr2_DeltaR < 0.03)";
/*  
TCut selcut;
  if (dobarrel) 
    selcut = "ph.genpt>25. && ph.isbarrel && ph.ispromptgen"; 
  else
    selcut = "ph.genpt>25. && !ph.isbarrel && ph.ispromptgen"; 
 */ 
  TCut selweight = "xsecweight(procidx)*puweight(numPU,procidx)";
  TCut prescale10 = "(Entry$%10==0)";
  TCut prescale10alt = "(Entry$%10==1)";
  TCut prescale25 = "(Entry$%25==0)";
  TCut prescale100 = "(Entry$%100==0)";  
  TCut prescale1000 = "(Entry$%1000==0)";  
  TCut evenevents = "(Entry$%2==0)";
  TCut oddevents = "(Entry$%2==1)";
  TCut prescale100alt = "(Entry$%100==1)";
  TCut prescale1000alt = "(Entry$%1000==1)";
  TCut prescale50alt = "(Entry$%50==1)";
  TCut Events3_4 = "(Entry$%4==3)";
  TCut Events1_4 = "(Entry$%4==1)";
  TCut Events2_4 = "(Entry$%4==2)";
  TCut Events0_4 = "(Entry$%4==0)";

  TCut Events01_4 = "(Entry$%4<2)";
  TCut Events23_4 = "(Entry$%4>1)";

  if (doele) 
    weightvar.SetTitle(prescale100alt*selcut);
  else
    weightvar.SetTitle(selcut);
  
  //make testing dataset
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",dtree,vars,weightvar);   

  if (doele) 
    weightvar.SetTitle(prescale1000alt*selcut);
  else
    weightvar.SetTitle(prescale10alt*selcut);
  //make reduced testing dataset for integration over conditional variables
  RooDataSet *hdatasmall = RooTreeConvert::CreateDataSet("hdatasmall",dtree,vars,weightvar);     
    
  //retrieve full pdf from workspace
  RooAbsPdf *sigpdf = ws->pdf("sigpdf");
  
  //input variable corresponding to sceta
  RooRealVar *scetavar = ws->var("var_1");
  RooRealVar *scphivar = ws->var("var_2");
  
 
  //regressed output functions
  RooAbsReal *sigmeanlim = ws->function("sigmeanlim");
  RooAbsReal *sigwidthlim = ws->function("sigwidthlim");
  RooAbsReal *signlim = ws->function("signlim");
  RooAbsReal *sign2lim = ws->function("sign2lim");

  RooAbsReal *sigalphalim = ws->function("sigalphalim");
  RooAbsReal *sigalpha2lim = ws->function("sigalpha2lim");


  //formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
  RooFormulaVar ecor("ecor","","1./(@0)*@1",RooArgList(*tgtvar,*sigmeanlim));
  RooRealVar *ecorvar = (RooRealVar*)hdata->addColumn(ecor);
  ecorvar->setRange(0.,2.);
  ecorvar->setBins(800);
  
  //formula for raw energy/true energy (1.0/(etrue/eraw))
  RooFormulaVar raw("raw","","1./@0",RooArgList(*tgtvar));
  RooRealVar *rawvar = (RooRealVar*)hdata->addColumn(raw);
  rawvar->setRange(0.,2.);
  rawvar->setBins(800);

  //clone data and add regression outputs for plotting
  RooDataSet *hdataclone = new RooDataSet(*hdata,"hdataclone");
  RooRealVar *meanvar = (RooRealVar*)hdataclone->addColumn(*sigmeanlim);
  RooRealVar *widthvar = (RooRealVar*)hdataclone->addColumn(*sigwidthlim);
  RooRealVar *nvar = (RooRealVar*)hdataclone->addColumn(*signlim);
  RooRealVar *n2var = (RooRealVar*)hdataclone->addColumn(*sign2lim);
 
  RooRealVar *alphavar = (RooRealVar*)hdataclone->addColumn(*sigalphalim);
  RooRealVar *alpha2var = (RooRealVar*)hdataclone->addColumn(*sigalpha2lim);
  
  
  //plot target variable and weighted regression prediction (using numerical integration over reduced testing dataset)
  TCanvas *craw = new TCanvas;
  //RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  RooPlot *plot = tgtvar->frame(0.6,2.0,100);
  hdata->plotOn(plot);
  sigpdf->plotOn(plot,ProjWData(*hdatasmall));
  plot->Draw();
  craw->SaveAs("RawE.pdf");
  craw->SaveAs("RawE.png");
  craw->SetLogy();
  plot->SetMinimum(0.1);
  craw->SaveAs("RawElog.pdf");
  craw->SaveAs("RawElog.png");
  
  //plot distribution of regressed functions over testing dataset
  TCanvas *cmean = new TCanvas;
  RooPlot *plotmean = meanvar->frame(0.8,2.0,100);
  hdataclone->plotOn(plotmean);
  plotmean->Draw();
  cmean->SaveAs("mean.pdf");
  cmean->SaveAs("mean.png");
  
  
  TCanvas *cwidth = new TCanvas;
  RooPlot *plotwidth = widthvar->frame(0.,0.05,100);
  hdataclone->plotOn(plotwidth);
  plotwidth->Draw();
  cwidth->SaveAs("width.pdf");
  cwidth->SaveAs("width.png");
  
  TCanvas *cn = new TCanvas;
  RooPlot *plotn = nvar->frame(0.,111.,200);
  hdataclone->plotOn(plotn);
  plotn->Draw();
  cn->SaveAs("n.pdf");
  cn->SaveAs("n.png");

  TCanvas *cn2 = new TCanvas;
  RooPlot *plotn2 = n2var->frame(0.,111.,100);
  hdataclone->plotOn(plotn2);
  plotn2->Draw();
  cn2->SaveAs("n2.pdf");
  cn2->SaveAs("n2.png");


  TCanvas *calpha = new TCanvas;
  RooPlot *plotalpha = alphavar->frame(0.,5.,200);
  hdataclone->plotOn(plotalpha);
  plotalpha->Draw();
  calpha->SaveAs("alpha.pdf");
  calpha->SaveAs("alpha.png");

  TCanvas *calpha2 = new TCanvas;
  RooPlot *plotalpha2 = alpha2var->frame(0.,5.,200);
  hdataclone->plotOn(plotalpha2);
  plotalpha2->Draw();
  calpha2->SaveAs("alpha2.pdf");
  calpha2->SaveAs("alpha2.png");

 
  TCanvas *ceta = new TCanvas;
  RooPlot *ploteta = scetavar->frame(-2.6,2.6,200);
  hdataclone->plotOn(ploteta);
  ploteta->Draw();      
  ceta->SaveAs("eta.pdf");  
  ceta->SaveAs("eta.png");  
  

  //create histograms for eraw/etrue and ecor/etrue to quantify regression performance
  TH1 *heraw;// = hdata->createHistogram("hraw",*rawvar,Binning(800,0.,2.));
  TH1 *hecor;// = hdata->createHistogram("hecor",*ecorvar);
  if (EEorEB == "EB")
  {
         heraw = hdata->createHistogram("hraw",*rawvar,Binning(800,0.8,1.1));
         hecor = hdata->createHistogram("hecor",*ecorvar, Binning(800,0.8,1.1));
  }
  else
  {
         heraw = hdata->createHistogram("hraw",*rawvar,Binning(200,0.,2.));
         hecor = hdata->createHistogram("hecor",*ecorvar, Binning(200,0.,2.));
  }

  
  
  //heold->SetLineColor(kRed);
  hecor->SetLineColor(kBlue);
  heraw->SetLineColor(kMagenta);
  
  hecor->GetYaxis()->SetRangeUser(1.0,1.3*hecor->GetMaximum());
  heraw->GetYaxis()->SetRangeUser(1.0,1.3*hecor->GetMaximum());

  hecor->GetXaxis()->SetRangeUser(0.0,1.5);
  heraw->GetXaxis()->SetRangeUser(0.0,1.5);
  
/*if(EEorEB == "EE")
{
  heraw->GetYaxis()->SetRangeUser(10.0,200.0);
  hecor->GetYaxis()->SetRangeUser(10.0,200.0);
}
*/ 
 
//heold->GetXaxis()->SetRangeUser(0.6,1.2);
  double effsigma_cor, effsigma_raw, fwhm_cor, fwhm_raw;

  if(EEorEB == "EB")
  {
  TH1 *hecorfine = hdata->createHistogram("hecorfine",*ecorvar,Binning(200,0.,2.));
  effsigma_cor = effSigma(hecorfine);
  fwhm_cor = FWHM(hecorfine);
  TH1 *herawfine = hdata->createHistogram("herawfine",*rawvar,Binning(200,0.,2.));
  effsigma_raw = effSigma(herawfine);
  fwhm_raw = FWHM(herawfine);
  }
  else
  {
  TH1 *hecorfine = hdata->createHistogram("hecorfine",*ecorvar,Binning(200,0.,2.));
  effsigma_cor = effSigma(hecorfine);
  fwhm_cor = FWHM(hecorfine);
  TH1 *herawfine = hdata->createHistogram("herawfine",*rawvar,Binning(200,0.,2.));
  effsigma_raw = effSigma(herawfine);
  fwhm_raw = FWHM(herawfine);
  }


  TCanvas *cresponse = new TCanvas;
  gStyle->SetOptStat(0); 
  gStyle->SetPalette(107);
  hecor->SetTitle("");
  heraw->SetTitle("");
  hecor->Draw("HIST");
  //heold->Draw("HISTSAME");
  heraw->Draw("HISTSAME");

  //show errSigma in the plot
  TLegend *leg = new TLegend(0.1, 0.75, 0.7, 0.9);
  leg->AddEntry(hecor,Form("E_{cor}/E_{true}, #sigma_{eff}=%4.3f, FWHM=%4.3f", effsigma_cor, fwhm_cor),"l");
  leg->AddEntry(heraw,Form("E_{raw}/E_{true}, #sigma_{eff}=%4.3f, FWHM=%4.3f", effsigma_raw, fwhm_raw),"l");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
 // leg->SetTextColor(kRed);
  leg->Draw();

  cresponse->SaveAs("response.pdf");
  cresponse->SaveAs("response.png");
  cresponse->SetLogy();
  cresponse->SaveAs("responselog.pdf");
  cresponse->SaveAs("responselog.png");
 

  // draw CCs vs eta and phi

  TCanvas *c_eta = new TCanvas;
  TH1 *h_eta = hdata->createHistogram("h_eta",*scetavar,Binning(100,-3.2,3.2));
  h_eta->Draw("HIST");
  c_eta->SaveAs("heta.pdf");
  c_eta->SaveAs("heta.png");

  TCanvas *c_phi = new TCanvas;
  TH1 *h_phi = hdata->createHistogram("h_phi",*scphivar,Binning(100,-3.2,3.2));
  h_phi->Draw("HIST");
  c_phi->SaveAs("hphi.pdf");
  c_phi->SaveAs("hphi.png");

  RooRealVar *scetaiXvar = ws->var("var_6");
  RooRealVar *scphiiYvar = ws->var("var_7");
 
   if(EEorEB=="EB")
   {
   scetaiXvar->setRange(-90,90);
   scetaiXvar->setBins(180);
   scphiiYvar->setRange(0,360);
   scphiiYvar->setBins(360);
   }
   else
   {
   scetaiXvar->setRange(0,50);
   scetaiXvar->setBins(50);
   scphiiYvar->setRange(0,50);
   scphiiYvar->setBins(50);
 
   }
   ecorvar->setRange(0.5,1.5);
   ecorvar->setBins(800);
   rawvar->setRange(0.5,1.5);
   rawvar->setBins(800);
  

  TCanvas *c_cor_eta = new TCanvas;
  TH2F *h_CC_eta = hdata->createHistogram(*scetaiXvar, *ecorvar, "","cor_vs_eta");
  if(EEorEB=="EB")
  {
  h_CC_eta->GetXaxis()->SetTitle("i#eta"); 
  }
  else
  {
  h_CC_eta->GetXaxis()->SetTitle("iX");
  }
  h_CC_eta->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_eta->Draw("COLZ");
  c_cor_eta->SaveAs("cor_vs_eta.pdf");
  c_cor_eta->SaveAs("cor_vs_eta.png");

  	
  TCanvas *c_cor_phi = new TCanvas;
  TH2F *h_CC_phi = hdata->createHistogram(*scphiiYvar, *ecorvar, "","cor_vs_phi"); 
  if(EEorEB=="EB")
  {
  h_CC_phi->GetXaxis()->SetTitle("i#phi"); 
  }
  else
  {
  h_CC_phi->GetXaxis()->SetTitle("iY");
  }

  h_CC_phi->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_phi->Draw("COLZ");
  c_cor_phi->SaveAs("cor_vs_phi.pdf");
  c_cor_phi->SaveAs("cor_vs_phi.png");
 
  TCanvas *c_raw_eta = new TCanvas;
  TH2F *h_RC_eta = hdata->createHistogram(*scetaiXvar, *rawvar, "","raw_vs_eta");
  if(EEorEB=="EB")
  {
  h_RC_eta->GetXaxis()->SetTitle("i#eta"); 
  }
  else
  {
  h_RC_eta->GetXaxis()->SetTitle("iX");
  }

  h_RC_eta->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_eta->Draw("COLZ");
  c_raw_eta->SaveAs("raw_vs_eta.pdf");
  c_raw_eta->SaveAs("raw_vs_eta.png");
	
  TCanvas *c_raw_phi = new TCanvas;
  TH2F *h_RC_phi = hdata->createHistogram(*scphiiYvar, *rawvar, "","raw_vs_phi"); 
  if(EEorEB=="EB")
  {
  h_RC_phi->GetXaxis()->SetTitle("i#phi"); 
  }
  else
  {
  h_RC_phi->GetXaxis()->SetTitle("iY");
  }

  h_RC_phi->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_phi->Draw("COLZ");
  c_raw_phi->SaveAs("raw_vs_phi.pdf");
  c_raw_phi->SaveAs("raw_vs_phi.png");


//on2,5,20, etc
if(EEorEB == "EB")
{

  TCanvas *myC_iCrystal_mod = new TCanvas;

  RooRealVar *iEtaOn5var = ws->var("var_8");
  iEtaOn5var->setRange(0,5);
  iEtaOn5var->setBins(5);
  TH2F *h_CC_iEtaOn5 = hdata->createHistogram(*iEtaOn5var, *ecorvar, "","cor_vs_iEtaOn5");
  h_CC_iEtaOn5->GetXaxis()->SetTitle("iEtaOn5"); 
  h_CC_iEtaOn5->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_iEtaOn5->Draw("COLZ");
  myC_iCrystal_mod->SaveAs("cor_vs_iEtaOn5.pdf");
  myC_iCrystal_mod->SaveAs("cor_vs_iEtaOn5.png");
  TH2F *h_RC_iEtaOn5 = hdata->createHistogram(*iEtaOn5var, *rawvar, "","raw_vs_iEtaOn5");
  h_RC_iEtaOn5->GetXaxis()->SetTitle("iEtaOn5"); 
  h_RC_iEtaOn5->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_iEtaOn5->Draw("COLZ");
  myC_iCrystal_mod->SaveAs("raw_vs_iEtaOn5.pdf");
  myC_iCrystal_mod->SaveAs("raw_vs_iEtaOn5.png");

  RooRealVar *iPhiOn2var = ws->var("var_9");
  iPhiOn2var->setRange(0,2);
  iPhiOn2var->setBins(2);
  TH2F *h_CC_iPhiOn2 = hdata->createHistogram(*iPhiOn2var, *ecorvar, "","cor_vs_iPhiOn2");
  h_CC_iPhiOn2->GetXaxis()->SetTitle("iPhiOn2"); 
  h_CC_iPhiOn2->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_iPhiOn2->Draw("COLZ");
  myC_iCrystal_mod->SaveAs("cor_vs_iPhiOn2.pdf");
  myC_iCrystal_mod->SaveAs("cor_vs_iPhiOn2.png");
  TH2F *h_RC_iPhiOn2 = hdata->createHistogram(*iPhiOn2var, *rawvar, "","raw_vs_iPhiOn2");
  h_RC_iPhiOn2->GetXaxis()->SetTitle("iPhiOn2"); 
  h_RC_iPhiOn2->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_iPhiOn2->Draw("COLZ");
  myC_iCrystal_mod->SaveAs("raw_vs_iPhiOn2.pdf");
  myC_iCrystal_mod->SaveAs("raw_vs_iPhiOn2.png");

  RooRealVar *iPhiOn20var = ws->var("var_10");
  iPhiOn20var->setRange(0,20);
  iPhiOn20var->setBins(20);
  TH2F *h_CC_iPhiOn20 = hdata->createHistogram(*iPhiOn20var, *ecorvar, "","cor_vs_iPhiOn20");
  h_CC_iPhiOn20->GetXaxis()->SetTitle("iPhiOn20"); 
  h_CC_iPhiOn20->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_iPhiOn20->Draw("COLZ");
  myC_iCrystal_mod->SaveAs("cor_vs_iPhiOn20.pdf");
  myC_iCrystal_mod->SaveAs("cor_vs_iPhiOn20.png");
  TH2F *h_RC_iPhiOn20 = hdata->createHistogram(*iPhiOn20var, *rawvar, "","raw_vs_iPhiOn20");
  h_RC_iPhiOn20->GetXaxis()->SetTitle("iPhiOn20"); 
  h_RC_iPhiOn20->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_iPhiOn20->Draw("COLZ");
  myC_iCrystal_mod->SaveAs("raw_vs_iPhiOn20.pdf");
  myC_iCrystal_mod->SaveAs("raw_vs_iPhiOn20.png");

  RooRealVar *iEtaOn2520var = ws->var("var_11");
  iEtaOn2520var->setRange(-25,25);
  iEtaOn2520var->setBins(50);
  TH2F *h_CC_iEtaOn2520 = hdata->createHistogram(*iEtaOn2520var, *ecorvar, "","cor_vs_iEtaOn2520");
  h_CC_iEtaOn2520->GetXaxis()->SetTitle("iEtaOn2520"); 
  h_CC_iEtaOn2520->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_iEtaOn2520->Draw("COLZ");
  myC_iCrystal_mod->SaveAs("cor_vs_iEtaOn2520.pdf");
  myC_iCrystal_mod->SaveAs("cor_vs_iEtaOn2520.png");
  TH2F *h_RC_iEtaOn2520 = hdata->createHistogram(*iEtaOn2520var, *rawvar, "","raw_vs_iEtaOn2520");
  h_RC_iEtaOn2520->GetXaxis()->SetTitle("iEtaOn2520"); 
  h_RC_iEtaOn2520->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_iEtaOn2520->Draw("COLZ");
  myC_iCrystal_mod->SaveAs("raw_vs_iEtaOn2520.pdf");
  myC_iCrystal_mod->SaveAs("raw_vs_iEtaOn2520.png");

}
	 

// other variables

  TCanvas *myC_variables = new TCanvas;

  RooRealVar *Nxtalvar = ws->var("var_3");
  Nxtalvar->setRange(0,10);
  Nxtalvar->setBins(10);
  TH2F *h_CC_Nxtal = hdata->createHistogram(*Nxtalvar, *ecorvar, "","cor_vs_Nxtal");
  h_CC_Nxtal->GetXaxis()->SetTitle("Nxtal"); 
  h_CC_Nxtal->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_Nxtal->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_Nxtal.pdf");
  myC_variables->SaveAs("cor_vs_Nxtal.png");
  TH2F *h_RC_Nxtal = hdata->createHistogram(*Nxtalvar, *rawvar, "","raw_vs_Nxtal");
  h_RC_Nxtal->GetXaxis()->SetTitle("Nxtal"); 
  h_RC_Nxtal->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_Nxtal->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_Nxtal.pdf");
  myC_variables->SaveAs("raw_vs_Nxtal.png");
	
  RooRealVar *S4S9var = ws->var("var_4");

  int Nbins_S4S9 = 100;
  double Low_S4S9 = 0.6;
  double High_S4S9 = 1.0; 
  S4S9var->setRange(Low_S4S9,High_S4S9);
  S4S9var->setBins(Nbins_S4S9);
 
  TH2F *h_CC_S4S9 = hdata->createHistogram(*S4S9var, *ecorvar, "","cor_vs_S4S9");
  h_CC_S4S9->GetXaxis()->SetTitle("S4S9"); 
  h_CC_S4S9->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_S4S9->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_S4S9.pdf");
  myC_variables->SaveAs("cor_vs_S4S9.png");
  TH2F *h_RC_S4S9 = hdata->createHistogram(*S4S9var, *rawvar, "","raw_vs_S4S9");
  h_RC_S4S9->GetXaxis()->SetTitle("S4S9"); 
  h_RC_S4S9->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_S4S9->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_S4S9.pdf");
  myC_variables->SaveAs("raw_vs_S4S9.png");
	
/* 
  RooRealVar *S1S9var = ws->var("var_5");
  S1S9var->setRange(0.3,1.0);
  S1S9var->setBins(100);
  TH2F *h_CC_S1S9 = hdata->createHistogram(*S1S9var, *ecorvar, "","cor_vs_S1S9");
  h_CC_S1S9->GetXaxis()->SetTitle("S1S9"); 
  h_CC_S1S9->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_S1S9->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_S1S9.pdf");
  TH2F *h_RC_S1S9 = hdata->createHistogram(*S1S9var, *rawvar, "","raw_vs_S1S9");
  h_RC_S1S9->GetXaxis()->SetTitle("S1S9"); 
  h_RC_S1S9->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_S1S9->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_S1S9.pdf");
 */

  RooRealVar *S2S9var = ws->var("var_5");
  int Nbins_S2S9 = 100;
  double Low_S2S9 = 0.5;
  double High_S2S9 = 1.0; 
  S2S9var->setRange(Low_S2S9,High_S2S9);
  S2S9var->setBins(Nbins_S2S9);
  TH2F *h_CC_S2S9 = hdata->createHistogram(*S2S9var, *ecorvar, "","cor_vs_S2S9");
  h_CC_S2S9->GetXaxis()->SetTitle("S2S9"); 
  h_CC_S2S9->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_S2S9->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_S2S9.pdf");
  myC_variables->SaveAs("cor_vs_S2S9.png");
  TH2F *h_RC_S2S9 = hdata->createHistogram(*S2S9var, *rawvar, "","raw_vs_S2S9");
  h_RC_S2S9->GetXaxis()->SetTitle("S2S9"); 
  h_RC_S2S9->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_S2S9->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_S2S9.pdf");
  myC_variables->SaveAs("raw_vs_S2S9.png");

  TH2F *h_S2S9_eta = hdata->createHistogram(*scetaiXvar, *S2S9var, "","S2S9_vs_eta");
  h_S2S9_eta->GetYaxis()->SetTitle("S2S9"); 
  if(EEorEB=="EB")
  {
  h_CC_eta->GetYaxis()->SetTitle("i#eta"); 
  }
  else
  {
  h_CC_eta->GetYaxis()->SetTitle("iX");
  }
  h_S2S9_eta->Draw("COLZ");
  myC_variables->SaveAs("S2S9_vs_eta.pdf");
  myC_variables->SaveAs("S2S9_vs_eta.png");
  
  TH2F *h_S4S9_eta = hdata->createHistogram(*scetaiXvar, *S4S9var, "","S4S9_vs_eta");
  h_S4S9_eta->GetYaxis()->SetTitle("S4S9"); 
  if(EEorEB=="EB")
  {
  h_CC_eta->GetYaxis()->SetTitle("i#eta"); 
  }
  else
  {
  h_CC_eta->GetYaxis()->SetTitle("iX");
  }
  h_S4S9_eta->Draw("COLZ");
  myC_variables->SaveAs("S4S9_vs_eta.pdf");
  myC_variables->SaveAs("S4S9_vs_eta.png");
  
  TH2F *h_S2S9_phi = hdata->createHistogram(*scphiiYvar, *S2S9var, "","S2S9_vs_phi");
  h_S2S9_phi->GetYaxis()->SetTitle("S2S9"); 
  if(EEorEB=="EB")
  {
  h_CC_phi->GetYaxis()->SetTitle("i#phi"); 
  }
  else
  {
  h_CC_phi->GetYaxis()->SetTitle("iY");
  }
  h_S2S9_phi->Draw("COLZ");
  myC_variables->SaveAs("S2S9_vs_phi.pdf");
  myC_variables->SaveAs("S2S9_vs_phi.png");
  
  TH2F *h_S4S9_phi = hdata->createHistogram(*scphiiYvar, *S4S9var, "","S4S9_vs_phi");
  h_S4S9_phi->GetYaxis()->SetTitle("S4S9"); 
  if(EEorEB=="EB")
  {
  h_CC_phi->GetYaxis()->SetTitle("i#phi"); 
  }
  else
  {
  h_CC_phi->GetYaxis()->SetTitle("iY");
  }
  h_S4S9_phi->Draw("COLZ");
  myC_variables->SaveAs("S4S9_vs_phi.pdf");
  myC_variables->SaveAs("S4S9_vs_phi.png");
  

 
/* 
  RooRealVar *DeltaRvar = ws->var("var_6");
  DeltaRvar->setRange(0.0,0.1);
  DeltaRvar->setBins(100);
  TH2F *h_CC_DeltaR = hdata->createHistogram(*DeltaRvar, *ecorvar, "","cor_vs_DeltaR");
  h_CC_DeltaR->GetXaxis()->SetTitle("#Delta R"); 
  h_CC_DeltaR->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_DeltaR->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_DeltaR.pdf");
  myC_variables->SaveAs("cor_vs_DeltaR.png");
  TH2F *h_RC_DeltaR = hdata->createHistogram(*DeltaRvar, *rawvar, "","raw_vs_DeltaR");
  h_RC_DeltaR->GetXaxis()->SetTitle("#Delta R"); 
  h_RC_DeltaR->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_DeltaR->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_DeltaR.pdf");
  myC_variables->SaveAs("raw_vs_DeltaR.png");
*/

  if(EEorEB=="EE")
{

/*  RooRealVar *Es_e1var = ws->var("var_9");
  Es_e1var->setRange(0.0,200.0);
  Es_e1var->setBins(1000);
  TH2F *h_CC_Es_e1 = hdata->createHistogram(*Es_e1var, *ecorvar, "","cor_vs_Es_e1");
  h_CC_Es_e1->GetXaxis()->SetTitle("Es_e1"); 
  h_CC_Es_e1->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_Es_e1->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_Es_e1.pdf");
  myC_variables->SaveAs("cor_vs_Es_e1.png");
  TH2F *h_RC_Es_e1 = hdata->createHistogram(*Es_e1var, *rawvar, "","raw_vs_Es_e1");
  h_RC_Es_e1->GetXaxis()->SetTitle("Es_e1"); 
  h_RC_Es_e1->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_Es_e1->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_Es_e1.pdf");
  myC_variables->SaveAs("raw_vs_Es_e1.png");

  RooRealVar *Es_e2var = ws->var("var_10");
  Es_e2var->setRange(0.0,200.0);
  Es_e2var->setBins(1000);
  TH2F *h_CC_Es_e2 = hdata->createHistogram(*Es_e2var, *ecorvar, "","cor_vs_Es_e2");
  h_CC_Es_e2->GetXaxis()->SetTitle("Es_e2"); 
  h_CC_Es_e2->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_Es_e2->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_Es_e2.pdf");
  myC_variables->SaveAs("cor_vs_Es_e2.png");
  TH2F *h_RC_Es_e2 = hdata->createHistogram(*Es_e2var, *rawvar, "","raw_vs_Es_e2");
  h_RC_Es_e2->GetXaxis()->SetTitle("Es_e2"); 
  h_RC_Es_e2->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_Es_e2->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_Es_e2.pdf");
  myC_variables->SaveAs("raw_vs_Es_e2.png");
*/
}
	
  TProfile *p_CC_eta = h_CC_eta->ProfileX("p_CC_eta",1,-1,"s");
  p_CC_eta->GetYaxis()->SetRangeUser(0.7,1.2);
  if(EEorEB == "EB")
  {
//   p_CC_eta->GetYaxis()->SetRangeUser(0.85,1.0);
//   p_CC_eta->GetXaxis()->SetRangeUser(-1.5,1.5);
  }
  p_CC_eta->GetYaxis()->SetTitle("E_{cor}/E_{true}");
  p_CC_eta->SetTitle("");
  p_CC_eta->Draw();
  myC_variables->SaveAs("profile_cor_vs_eta.pdf"); 
  myC_variables->SaveAs("profile_cor_vs_eta.png"); 
  
  TProfile *p_RC_eta = h_RC_eta->ProfileX("p_RC_eta",1,-1,"s");
  p_RC_eta->GetYaxis()->SetRangeUser(0.7,1.2);
  if(EEorEB=="EB")
  {
//   p_RC_eta->GetYaxis()->SetRangeUser(0.80,0.95);
  // p_RC_eta->GetXaxis()->SetRangeUser(-1.5,1.5);
  }
  p_RC_eta->GetYaxis()->SetTitle("E_{raw}/E_{true}");
  p_RC_eta->SetTitle("");
  p_RC_eta->Draw();
  myC_variables->SaveAs("profile_raw_vs_eta.pdf"); 
  myC_variables->SaveAs("profile_raw_vs_eta.png"); 

  int Nbins_iEta = EEorEB=="EB" ? 180 : 50;
  int nLow_iEta  = EEorEB=="EB" ? -90 : 0;
  int nHigh_iEta = EEorEB=="EB" ? 90 : 50;
  
  TH1F *h1_RC_eta = new TH1F("h1_RC_eta","h1_RC_eta",Nbins_iEta,nLow_iEta,nHigh_iEta);
  for(int i=1;i<=Nbins_iEta;i++)
  {
    h1_RC_eta->SetBinContent(i,p_RC_eta->GetBinError(i)); 
  } 
  h1_RC_eta->GetXaxis()->SetTitle("i#eta");
  h1_RC_eta->GetYaxis()->SetTitle("#sigma_{E_{raw}/E_{true}}");
  h1_RC_eta->SetTitle("");
  h1_RC_eta->Draw();
  myC_variables->SaveAs("sigma_Eraw_Etrue_vs_eta.pdf");
  myC_variables->SaveAs("sigma_Eraw_Etrue_vs_eta.png");
 
  TH1F *h1_CC_eta = new TH1F("h1_CC_eta","h1_CC_eta",Nbins_iEta,nLow_iEta,nHigh_iEta);
  for(int i=1;i<=Nbins_iEta;i++)
  {
    h1_CC_eta->SetBinContent(i,p_CC_eta->GetBinError(i)); 
  } 
  h1_CC_eta->GetXaxis()->SetTitle("i#eta");
  h1_CC_eta->GetYaxis()->SetTitle("#sigma_{E_{cor}/E_{true}}");
  h1_CC_eta->SetTitle("");
  h1_CC_eta->Draw();
  myC_variables->SaveAs("sigma_Ecor_Etrue_vs_eta.pdf");
  myC_variables->SaveAs("sigma_Ecor_Etrue_vs_eta.png");
 
  TProfile *p_CC_phi = h_CC_phi->ProfileX("p_CC_phi",1,-1,"s");
  p_CC_phi->GetYaxis()->SetRangeUser(0.7,1.2);
  if(EEorEB == "EB")
  {
//   p_CC_phi->GetYaxis()->SetRangeUser(0.94,1.00);
  }
  p_CC_phi->GetYaxis()->SetTitle("E_{cor}/E_{true}");
  p_CC_phi->SetTitle("");
  p_CC_phi->Draw();
  myC_variables->SaveAs("profile_cor_vs_phi.pdf"); 
  myC_variables->SaveAs("profile_cor_vs_phi.png"); 
  
  TProfile *p_RC_phi = h_RC_phi->ProfileX("p_RC_phi",1,-1,"s");
  p_RC_phi->GetYaxis()->SetRangeUser(0.7,1.2);
  if(EEorEB=="EB")
  {
 //  p_RC_phi->GetYaxis()->SetRangeUser(0.89,0.95);
  }
  p_RC_phi->GetYaxis()->SetTitle("E_{raw}/E_{true}");
  p_RC_phi->SetTitle("");
  p_RC_phi->Draw();
  myC_variables->SaveAs("profile_raw_vs_phi.pdf"); 
  myC_variables->SaveAs("profile_raw_vs_phi.png"); 

  int Nbins_iPhi = EEorEB=="EB" ? 360 : 50;
  int nLow_iPhi  = EEorEB=="EB" ? 0 : 0;
  int nHigh_iPhi = EEorEB=="EB" ? 360 : 50;
  
  TH1F *h1_RC_phi = new TH1F("h1_RC_phi","h1_RC_phi",Nbins_iPhi,nLow_iPhi,nHigh_iPhi);
  for(int i=1;i<=Nbins_iPhi;i++)
  {
    h1_RC_phi->SetBinContent(i,p_RC_phi->GetBinError(i)); 
  } 
  h1_RC_phi->GetXaxis()->SetTitle("i#phi");
  h1_RC_phi->GetYaxis()->SetTitle("#sigma_{E_{raw}/E_{true}}");
  h1_RC_phi->SetTitle("");
  h1_RC_phi->Draw();
  myC_variables->SaveAs("sigma_Eraw_Etrue_vs_phi.pdf");
  myC_variables->SaveAs("sigma_Eraw_Etrue_vs_phi.png");
 
  TH1F *h1_CC_phi = new TH1F("h1_CC_phi","h1_CC_phi",Nbins_iPhi,nLow_iPhi,nHigh_iPhi);
  for(int i=1;i<=Nbins_iPhi;i++)
  {
    h1_CC_phi->SetBinContent(i,p_CC_phi->GetBinError(i)); 
  } 
  h1_CC_phi->GetXaxis()->SetTitle("i#phi");
  h1_CC_phi->GetYaxis()->SetTitle("#sigma_{E_{cor}/E_{true}}");
  h1_CC_phi->SetTitle("");
  h1_CC_phi->Draw();
  myC_variables->SaveAs("sigma_Ecor_Etrue_vs_phi.pdf");
  myC_variables->SaveAs("sigma_Ecor_Etrue_vs_phi.png");


// FWHM over sigma_eff vs. eta/phi
   
  TH1F *h1_FoverS_RC_phi = new TH1F("h1_FoverS_RC_phi","h1_FoverS_RC_phi",Nbins_iPhi,nLow_iPhi,nHigh_iPhi);
  TH1F *h1_FoverS_CC_phi = new TH1F("h1_FoverS_CC_phi","h1_FoverS_CC_phi",Nbins_iPhi,nLow_iPhi,nHigh_iPhi);
  TH1F *h1_FoverS_RC_eta = new TH1F("h1_FoverS_RC_eta","h1_FoverS_RC_eta",Nbins_iEta,nLow_iEta,nHigh_iEta);
  TH1F *h1_FoverS_CC_eta = new TH1F("h1_FoverS_CC_eta","h1_FoverS_CC_eta",Nbins_iEta,nLow_iEta,nHigh_iEta);
  TH1F *h1_FoverS_CC_S2S9 = new TH1F("h1_FoverS_CC_S2S9","h1_FoverS_CC_S2S9",Nbins_S2S9,Low_S2S9,High_S2S9);
  TH1F *h1_FoverS_RC_S2S9 = new TH1F("h1_FoverS_RC_S2S9","h1_FoverS_RC_S2S9",Nbins_S2S9,Low_S2S9,High_S2S9);
  TH1F *h1_FoverS_CC_S4S9 = new TH1F("h1_FoverS_CC_S4S9","h1_FoverS_CC_S4S9",Nbins_S4S9,Low_S4S9,High_S4S9);
  TH1F *h1_FoverS_RC_S4S9 = new TH1F("h1_FoverS_RC_S4S9","h1_FoverS_RC_S4S9",Nbins_S4S9,Low_S4S9,High_S4S9);

  float FWHMoverSigmaEff = 0.0;  
  TH1F *h_tmp_rawvar = new TH1F("tmp_rawvar","tmp_rawvar",800,0.5,1.5);
  TH1F *h_tmp_corvar = new TH1F("tmp_corvar","tmp_corvar",800,0.5,1.5);

  for(int i=1;i<=Nbins_iPhi;i++)
  {
    float FWHM_tmp = 0.0;
    float effSigma_tmp = 0.0;
    for(int j=1;j<=800;j++) 
    {
	h_tmp_rawvar->SetBinContent(j,h_RC_phi->GetBinContent(i,j));
	h_tmp_corvar->SetBinContent(j,h_CC_phi->GetBinContent(i,j));
    }

    FWHMoverSigmaEff = 0.0;
    FWHM_tmp= FWHM(h_tmp_rawvar);
    effSigma_tmp = effSigma(h_tmp_rawvar);
    if(effSigma_tmp>0.000001)  FWHMoverSigmaEff = FWHM_tmp/effSigma_tmp;
    h1_FoverS_RC_phi->SetBinContent(i, FWHMoverSigmaEff); 

    FWHMoverSigmaEff = 0.0;
    FWHM_tmp= FWHM(h_tmp_corvar);
    effSigma_tmp = effSigma(h_tmp_corvar);
    if(effSigma_tmp>0.000001)  FWHMoverSigmaEff = FWHM_tmp/effSigma_tmp;
    h1_FoverS_CC_phi->SetBinContent(i, FWHMoverSigmaEff); 
  }
  
  h1_FoverS_CC_phi->GetXaxis()->SetTitle("i#phi");
  h1_FoverS_CC_phi->GetYaxis()->SetTitle("FWHM/#sigma_{eff} of E_{cor}/E_{true}");
  h1_FoverS_CC_phi->SetTitle("");
  h1_FoverS_CC_phi->Draw();
  myC_variables->SaveAs("FoverS_Ecor_Etrue_vs_phi.pdf");
  myC_variables->SaveAs("FoverS_Ecor_Etrue_vs_phi.png");

  h1_FoverS_RC_phi->GetXaxis()->SetTitle("i#phi");
  h1_FoverS_RC_phi->GetYaxis()->SetTitle("FWHM/#sigma_{eff} of E_{raw}/E_{true}");
  h1_FoverS_RC_phi->SetTitle("");
  h1_FoverS_RC_phi->Draw();
  myC_variables->SaveAs("FoverS_Eraw_Etrue_vs_phi.pdf");
  myC_variables->SaveAs("FoverS_Eraw_Etrue_vs_phi.png");


  for(int i=1;i<=Nbins_iEta;i++)
  {
    float FWHM_tmp = 0.0;
    float effSigma_tmp = 0.0;
    for(int j=1;j<=800;j++) 
    {
	h_tmp_rawvar->SetBinContent(j,h_RC_eta->GetBinContent(i,j));
	h_tmp_corvar->SetBinContent(j,h_CC_eta->GetBinContent(i,j));
    }

    FWHMoverSigmaEff = 0.0;
    FWHM_tmp= FWHM(h_tmp_rawvar);
    effSigma_tmp = effSigma(h_tmp_rawvar);
    if(effSigma_tmp>0.000001)  FWHMoverSigmaEff = FWHM_tmp/effSigma_tmp;
    h1_FoverS_RC_eta->SetBinContent(i, FWHMoverSigmaEff); 

    FWHMoverSigmaEff = 0.0;
    FWHM_tmp= FWHM(h_tmp_corvar);
    effSigma_tmp = effSigma(h_tmp_corvar);
    if(effSigma_tmp>0.000001)  FWHMoverSigmaEff = FWHM_tmp/effSigma_tmp;
    h1_FoverS_CC_eta->SetBinContent(i, FWHMoverSigmaEff); 
  }
  
  h1_FoverS_CC_eta->GetXaxis()->SetTitle("i#eta");
  h1_FoverS_CC_eta->GetYaxis()->SetTitle("FWHM/#sigma_{eff} of E_{cor}/E_{true}");
  h1_FoverS_CC_eta->SetTitle("");
  h1_FoverS_CC_eta->Draw();
  myC_variables->SaveAs("FoverS_Ecor_Etrue_vs_eta.pdf");
  myC_variables->SaveAs("FoverS_Ecor_Etrue_vs_eta.png");

  h1_FoverS_RC_eta->GetXaxis()->SetTitle("i#eta");
  h1_FoverS_RC_eta->GetYaxis()->SetTitle("FWHM/#sigma_{eff} of E_{raw}/E_{true}");
  h1_FoverS_RC_eta->SetTitle("");
  h1_FoverS_RC_eta->Draw();
  myC_variables->SaveAs("FoverS_Eraw_Etrue_vs_eta.pdf");
  myC_variables->SaveAs("FoverS_Eraw_Etrue_vs_eta.png");


  for(int i=1;i<=Nbins_S2S9;i++)
  {
    float FWHM_tmp = 0.0;
    float effSigma_tmp = 0.0;
    for(int j=1;j<=800;j++) 
    {
	h_tmp_rawvar->SetBinContent(j,h_RC_S2S9->GetBinContent(i,j));
	h_tmp_corvar->SetBinContent(j,h_CC_S2S9->GetBinContent(i,j));
    }

    FWHMoverSigmaEff = 0.0;
    FWHM_tmp= FWHM(h_tmp_rawvar);
    effSigma_tmp = effSigma(h_tmp_rawvar);
    if(effSigma_tmp>0.000001)  FWHMoverSigmaEff = FWHM_tmp/effSigma_tmp;
    h1_FoverS_RC_S2S9->SetBinContent(i, FWHMoverSigmaEff); 

    FWHMoverSigmaEff = 0.0;
    FWHM_tmp= FWHM(h_tmp_corvar);
    effSigma_tmp = effSigma(h_tmp_corvar);
    if(effSigma_tmp>0.000001)  FWHMoverSigmaEff = FWHM_tmp/effSigma_tmp;
    h1_FoverS_CC_S2S9->SetBinContent(i, FWHMoverSigmaEff); 
  }
  
  h1_FoverS_CC_S2S9->GetXaxis()->SetTitle("S2S9");
  h1_FoverS_CC_S2S9->GetYaxis()->SetTitle("FWHM/#sigma_{eff} of E_{cor}/E_{true}");
  h1_FoverS_CC_S2S9->GetYaxis()->SetRangeUser(0.0,1.0);
  h1_FoverS_CC_S2S9->SetTitle("");
  h1_FoverS_CC_S2S9->Draw();
  myC_variables->SaveAs("FoverS_Ecor_Etrue_vs_S2S9.pdf");
  myC_variables->SaveAs("FoverS_Ecor_Etrue_vs_S2S9.png");

  h1_FoverS_RC_S2S9->GetXaxis()->SetTitle("S2S9");
  h1_FoverS_RC_S2S9->GetYaxis()->SetTitle("FWHM/#sigma_{eff} of E_{raw}/E_{true}");
  h1_FoverS_RC_S2S9->GetYaxis()->SetRangeUser(0.0,2.0);
  h1_FoverS_RC_S2S9->SetTitle("");
  h1_FoverS_RC_S2S9->Draw();
  myC_variables->SaveAs("FoverS_Eraw_Etrue_vs_S2S9.pdf");
  myC_variables->SaveAs("FoverS_Eraw_Etrue_vs_S2S9.png");


  for(int i=1;i<=Nbins_S4S9;i++)
  {
    float FWHM_tmp = 0.0;
    float effSigma_tmp = 0.0;
    for(int j=1;j<=800;j++) 
    {
	h_tmp_rawvar->SetBinContent(j,h_RC_S4S9->GetBinContent(i,j));
	h_tmp_corvar->SetBinContent(j,h_CC_S4S9->GetBinContent(i,j));
    }

    FWHMoverSigmaEff = 0.0;
    FWHM_tmp= FWHM(h_tmp_rawvar);
    effSigma_tmp = effSigma(h_tmp_rawvar);
    if(effSigma_tmp>0.000001)  FWHMoverSigmaEff = FWHM_tmp/effSigma_tmp;
    h1_FoverS_RC_S4S9->SetBinContent(i, FWHMoverSigmaEff); 

    FWHMoverSigmaEff = 0.0;
    FWHM_tmp= FWHM(h_tmp_corvar);
    effSigma_tmp = effSigma(h_tmp_corvar);
    if(effSigma_tmp>0.000001)  FWHMoverSigmaEff = FWHM_tmp/effSigma_tmp;
    h1_FoverS_CC_S4S9->SetBinContent(i, FWHMoverSigmaEff); 
  }
  
  h1_FoverS_CC_S4S9->GetXaxis()->SetTitle("S4S9");
  h1_FoverS_CC_S4S9->GetYaxis()->SetTitle("FWHM/#sigma_{eff} of E_{cor}/E_{true}");
  h1_FoverS_CC_S4S9->GetYaxis()->SetRangeUser(0.0,1.0);
  h1_FoverS_CC_S4S9->SetTitle("");
  h1_FoverS_CC_S4S9->Draw();
  myC_variables->SaveAs("FoverS_Ecor_Etrue_vs_S4S9.pdf");
  myC_variables->SaveAs("FoverS_Ecor_Etrue_vs_S4S9.png");

  h1_FoverS_RC_S4S9->GetXaxis()->SetTitle("S4S9");
  h1_FoverS_RC_S4S9->GetYaxis()->SetTitle("FWHM/#sigma_{eff} of E_{raw}/E_{true}");
  h1_FoverS_RC_S4S9->GetYaxis()->SetRangeUser(0.0,2.0);
  h1_FoverS_RC_S4S9->SetTitle("");
  h1_FoverS_RC_S4S9->Draw();
  myC_variables->SaveAs("FoverS_Eraw_Etrue_vs_S4S9.pdf");
  myC_variables->SaveAs("FoverS_Eraw_Etrue_vs_S4S9.png");




  printf("calc effsigma\n");
  std::cout<<"_"<<EEorEB<<std::endl;
  printf("corrected curve effSigma= %5f, FWHM=%5f \n",effsigma_cor, fwhm_cor);
  printf("raw curve effSigma= %5f FWHM=%5f \n",effsigma_raw, fwhm_raw);

  
/*  new TCanvas;
  RooPlot *ploteold = testvar.frame(0.6,1.2,100);
  hdatasigtest->plotOn(ploteold);
  ploteold->Draw();    
  
  new TCanvas;
  RooPlot *plotecor = ecorvar->frame(0.6,1.2,100);
  hdatasig->plotOn(plotecor);
  plotecor->Draw(); */   
  
  
}
