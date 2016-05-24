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

void eregtesting_13TeV_Pi0_lessP3(bool dobarrel=true, bool doele=false,int gammaID=0) {
  
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
  TString dirname = TString::Format("ereg_test_plots_PU_lessP_3/%s_%s",gammaDir.Data(),EEorEB.Data());
  
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
  
  TString infile = TString::Format("../../ereg_ws_PU_lessP_3/%s/%s",gammaDir.Data(),fname.Data());
  
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
    TFile *fdin = TFile::Open("/afs/cern.ch/work/z/zhicaiz/public/ECALpro_MC_TreeForRegression/sum_Pi0Gun_Flat0to50bx25_EB_combine.root");//("root://eoscms.cern.ch///eos/cms/store/cmst3/user/bendavid/idTreesAug1/hgg-2013Final8TeV_ID_s12-h124gg-gf-v7n_noskim.root");
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
  TFile *fdin = TFile::Open("/afs/cern.ch/work/z/zhicaiz/public/ECALpro_MC_TreeForRegression/sum_Pi0Gun_Flat0to50bx25_EE_combine.root");//("root://eoscms.cern.ch///eos/cms/store/cmst3/user/bendavid/idTreesAug1/hgg-2013Final8TeV_ID_s12-h124gg-gf-v7n_noskim.root");
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
  TCut selcut = "(STr2_enG_nocor/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.75) && (STr2_isMerging < 2) && (STr2_DeltaR < 0.05)";
//  TCut selcut = "(STr2_enG_nocor/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.75) && (STr2_S4S9 < 0.999) && (STr2_S2S9 < 0.999)";
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
    weightvar.SetTitle(Events01_4*selcut);
  
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
////  RooRealVar *scetavar = ws->var("var_1");
////  RooRealVar *scphivar = ws->var("var_2");
  
 
  //regressed output functions
  RooAbsReal *sigmeanlim = ws->function("sigmeanlim");
  RooAbsReal *sigwidthlim = ws->function("sigwidthlim");
  RooAbsReal *signlim = ws->function("signlim");
  RooAbsReal *sign2lim = ws->function("sign2lim");

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
  
  
  //plot target variable and weighted regression prediction (using numerical integration over reduced testing dataset)
  TCanvas *craw = new TCanvas;
  //RooPlot *plot = tgtvar->frame(0.6,1.2,100);
  RooPlot *plot = tgtvar->frame(0.6,2.0,100);
  hdata->plotOn(plot);
  sigpdf->plotOn(plot,ProjWData(*hdatasmall));
  plot->Draw();
  craw->SaveAs("RawE.eps");
  craw->SetLogy();
  plot->SetMinimum(0.1);
  craw->SaveAs("RawElog.eps");
  
  //plot distribution of regressed functions over testing dataset
  TCanvas *cmean = new TCanvas;
  RooPlot *plotmean = meanvar->frame(0.8,2.0,100);
  hdataclone->plotOn(plotmean);
  plotmean->Draw();
  cmean->SaveAs("mean.eps");
  
  
  TCanvas *cwidth = new TCanvas;
  RooPlot *plotwidth = widthvar->frame(0.,0.05,100);
  hdataclone->plotOn(plotwidth);
  plotwidth->Draw();
  cwidth->SaveAs("width.eps");
  
  TCanvas *cn = new TCanvas;
  RooPlot *plotn = nvar->frame(0.,111.,200);
  hdataclone->plotOn(plotn);
  plotn->Draw();
  cn->SaveAs("n.eps");

  TCanvas *cn2 = new TCanvas;
  RooPlot *plotn2 = n2var->frame(0.,111.,100);
  hdataclone->plotOn(plotn2);
  plotn2->Draw();
  cn2->SaveAs("n2.eps");
 
////
/* 
  TCanvas *ceta = new TCanvas;
  RooPlot *ploteta = scetavar->frame(-2.6,2.6,200);
  hdataclone->plotOn(ploteta);
  ploteta->Draw();      
  ceta->SaveAs("eta.eps");  
  
*/

  //create histograms for eraw/etrue and ecor/etrue to quantify regression performance
  TH1 *heraw;// = hdata->createHistogram("hraw",*rawvar,Binning(800,0.,2.));
  TH1 *hecor;// = hdata->createHistogram("hecor",*ecorvar);
  if (EEorEB == "EB")
  {
         heraw = hdata->createHistogram("hraw",*rawvar,Binning(800,0.,2.));
         hecor = hdata->createHistogram("hecor",*ecorvar, Binning(800,0.,2.));
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
  TH1 *hecorfine = hdata->createHistogram("hecorfine",*ecorvar,Binning(800,0.,2.));
  effsigma_cor = effSigma(hecorfine);
  fwhm_cor = FWHM(hecorfine);
  TH1 *herawfine = hdata->createHistogram("herawfine",*rawvar,Binning(800,0.,2.));
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

  cresponse->SaveAs("response.eps");
  cresponse->SetLogy();
  cresponse->SaveAs("responselog.eps");
 

  // draw CCs vs eta and phi

//////
/*
  TCanvas *c_eta = new TCanvas;
  TH1 *h_eta = hdata->createHistogram("h_eta",*scetavar,Binning(100,-3.2,3.2));
  h_eta->Draw("HIST");
  c_eta->SaveAs("heta.eps");

  TCanvas *c_phi = new TCanvas;
  TH1 *h_phi = hdata->createHistogram("h_phi",*scphivar,Binning(100,-3.2,3.2));
  h_phi->Draw("HIST");
  c_phi->SaveAs("hphi.eps");
*/
/*
  RooRealVar *scetaiXvar = ws->var("var_1");
  RooRealVar *scphiiYvar = ws->var("var_2");
 
   if(EEorEB=="EB")
   {
   scetaiXvar->setRange(-90,90);
   scetaiXvar->setBins(1800);
   scphiiYvar->setRange(0,360);
   scphiiYvar->setBins(3600);
   }
   else
   {
   scetaiXvar->setRange(0,50);
   scetaiXvar->setBins(500);
   scphiiYvar->setRange(0,50);
   scphiiYvar->setBins(500);
 
   }
   ecorvar->setRange(0.5,1.5);
   ecorvar->setBins(100);
   rawvar->setRange(0.5,1.5);
   rawvar->setBins(100);
  

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
  c_cor_eta->SaveAs("cor_vs_eta.eps");

  	
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
  c_cor_phi->SaveAs("cor_vs_phi.eps");
 
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
  c_raw_eta->SaveAs("raw_vs_eta.eps");
	
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
  c_raw_phi->SaveAs("raw_vs_phi.eps");

// other variables

  TCanvas *myC_variables = new TCanvas;
*/

/////
/*  RooRealVar *Nxtalvar = ws->var("var_3");
  Nxtalvar->setRange(0,10);
  Nxtalvar->setBins(10);
  TH2F *h_CC_Nxtal = hdata->createHistogram(*Nxtalvar, *ecorvar, "","cor_vs_Nxtal");
  h_CC_Nxtal->GetXaxis()->SetTitle("Nxtal"); 
  h_CC_Nxtal->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_Nxtal->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_Nxtal.eps");
  TH2F *h_RC_Nxtal = hdata->createHistogram(*Nxtalvar, *rawvar, "","raw_vs_Nxtal");
  h_RC_Nxtal->GetXaxis()->SetTitle("Nxtal"); 
  h_RC_Nxtal->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_Nxtal->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_Nxtal.eps");
	
  RooRealVar *S4S9var = ws->var("var_4");
  S4S9var->setRange(0.6,1.0);
  S4S9var->setBins(100);
  TH2F *h_CC_S4S9 = hdata->createHistogram(*S4S9var, *ecorvar, "","cor_vs_S4S9");
  h_CC_S4S9->GetXaxis()->SetTitle("S4S9"); 
  h_CC_S4S9->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_S4S9->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_S4S9.eps");
  TH2F *h_RC_S4S9 = hdata->createHistogram(*S4S9var, *rawvar, "","raw_vs_S4S9");
  h_RC_S4S9->GetXaxis()->SetTitle("S4S9"); 
  h_RC_S4S9->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_S4S9->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_S4S9.eps");
*/	
//////
/* 
  RooRealVar *S1S9var = ws->var("var_5");
  S1S9var->setRange(0.3,1.0);
  S1S9var->setBins(100);
  TH2F *h_CC_S1S9 = hdata->createHistogram(*S1S9var, *ecorvar, "","cor_vs_S1S9");
  h_CC_S1S9->GetXaxis()->SetTitle("S1S9"); 
  h_CC_S1S9->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_S1S9->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_S1S9.eps");
  TH2F *h_RC_S1S9 = hdata->createHistogram(*S1S9var, *rawvar, "","raw_vs_S1S9");
  h_RC_S1S9->GetXaxis()->SetTitle("S1S9"); 
  h_RC_S1S9->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_S1S9->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_S1S9.eps");
 */
//////
/*
  RooRealVar *S2S9var = ws->var("var_5");
  S2S9var->setRange(0.5,1.0);
  S2S9var->setBins(100);
  TH2F *h_CC_S2S9 = hdata->createHistogram(*S2S9var, *ecorvar, "","cor_vs_S2S9");
  h_CC_S2S9->GetXaxis()->SetTitle("S2S9"); 
  h_CC_S2S9->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_S2S9->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_S2S9.eps");
  TH2F *h_RC_S2S9 = hdata->createHistogram(*S2S9var, *rawvar, "","raw_vs_S2S9");
  h_RC_S2S9->GetXaxis()->SetTitle("S2S9"); 
  h_RC_S2S9->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_S2S9->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_S2S9.eps");
  
  RooRealVar *DeltaRvar = ws->var("var_6");
  DeltaRvar->setRange(0.0,0.1);
  DeltaRvar->setBins(100);
  TH2F *h_CC_DeltaR = hdata->createHistogram(*DeltaRvar, *ecorvar, "","cor_vs_DeltaR");
  h_CC_DeltaR->GetXaxis()->SetTitle("#Delta R"); 
  h_CC_DeltaR->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_DeltaR->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_DeltaR.eps");
  TH2F *h_RC_DeltaR = hdata->createHistogram(*DeltaRvar, *rawvar, "","raw_vs_DeltaR");
  h_RC_DeltaR->GetXaxis()->SetTitle("#Delta R"); 
  h_RC_DeltaR->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_DeltaR->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_DeltaR.eps");

  if(EEorEB=="EE")
{
  RooRealVar *Es_e1var = ws->var("var_9");
  Es_e1var->setRange(0.0,200.0);
  Es_e1var->setBins(1000);
  TH2F *h_CC_Es_e1 = hdata->createHistogram(*Es_e1var, *ecorvar, "","cor_vs_Es_e1");
  h_CC_Es_e1->GetXaxis()->SetTitle("Es_e1"); 
  h_CC_Es_e1->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_Es_e1->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_Es_e1.eps");
  TH2F *h_RC_Es_e1 = hdata->createHistogram(*Es_e1var, *rawvar, "","raw_vs_Es_e1");
  h_RC_Es_e1->GetXaxis()->SetTitle("Es_e1"); 
  h_RC_Es_e1->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_Es_e1->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_Es_e1.eps");

  RooRealVar *Es_e2var = ws->var("var_10");
  Es_e2var->setRange(0.0,200.0);
  Es_e2var->setBins(1000);
  TH2F *h_CC_Es_e2 = hdata->createHistogram(*Es_e2var, *ecorvar, "","cor_vs_Es_e2");
  h_CC_Es_e2->GetXaxis()->SetTitle("Es_e2"); 
  h_CC_Es_e2->GetYaxis()->SetTitle("E_{cor}/E_{true}"); 
  h_CC_Es_e2->Draw("COLZ");
  myC_variables->SaveAs("cor_vs_Es_e2.eps");
  TH2F *h_RC_Es_e2 = hdata->createHistogram(*Es_e2var, *rawvar, "","raw_vs_Es_e2");
  h_RC_Es_e2->GetXaxis()->SetTitle("Es_e2"); 
  h_RC_Es_e2->GetYaxis()->SetTitle("E_{raw}/E_{true}"); 
  h_RC_Es_e2->Draw("COLZ");
  myC_variables->SaveAs("raw_vs_Es_e2.eps");

}

*/

/*	
  TProfile *p_CC_eta = h_CC_eta->ProfileX();
  p_CC_eta->GetYaxis()->SetRangeUser(0.5,1.5);
  if(EEorEB == "EB")
  {
   p_CC_eta->GetYaxis()->SetRangeUser(0.85,1.0);
//   p_CC_eta->GetXaxis()->SetRangeUser(-1.5,1.5);
  }
  p_CC_eta->GetYaxis()->SetTitle("E_{cor}/E_{true}");
  p_CC_eta->SetTitle("");
  p_CC_eta->Draw();
  myC_variables->SaveAs("profile_cor_vs_eta.eps"); 
  
  TProfile *p_RC_eta = h_RC_eta->ProfileX();
  p_RC_eta->GetYaxis()->SetRangeUser(0.5,1.5);
  if(EEorEB=="EB")
  {
   p_RC_eta->GetYaxis()->SetRangeUser(0.80,0.95);
  // p_RC_eta->GetXaxis()->SetRangeUser(-1.5,1.5);
  }
  p_RC_eta->GetYaxis()->SetTitle("E_{raw}/E_{true}");
  p_RC_eta->SetTitle("");
  p_RC_eta->Draw();
  myC_variables->SaveAs("profile_raw_vs_eta.eps"); 

  TProfile *p_CC_phi = h_CC_phi->ProfileX();
  p_CC_phi->GetYaxis()->SetRangeUser(0.84,1.06);
  if(EEorEB == "EB")
  {
   p_CC_phi->GetYaxis()->SetRangeUser(0.91,1.00);
  }
  p_CC_phi->GetYaxis()->SetTitle("E_{cor}/E_{true}");
  p_CC_phi->SetTitle("");
  p_CC_phi->Draw();
  myC_variables->SaveAs("profile_cor_vs_phi.eps"); 
  
  TProfile *p_RC_phi = h_RC_phi->ProfileX();
  p_RC_phi->GetYaxis()->SetRangeUser(0.82,1.04);
  if(EEorEB=="EB")
  {
   p_RC_phi->GetYaxis()->SetRangeUser(0.86,0.95);
  }
  p_RC_phi->GetYaxis()->SetTitle("E_{raw}/E_{true}");
  p_RC_phi->SetTitle("");
  p_RC_phi->Draw();
  myC_variables->SaveAs("profile_raw_vs_phi.eps"); 



  printf("calc effsigma\n");
  std::cout<<"_"<<EEorEB<<std::endl;
  printf("corrected curve effSigma= %5f, FWHM=%5f \n",effsigma_cor, fwhm_cor);
  printf("raw curve effSigma= %5f FWHM=%5f \n",effsigma_raw, fwhm_raw);
*/
  
/*  new TCanvas;
  RooPlot *ploteold = testvar.frame(0.6,1.2,100);
  hdatasigtest->plotOn(ploteold);
  ploteold->Draw();    
  
  new TCanvas;
  RooPlot *plotecor = ecorvar->frame(0.6,1.2,100);
  hdatasig->plotOn(plotecor);
  plotecor->Draw(); */   
  
  
}
