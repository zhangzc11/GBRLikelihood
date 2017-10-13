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
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
#include "RooLinearVar.h"

//#include "Cintex/Cintex.h"

#include "HybridGBRForestFlex.h"

#ifndef __CINT__
  #include "CondFormats/EgammaObjects/interface/GBRForestD.h"
#endif

using namespace RooFit;

void regcmsformat_13TeV_Pi0_gamma12(bool dobarrel=true, bool doele=false) {
  
  //output dir
   TString dirname = "ereg_cms_Pi0/EB";
   if(!dobarrel)
	dirname = "ereg_cms_Pi0/EE"; 
   gSystem->mkdir(dirname,true);
   gSystem->cd(dirname);    
  
  //read workspace from training
  TString fname_gamma1;
  TString fname_gamma2;
  if (doele && dobarrel) 
  {  fname_gamma1 = "../../ereg_ws_Pi0/gamma1/wereg_ele_eb.root";
    fname_gamma2 = "../../ereg_ws_Pi0/gamma2/wereg_ele_eb.root";
  }else if (doele && !dobarrel) 
  { fname_gamma1 = "../../ereg_ws_Pi0/gamma1/wereg_ele_ee.root";
    fname_gamma2 = "../../ereg_ws_Pi0/gamma2/wereg_ele_ee.root";
  }else if (!doele && dobarrel) 
  {  fname_gamma1 = "../../ereg_ws_Pi0/gamma1/wereg_ph_eb.root";
    fname_gamma2 = "../../ereg_ws_Pi0/gamma2/wereg_ph_eb.root";
  }else if (!doele && !dobarrel) 
  {  fname_gamma1 = "../../ereg_ws_Pi0/gamma1/wereg_ph_ee.root";
    fname_gamma2 = "../../ereg_ws_Pi0/gamma2/wereg_ph_ee.root";
  }

  TString infile_gamma1 = TString::Format("%s",fname_gamma1.Data());
  TString infile_gamma2 = TString::Format("%s",fname_gamma2.Data());
  
  TFile *fws_gamma1 = TFile::Open(infile_gamma1); 
  RooWorkspace *ws_gamma1 = (RooWorkspace*)fws_gamma1->Get("wereg");
  //read variables from workspace
  RooGBRTargetFlex *meantgt_gamma1 = static_cast<RooGBRTargetFlex*>(ws_gamma1->arg("sigmeant"));  
  RooRealVar *tgtvar_gamma1 = ws_gamma1->var("tgtvar");
  
  RooArgList vars_gamma1;
  vars_gamma1.add(meantgt_gamma1->FuncVars());
  vars_gamma1.add(*tgtvar_gamma1);
  
//  ROOT::Cintex::Cintex::Enable();

  TFile *fout_gamma1 = new TFile("wereg_cms_gamma1.root","RECREATE");
  GBRForestD *forestout_gamma1 = new GBRForestD(*meantgt_gamma1->Forest());
  fout_gamma1->WriteObject(forestout_gamma1,"Correction");
  fout_gamma1->Close();
  
   TFile *fws_gamma2 = TFile::Open(infile_gamma2); 
  RooWorkspace *ws_gamma2 = (RooWorkspace*)fws_gamma2->Get("wereg");
  //read variables from workspace
  RooGBRTargetFlex *meantgt_gamma2 = static_cast<RooGBRTargetFlex*>(ws_gamma2->arg("sigmeant"));  
  RooRealVar *tgtvar_gamma2 = ws_gamma2->var("tgtvar");
  
  RooArgList vars_gamma2;
  vars_gamma2.add(meantgt_gamma2->FuncVars());
  vars_gamma2.add(*tgtvar_gamma2);
  
//  ROOT::Cintex::Cintex::Enable();

  TFile *fout_gamma2 = new TFile("wereg_cms_gamma2.root","RECREATE");
  GBRForestD *forestout_gamma2 = new GBRForestD(*meantgt_gamma2->Forest());
  fout_gamma2->WriteObject(forestout_gamma2,"Correction");
  fout_gamma2->Close();
  
}
