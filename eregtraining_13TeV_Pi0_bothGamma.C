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
#include "RooCBExp.h"
#include "RooCBFast.h"
#include "RooGaussianFast.h"

 
using namespace RooFit;
  
double getweight(TFile *file, double xsec) {
 
  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  TH1D *hallevts = (TH1D*)dir->Get("hDAllEvents");
  
  return xsec/hallevts->GetSumOfWeights();
  
}

float xsecweights[50];
float xsecweight(int procidx=0) {
  return xsecweights[procidx];
}

void initweights(TChain *chain, float *xsecs, float lumi) {
 
  TObjArray *files = chain->GetListOfFiles();
  for (int i=0; i<files->GetEntries(); ++i) {    
    TFile *file = TFile::Open(files->At(i)->GetTitle(),"READ");
    
    xsecweights[i] = getweight(file,lumi*xsecs[i]);
    
    file->Close();    
  } 
  
  chain->SetAlias("procidx","This->GetTreeNumber()");
  
}

void eregtraining_13TeV_Pi0_bothGamma(bool dobarrel=true, bool doele=false) {
//  doele=false;
//  dobarrel=true;
 
  //output dir
  //TString dirname = "/data/bendavid/eregexampletest/"; 
  TString dirname = "ereg_ws/";
  gSystem->mkdir(dirname,true);
  gSystem->cd(dirname);  
  
  
	cout<<"DEBUG 001..."<<endl; 
  //build vectors with list of input variables
  std::vector<std::string> *varsf = new std::vector<std::string>;
 ///////////////zzc, adjust variables for pi0 
//  varsf->push_back("ph.scrawe");
//  varsf->push_back("ph.sceta");
//  varsf->push_back("ph.scphi");
//  varsf->push_back("ph.r9");  
//  varsf->push_back("ph.scetawidth");
//  varsf->push_back("ph.scphiwidth");  
//  varsf->push_back("ph.scnclusters");
//  varsf->push_back("ph.hoveretower");
//  varsf->push_back("rho");
//  varsf->push_back("nVtx");  
 
//  varsf->push_back("ph.etaseed-ph.sceta");
//  varsf->push_back("atan2(sin(ph.phiseed-ph.scphi),cos(ph.phiseed-ph.scphi))");
//  varsf->push_back("ph.eseed/ph.scrawe");
  
//  varsf->push_back("ph.e3x3seed/ph.e5x5seed");
//  varsf->push_back("ph.sigietaietaseed");   
//  varsf->push_back("ph.sigiphiphiseed");   
//  varsf->push_back("ph.covietaiphiseed");
//  varsf->push_back("ph.emaxseed/ph.e5x5seed");
//  varsf->push_back("ph.e2ndseed/ph.e5x5seed");
//  varsf->push_back("ph.etopseed/ph.e5x5seed");
//  varsf->push_back("ph.ebottomseed/ph.e5x5seed");
//  varsf->push_back("ph.eleftseed/ph.e5x5seed");
//  varsf->push_back("ph.erightseed/ph.e5x5seed");
//  varsf->push_back("ph.e2x5maxseed/ph.e5x5seed");
//  varsf->push_back("ph.e2x5topseed/ph.e5x5seed");
//  varsf->push_back("ph.e2x5bottomseed/ph.e5x5seed");
//  varsf->push_back("ph.e2x5leftseed/ph.e5x5seed");
//  varsf->push_back("ph.e2x5rightseed/ph.e5x5seed");


///////////////////////////////////////////////////////////////
   	std::vector<std::string> *varseb;
	std::vector<std::string> *varsee;
   	
	//common for EE and EB
    	varsf->push_back("STr2_enG_nocor");// /cosh(STr2_Eta_1)");
    	varsf->push_back("STr2_Nxtal");
    	
    	varsf->push_back("STr2_S4S9");
   	varsf->push_back("STr2_S1S9");
   	varsf->push_back("STr2_S2S9");
    	varsf->push_back("STr2_DeltaR");

        varseb = new std::vector<std::string>(*varsf);
        varsee = new std::vector<std::string>(*varsf);
        //EE
	varsee->push_back("STr2_iEtaiX");
   	varsee->push_back("STr2_iPhiiY");
	varsee->push_back("STr2_Eta");
	varsee->push_back("STr2_Es_e1");
    	varsee->push_back("STr2_Es_e2");
   	//EB
	varseb->push_back("STr2_Eta");
   	varseb->push_back("STr2_Phi");
	varseb->push_back("STr2_iEta_on5");
   	varseb->push_back("STr2_iPhi_on2");
   	varseb->push_back("STr2_iPhi_on20");
   	varseb->push_back("STr2_iEta_on2520");
  
//  varseb->push_back("ph.e5x5seed/ph.eseed");
 
/* 
  varseb->push_back("STr2_iEtaiX_1");//("ph.ietaseed");
  varseb->push_back("STr2_iPhiiY_1");//("ph.iphiseed");
  varseb->push_back("(STr2_iEtaiX_1-1*abs(STr2_iEtaiX_1)/STr2_iEtaiX_1)%5");
  varseb->push_back("(STr2_iPhiiY_1-1)%2");       
  varseb->push_back("(abs(STr2_iEtaiX_1)<=25)*((STr2_iEtaiX_1-1*abs(STr2_iEtaiX_1)/STr2_iEtaiX_1)%25) + (abs(STr2_iEtaiX_1)>25)*((STr2_iEtaiX_1-26*abs(STr2_iEtaiX_1)/STr2_iEtaiX_1)%20)");
  varseb->push_back("(STr2_iPhiiY_1-1)%20"); 
  varseb->push_back("STr2_Eta_1");//("ph.etacryseed");
  varseb->push_back("STr2_Phi_1");//("ph.phicryseed");
*/
//  varsee->push_back("ph.scpse/ph.scrawe");
    
	cout<<"DEBUG 002..."<<endl; 
  //select appropriate input list for barrel or endcap
  std::vector<std::string> *varslist;
  if (dobarrel) varslist = varseb;
  else varslist = varsee;
  
	cout<<"DEBUG 003..."<<endl; 
  //create RooRealVars for each input variable
  RooArgList vars;
  for (unsigned int ivar=0; ivar<varslist->size(); ++ivar) {
    RooRealVar *var = new RooRealVar(TString::Format("var_%i",ivar),varslist->at(ivar).c_str(),0.);
    vars.addOwned(*var);
  }
  
	cout<<"DEBUG 004..."<<endl; 
  //make list of input variable RooRealVars
  RooArgList condvars(vars);
  
  //create RooRealVar for target
  RooRealVar *tgtvar = new RooRealVar("tgtvar","STr2_enG_true/STr2_enG_nocor",1.);//("tgtvar","ph.gene/ph.scrawe",1.);
 // if (!dobarrel) tgtvar->SetTitle("ph.gene/(ph.scrawe + ph.scpse)");  
  
  //add target to full list
  vars.addOwned(*tgtvar);
    
  //RooRealVar for event weight 
  RooRealVar weightvar("weightvar","",1.);

	cout<<"DEBUG 005..."<<endl; 
  //Initialize TChains with event weights if needed
  TString treeloc;
  if (doele) {
    treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPreselInvert/PhotonTreeWriterSingleInvert/hPhotonTreeSingle";
  }
  else {
    treeloc = "RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPresel/PhotonTreeWriterSingle/hPhotonTreeSingle";
  }

  TChain *tree;
  float xsecs[50];

      
  if (doele) {
    tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPreselInvert/PhotonTreeWriterSingleInvert/hPhotonTreeSingle");
    //tree->Add("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
    tree->Add("/data/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-zllm50-v7n_noskim.root");
    
    xsecs[0] = 1.;
    initweights(tree,xsecs,1.);      
    
    xsecweights[0] = 1.0;
    
  }
  else {
//    tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonIDModPresel/PhotonTreeWriterSingle/hPhotonTreeSingle");
//     tree->Add("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-pj20_40-2em-v7n_noskim.root");
//     tree->Add("root://eoscms.cern.ch//eos/cms/store/cmst3/user/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-pj40-2em-v7n_noskim.root");
//    tree->Add("/data/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-pj20_40-2em-v7n_noskim.root");
//    tree->Add("/data/bendavid/regTreesAug1/hgg-2013Final8TeV_reg_s12-pj40-2em-v7n_noskim.root");    
    tree = new TChain("Tree_Optim_gamma");
    if(dobarrel)
    {
    tree->Add("/afs/cern.ch/work/z/zhicaiz/public/ECALpro_MC_TreeForRegression/sum_Pi0Gun_Flat0to50bx25_EB_combine.root");      
    }
    else
     {
    tree->Add("/afs/cern.ch/work/z/zhicaiz/public/ECALpro_MC_TreeForRegression/sum_Pi0Gun_Flat0to50bx25_EE_combine.root");      
     }
	cout<<"DEBUG 006..."<<endl; 
    xsecs[0] = 0.001835*81930.0;
    xsecs[1] = 0.05387*8884.0;    
    //initweights(tree,xsecs,1.);
 
	cout<<"DEBUG 007..."<<endl; 
    
    //double weightscale = xsecweights[1];
    tree->SetAlias("procidx","This->GetTreeNumber()");
    xsecweights[0] = 1.0;///= weightscale;
    //xsecweights[1] /= weightscale;
	cout<<"DEBUG 008..."<<endl; 
  }
  
  //training selection cut
  ///////////////////////////////zzc, photon gen pt?
  TCut selcut = "(STr2_enG_true/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.75)";
/*
  if (dobarrel) {
    selcut = "ph.genpt>0.5 && ph.isbarrel && ph.ispromptgen"; 
  }
  else {
    selcut = "ph.genpt>0.5 && !ph.isbarrel && ph.ispromptgen";     
  }
  */

	cout<<"DEBUG 009..."<<endl; 
  
  TCut selweight = "xsecweight(procidx)";
  TCut prescale10 = "(Entry$%10==0)";
  TCut prescale20 = "(Entry$%20==0)";
  TCut prescale25 = "(Entry$%25==0)";
  TCut prescale50 = "(Entry$%50==0)";
  TCut prescale100 = "(Entry$%100==0)";  
  TCut prescale1000 = "(Entry$%1000==0)";  
  TCut evenevents = "(Entry$%2==0)";
  TCut oddevents = "(Entry$%2==1)";  

  TCut Events3_4 = "((Entry$%4==0)||(Entry$%4==1)||(Entry$%4==2))";//75% of events
  TCut Events4_5 = "((Entry$%5==0)||(Entry$%5==1)||(Entry$%5==2)||(Entry$%5==3))";//80% of events
  
  weightvar.SetTitle(Events4_5*selcut);

  //weightvar title used for per-event weights and selection cuts
/////////////////////////////////zzc, no evt in current tree//////////
/* 
 if (doele) {
    weightvar.SetTitle(prescale100*evenevents*selcut);
  }
  else {
    weightvar.SetTitle(prescale100*selweight*selcut);
  }
*/
	cout<<"DEBUG 010..."<<endl; 
  //create RooDataSet from TChain
  RooDataSet *hdata = RooTreeConvert::CreateDataSet("hdata",tree,vars,weightvar);   
  
  
	cout<<"DEBUG 011..."<<endl; 
  //RooRealVars corresponding to regressed parameters (sigma, mean, left tail parameter, right tail parameter)
  RooRealVar sigwidthtvar("sigwidthtvar","",0.01);
  sigwidthtvar.setConstant(false);
  
  RooRealVar sigmeantvar("sigmeantvar","",1.);
  sigmeantvar.setConstant(false); 
  
  RooRealVar signvar("signvar","",3.);
  signvar.setConstant(false);       
  
  RooRealVar sign2var("sign2var","",3.);
  sign2var.setConstant(false);     

  //define non-parametric functions for each regressed parameter
  RooGBRFunctionFlex *sigwidthtfunc = new RooGBRFunctionFlex("sigwidthtfunc","");
  RooGBRFunctionFlex *sigmeantfunc = new RooGBRFunctionFlex("sigmeantfunc","");
  RooGBRFunctionFlex *signfunc = new RooGBRFunctionFlex("signfunc","");
  RooGBRFunctionFlex *sign2func = new RooGBRFunctionFlex("sign2func","");

  //define mapping of input variables to non-parametric functions (in this case trivial since all 4 functions depend on the same inputs, but this is not a requirement)
  RooGBRTargetFlex *sigwidtht = new RooGBRTargetFlex("sigwidtht","",*sigwidthtfunc,sigwidthtvar,condvars);  
  RooGBRTargetFlex *sigmeant = new RooGBRTargetFlex("sigmeant","",*sigmeantfunc,sigmeantvar,condvars);  
  RooGBRTargetFlex *signt = new RooGBRTargetFlex("signt","",*signfunc,signvar,condvars);  
  RooGBRTargetFlex *sign2t = new RooGBRTargetFlex("sign2t","",*sign2func,sign2var,condvars);  

  //define list of mapped functions to regress
  RooArgList tgts;
  tgts.add(*sigwidtht);
  tgts.add(*sigmeant);
  tgts.add(*signt);
  tgts.add(*sign2t);  
  
  //define transformations corresponding to parameter bounds for non-parametric outputs  
  RooRealConstraint sigwidthlim("sigwidthlim","",*sigwidtht,0.0002,0.5);
  RooRealConstraint sigmeanlim("sigmeanlim","",*sigmeant,0.2,2.0);
  RooRealConstraint signlim("signlim","",*signt,1.01,5000.); 
  RooRealConstraint sign2lim("sign2lim","",*sign2t,1.01,5000.); 

  //define pdf, which depends on transformed outputs (and is intended to be treated as a conditional pdf over the
  //regression inputs in this case)
  //The actual pdf below is a double crystal ball, with crossover points alpha_1 and alpha_2 set constant, but all other
  //parameters regressed
  RooDoubleCBFast sigpdf("sigpdf","",*tgtvar,sigmeanlim,sigwidthlim,RooConst(2.),signlim,RooConst(1.),sign2lim);
  
  //dummy variable
  RooConstVar etermconst("etermconst","",0.);  
   
  //dummy variable
  RooRealVar r("r","",1.);
  r.setConstant();

  //define list of pdfs
  std::vector<RooAbsReal*> vpdf;
  vpdf.push_back(&sigpdf);  

  //define list of training datasets
  std::vector<RooAbsData*> vdata;
  vdata.push_back(hdata);     
  
  //define minimum event weight per tree node
  double minweight = 200;
  std::vector<double> minweights;
  minweights.push_back(minweight);
  
  //temp output file
  TFile *fres = new TFile("fres.root","RECREATE");

	cout<<"DEBUG 012..."<<endl; 
  //run training
  if (1) {
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff","",tgts,etermconst,r,vdata,vpdf);
	cout<<"DEBUG 012.1..."<<endl; 
    bdtpdfdiff.SetMinCutSignificance(2.);
	cout<<"DEBUG 012.2..."<<endl; 
    //bdtpdfdiff.SetPrescaleInit(100);
    bdtpdfdiff.SetShrinkage(0.1);
	cout<<"DEBUG 012.3..."<<endl; 
    bdtpdfdiff.SetMinWeights(minweights);
	cout<<"DEBUG 012.4..."<<endl; 
    bdtpdfdiff.SetMaxNodes(750);
	cout<<"DEBUG 012.5..."<<endl; 
    bdtpdfdiff.TrainForest(1e6);//(1e6);   
	cout<<"DEBUG 012.6..."<<endl; 
  }
     
	cout<<"DEBUG 013..."<<endl; 
  //create workspace and output to file
  RooWorkspace *wereg = new RooWorkspace("wereg");
  wereg->import(sigpdf);
  
  if (doele && dobarrel)
    wereg->writeToFile("wereg_ele_eb.root");    
  else if (doele && !dobarrel) 
    wereg->writeToFile("wereg_ele_ee.root");    
  else if (!doele && dobarrel)
    wereg->writeToFile("wereg_ph_eb.root");    
  else if (!doele && !dobarrel)
    wereg->writeToFile("wereg_ph_ee.root");    
  
  
  return;
  
  
}
