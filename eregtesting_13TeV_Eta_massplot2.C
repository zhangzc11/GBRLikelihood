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

void eregtesting_13TeV_Eta_massplot2(bool dobarrel=true, bool doele=false) {
  
  cout<<"DEBUG ----0000.1"<<endl;
  //output dir
  TString EEorEB = "EE";
  if(dobarrel)
	{
	EEorEB = "EB";
	}
  TString gammaDir = "bothGammas";

  TString gamma1Dir = "gamma1";
  TString gamma2Dir = "gamma2";

  TString dirname = TString::Format("ereg_test_plots_trainetatesteta/%s_%s",gammaDir.Data(),EEorEB.Data());

  
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
  
  TString infile_gamma1 = TString::Format("../../ereg_ws_Eta/%s/%s",gamma1Dir.Data(),fname.Data());
  TString infile_gamma2 = TString::Format("../../ereg_ws_Eta/%s/%s",gamma2Dir.Data(),fname.Data());
  
  TFile *fws_gamma1 = TFile::Open(infile_gamma1); 
  TFile *fws_gamma2 = TFile::Open(infile_gamma2); 
  RooWorkspace *ws_gamma1 = (RooWorkspace*)fws_gamma1->Get("wereg");
  RooWorkspace *ws_gamma2 = (RooWorkspace*)fws_gamma2->Get("wereg");
  cout<<"DEBUG ----0000.2"<<endl;
 
   
  //read variables from workspace
  RooGBRTargetFlex *meantgt_gamma1 = static_cast<RooGBRTargetFlex*>(ws_gamma1->arg("sigmeant"));  
  RooRealVar *tgtvar_gamma1 = ws_gamma1->var("tgtvar");
  RooRealVar *scetavar_gamma1 = new RooRealVar("Eta_gamma1","STr2_Eta",0.); 
  RooRealVar *scphivar_gamma1 = new RooRealVar("Phi_gamma1","STr2_phi",0.); 
  RooRealVar *scNxtalvar_gamma1 = new RooRealVar("Nxtal_gamma1","STr2_Nxtal",0.); 
  RooRealVar *scptPi0var_gamma1 = new RooRealVar("ptPi0_gamma1","STr2_ptPi0_nocor",0.); 
//  RooRealVar *EOverEOthervar = new RooRealVar("EOverEOther","STr2_EOverEOther",0.); 
  
  RooArgList vars_gamma1;
  vars_gamma1.add(meantgt_gamma1->FuncVars());
  vars_gamma1.add(*tgtvar_gamma1);
  vars_gamma1.add(*scetavar_gamma1);
  vars_gamma1.add(*scphivar_gamma1);
  vars_gamma1.add(*scNxtalvar_gamma1);
  vars_gamma1.add(*scptPi0var_gamma1);
//  vars.add(*EOverEOthervar);
   
  RooGBRTargetFlex *meantgt_gamma2 = static_cast<RooGBRTargetFlex*>(ws_gamma2->arg("sigmeant"));  
  RooRealVar *tgtvar_gamma2 = ws_gamma2->var("tgtvar");
  RooRealVar *scetavar_gamma2 = new RooRealVar("Eta_gamma2","STr2_Eta",0.); 
  RooRealVar *scphivar_gamma2 = new RooRealVar("Phi_gamma2","STr2_phi",0.); 
  RooRealVar *scNxtalvar_gamma2 = new RooRealVar("Nxtal_gamma2","STr2_Nxtal",0.); 
  RooRealVar *scptPi0var_gamma2 = new RooRealVar("ptPi0_gamma2","STr2_ptPi0_nocor",0.); 
//  RooRealVar *EOverEOthervar = new RooRealVar("EOverEOther","STr2_EOverEOther",0.); 
  
  RooArgList vars_gamma2;
  vars_gamma2.add(meantgt_gamma2->FuncVars());
  vars_gamma2.add(*tgtvar_gamma2);
  vars_gamma2.add(*scetavar_gamma2);
  vars_gamma2.add(*scphivar_gamma2);
  vars_gamma2.add(*scNxtalvar_gamma2);
  vars_gamma2.add(*scptPi0var_gamma2);
   
  //read testing dataset from TTree
  RooRealVar weightvar_gamma1("weightvar_gamma1","",1.);
  RooRealVar weightvar_gamma2("weightvar_gamma2","",1.);

  TTree *dtree;
  
  //TFile *fdin = TFile::Open("/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/zhicaiz/Gun_MultiPion_FlatPt-1To15/Gun_FlatPt1to15_MultiPion_withPhotonPtFilter_pythia8/photons_0_half2.root"); 
  //TFile *fdin = TFile::Open("/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/zhicaiz/Gun_MultiEta_FlatPt-1To15/Gun_FlatPt1to15_MultiEta_withPhotonPtFilter_pythia8/photons_22Aug2017_V3_half2.root"); 
  TFile *fdin = TFile::Open("/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/zhicaiz/Gun_MultiEta_FlatPt-1To15/Gun_FlatPt1to15_MultiEtaToGG_withPhotonPtFilter_pythia8/photons_20171008_half2.root");
	dtree = (TTree*)fdin->Get("Tree_Optim_gamma");
 

  cout<<"DEBUG ----0001"<<endl;

 
  //selection cuts for testing
//  TCut selcut = "(STr2_enG1_true/cosh(STr2_Eta_1)>1.0) && (STr2_S4S9_1>0.75)";
  //TCut selcut = "(STr2_enG_nocor/cosh(STr2_Eta)>1.0) && (STr2_S4S9 > 0.75) && (STr2_isMerging < 2)";

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

   
  TCut selcut_gamma1 = "(!STr2_fromPi0)";
  TCut selcut_gamma2 = "(!STr2_fromPi0)";
 
    weightvar_gamma1.SetTitle(selcut_gamma1);
    weightvar_gamma2.SetTitle(selcut_gamma2);
  
  cout<<"DEBUG ----0002"<<endl;
  cout<<"dtree: "<<dtree->GetEntries()<<endl;
  //make testing dataset

  RooDataSet *hdata_gamma1 = RooTreeConvert::CreateDataSet("hdata_gamma1",dtree,vars_gamma1,weightvar_gamma1);   
  RooDataSet *hdata_gamma2 = RooTreeConvert::CreateDataSet("hdata_gamma2",dtree,vars_gamma2,weightvar_gamma2);   
  cout<<"DEBUG ----0003"<<endl;

  //retrieve full pdf from workspace
  RooAbsPdf *sigpdf_gamma1 = ws_gamma1->pdf("sigpdf");
  
  //input variable corresponding to sceta
  //
  RooRealVar *abs_erawvar_gamma1 = ws_gamma1->var("var_0");
  RooRealVar *S4S9var_gamma1 = ws_gamma1->var("var_2");
  
  
  RooRealVar *scetaiXvar_gamma1 = ws_gamma1->var("var_4");
  RooRealVar *scphiiYvar_gamma1 = ws_gamma1->var("var_5");
//  RooRealVar *EOverEOthervar = ws->var("var_10");
  
  //regressed output functions
  RooAbsReal *sigmeanlim_gamma1 = ws_gamma1->function("sigmeanlim");
  RooAbsReal *sigwidthlim_gamma1 = ws_gamma1->function("sigwidthlim");
  RooAbsReal *signlim_gamma1 = ws_gamma1->function("signlim");
  RooAbsReal *sign2lim_gamma1 = ws_gamma1->function("sign2lim");
  //cout<<"Function of sigmeanlim:  "<<sigmeanlim->getTitle()<<endl;

  //formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
  RooFormulaVar ecor_gamma1("ecor_gamma1","","1./(@0)*@1",RooArgList(*tgtvar_gamma1,*sigmeanlim_gamma1));
  RooRealVar *ecorvar_gamma1 = (RooRealVar*)hdata_gamma1->addColumn(ecor_gamma1);
  ecorvar_gamma1->setRange(0.,2.);
  ecorvar_gamma1->setBins(800);
  
  RooFormulaVar abs_ecor_gamma1("abs_ecor_gamma1","","1.*@0*@1",RooArgList(*abs_erawvar_gamma1,*sigmeanlim_gamma1));
  RooRealVar *abs_ecorvar_gamma1 = (RooRealVar*)hdata_gamma1->addColumn(abs_ecor_gamma1);

  //formula for raw energy/true energy (1.0/(etrue/eraw))
  RooFormulaVar raw_gamma1("raw_gamma1","","1./@0",RooArgList(*tgtvar_gamma1));
  RooRealVar *rawvar_gamma1 = (RooRealVar*)hdata_gamma1->addColumn(raw_gamma1);
  rawvar_gamma1->setRange(0.,2.);
  rawvar_gamma1->setBins(800);

   //retrieve full pdf from workspace
  RooAbsPdf *sigpdf_gamma2 = ws_gamma2->pdf("sigpdf");
  
  //input variable corresponding to sceta
  //
  RooRealVar *abs_erawvar_gamma2 = ws_gamma2->var("var_0");
  RooRealVar *S4S9var_gamma2 = ws_gamma2->var("var_2");
  
  
  RooRealVar *scetaiXvar_gamma2 = ws_gamma2->var("var_4");
  RooRealVar *scphiiYvar_gamma2 = ws_gamma2->var("var_5");
//  RooRealVar *EOverEOthervar = ws->var("var_10");
  
  //regressed output functions
  RooAbsReal *sigmeanlim_gamma2 = ws_gamma2->function("sigmeanlim");
  RooAbsReal *sigwidthlim_gamma2 = ws_gamma2->function("sigwidthlim");
  RooAbsReal *signlim_gamma2 = ws_gamma2->function("signlim");
  RooAbsReal *sign2lim_gamma2 = ws_gamma2->function("sign2lim");
  //cout<<"Function of sigmeanlim:  "<<sigmeanlim->getTitle()<<endl;


  //formula for corrected energy/true energy ( 1.0/(etrue/eraw) * regression mean)
  RooFormulaVar ecor_gamma2("ecor_gamma2","","1./(@0)*@1",RooArgList(*tgtvar_gamma2,*sigmeanlim_gamma2));
  RooRealVar *ecorvar_gamma2 = (RooRealVar*)hdata_gamma2->addColumn(ecor_gamma2);
  ecorvar_gamma2->setRange(0.,2.);
  ecorvar_gamma2->setBins(800);
  
  RooFormulaVar abs_ecor_gamma2("abs_ecor_gamma2","","1.*@0*@1",RooArgList(*abs_erawvar_gamma2,*sigmeanlim_gamma2));
  RooRealVar *abs_ecorvar_gamma2 = (RooRealVar*)hdata_gamma2->addColumn(abs_ecor_gamma2);

  //formula for raw energy/true energy (1.0/(etrue/eraw))
  RooFormulaVar raw_gamma2("raw_gamma2","","1./@0",RooArgList(*tgtvar_gamma2));
  RooRealVar *rawvar_gamma2 = (RooRealVar*)hdata_gamma2->addColumn(raw_gamma2);
  rawvar_gamma2->setRange(0.,2.);
  rawvar_gamma2->setBins(800);


  TH2F *h2_g1g2_raw = new TH2F("h2_g1g2_raw","h2_g1g2_raw",100,0.75,1.15,100,0.75,1.15); 
  TH2F *h2_g1g2_cor = new TH2F("h2_g1g2_cor","h2_g1g2_cor",100,0.8,1.2,100,0.8,1.2); 

//get pi0 mass peak from the data
	TH1F *hraw_pi0_mass = new TH1F("raw_m_pi0", "raw_m_pi0", 800,0.0,1.0);
	TH1F *hcor_pi0_mass = new TH1F("cor_m_pi0", "cor_m_pi0", 800,0.0,1.0);
	const RooArgSet* set_gamma1;
	const RooArgSet* set_gamma2;
    	int entries=hdata_gamma1->numEntries();
	double eta_val[2], ErawOverEtrue[2], EcorOverEtrue[2],ieta_val[2], phi_val[2], eraw_val[2], ecor_val[2], S4S9_val[2], DeltaR_val[2], massraw_pi0, masscor_pi0, EOverEOther_val[2], Nxtal_val[2], ptPi0_val[2];
	cout<<"Total number of photons: "<<entries<<endl;
	bool passCut[2] = {false,false};
 	for(int i=0;i<entries;i++)
	{
		set_gamma1 = hdata_gamma1->get(i);
		S4S9var_gamma1 = (RooRealVar*)set_gamma1->find(S4S9var_gamma1->GetName());
		scetavar_gamma1 = (RooRealVar*)set_gamma1->find(scetavar_gamma1->GetName());
		scptPi0var_gamma1 = (RooRealVar*)set_gamma1->find(scptPi0var_gamma1->GetName());
		scNxtalvar_gamma1 = (RooRealVar*)set_gamma1->find(scNxtalvar_gamma1->GetName());
		scetaiXvar_gamma1 = (RooRealVar*)set_gamma1->find(scetaiXvar_gamma1->GetName());
		scphivar_gamma1 = (RooRealVar*)set_gamma1->find(scphivar_gamma1->GetName());
		abs_erawvar_gamma1 = (RooRealVar*)set_gamma1->find(abs_erawvar_gamma1->GetName());
		abs_ecorvar_gamma1 = (RooRealVar*)set_gamma1->find(abs_ecorvar_gamma1->GetName());
		ecorvar_gamma1 = (RooRealVar*)set_gamma1->find(ecorvar_gamma1->GetName());
		tgtvar_gamma1 = (RooRealVar*)set_gamma1->find(tgtvar_gamma1->GetName());

		set_gamma2 = hdata_gamma2->get(i);
		S4S9var_gamma2 = (RooRealVar*)set_gamma2->find(S4S9var_gamma2->GetName());
		scetavar_gamma2 = (RooRealVar*)set_gamma2->find(scetavar_gamma2->GetName());
		scptPi0var_gamma2 = (RooRealVar*)set_gamma2->find(scptPi0var_gamma2->GetName());
		scNxtalvar_gamma2 = (RooRealVar*)set_gamma2->find(scNxtalvar_gamma2->GetName());
		scetaiXvar_gamma2 = (RooRealVar*)set_gamma2->find(scetaiXvar_gamma2->GetName());
		scphivar_gamma2 = (RooRealVar*)set_gamma2->find(scphivar_gamma2->GetName());
		abs_erawvar_gamma2 = (RooRealVar*)set_gamma2->find(abs_erawvar_gamma2->GetName());
		abs_ecorvar_gamma2 = (RooRealVar*)set_gamma2->find(abs_ecorvar_gamma2->GetName());
		ecorvar_gamma2 = (RooRealVar*)set_gamma2->find(ecorvar_gamma2->GetName());
		tgtvar_gamma2 = (RooRealVar*)set_gamma2->find(tgtvar_gamma2->GetName());

		if(i%2==0)
		{
			if (tgtvar_gamma1->getVal()>0.0001)
			{
				ErawOverEtrue[i%2] = 1./(tgtvar_gamma1->getVal());
			}
			else
			{
				ErawOverEtrue[i%2] = 999.;
			}

			EcorOverEtrue[i%2] = ecorvar_gamma1->getVal();
			S4S9_val[i%2] = S4S9var_gamma1->getVal();
			eta_val[i%2] = scetavar_gamma1->getVal();
			Nxtal_val[i%2] = scNxtalvar_gamma1->getVal();
			ptPi0_val[i%2] = scptPi0var_gamma1->getVal();
			ieta_val[i%2] = scetaiXvar_gamma1->getVal();
			phi_val[i%2] = scphivar_gamma1->getVal();
			eraw_val[i%2] = abs_erawvar_gamma1->getVal();
			ecor_val[i%2] = abs_ecorvar_gamma1->getVal();
		}
		else
		{
			if (tgtvar_gamma2->getVal()>0.0001)
			{
				ErawOverEtrue[i%2] = 1./(tgtvar_gamma2->getVal());
			}
			else
			{
				ErawOverEtrue[i%2] = 999.;
			}

			EcorOverEtrue[i%2] = ecorvar_gamma2->getVal();
			S4S9_val[i%2] = S4S9var_gamma2->getVal();
			eta_val[i%2] = scetavar_gamma2->getVal();
			Nxtal_val[i%2] = scNxtalvar_gamma2->getVal();
			ptPi0_val[i%2] = scptPi0var_gamma2->getVal();
			ieta_val[i%2] = scetaiXvar_gamma2->getVal();
			phi_val[i%2] = scphivar_gamma2->getVal();
			eraw_val[i%2] = abs_erawvar_gamma2->getVal();
			ecor_val[i%2] = abs_ecorvar_gamma2->getVal();
		}

		if((eraw_val[i%2]/cosh(eta_val[i%2])>1.0)&&(S4S9_val[i%2]>0.75) && Nxtal_val[i%2]>6.5 && ptPi0_val[i%2]>2.0)
		{
			passCut[i%2] = true;
			if(dobarrel && abs(eta_val[i%2])>1.479) passCut[i%2] = false;
			if((!dobarrel) && abs(eta_val[i%2])<1.479) passCut[i%2] = false;
		}
		else
		{
			passCut[i%2] = false;
		}
		if((i%2==1)&&(i>0)&&passCut[0]&&passCut[1])
		{
			massraw_pi0 = sqrt(2*eraw_val[0]*eraw_val[1]*(cosh(eta_val[0]-eta_val[1])-cos(phi_val[0]-phi_val[1]))/(cosh(eta_val[0])*cosh(eta_val[1])));
			masscor_pi0 = sqrt(2*ecor_val[0]*ecor_val[1]*(cosh(eta_val[0]-eta_val[1])-cos(phi_val[0]-phi_val[1]))/(cosh(eta_val[0])*cosh(eta_val[1])));
		//	cout<<i<<"  eraw: "<<eraw_val[0]<<"  "<<eraw_val[1]<<"   ecor: "<<ecor_val[0]<<"  "<<ecor_val[1]<<"   eta: "<<eta_val[0]<<"   "<<eta_val[1]<<"   phi: "<<phi_val[0]<<"   "<<phi_val[1]    <<"   massraw: "<<massraw_pi0<<"   masscor: "<<masscor_pi0<<endl;
			hraw_pi0_mass->Fill(massraw_pi0);
			hcor_pi0_mass->Fill(masscor_pi0);
			h2_g1g2_raw->Fill(ErawOverEtrue[0],ErawOverEtrue[1]);
			h2_g1g2_cor->Fill(EcorOverEtrue[0],EcorOverEtrue[1]);

		}	
	}
  double effsigma_mpi0_cor, effsigma_mpi0_raw, fwhm_mpi0_cor, fwhm_mpi0_raw;

  if(EEorEB == "EE")
  {
	hraw_pi0_mass->SetBins(200,0.0,1.0);
	hcor_pi0_mass->SetBins(200,0.0,1.0);
  }

  effsigma_mpi0_cor = effSigma(hcor_pi0_mass);
  fwhm_mpi0_cor = FWHM(hcor_pi0_mass);
  effsigma_mpi0_raw = effSigma(hraw_pi0_mass);
  fwhm_mpi0_raw = FWHM(hraw_pi0_mass);


//draw pi0 mass peak
	TH1F *h_theoretical = new TH1F("theo","theo",1000000,0.0,1.0);
	h_theoretical->SetBinContent(547862,1.05*hcor_pi0_mass->GetMaximum());

  hcor_pi0_mass->SetLineColor(kBlue);
  hraw_pi0_mass->SetLineColor(kMagenta);
  h_theoretical->SetLineColor(kBlack);
  
  hcor_pi0_mass->GetXaxis()->SetTitle("m_{#eta}/GeV");
  hcor_pi0_mass->GetYaxis()->SetRangeUser(1.0,1.4*hcor_pi0_mass->GetMaximum());
  hraw_pi0_mass->GetYaxis()->SetRangeUser(1.0,1.4*hcor_pi0_mass->GetMaximum());

  TCanvas *cpi0mass = new TCanvas;
  gStyle->SetOptStat(0); 
  hcor_pi0_mass->SetTitle("");
  hraw_pi0_mass->SetTitle("");
  hcor_pi0_mass->Draw("HIST");
  hraw_pi0_mass->Draw("HISTSAME");
  h_theoretical->Draw("HISTSAME");

  //show errSigma in the plot
  TLegend *leg_mpi0 = new TLegend(0.1, 0.74, 0.7, 0.9);
  leg_mpi0->AddEntry(hcor_pi0_mass,Form("m_{#pi_{0}, cor}, #sigma_{eff}=%4.3f, FWHM=%4.3f", effsigma_mpi0_cor, fwhm_mpi0_cor),"l");
  leg_mpi0->AddEntry(hraw_pi0_mass,Form("m_{#pi_{0}, raw}, #sigma_{eff}=%4.3f, FWHM=%4.3f", effsigma_mpi0_raw, fwhm_mpi0_raw),"l");
  leg_mpi0->AddEntry(h_theoretical,"m_{#eta} from PDG (0.548 GeV)","l");

  leg_mpi0->SetFillStyle(0);
  leg_mpi0->SetBorderSize(0);
  leg_mpi0->Draw();

  cpi0mass->SaveAs("mpi0.pdf");
  cpi0mass->SaveAs("mpi0.png");
  cpi0mass->SetLogy();
  cpi0mass->SaveAs("mpi0log.pdf");
  cpi0mass->SaveAs("mpi0log.png");

  TCanvas *myC_EoverEtrue = new TCanvas;
  gStyle->SetOptStat(0);
  h2_g1g2_raw->GetXaxis()->SetTitle("E_{raw}/E_{true} of #gamma_{1}");
  h2_g1g2_raw->GetYaxis()->SetTitle("E_{raw}/E_{true} of #gamma_{2}");
  h2_g1g2_raw->GetYaxis()->SetRangeUser(0.75,1.15);
  h2_g1g2_raw->GetXaxis()->SetRangeUser(0.75,1.15);
  h2_g1g2_raw->SetTitle("");
  h2_g1g2_raw->Draw("COLZ");
  myC_EoverEtrue->SaveAs("ErawOverEtrue_gamma1_vs_gamma2.pdf");
  myC_EoverEtrue->SaveAs("ErawOverEtrue_gamma1_vs_gamma2.png");
  myC_EoverEtrue->SaveAs("ErawOverEtrue_gamma1_vs_gamma2.C");

  h2_g1g2_cor->GetXaxis()->SetTitle("E_{cor}/E_{true} of #gamma_{1}");
  h2_g1g2_cor->GetYaxis()->SetTitle("E_{cor}/E_{true} of #gamma_{2}");
  h2_g1g2_cor->GetYaxis()->SetRangeUser(0.8,1.2);
  h2_g1g2_cor->GetXaxis()->SetRangeUser(0.8,1.2);
  h2_g1g2_cor->SetTitle("");
  h2_g1g2_cor->Draw("COLZ");
  myC_EoverEtrue->SaveAs("EcorOverEtrue_gamma1_vs_gamma2.pdf");
  myC_EoverEtrue->SaveAs("EcorOverEtrue_gamma1_vs_gamma2.png");
  myC_EoverEtrue->SaveAs("EcorOverEtrue_gamma1_vs_gamma2.C");


/* 
//((RooTreeDataStore*)(hdata->store())->tree())->Write(); 
  
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
  
  TCanvas *ceta = new TCanvas;
  RooPlot *ploteta = scetavar->frame(-2.6,2.6,200);
  hdataclone->plotOn(ploteta);
  ploteta->Draw();      
  ceta->SaveAs("eta.eps");  
  

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

  TCanvas *c_eta = new TCanvas;
  TH1 *h_eta = hdata->createHistogram("h_eta",*scetavar,Binning(100,-3.2,3.2));
  h_eta->Draw("HIST");
  c_eta->SaveAs("heta.eps");

  TCanvas *c_phi = new TCanvas;
  TH1 *h_phi = hdata->createHistogram("h_phi",*scphivar,Binning(100,-3.2,3.2));
  h_phi->Draw("HIST");
  c_phi->SaveAs("hphi.eps");


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
   ecorvar->setRange(1.2,2.0);
   ecorvar->setBins(100);
   rawvar->setRange(1.2,2.0);
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

  RooRealVar *Nxtalvar = ws->var("var_3");
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
	
//  RooRealVar *S4S9var = ws->var("var_4");
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
 
  RooRealVar *S2S9var = ws->var("var_6");
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
  
  RooRealVar *DeltaRvar = ws->var("var_7");
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
  RooRealVar *Es_e1var = ws->var("var_10");
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

  RooRealVar *Es_e2var = ws->var("var_11");
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
	
  TProfile *p_CC_eta = h_CC_eta->ProfileX();
  p_CC_eta->GetYaxis()->SetRangeUser(1.2,2.0);
  if(EEorEB == "EB")
  {
   p_CC_eta->GetYaxis()->SetRangeUser(1.2,2.0);
//   p_CC_eta->GetXaxis()->SetRangeUser(-1.5,1.5);
  }
  p_CC_eta->GetYaxis()->SetTitle("E_{cor}/E_{true}");
  p_CC_eta->SetTitle("");
  p_CC_eta->Draw();
  myC_variables->SaveAs("profile_cor_vs_eta.eps"); 
  
  TProfile *p_RC_eta = h_RC_eta->ProfileX();
  p_RC_eta->GetYaxis()->SetRangeUser(1.2,2.0);
  if(EEorEB=="EB")
  {
   p_RC_eta->GetYaxis()->SetRangeUser(1.2,2.0);
  // p_RC_eta->GetXaxis()->SetRangeUser(-1.5,1.5);
  }
  p_RC_eta->GetYaxis()->SetTitle("E_{raw}/E_{true}");
  p_RC_eta->SetTitle("");
  p_RC_eta->Draw();
  myC_variables->SaveAs("profile_raw_vs_eta.eps"); 

  TProfile *p_CC_phi = h_CC_phi->ProfileX();
  p_CC_phi->GetYaxis()->SetRangeUser(0.94,1.06);
  if(EEorEB == "EB")
  {
   p_CC_phi->GetYaxis()->SetRangeUser(0.96,1.02);
  }
  p_CC_phi->GetYaxis()->SetTitle("E_{cor}/E_{true}");
  p_CC_phi->SetTitle("");
  p_CC_phi->Draw();
  myC_variables->SaveAs("profile_cor_vs_phi.eps"); 
  
  TProfile *p_RC_phi = h_RC_phi->ProfileX();
  p_RC_phi->GetYaxis()->SetRangeUser(0.92,1.04);
  if(EEorEB=="EB")
  {
   p_RC_phi->GetYaxis()->SetRangeUser(0.92,0.98);
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
