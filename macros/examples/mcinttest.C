//Example macro for training a simple BDT classifier using the GBRLikelihood machinery
//J.Bendavid

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
#include "RooGenericPdf.h"
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
#include "MCGBRIntegrator.h"
#include "TF1.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TStyle.h"
  
void mcinttest() {
     
  
  

  MCGBRIntegrator mcint("mcint","",1000);
//   mcint.SetMinCutSignificance(std::numeric_limits<double>::min());
  mcint.SetMinCutSignificance(10.);
//   mcint.SetMinCutSignificance(1e-9);
//     mcint.SetMinCutSignificance(10.);
//   mcint.SetMaxDepth(3);
  mcint.SetMinEvents(10);
//   mcint.SetMinEvents(20);
  mcint.SetShrinkage(0.01);
//   mcint.SetShrinkage(0.01);
//   mcint.SetMaxNodes(100);    
  mcint.TrainForest(1500);
  
  
//   const MCGBRTreeD *gentree = mcint.GenTree();
  const MCGBRForest *gentree = mcint.Forest();
  const MCGBRTreeD *lasttree = &mcint.Forest()->Trees().back();
  
  int nbins = 501;
  double low = 0.;
  double high = 1.;
  double step = (high-low)/double(nbins);
  
  double rlow = 0.;
  double rhigh = 2.0;
  double rstep = (rhigh-rlow)/double(nbins);
  
//   gStyle->SetOptStat(0);
  
  TGraph *hfunc = new TGraph(nbins);
  TGraph *hfunc2 = new TGraph(nbins);
  TGraph *hfunc3 = new TGraph(nbins);
//   TGraph2D *hfunc2d = new TGraph2D(nbins*nbins);
  TH2D *hfunc2d = new TH2D("hfunc2d","",nbins,-5.,5.,nbins,-10.,log(2.));
  TH2D *hfunc2d2 = new TH2D("hfunc2d2","",nbins,-5.,5.,nbins,0.,2.);
  TH2D *hfunc2d3 = new TH2D("hfunc2d3","",nbins,-5.,5.,nbins,0.,2.);
  TH2D *hfunc2d4 = new TH2D("hfunc2d4","",nbins,-5.,5.,nbins,-10.,log(2.));

  std::vector<float> vals(10,0.);
  std::vector<double> valsd(10,0.);
  
  const unsigned int nvars = 4;
  
//   int ipoint2d=0;
 for (int ibin=0; ibin<nbins; ++ibin) {
    float x = low + ibin*step;
//     float x = hfunc2d->GetXaxis()->GetBinCenter(ibin+1);
    
    for (unsigned int ivar=0; ivar<nvars; ++ivar) {
      vals[ivar] = x;
      valsd[ivar] = x;
//       printf("ivar = %i, x = %5f\n",ivar, x);
    }
    
    
    double targetmin;
    double target3;
    double yval = gentree->GetResponse(vals.data(),targetmin,target3);      
//           
    double yvaltrue = mcint.Camel(nvars,valsd.data());
    
//     printf("yval = %5e, yvaltrue = %5e\n",yval,yvaltrue);
    
//     for (int jbin=0; jbin<nbins; ++jbin) {
// //       float r = rlow + jbin*rstep;
//       float r = hfunc2d2->GetYaxis()->GetBinCenter(jbin+1);
//       vals[5] = r;
//       
//       double targetmin;
//       double target3;
//       double yval = gentree->GetResponse(vals.data(),targetmin,target3);      
// //       
// //       double ssb = 1./(1.+exp(-targetmin));
//       double ssb = exp(targetmin);
//       
// //       printf("x = %5f, r = %5f, targetmin = %5f, target = %5f,ssb = %5f\n",x,r,targetmin,yval,ssb);
//       
// //       hfunc2d->Fill(x,r,ssb);
//       hfunc2d2->Fill(x,r,ssb);
//       hfunc2d3->Fill(x,r,yval);
// //       hfunc2d3->Fill(x,r,targetmin);
// //       hfunc2d->SetPoint(ipoint2d,x,r,ssb);
// //       ++ipoint2d;
//     }
//     
//     for (int jbin=0; jbin<nbins; ++jbin) {
// //       float r = rlow + jbin*rstep;
//       double logr = hfunc2d->GetYaxis()->GetBinCenter(jbin+1);
//       float r = exp(logr);
//       vals[5] = r;
//       
//       double targetmin;
//       double targetabs;
//       double yval = gentree->GetResponse(vals.data(),targetmin,targetabs);      
//       
// //       double ssb = 1./(1.+exp(-targetmin));
//       double ssb = exp(targetmin);
//       
// //       double ssb = exp(targetmin)/(1+exp(targetmin));
//       
// //       printf("x = %5f, r = %5f, targetmin = %5f, target = %5f,ssb = %5f\n",x,r,targetmin,yval,ssb);
//       
//       hfunc2d4->Fill(x,logr,ssb);
//       hfunc2d->Fill(x,logr,yval);
// //       hfunc2d3->Fill(x,r,targetmin);
// //       hfunc2d->SetPoint(ipoint2d,x,r,ssb);
// //       ++ipoint2d;
//     }
    
//     for (int ivar=1; ivar<10; ++ivar) {
//       vals[ivar] = gRandom->Uniform(-5.,5.);
//     }
// //     double targetmin;
// //     double targetabs;
//     double yval = gentree->GetResponse(vals.data(),targetmin,targetabs);
//     double lastval = lasttree->GetResponse(vals.data());
//     printf("ibin = %i, x = %5f, yval = %5f\n",ibin,x,yval);
    hfunc->SetPoint(ibin,x,yval);
    hfunc2->SetPoint(ibin,x,yvaltrue);
    
//     hfunc2->SetPoint(ibin,x,exp(targetmin));
//     hfunc3->SetPoint(ibin,x,target3);
// //     hfunc2->SetPoint(ibin,x,exp(-targetmin));
// //     hfunc2->SetPoint(ibin,x,yval-lastval);
//     hfunc2->SetPoint(ibin,x,yval);
//     hfunc3->SetPoint(ibin,x,targetmin);
//     printf("ibin = %i, x = %5f, yval = %5f, targetmin = %5f, sigma = %5f\n",ibin,x,yval,targetmin,exp(-targetmin));
  }
  
  hfunc2->SetLineColor(kRed);
  hfunc2->SetMarkerColor(kRed);
  
  TF1 *gaus = new TF1("gaus","exp(-0.5*x*x)",low,high);
  
  
  TCanvas *cfunc = new TCanvas;
//   gaus->Draw();
  hfunc2->Draw("APL");
  hfunc->Draw("LPSAME");
  
  hfunc2->GetXaxis()->SetTitle("d, distance along multidimensional diagonal (a.u.)");
  hfunc2->GetYaxis()->SetTitle("function value (a.u.)");
  
  cfunc->SaveAs("func.pdf");
  cfunc->SetLogy();
  cfunc->SaveAs("funclog.pdf");
  
//   
/*  new TCanvas;
  gaus->Draw();
  hfunc2->Draw("LPSAME");
  
  new TCanvas;
  gaus->Draw();
  hfunc3->Draw("LPSAME"); */ 
//     
  
  new TCanvas;
  hfunc->Draw("ALP");
//   
//   new TCanvas;
//   hfunc2->Draw("ALP");  
//   
//   new TCanvas;
//   hfunc3->Draw("ALP");    

//   new TCanvas;
//   hfunc2d->Draw("COLZ");
// 
//   new TCanvas;
//   hfunc2d2->Draw("COLZ");
// 
//   new TCanvas;
//   hfunc2d3->Draw("COLZ");  
//  
//   new TCanvas;
//   hfunc2d4->Draw("COLZ");  
  
/*  new TCanvas;
  hfunc3->Draw("ALP"); */   
  
}