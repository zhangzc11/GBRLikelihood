#include <iostream>

void treeFormulaTest()
{
	string file_in_Name = "/afs/cern.ch/work/z/zhicaiz/public/ECALpro_MC_TreeForRegression/sum_Pi0Gun_Flat0to50bx25_EB.root";
	string file_out_Name = "/afs/cern.ch/work/z/zhicaiz/public/ECALpro_MC_TreeForRegression/sum_Pi0Gun_Flat0to50bx25_EB_combine.root";
      
      Float_t Op_enG1_rec;
      Float_t Op_enG2_rec;
      //Float_t Op_DeltaRG1G2;
      Float_t Op_Es_e1_1;
      Float_t Op_Es_e1_2;
      Float_t Op_Es_e2_1;
      Float_t Op_Es_e2_2;
      Float_t Op_S4S9_1;
      Float_t Op_S4S9_2;
      Float_t Op_S1S9_1;
      Float_t Op_S1S9_2;
      Float_t Op_S2S9_1;
      Float_t Op_S2S9_2;
      Float_t Op_Eta_1;
      Float_t Op_Eta_2;
      Float_t Op_Phi_1;
      Float_t Op_Phi_2;
      Float_t Op_Time_1;
      Float_t Op_Time_2;
      Float_t Op_DeltaR_1;
      Float_t Op_DeltaR_2;
      Float_t Op_enG1_nocor;
      Float_t Op_enG2_nocor;
      Float_t Op_enG1_true;
      Float_t Op_enG2_true;
      Int_t Op_Nxtal_1;
      Int_t Op_Nxtal_2;
      Int_t Op_iEtaiX_1;
      Int_t Op_iEtaiX_2;
      Int_t Op_iPhiiY_1;
      Int_t Op_iPhiiY_2;
      Int_t Op_iEta_1on5;
      Int_t Op_iEta_2on5;
      Int_t Op_iPhi_1on2;
      Int_t Op_iPhi_2on2;
      Int_t Op_iEta_1on2520;
      Int_t Op_iEta_2on2520;
      Int_t Op_iPhi_1on20;
      Int_t Op_iPhi_2on20;


      Float_t Op_enG_rec;
     // Float_t Op_DeltaRG1G2_common;
      Float_t Op_Es_e1;
      Float_t Op_Es_e2;
      Float_t Op_S4S9;
      Float_t Op_S1S9;
      Float_t Op_S2S9;
      Float_t Op_Eta;
      Float_t Op_Phi;
      Float_t Op_Time;
      Float_t Op_DeltaR;
      Float_t Op_enG_nocor;
      Float_t Op_enG_true;
      Int_t Op_Nxtal;
      Int_t Op_iEtaiX;
      Int_t Op_iPhiiY;
      Int_t Op_iEta_on5;
      Int_t Op_iPhi_on2;
      Int_t Op_iEta_on2520;
      Int_t Op_iPhi_on20;


      double Op_enG12_rec;

	TFile *f_in = new TFile(file_in_Name.c_str(),"READ"); 
	TTree *tree_in = (TTree*)f_in->Get("Tree_Optim");

        tree_in->SetBranchAddress( "STr2_enG1_rec",      &Op_enG1_rec);//,         "STr2_enG1_rec[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_enG1_rec/STr2_enG2_rec",      &Op_enG12_rec);//,         "STr2_enG1_rec[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_enG2_rec",      &Op_enG2_rec);//,         "STr2_enG2_rec[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_enG1_nocor",    &Op_enG1_nocor);//,       "STr2_enG1_nocor[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_enG2_nocor",    &Op_enG2_nocor);//,       "STr2_enG2_nocor[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_enG1_true",     &Op_enG1_true);//,        "STr2_enG1_true[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_enG2_true",     &Op_enG2_true);//,        "STr2_enG2_true[STr2_NPi0_rec]/F");
      //  tree_in->SetBranchAddress( "STr2_DeltaRG1G2",    &Op_DeltaRG1G2);//,       "STr2_DeltaRG1G2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Es_e1_1",       &Op_Es_e1_1);//,          "STr2_Es_e1_1[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Es_e1_2",       &Op_Es_e1_2);//,          "STr2_Es_e1_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Es_e2_1",       &Op_Es_e2_1);//,          "STr2_Es_e2_1[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Es_e2_2",       &Op_Es_e2_2);//,          "STr2_Es_e2_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_S4S9_1",        &Op_S4S9_1);//,           "STr2_S4S9_1[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_S4S9_2",        &Op_S4S9_2);//,           "STr2_S4S9_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_S2S9_1",        &Op_S2S9_1);//,           "STr2_S2S9_1[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_S2S9_2",        &Op_S2S9_2);//,           "STr2_S2S9_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_S1S9_1",        &Op_S1S9_1);//,           "STr2_S1S9_1[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_S1S9_2",        &Op_S1S9_2);//,           "STr2_S1S9_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Eta_1",         &Op_Eta_1);//,            "STr2_Eta_1[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Eta_2",         &Op_Eta_2);//,            "STr2_Eta_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Phi_1",         &Op_Phi_1);//,            "STr2_Phi_1[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Phi_2",         &Op_Phi_2);//,            "STr2_Phi_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Time_1",        &Op_Time_1);//,           "STr2_Time_1[STr2_NPi0_rec]/F");
	tree_in->SetBranchAddress( "STr2_Time_2",        &Op_Time_2);//,           "STr2_Time_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_DeltaR_1",      &Op_DeltaR_1);//,         "STr2_DeltaR_1[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_DeltaR_2",      &Op_DeltaR_2);//,         "STr2_DeltaR_2[STr2_NPi0_rec]/F");
        tree_in->SetBranchAddress( "STr2_Nxtal_1",       &Op_Nxtal_1);//,          "STr2_Nxtal_1[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_Nxtal_2",       &Op_Nxtal_2);//,          "STr2_Nxtal_2[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iEtaiX_1",      &Op_iEtaiX_1);//,         "STr2_iEtaiX_1[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iEtaiX_2",      &Op_iEtaiX_2);//,         "STr2_iEtaiX_2[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iPhiiY_1",      &Op_iPhiiY_1);//,         "STr2_iPhiiY_1[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iPhiiY_2",      &Op_iPhiiY_2);//,         "STr2_iPhiiY_2[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iEta_1on5",     &Op_iEta_1on5);//,        "STr2_iEta_1on5[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iEta_2on5",     &Op_iEta_2on5);//,        "STr2_iEta_2on5[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iPhi_1on2",     &Op_iPhi_1on2);//,        "STr2_iPhi_1on2[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iPhi_2on2",     &Op_iPhi_2on2);//,        "STr2_iPhi_2on2[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iEta_1on2520",  &Op_iEta_1on2520);//,     "STr2_iEta_1on2520[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iEta_2on2520",  &Op_iEta_2on2520);//,     "STr2_iEta_2on2520[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iPhi_1on20",    &Op_iPhi_1on20);//,       "STr2_iPhi_1on20[STr2_NPi0_rec]/I");
        tree_in->SetBranchAddress( "STr2_iPhi_2on20",    &Op_iPhi_2on20);//,       "STr2_iPhi_2on20[STr2_NPi0_rec]/I");

	
	int N_Entries_Pi0 = tree_in->GetEntries();
/*
	TFile *f_out = new TFile(file_out_Name.c_str(),"RECREATE");
	Tree_Optim_gamma = new TTree("Tree_Optim_gamma","Output TTree gamma");
	Tree_Optim_gamma->Branch( "STr2_enG_rec",      &Op_enG_rec,         "STr2_enG_rec/F");
        Tree_Optim_gamma->Branch( "STr2_enG_nocor",    &Op_enG_nocor,       "STr2_enG_nocor/F");
        Tree_Optim_gamma->Branch( "STr2_enG_true",     &Op_enG_true,        "STr2_enG_true/F");
       // Tree_Optim_gamma->Branch( "STr2_DeltaRG1G2",    &Op_DeltaRG1G2_common,       "STr2_DeltaRG1G2/F");
        Tree_Optim_gamma->Branch( "STr2_Es_e1",       &Op_Es_e1,          "STr2_Es_e1/F");
        Tree_Optim_gamma->Branch( "STr2_Es_e2",       &Op_Es_e2,          "STr2_Es_e2/F");
        Tree_Optim_gamma->Branch( "STr2_S4S9",        &Op_S4S9,           "STr2_S4S9/F");
        Tree_Optim_gamma->Branch( "STr2_S2S9",        &Op_S2S9,           "STr2_S2S9/F");
        Tree_Optim_gamma->Branch( "STr2_S1S9",        &Op_S1S9,           "STr2_S1S9/F");
        Tree_Optim_gamma->Branch( "STr2_Eta",         &Op_Eta,            "STr2_Eta/F");
        Tree_Optim_gamma->Branch( "STr2_Phi",         &Op_Phi,            "STr2_Phi/F");
        Tree_Optim_gamma->Branch( "STr2_Time",        &Op_Time,           "STr2_Time/F");
        Tree_Optim_gamma->Branch( "STr2_DeltaR",      &Op_DeltaR,         "STr2_DeltaR/F");
        Tree_Optim_gamma->Branch( "STr2_Nxtal",       &Op_Nxtal,          "STr2_Nxtal/I");
        Tree_Optim_gamma->Branch( "STr2_iEtaiX",      &Op_iEtaiX,         "STr2_iEtaiX/I");
        Tree_Optim_gamma->Branch( "STr2_iPhiiY",      &Op_iPhiiY,         "STr2_iPhiiY/I");
        Tree_Optim_gamma->Branch( "STr2_iEta_on5",     &Op_iEta_on5,        "STr2_iEta_on5/I");
        Tree_Optim_gamma->Branch( "STr2_iPhi_on2",     &Op_iPhi_on2,        "STr2_iPhi_on2/I");
        Tree_Optim_gamma->Branch( "STr2_iEta_on2520",  &Op_iEta_on2520,     "STr2_iEta_on2520/I");
        Tree_Optim_gamma->Branch( "STr2_iPhi_on20",    &Op_iPhi_on20,       "STr2_iPhi_on20/I");

	Tree_Optim_gamma1 = new TTree("Tree_Optim_gamma1","Output TTree gamma1");
	Tree_Optim_gamma1->Branch( "STr2_enG_rec",      &Op_enG_rec,         "STr2_enG_rec/F");
        Tree_Optim_gamma1->Branch( "STr2_enG_nocor",    &Op_enG_nocor,       "STr2_enG_nocor/F");
        Tree_Optim_gamma1->Branch( "STr2_enG_true",     &Op_enG_true,        "STr2_enG_true/F");
       // Tree_Optim_gamma1->Branch( "STr2_DeltaRG1G2",    &Op_DeltaRG1G2_common,       "STr2_DeltaRG1G2/F");
        Tree_Optim_gamma1->Branch( "STr2_Es_e1",       &Op_Es_e1,          "STr2_Es_e1/F");
        Tree_Optim_gamma1->Branch( "STr2_Es_e2",       &Op_Es_e2,          "STr2_Es_e2/F");
        Tree_Optim_gamma1->Branch( "STr2_S4S9",        &Op_S4S9,           "STr2_S4S9/F");
        Tree_Optim_gamma1->Branch( "STr2_S2S9",        &Op_S2S9,           "STr2_S2S9/F");
        Tree_Optim_gamma1->Branch( "STr2_S1S9",        &Op_S1S9,           "STr2_S1S9/F");
        Tree_Optim_gamma1->Branch( "STr2_Eta",         &Op_Eta,            "STr2_Eta/F");
        Tree_Optim_gamma1->Branch( "STr2_Phi",         &Op_Phi,            "STr2_Phi/F");
        Tree_Optim_gamma1->Branch( "STr2_Time",        &Op_Time,           "STr2_Time/F");
        Tree_Optim_gamma1->Branch( "STr2_DeltaR",      &Op_DeltaR,         "STr2_DeltaR/F");
        Tree_Optim_gamma1->Branch( "STr2_Nxtal",       &Op_Nxtal,          "STr2_Nxtal/I");
        Tree_Optim_gamma1->Branch( "STr2_iEtaiX",      &Op_iEtaiX,         "STr2_iEtaiX/I");
        Tree_Optim_gamma1->Branch( "STr2_iPhiiY",      &Op_iPhiiY,         "STr2_iPhiiY/I");
        Tree_Optim_gamma1->Branch( "STr2_iEta_on5",     &Op_iEta_on5,        "STr2_iEta_on5/I");
        Tree_Optim_gamma1->Branch( "STr2_iPhi_on2",     &Op_iPhi_on2,        "STr2_iPhi_on2/I");
        Tree_Optim_gamma1->Branch( "STr2_iEta_on2520",  &Op_iEta_on2520,     "STr2_iEta_on2520/I");
        Tree_Optim_gamma1->Branch( "STr2_iPhi_on20",    &Op_iPhi_on20,       "STr2_iPhi_on20/I");

	Tree_Optim_gamma2 = new TTree("Tree_Optim_gamma2","Output TTree gamma2");
	Tree_Optim_gamma2->Branch( "STr2_enG_rec",      &Op_enG_rec,         "STr2_enG_rec/F");
        Tree_Optim_gamma2->Branch( "STr2_enG_nocor",    &Op_enG_nocor,       "STr2_enG_nocor/F");
        Tree_Optim_gamma2->Branch( "STr2_enG_true",     &Op_enG_true,        "STr2_enG_true/F");
       // Tree_Optim_gamma2->Branch( "STr2_DeltaRG1G2",    &Op_DeltaRG1G2_common,       "STr2_DeltaRG1G2/F");
        Tree_Optim_gamma2->Branch( "STr2_Es_e1",       &Op_Es_e1,          "STr2_Es_e1/F");
        Tree_Optim_gamma2->Branch( "STr2_Es_e2",       &Op_Es_e2,          "STr2_Es_e2/F");
        Tree_Optim_gamma2->Branch( "STr2_S4S9",        &Op_S4S9,           "STr2_S4S9/F");
        Tree_Optim_gamma2->Branch( "STr2_S2S9",        &Op_S2S9,           "STr2_S2S9/F");
        Tree_Optim_gamma2->Branch( "STr2_S1S9",        &Op_S1S9,           "STr2_S1S9/F");
        Tree_Optim_gamma2->Branch( "STr2_Eta",         &Op_Eta,            "STr2_Eta/F");
        Tree_Optim_gamma2->Branch( "STr2_Phi",         &Op_Phi,            "STr2_Phi/F");
        Tree_Optim_gamma2->Branch( "STr2_Time",        &Op_Time,           "STr2_Time/F");
        Tree_Optim_gamma2->Branch( "STr2_DeltaR",      &Op_DeltaR,         "STr2_DeltaR/F");
        Tree_Optim_gamma2->Branch( "STr2_Nxtal",       &Op_Nxtal,          "STr2_Nxtal/I");
        Tree_Optim_gamma2->Branch( "STr2_iEtaiX",      &Op_iEtaiX,         "STr2_iEtaiX/I");
        Tree_Optim_gamma2->Branch( "STr2_iPhiiY",      &Op_iPhiiY,         "STr2_iPhiiY/I");
        Tree_Optim_gamma2->Branch( "STr2_iEta_on5",     &Op_iEta_on5,        "STr2_iEta_on5/I");
        Tree_Optim_gamma2->Branch( "STr2_iPhi_on2",     &Op_iPhi_on2,        "STr2_iPhi_on2/I");
        Tree_Optim_gamma2->Branch( "STr2_iEta_on2520",  &Op_iEta_on2520,     "STr2_iEta_on2520/I");
        Tree_Optim_gamma2->Branch( "STr2_iPhi_on20",    &Op_iPhi_on20,       "STr2_iPhi_on20/I");
*/

	for(int i=0;i<N_Entries_Pi0;i++)
	{
	      	tree_in->GetEntry(i);
		cout<<Op_enG1_rec<<"      "<<Op_enG2_rec<<"   "<<Op_enG12_rec<<endl;


	}

}
