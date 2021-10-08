#define photonTreeProducer_cxx
#include "photonTreeProducer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <vector>
#include <TVector3.h>
#include "TStopwatch.h"

using namespace std;

bool debug = false;

void photonTreeProducer::Loop(string outputfile, int ntotEvent, int nPrintEvent, float eventWeight)
{

  TStopwatch sw;
   sw.Start();
   
  if (fChain == 0) return;

   cout<<"Output file : " << outputfile << endl;
   char *outfilename = const_cast<char*>(outputfile.c_str());
   TFile *fileout = new TFile(outfilename,"RECREATE");
   TTree *tree = new TTree("phoTree","Photon Tree");
   TTree *tree_sig = new TTree("sigPhoTree","Signal Photon Tree");
   TTree *tree_bkg = new TTree("bkgPhoTree","Backgroud Photon Tree");
   
   Long64_t nentries = fChain->GetEntriesFast();

   std::cout << "Total entries: " << nentries << std::endl;
   if (ntotEvent >= 0) nentries = ntotEvent; // number of events you want to process..
   std::cout << "Running on Total entries: " << nentries << std::endl;



   float  evtWeightXS_;
   int    evt_nVtx_;
   int    evt_nGoodVtx_;
   float  evt_vtx_;
   float  evt_vty_;
   float  evt_vtz_;
   float  evt_rho_;

   bool   isTruePho_;
   float  pho_genPt_;
   float  pho_E_;
   float  pho_Pt_;
   float  pho_Eta_;
   float  pho_Phi_;
   float  pho_SCE_;
   float  pho_SCRawE_;
   float  pho_PreShE_;
   float  pho_PreShEbySCRawE_;
   float  pho_SCEta_;
   float  pho_SCPhi_;
   float  pho_SCEtaWidth_;
   float  pho_SCPhiWidth_;
   int    pho_hasPixelSeed_;
   int    pho_EleVeto_;
   float  pho_R9_;
   float  pho_HoverE_;
   float  pho_ESEffSigmaRR_;
   float  pho_SigmaIEtaIEtaFull5x5_;
   float  pho_SigmaIEtaIPhiFull5x5_;
   float  pho_SigmaIPhiIPhiFull5x5_;
   float  pho_E2x2Full5x5_;
   float  pho_E5x5Full5x5_;
   float  pho_R9Full5x5_;
   float  pho_PFChIso_;
   float  pho_PFChPVIso_;
   float  pho_PFPhoIso_;
   float  pho_PFNeuIso_;
   float  pho_PFChWorstVetoIso_;
   float  pho_PFChWorstIso_;
   float  pho_EcalPFClusterIso_;
   float  pho_HcalPFClusterIso_;
   float  pho_SeedTime_;
   float  pho_SeedEnergy_;
   int    pho_SeediEta_;
   int    pho_SeediPhi_;
   vector<float> pho_EnergyMatrix5x5_;
   vector<float> pho_EnergyMatrix7x7_;
   vector<float> pho_EnergyMatrix9x9_;
   vector<float> pho_EnergyMatrix11x11_;
   vector<float> pho_EnergyMatrix15x15_;
   float  mc_Pt_;
   float  mc_PID_;
   float  mc_MomPID_;

   tree->Branch("evt_nVtx",                   &evt_nVtx_);
   tree->Branch("evt_nGoodVtx",               &evt_nGoodVtx_);
   tree->Branch("evt_vtx",                    &evt_vtx_);
   tree->Branch("evt_vty",                    &evt_vty_);
   tree->Branch("evt_vtz",                    &evt_vtz_);
   tree->Branch("evt_rho",                    &evt_rho_);

   tree->Branch("mc_Pt",                    &mc_Pt_);
   tree->Branch("mc_PID",                   &mc_PID_);
   tree->Branch("mc_MomPID",                &mc_MomPID_);
   tree->Branch("isTruePho",                &isTruePho_);
   tree->Branch("evtWeightXS",              &evtWeightXS_);
   tree->Branch("pho_genPt",                &pho_genPt_);
   tree->Branch("pho_E",                    &pho_E_);
   tree->Branch("pho_Pt",                   &pho_Pt_);
   tree->Branch("pho_Eta",                  &pho_Eta_);
   tree->Branch("pho_Phi",                  &pho_Phi_);
   tree->Branch("pho_SCE",                  &pho_SCE_);
   tree->Branch("pho_PreShE",               &pho_PreShE_);
   tree->Branch("pho_SCRawE",               &pho_SCRawE_);
   tree->Branch("pho_PreShEbySCRawE",       &pho_PreShEbySCRawE_);
   tree->Branch("pho_SCEta",                &pho_SCEta_);
   tree->Branch("pho_SCPhi",                &pho_SCPhi_);
   tree->Branch("pho_SCEtaWidth",           &pho_SCEtaWidth_);
   tree->Branch("pho_SCPhiWidth",           &pho_SCPhiWidth_);
   tree->Branch("pho_hasPixelSeed",         &pho_hasPixelSeed_);
   tree->Branch("pho_EleVeto",              &pho_EleVeto_);
   tree->Branch("pho_R9",                   &pho_R9_);
   tree->Branch("pho_HoverE",               &pho_HoverE_);
   tree->Branch("pho_ESEffSigmaRR",         &pho_ESEffSigmaRR_);
   tree->Branch("pho_SigmaIEtaIEtaFull5x5", &pho_SigmaIEtaIEtaFull5x5_);
   tree->Branch("pho_SigmaIEtaIPhiFull5x5", &pho_SigmaIEtaIPhiFull5x5_);
   tree->Branch("pho_SigmaIPhiIPhiFull5x5", &pho_SigmaIPhiIPhiFull5x5_);
   tree->Branch("pho_E2x2Full5x5",          &pho_E2x2Full5x5_);
   tree->Branch("pho_E5x5Full5x5",          &pho_E5x5Full5x5_);
   tree->Branch("pho_R9Full5x5",            &pho_R9Full5x5_);
   tree->Branch("pho_PFChIso",              &pho_PFChIso_);
   tree->Branch("pho_PFChPVIso",            &pho_PFChPVIso_);
   tree->Branch("pho_PFPhoIso",             &pho_PFPhoIso_);
   tree->Branch("pho_PFNeuIso",             &pho_PFNeuIso_);
   tree->Branch("pho_PFChWorstIso",         &pho_PFChWorstIso_);
   tree->Branch("pho_PFChWorstVetoIso",     &pho_PFChWorstVetoIso_);
   tree->Branch("pho_EcalPFClusterIso",     &pho_EcalPFClusterIso_);
   tree->Branch("pho_HcalPFClusterIso",     &pho_HcalPFClusterIso_);
   //tree->Branch("pho_SeedTime",             &pho_SeedTime_);
   //tree->Branch("pho_SeedEnergy",           &pho_SeedEnergy_);
   //tree->Branch("pho_SeediEta",             &pho_SeediEta_);
   //tree->Branch("pho_SeediPhi",             &pho_SeediPhi_);
   tree->Branch("pho_EnergyMatrix5x5",      &pho_EnergyMatrix5x5_);
   //tree->Branch("pho_EnergyMatrix7x7",      &pho_EnergyMatrix7x7_);
   tree->Branch("pho_EnergyMatrix9x9",      &pho_EnergyMatrix9x9_);
   tree->Branch("pho_EnergyMatrix11x11",    &pho_EnergyMatrix11x11_);
   tree->Branch("pho_EnergyMatrix15x15",    &pho_EnergyMatrix15x15_);

   tree_sig->Branch("evt_nVtx",                   &evt_nVtx_);
   tree_sig->Branch("evt_nGoodVtx",               &evt_nGoodVtx_);
   tree_sig->Branch("evt_vtx",                    &evt_vtx_);
   tree_sig->Branch("evt_vty",                    &evt_vty_);
   tree_sig->Branch("evt_vtz",                    &evt_vtz_);
   tree_sig->Branch("evt_rho",                    &evt_rho_);   
   tree_sig->Branch("mc_Pt",                    &mc_Pt_);
   tree_sig->Branch("mc_PID",                   &mc_PID_);
   tree_sig->Branch("mc_MomPID",                &mc_MomPID_);
   tree_sig->Branch("isTruePho",                &isTruePho_);
   tree_sig->Branch("evtWeightXS",              &evtWeightXS_);
   tree_sig->Branch("pho_genPt",                &pho_genPt_);
   tree_sig->Branch("pho_E",                    &pho_E_);
   tree_sig->Branch("pho_Pt",                   &pho_Pt_);
   tree_sig->Branch("pho_Eta",                  &pho_Eta_);
   tree_sig->Branch("pho_Phi",                  &pho_Phi_);
   tree_sig->Branch("pho_SCE",                  &pho_SCE_);
   tree_sig->Branch("pho_PreShE",               &pho_PreShE_);
   tree_sig->Branch("pho_SCRawE",               &pho_SCRawE_);
   tree_sig->Branch("pho_PreShEbySCRawE",       &pho_PreShEbySCRawE_);
   tree_sig->Branch("pho_SCEta",                &pho_SCEta_);
   tree_sig->Branch("pho_SCPhi",                &pho_SCPhi_);
   tree_sig->Branch("pho_SCEtaWidth",           &pho_SCEtaWidth_);
   tree_sig->Branch("pho_SCPhiWidth",           &pho_SCPhiWidth_);
   tree_sig->Branch("pho_hasPixelSeed",         &pho_hasPixelSeed_);
   tree_sig->Branch("pho_EleVeto",              &pho_EleVeto_);
   tree_sig->Branch("pho_R9",                   &pho_R9_);
   tree_sig->Branch("pho_HoverE",               &pho_HoverE_);
   tree_sig->Branch("pho_ESEffSigmaRR",         &pho_ESEffSigmaRR_);
   tree_sig->Branch("pho_SigmaIEtaIEtaFull5x5", &pho_SigmaIEtaIEtaFull5x5_);
   tree_sig->Branch("pho_SigmaIEtaIPhiFull5x5", &pho_SigmaIEtaIPhiFull5x5_);
   tree_sig->Branch("pho_SigmaIPhiIPhiFull5x5", &pho_SigmaIPhiIPhiFull5x5_);
   tree_sig->Branch("pho_E2x2Full5x5",          &pho_E2x2Full5x5_);
   tree_sig->Branch("pho_E5x5Full5x5",          &pho_E5x5Full5x5_);
   tree_sig->Branch("pho_R9Full5x5",            &pho_R9Full5x5_);
   tree_sig->Branch("pho_PFChIso",              &pho_PFChIso_);
   tree_sig->Branch("pho_PFChPVIso",            &pho_PFChPVIso_);
   tree_sig->Branch("pho_PFPhoIso",             &pho_PFPhoIso_);
   tree_sig->Branch("pho_PFNeuIso",             &pho_PFNeuIso_);
   tree_sig->Branch("pho_PFChWorstIso",         &pho_PFChWorstIso_);
   tree_sig->Branch("pho_PFChWorstVetoIso",     &pho_PFChWorstVetoIso_);
   tree_sig->Branch("pho_EcalPFClusterIso",     &pho_EcalPFClusterIso_);
   tree_sig->Branch("pho_HcalPFClusterIso",     &pho_HcalPFClusterIso_);
   //tree_sig->Branch("pho_SeedTime",             &pho_SeedTime_);
   //tree_sig->Branch("pho_SeedEnergy",           &pho_SeedEnergy_);
   //tree_sig->Branch("pho_SeediEta",             &pho_SeediEta_);
   //tree_sig->Branch("pho_SeediPhi",             &pho_SeediPhi_);
   tree_sig->Branch("pho_EnergyMatrix5x5",      &pho_EnergyMatrix5x5_);
   //tree_sig->Branch("pho_EnergyMatrix7x7",      &pho_EnergyMatrix7x7_);
   tree_sig->Branch("pho_EnergyMatrix9x9",      &pho_EnergyMatrix9x9_);
   tree_sig->Branch("pho_EnergyMatrix11x11",    &pho_EnergyMatrix11x11_);
   tree_sig->Branch("pho_EnergyMatrix15x15",    &pho_EnergyMatrix15x15_);

   tree_bkg->Branch("evt_nVtx",                   &evt_nVtx_);
   tree_bkg->Branch("evt_nGoodVtx",               &evt_nGoodVtx_);
   tree_bkg->Branch("evt_vtx",                    &evt_vtx_);
   tree_bkg->Branch("evt_vty",                    &evt_vty_);
   tree_bkg->Branch("evt_vtz",                    &evt_vtz_);
   tree_bkg->Branch("evt_rho",                    &evt_rho_);
   tree_bkg->Branch("mc_Pt",                    &mc_Pt_);
   tree_bkg->Branch("mc_PID",                   &mc_PID_);
   tree_bkg->Branch("mc_MomPID",                &mc_MomPID_);
   tree_bkg->Branch("isTruePho",                &isTruePho_);
   tree_bkg->Branch("evtWeightXS",              &evtWeightXS_);
   tree_bkg->Branch("pho_genPt",                &pho_genPt_);
   tree_bkg->Branch("pho_E",                    &pho_E_);
   tree_bkg->Branch("pho_Pt",                   &pho_Pt_);
   tree_bkg->Branch("pho_Eta",                  &pho_Eta_);
   tree_bkg->Branch("pho_Phi",                  &pho_Phi_);
   tree_bkg->Branch("pho_SCE",                  &pho_SCE_);
   tree_bkg->Branch("pho_PreShE",               &pho_PreShE_);
   tree_bkg->Branch("pho_SCRawE",               &pho_SCRawE_);
   tree_bkg->Branch("pho_PreShEbySCRawE",       &pho_PreShEbySCRawE_);
   tree_bkg->Branch("pho_SCEta",                &pho_SCEta_);
   tree_bkg->Branch("pho_SCPhi",                &pho_SCPhi_);
   tree_bkg->Branch("pho_SCEtaWidth",           &pho_SCEtaWidth_);
   tree_bkg->Branch("pho_SCPhiWidth",           &pho_SCPhiWidth_);
   tree_bkg->Branch("pho_hasPixelSeed",         &pho_hasPixelSeed_);
   tree_bkg->Branch("pho_EleVeto",              &pho_EleVeto_);
   tree_bkg->Branch("pho_R9",                   &pho_R9_);
   tree_bkg->Branch("pho_HoverE",               &pho_HoverE_);
   tree_bkg->Branch("pho_ESEffSigmaRR",         &pho_ESEffSigmaRR_);
   tree_bkg->Branch("pho_SigmaIEtaIEtaFull5x5", &pho_SigmaIEtaIEtaFull5x5_);
   tree_bkg->Branch("pho_SigmaIEtaIPhiFull5x5", &pho_SigmaIEtaIPhiFull5x5_);
   tree_bkg->Branch("pho_SigmaIPhiIPhiFull5x5", &pho_SigmaIPhiIPhiFull5x5_);
   tree_bkg->Branch("pho_E2x2Full5x5",          &pho_E2x2Full5x5_);
   tree_bkg->Branch("pho_E5x5Full5x5",          &pho_E5x5Full5x5_);
   tree_bkg->Branch("pho_R9Full5x5",            &pho_R9Full5x5_);
   tree_bkg->Branch("pho_PFChIso",              &pho_PFChIso_);
   tree_bkg->Branch("pho_PFChPVIso",            &pho_PFChPVIso_);
   tree_bkg->Branch("pho_PFPhoIso",             &pho_PFPhoIso_);
   tree_bkg->Branch("pho_PFNeuIso",             &pho_PFNeuIso_);
   tree_bkg->Branch("pho_PFChWorstIso",         &pho_PFChWorstIso_);
   tree_bkg->Branch("pho_PFChWorstVetoIso",     &pho_PFChWorstVetoIso_);
   tree_bkg->Branch("pho_EcalPFClusterIso",     &pho_EcalPFClusterIso_);
   tree_bkg->Branch("pho_HcalPFClusterIso",     &pho_HcalPFClusterIso_);
   //tree_bkg->Branch("pho_SeedTime",             &pho_SeedTime_);
   //tree_bkg->Branch("pho_SeedEnergy",           &pho_SeedEnergy_);
   //tree_bkg->Branch("pho_SeediEta",             &pho_SeediEta_);
   //tree_bkg->Branch("pho_SeediPhi",             &pho_SeediPhi_);
   tree_bkg->Branch("pho_EnergyMatrix5x5",      &pho_EnergyMatrix5x5_);
   //tree_bkg->Branch("pho_EnergyMatrix7x7",      &pho_EnergyMatrix7x7_);
   tree_bkg->Branch("pho_EnergyMatrix9x9",      &pho_EnergyMatrix9x9_);
   tree_bkg->Branch("pho_EnergyMatrix11x11",    &pho_EnergyMatrix11x11_);
   tree_bkg->Branch("pho_EnergyMatrix15x15",    &pho_EnergyMatrix15x15_);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry % nPrintEvent == 0) std::cout << "  " << jentry  << "  Events Processed... " << std::endl;

      evtWeightXS_ = eventWeight;
      evt_nVtx_ = nVtx;
      evt_nGoodVtx_ = nGoodVtx;      
      evt_vtx_ = vtx;
      evt_vty_ = vty;
      evt_vtz_ = vtz;
      evt_rho_ = rho;

      if (nPho == 0) continue;
      for (int ipho=0; ipho<nPho; ipho++){

	double pEta = (*phoEta)[ipho]; 
	double pPhi = (*phoPhi)[ipho];

	pho_EnergyMatrix5x5_.clear();
	pho_EnergyMatrix7x7_.clear();
	pho_EnergyMatrix9x9_.clear();
	pho_EnergyMatrix11x11_.clear();
	pho_EnergyMatrix15x15_.clear();

	pho_E_ = (*phoE)[ipho];
	pho_Pt_ = (*phoEt)[ipho];
	pho_Eta_ = (*phoEta)[ipho];
	pho_Phi_ = (*phoPhi)[ipho];
	pho_SCE_ = (*phoSCE)[ipho];
	pho_SCRawE_ = (*phoSCRawE)[ipho];
	pho_PreShE_ = (*phoESEn)[ipho];
	pho_PreShEbySCRawE_ = ((*phoESEn)[ipho])/((*phoSCRawE)[ipho]);
	pho_SCEta_ = (*phoSCEta)[ipho];
	pho_SCPhi_ = (*phoSCPhi)[ipho];
	pho_SCEtaWidth_ = (*phoSCEtaWidth)[ipho];
	pho_SCPhiWidth_ = (*phoSCPhiWidth)[ipho];
	pho_hasPixelSeed_ = (*phohasPixelSeed)[ipho];
	pho_EleVeto_ = (*phoEleVeto)[ipho];
	pho_R9_ = (*phoR9)[ipho];
	pho_HoverE_ = (*phoHoverE)[ipho];
	pho_ESEffSigmaRR_ = (*phoESEffSigmaRR)[ipho];
	pho_SigmaIEtaIEtaFull5x5_ = (*phoSigmaIEtaIEtaFull5x5)[ipho];
	pho_SigmaIEtaIPhiFull5x5_ = (*phoSigmaIEtaIPhiFull5x5)[ipho];
	pho_SigmaIPhiIPhiFull5x5_ = (*phoSigmaIPhiIPhiFull5x5)[ipho];
	pho_E2x2Full5x5_ = (*phoE2x2Full5x5)[ipho];
	pho_E5x5Full5x5_ = (*phoE5x5Full5x5)[ipho];
	pho_R9Full5x5_ = (*phoR9Full5x5)[ipho];
	pho_PFChIso_ = (*phoPFChIso)[ipho];
	pho_PFChPVIso_ = (*phoPFChPVIso)[ipho];
	pho_PFPhoIso_ = (*phoPFPhoIso)[ipho];
	pho_PFNeuIso_ = (*phoPFNeuIso)[ipho];
	pho_PFChWorstVetoIso_ = (*phoPFChWorstVetoIso)[ipho];
	pho_PFChWorstIso_ = (*phoPFChWorstIso)[ipho];
	pho_EcalPFClusterIso_ = (*phoEcalPFClusterIso)[ipho];
	pho_HcalPFClusterIso_ = (*phoHcalPFClusterIso)[ipho];
	pho_SeedTime_  = (*phoSeedTime)[ipho];   
	pho_SeedEnergy_ = (*phoSeedEnergy)[ipho];
	pho_SeediEta_  = (*phoSeediEta)[ipho];
	pho_SeediPhi_  = (*phoSeediPhi)[ipho];

	pho_EnergyMatrix5x5_    = phoEnergyMatrix5x5->at(ipho);
	pho_EnergyMatrix7x7_    = phoEnergyMatrix7x7->at(ipho);
	pho_EnergyMatrix9x9_    = phoEnergyMatrix9x9->at(ipho); 
	pho_EnergyMatrix11x11_  = phoEnergyMatrix11x11->at(ipho);
	pho_EnergyMatrix15x15_  = phoEnergyMatrix15x15->at(ipho);

	mc_Pt_       = (*mcPt)[0]; 
	mc_PID_      = (*mcPID)[0]; 
	mc_MomPID_   = (*mcMomPID)[0]; 
	
	int pass = 0;
	bool isPhoPrompt = false; 
	float genPt = -1;
	if (debug) cout << " (*mcPID).size(): " << (*mcPID).size() << endl;
	for(unsigned int imc = 0; imc < (*mcPID).size(); imc++){
	  if((*mcPt)[imc] < 10 ) continue;	  
	  if((*mcStatus)[imc] != 1)continue; 
	  if((*mcPID)[imc] != 22)continue;
	  if(TMath::Abs((*mcMomPID)[imc]) >  21 )continue;
	  
	  double meta = (*mcEta)[imc];
	  double mphi = (*mcPhi)[imc];	  
	  TVector3 mcPhoton;
	  TVector3 recoPhoton;
	  mcPhoton.SetPtEtaPhi(1.0,meta,mphi);
	  recoPhoton.SetPtEtaPhi(1.0,pEta,pPhi);			       
	  double DR = mcPhoton.DrEtaPhi(recoPhoton);
	  double dp = fabs((*mcPt)[imc] - (*phoEt)[ipho] )/(*mcPt)[imc];
	  if (debug) cout << "DR: " << DR << "    dp: " << dp << endl;
	  if(DR < 0.2 && dp < 0.2  ){
	    if (debug) cout<< "pass : DR < 0.2 && dp < 0.1 " << endl;
	    pass++;
	    if (pass!=1) cout<< "*************************** Warning **************************   pass value: " << pass << endl;
	    if(pass == 1 ){
	      genPt = (*mcPt)[imc];
	      isPhoPrompt = true;
	    }
	  }
	}//EOF MC Particles loop
	
	if ( debug && genPt<0) cout << "genPt: " << genPt << endl;
	isTruePho_ = isPhoPrompt;
	pho_genPt_ = genPt;

	if (isPhoPrompt) tree_sig->Fill();
	else tree_bkg->Fill();
	
        if (debug) cout << "ipho: " << ipho << "  EnergyMatrix5x5 size: " << phoEnergyMatrix5x5->at(ipho).size() << endl;
        for (unsigned int ic=0; ic < phoEnergyMatrix5x5->at(ipho).size(); ic++){
          //cout << "phoEnergyMatix5x5 : " <<  phoEnergyMatrix5x5->at(ipho)[ic] << endl;
        }

	tree->Fill();
      } // End of iPho Loop 

   } // End of Event Loop
   
   fileout->Write();
   fileout->Close();
}
