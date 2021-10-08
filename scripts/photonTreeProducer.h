#ifndef photonTreeProducer_h
#define photonTreeProducer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

//Included by User
#include<iostream> // for std // cout etc..

using namespace std;

class photonTreeProducer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Int_t           nGoodVtx;
   Int_t           nTrksPV;
   Bool_t          isPVGood;
   Float_t         vtx;
   Float_t         vty;
   Float_t         vtz;
   Float_t         rho;
   Float_t         rhoCentral;
   vector<float>   *pdf;
   Float_t         pthat;
   Float_t         processID;
   Float_t         genWeight;
   Float_t         genHT;
   Float_t         genPho1;
   Float_t         genPho2;
   TString         *EventTag;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue;
   Int_t           nLHE;
   vector<int>     *lhePID;
   vector<float>   *lhePx;
   vector<float>   *lhePy;
   vector<float>   *lhePz;
   vector<float>   *lheE;
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx;
   vector<float>   *mcVty;
   vector<float>   *mcVtz;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcGMomPID;
   vector<int>     *mcMomPID;
   vector<float>   *mcMomPt;
   vector<float>   *mcMomMass;
   vector<float>   *mcMomEta;
   vector<float>   *mcMomPhi;
   vector<unsigned short> *mcStatusFlag;
   vector<int>     *mcParentage;
   vector<int>     *mcStatus;
   vector<float>   *mcCalIsoDR03;
   vector<float>   *mcTrkIsoDR03;
   vector<float>   *mcCalIsoDR04;
   vector<float>   *mcTrkIsoDR04;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoESEn;
   vector<float>   *phoESEnP1;
   vector<float>   *phoESEnP2;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9;
   vector<float>   *phoHoverE;
   vector<float>   *phoESEffSigmaRR;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIPhiFull5x5;
   vector<float>   *phoSigmaIPhiIPhiFull5x5;
   vector<float>   *phoE2x2Full5x5;
   vector<float>   *phoE5x5Full5x5;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFChPVIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoPFChWorstIso;
   vector<float>   *phoPFChWorstVetoIso;
   vector<float>   *phoEcalPFClusterIso;
   vector<float>   *phoHcalPFClusterIso;
   vector<float>   *phoSeedTime;
   vector<float>   *phoSeedEnergy;
   vector<int>     *phoSeediEta;
   vector<int>     *phoSeediPhi;
   vector<vector<float> > *phoEnergyMatrix5x5;
   vector<vector<float> > *phoEnergyMatrix7x7;
   vector<vector<float> > *phoEnergyMatrix9x9;
   vector<vector<float> > *phoEnergyMatrix11x11;
   vector<vector<float> > *phoEnergyMatrix15x15;
   Int_t           nEle;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleEn;
   vector<float>   *eleSCEn;
   vector<float>   *eleEcalEn;
   vector<float>   *eleESEnP1;
   vector<float>   *eleESEnP2;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleSIP;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *elePFClusEcalIso;
   vector<float>   *elePFClusHcalIso;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *eleR9Full5x5;
   vector<int>     *eleEcalDrivenSeed;
   vector<vector<float> > *eleEnergyMatrix5x5;
   vector<vector<float> > *eleEnergyMatrix7x7;
   vector<vector<float> > *eleEnergyMatrix9x9;
   vector<vector<float> > *eleEnergyMatrix11x11;
   vector<vector<float> > *eleEnergyMatrix15x15;
   Int_t           nebRechit;
   vector<float>   *ebRechitE;
   vector<int>     *ebRechitiEta;
   vector<int>     *ebRechitiPhi;
   vector<int>     *ebRechitZside;
   vector<float>   *ebRechitEta;
   vector<float>   *ebRechitPhi;
   vector<float>   *ebRechitX;
   vector<float>   *ebRechitY;
   vector<float>   *ebRechitZ;
   vector<int>     *ebRechitDetID;
   Int_t           neeRechit;
   vector<float>   *eeRechitE;
   vector<int>     *eeRechitiEta;
   vector<int>     *eeRechitiPhi;
   vector<int>     *eeRechitZside;
   vector<float>   *eeRechitEta;
   vector<float>   *eeRechitPhi;
   vector<float>   *eeRechitX;
   vector<float>   *eeRechitY;
   vector<float>   *eeRechitZ;
   vector<int>     *eeRechitDetID;
   Int_t           nesRechit;
   vector<float>   *esRechitE;
   vector<int>     *esRechitiEta;
   vector<int>     *esRechitiPhi;
   vector<int>     *esRechitZside;
   vector<int>     *esRechitPlane;
   vector<int>     *esRechitStrip;
   vector<float>   *esRechitEta;
   vector<float>   *esRechitPhi;
   vector<float>   *esRechitX;
   vector<float>   *esRechitY;
   vector<float>   *esRechitZ;

   // List of branches
   TBranch        *b_run;   //!                                                                                                                  
   TBranch        *b_event;   //!                                                                                                                
   TBranch        *b_lumis;   //!                                                                                                                
   TBranch        *b_isData;   //!                                                                                                               
   TBranch        *b_nVtx;   //!                                                                                                                 
   TBranch        *b_nGoodVtx;   //!                                                                                                             
   TBranch        *b_nTrksPV;   //!                                                                                                              
   TBranch        *b_isPVGood;   //!                                                                                                             
   TBranch        *b_vtx;   //!                                                                                                                  
   TBranch        *b_vty;   //!                                                                                                                  
   TBranch        *b_vtz;   //!                                                                                                                  
   TBranch        *b_rho;   //!                                                                                                                   
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genHT;   //!
   TBranch        *b_genPho1;   //!
   TBranch        *b_genPho2;   //!
   TBranch        *b_EventTag;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_nLHE;   //!
   TBranch        *b_lhePID;   //!
   TBranch        *b_lhePx;   //!
   TBranch        *b_lhePy;   //!
   TBranch        *b_lhePz;   //!
   TBranch        *b_lheE;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcVty;   //!
   TBranch        *b_mcVtz;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcStatusFlag;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoESEnP1;   //!
   TBranch        *b_phoESEnP2;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_phoE2x2Full5x5;   //!
   TBranch        *b_phoE5x5Full5x5;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFChPVIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   TBranch        *b_phoPFChWorstVetoIso;   //!
   TBranch        *b_phoEcalPFClusterIso;   //!
   TBranch        *b_phoHcalPFClusterIso;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoSeediEta;   //!
   TBranch        *b_phoSeediPhi;   //!
   TBranch        *b_phoEnergyMatrix5x5;   //!
   TBranch        *b_phoEnergyMatrix7x7;   //!
   TBranch        *b_phoEnergyMatrix9x9;   //!
   TBranch        *b_phoEnergyMatrix11x11;   //!
   TBranch        *b_phoEnergyMatrix15x15;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_eleESEnP1;   //!
   TBranch        *b_eleESEnP2;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleSIP;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_elePFClusEcalIso;   //!
   TBranch        *b_elePFClusHcalIso;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleEcalDrivenSeed;   //!
   TBranch        *b_eleEnergyMatrix5x5;   //!
   TBranch        *b_eleEnergyMatrix7x7;   //!
   TBranch        *b_eleEnergyMatrix9x9;   //!
   TBranch        *b_eleEnergyMatrix11x11;   //!
   TBranch        *b_eleEnergyMatrix15x15;   //!
   TBranch        *b_nebRechit;   //!
   TBranch        *b_ebRechitE;   //!
   TBranch        *b_ebRechitiEta;   //!
   TBranch        *b_ebRechitiPhi;   //!
   TBranch        *b_ebRechitZside;   //!
   TBranch        *b_ebRechitEta;   //!
   TBranch        *b_ebRechitPhi;   //!
   TBranch        *b_ebRechitX;   //!
   TBranch        *b_ebRechitY;   //!
   TBranch        *b_ebRechitZ;   //!
   TBranch        *b_ebRechitDetID;   //!
   TBranch        *b_neeRechit;   //!
   TBranch        *b_eeRechitE;   //!
   TBranch        *b_eeRechitiEta;   //!
   TBranch        *b_eeRechitiPhi;   //!
   TBranch        *b_eeRechitZside;   //!
   TBranch        *b_eeRechitEta;   //!
   TBranch        *b_eeRechitPhi;   //!
   TBranch        *b_eeRechitX;   //!
   TBranch        *b_eeRechitY;   //!
   TBranch        *b_eeRechitZ;   //!
   TBranch        *b_eeRechitDetID;   //!
   TBranch        *b_nesRechit;   //!
   TBranch        *b_esRechitE;   //!
   TBranch        *b_esRechitiEta;   //!
   TBranch        *b_esRechitiPhi;   //!
   TBranch        *b_esRechitZside;   //!
   TBranch        *b_esRechitPlane;   //!
   TBranch        *b_esRechitStrip;   //!
   TBranch        *b_esRechitEta;   //!
   TBranch        *b_esRechitPhi;   //!
   TBranch        *b_esRechitX;   //!
   TBranch        *b_esRechitY;   //!
   TBranch        *b_esRechitZ;   //!

  //photonTreeProducer(TTree *tree=0);
   photonTreeProducer(string infile, TTree *tree=0);
   virtual ~photonTreeProducer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
  virtual void     Loop(string outputfile, int nTotEvent, int nPrintEvent, float eventWeight);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef photonTreeProducer_cxx
//photonTreeProducer::photonTreeProducer(TTree *tree) : fChain(0)
photonTreeProducer::photonTreeProducer(string infile, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  /*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("egmTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("egmTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("egmTree.root:/egmNtuplizer");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);
  */

  
   TString inputfile = infile;
   cout<<"Input file : " << inputfile << endl;

   if (tree == 0) {
     TFile *f = TFile::Open(inputfile);

     TDirectory *dir = dynamic_cast<TDirectory*>(f->Get("egmNtuplizer"));
     dir->GetObject("EventTree",tree);
   }
   Init(tree);

}

photonTreeProducer::~photonTreeProducer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t photonTreeProducer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t photonTreeProducer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void photonTreeProducer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pdf = 0;
   EventTag = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   lhePID = 0;
   lhePx = 0;
   lhePy = 0;
   lhePz = 0;
   lheE = 0;
   mcPID = 0;
   mcVtx = 0;
   mcVty = 0;
   mcVtz = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcStatusFlag = 0;
   mcParentage = 0;
   mcStatus = 0;
   mcCalIsoDR03 = 0;
   mcTrkIsoDR03 = 0;
   mcCalIsoDR04 = 0;
   mcTrkIsoDR04 = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoESEn = 0;
   phoESEnP1 = 0;
   phoESEnP2 = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoHoverE = 0;
   phoESEffSigmaRR = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   phoE2x2Full5x5 = 0;
   phoE5x5Full5x5 = 0;
   phoR9Full5x5 = 0;
   phoPFChIso = 0;
   phoPFChPVIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoPFChWorstIso = 0;
   phoPFChWorstVetoIso = 0;
   phoEcalPFClusterIso = 0;
   phoHcalPFClusterIso = 0;
   phoSeedTime = 0;
   phoSeedEnergy = 0;
   phoSeediEta = 0;
   phoSeediPhi = 0;
   phoEnergyMatrix5x5 = 0;
   phoEnergyMatrix7x7 = 0;
   phoEnergyMatrix9x9 = 0;
   phoEnergyMatrix11x11 = 0;
   phoEnergyMatrix15x15 = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleSCEn = 0;
   eleEcalEn = 0;
   eleESEnP1 = 0;
   eleESEnP2 = 0;
   eleD0 = 0;
   eleDz = 0;
   eleSIP = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCRawEn = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   eleR9Full5x5 = 0;
   eleEcalDrivenSeed = 0;
   eleEnergyMatrix5x5 = 0;
   eleEnergyMatrix7x7 = 0;
   eleEnergyMatrix9x9 = 0;
   eleEnergyMatrix11x11 = 0;
   eleEnergyMatrix15x15 = 0;
   ebRechitE = 0;
   ebRechitiEta = 0;
   ebRechitiPhi = 0;
   ebRechitZside = 0;
   ebRechitEta = 0;
   ebRechitPhi = 0;
   ebRechitX = 0;
   ebRechitY = 0;
   ebRechitZ = 0;
   ebRechitDetID = 0;
   eeRechitE = 0;
   eeRechitiEta = 0;
   eeRechitiPhi = 0;
   eeRechitZside = 0;
   eeRechitEta = 0;
   eeRechitPhi = 0;
   eeRechitX = 0;
   eeRechitY = 0;
   eeRechitZ = 0;
   eeRechitDetID = 0;
   esRechitE = 0;
   esRechitiEta = 0;
   esRechitiPhi = 0;
   esRechitZside = 0;
   esRechitPlane = 0;
   esRechitStrip = 0;
   esRechitEta = 0;
   esRechitPhi = 0;
   esRechitX = 0;
   esRechitY = 0;
   esRechitZ = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("genPho1", &genPho1, &b_genPho1);
   fChain->SetBranchAddress("genPho2", &genPho2, &b_genPho2);
   fChain->SetBranchAddress("EventTag", &EventTag, &b_EventTag);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nLHE", &nLHE, &b_nLHE);
   fChain->SetBranchAddress("lhePID", &lhePID, &b_lhePID);
   fChain->SetBranchAddress("lhePx", &lhePx, &b_lhePx);
   fChain->SetBranchAddress("lhePy", &lhePy, &b_lhePy);
   fChain->SetBranchAddress("lhePz", &lhePz, &b_lhePz);
   fChain->SetBranchAddress("lheE", &lheE, &b_lheE);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);
   fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03, &b_mcCalIsoDR03);
   fChain->SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03, &b_mcTrkIsoDR03);
   fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04);
   fChain->SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04, &b_mcTrkIsoDR04);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   fChain->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFChPVIso", &phoPFChPVIso, &b_phoPFChPVIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   fChain->SetBranchAddress("phoPFChWorstVetoIso", &phoPFChWorstVetoIso, &b_phoPFChWorstVetoIso);
   fChain->SetBranchAddress("phoEcalPFClusterIso", &phoEcalPFClusterIso, &b_phoEcalPFClusterIso);
   fChain->SetBranchAddress("phoHcalPFClusterIso", &phoHcalPFClusterIso, &b_phoHcalPFClusterIso);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoSeediEta", &phoSeediEta, &b_phoSeediEta);
   fChain->SetBranchAddress("phoSeediPhi", &phoSeediPhi, &b_phoSeediPhi);
   fChain->SetBranchAddress("phoEnergyMatrix5x5", &phoEnergyMatrix5x5, &b_phoEnergyMatrix5x5);
   fChain->SetBranchAddress("phoEnergyMatrix7x7", &phoEnergyMatrix7x7, &b_phoEnergyMatrix7x7);
   fChain->SetBranchAddress("phoEnergyMatrix9x9", &phoEnergyMatrix9x9, &b_phoEnergyMatrix9x9);
   fChain->SetBranchAddress("phoEnergyMatrix11x11", &phoEnergyMatrix11x11, &b_phoEnergyMatrix11x11);
   fChain->SetBranchAddress("phoEnergyMatrix15x15", &phoEnergyMatrix15x15, &b_phoEnergyMatrix15x15);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleEcalEn", &eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("eleESEnP1", &eleESEnP1, &b_eleESEnP1);
   fChain->SetBranchAddress("eleESEnP2", &eleESEnP2, &b_eleESEnP2);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleSIP", &eleSIP, &b_eleSIP);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   fChain->SetBranchAddress("eleEnergyMatrix5x5", &eleEnergyMatrix5x5, &b_eleEnergyMatrix5x5);
   fChain->SetBranchAddress("eleEnergyMatrix7x7", &eleEnergyMatrix7x7, &b_eleEnergyMatrix7x7);
   fChain->SetBranchAddress("eleEnergyMatrix9x9", &eleEnergyMatrix9x9, &b_eleEnergyMatrix9x9);
   fChain->SetBranchAddress("eleEnergyMatrix11x11", &eleEnergyMatrix11x11, &b_eleEnergyMatrix11x11);
   fChain->SetBranchAddress("eleEnergyMatrix15x15", &eleEnergyMatrix15x15, &b_eleEnergyMatrix15x15);
   fChain->SetBranchAddress("nebRechit", &nebRechit, &b_nebRechit);
   fChain->SetBranchAddress("ebRechitE", &ebRechitE, &b_ebRechitE);
   fChain->SetBranchAddress("ebRechitiEta", &ebRechitiEta, &b_ebRechitiEta);
   fChain->SetBranchAddress("ebRechitiPhi", &ebRechitiPhi, &b_ebRechitiPhi);
   fChain->SetBranchAddress("ebRechitZside", &ebRechitZside, &b_ebRechitZside);
   fChain->SetBranchAddress("ebRechitEta", &ebRechitEta, &b_ebRechitEta);
   fChain->SetBranchAddress("ebRechitPhi", &ebRechitPhi, &b_ebRechitPhi);
   fChain->SetBranchAddress("ebRechitX", &ebRechitX, &b_ebRechitX);
   fChain->SetBranchAddress("ebRechitY", &ebRechitY, &b_ebRechitY);
   fChain->SetBranchAddress("ebRechitZ", &ebRechitZ, &b_ebRechitZ);
   fChain->SetBranchAddress("ebRechitDetID", &ebRechitDetID, &b_ebRechitDetID);
   fChain->SetBranchAddress("neeRechit", &neeRechit, &b_neeRechit);
   fChain->SetBranchAddress("eeRechitE", &eeRechitE, &b_eeRechitE);
   fChain->SetBranchAddress("eeRechitiEta", &eeRechitiEta, &b_eeRechitiEta);
   fChain->SetBranchAddress("eeRechitiPhi", &eeRechitiPhi, &b_eeRechitiPhi);
   fChain->SetBranchAddress("eeRechitZside", &eeRechitZside, &b_eeRechitZside);
   fChain->SetBranchAddress("eeRechitEta", &eeRechitEta, &b_eeRechitEta);
   fChain->SetBranchAddress("eeRechitPhi", &eeRechitPhi, &b_eeRechitPhi);
   fChain->SetBranchAddress("eeRechitX", &eeRechitX, &b_eeRechitX);
   fChain->SetBranchAddress("eeRechitY", &eeRechitY, &b_eeRechitY);
   fChain->SetBranchAddress("eeRechitZ", &eeRechitZ, &b_eeRechitZ);
   fChain->SetBranchAddress("eeRechitDetID", &eeRechitDetID, &b_eeRechitDetID);
   fChain->SetBranchAddress("nesRechit", &nesRechit, &b_nesRechit);
   fChain->SetBranchAddress("esRechitE", &esRechitE, &b_esRechitE);
   fChain->SetBranchAddress("esRechitiEta", &esRechitiEta, &b_esRechitiEta);
   fChain->SetBranchAddress("esRechitiPhi", &esRechitiPhi, &b_esRechitiPhi);
   fChain->SetBranchAddress("esRechitZside", &esRechitZside, &b_esRechitZside);
   fChain->SetBranchAddress("esRechitPlane", &esRechitPlane, &b_esRechitPlane);
   fChain->SetBranchAddress("esRechitStrip", &esRechitStrip, &b_esRechitStrip);
   fChain->SetBranchAddress("esRechitEta", &esRechitEta, &b_esRechitEta);
   fChain->SetBranchAddress("esRechitPhi", &esRechitPhi, &b_esRechitPhi);
   fChain->SetBranchAddress("esRechitX", &esRechitX, &b_esRechitX);
   fChain->SetBranchAddress("esRechitY", &esRechitY, &b_esRechitY);
   fChain->SetBranchAddress("esRechitZ", &esRechitZ, &b_esRechitZ);
   Notify();
}

Bool_t photonTreeProducer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void photonTreeProducer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t photonTreeProducer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef photonTreeProducer_cxx
