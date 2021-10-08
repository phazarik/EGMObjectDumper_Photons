// -*- C++ -*-
//
// Package:    EGMObjectDumper/egmNtuplizer
// Class:      egmNtuplizer
//
/**\class egmNtuplizer egmNtuplizer.cc EGMObjectDumper/egmNtuplizer/plugins/egmNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ashim Roy
//         Created:  Sat, 13 Jun 2020 20:37:58 GMT
//
//

#include "../interface/egmNtuplizer.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
using namespace std;
using namespace edm;

void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}

egmNtuplizer::egmNtuplizer(const edm::ParameterSet& iConfig) 
//egmNtuplizer::egmNtuplizer(const edm::ParameterSet& iConfig) :
  //electronCollection_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons")))
  //hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif

  doGenParticles_            =                                          iConfig.getParameter<bool>("doGenParticles");  
  runOnParticleGun_          =                                          iConfig.getParameter<bool>("runOnParticleGun");
  runOnSherpa_               =                                          iConfig.getParameter<bool>("runOnSherpa");
  dumpPDFSystWeight_         =                                          iConfig.getParameter<bool>("dumpPDFSystWeight");
  dumpGenScaleSystWeights_   =                                          iConfig.getParameter<bool>("dumpGenScaleSystWeights");
  newparticles_              =                                          iConfig.getParameter< vector<int > >("newParticles");

  vtxLabel_                  = consumes<reco::VertexCollection>        (iConfig.getParameter<InputTag>("VtxLabel"));
  //vtxBSLabel_                = consumes<reco::VertexCollection>        (iConfig.getParameter<InputTag>("VtxBSLabel"));
  generatorLabel_            = consumes<GenEventInfoProduct>           (iConfig.getParameter<InputTag>("generatorLabel"));
  lheEventLabel_             = consumes<LHEEventProduct>               (iConfig.getParameter<InputTag>("LHEEventLabel"));
  puCollection_              = consumes<vector<PileupSummaryInfo> >    (iConfig.getParameter<InputTag>("pileupCollection"));
  genParticlesCollection_    = consumes<vector<reco::GenParticle> >    (iConfig.getParameter<InputTag>("genParticleSrc"));

  rhoLabel_                  = consumes<double>                        (iConfig.getParameter<InputTag>("rhoLabel"));
  rhoCentralLabel_           = consumes<double>                        (iConfig.getParameter<InputTag>("rhoCentralLabel"));
  //trgEventLabel_             = consumes<trigger::TriggerEvent>         (iConfig.getParameter<InputTag>("triggerEvent"));
  //triggerObjectsLabel_       = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerEvent"));
  trgResultsLabel_           = consumes<edm::TriggerResults>           (iConfig.getParameter<InputTag>("triggerResults"));
  //patTrgResultsLabel_        = consumes<edm::TriggerResults>           (iConfig.getParameter<InputTag>("patTriggerResults"));
  trgResultsProcess_         =                                          iConfig.getParameter<InputTag>("triggerResults").process();

  electronCollection_         = consumes<View<pat::Electron> > (iConfig.getParameter<InputTag>("electronSrc"));
  photonCollection_           = consumes<View<pat::Photon> >   (iConfig.getParameter<InputTag>("photonSrc"));
  ebReducedRecHitCollection_  = consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("ebReducedRecHitCollection"));
  eeReducedRecHitCollection_  = consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("eeReducedRecHitCollection"));
  esReducedRecHitCollection_  = consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("esReducedRecHitCollection"));

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("EventTree", "EventInfo");

  branchesGlobalEvent(tree);
  if (doGenParticles_) {
    branchesGenInfo(tree, fs);
    branchesGenPart(tree);
  }

  branchesPhotons(tree);
  branchesElectrons(tree);
  branchesRechits(tree);

}

egmNtuplizer::~egmNtuplizer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void egmNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //AR  for (const auto& track : iEvent.get(tracksToken_)) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
  //AR }

  fillGlobalEvent(iEvent, iSetup);
  if (!iEvent.isRealData()) {
    fillGenInfo(iEvent);
    if (doGenParticles_)
      fillGenPart(iEvent);
  }

  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(vtxLabel_, vtxHandle);

  reco::Vertex vtx;

  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0);
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    // replace isFake() for miniAOD since it requires tracks while miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = (v->chi2() == 0 && v->ndof() == 0);

    if (!isFake) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = *v;
      break;
    }
  }

  fillPhotons(iEvent, iSetup); 
  fillElectrons(iEvent, iSetup, pv); 
  fillRechits(iEvent, iSetup); 

  tree->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

/*
// ------------ method called once each job just before starting event loop  ------------
void egmNtuplizer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void egmNtuplizer::endJob() {
  // please remove this method if not needed
}
*/
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void egmNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(egmNtuplizer);
