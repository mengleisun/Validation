// -*- C++ -*-
//
// Package:    EcalValidation
// Class:      EcalValidation

// 
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/EventBase.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h" 
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "TLorentzVector.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalCleaningAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoEgamma/EgammaTools/interface/ECALPositionCalculator.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "Validation/EcalValidation/interface/EcalValidationDigiless_ES.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"


#include "CondFormats/DataRecord/interface/ESChannelStatusRcd.h"
#include "CondFormats/ESObjects/interface/ESChannelStatus.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"



#include "TVector3.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <map>

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;

//
// constructors and destructor
//
EcalValidationDigiless_ES::EcalValidationDigiless_ES(const edm::ParameterSet& ps)
{
  //  std::cout << " >>> in EcalValidationDigiless_ES::EcalValidationDigiless_ES " << std::endl;

  //now do what ever initialization is needed
  PVTag_                      = consumes<reco::VertexCollection>(ps.getParameter<edm::InputTag>("PVTag"));
  recHitCollection_EB_       = ps.getParameter<edm::InputTag>("recHitCollection_EB");
  recHitCollection_EE_       = ps.getParameter<edm::InputTag>("recHitCollection_EE");
  redRecHitCollection_EB_    = ps.getParameter<edm::InputTag>("redRecHitCollection_EB");
  redRecHitCollection_EE_    = ps.getParameter<edm::InputTag>("redRecHitCollection_EE");
  basicClusterCollection_EB_ = ps.getParameter<edm::InputTag>("basicClusterCollection_EB");
  basicClusterCollection_EE_ = ps.getParameter<edm::InputTag>("basicClusterCollection_EE");
  superClusterCollection_EB_ = ps.getParameter<edm::InputTag>("superClusterCollection_EB");


  //questo e' con pf ES
  superClusterCollection_EE_ = consumes<reco::SuperClusterCollection>(ps.getParameter<edm::InputTag>("superClusterCollection_EE"));

  //reducedRecHit
  esRecHitCollection_        = ps.getParameter<edm::InputTag>("recHitCollection_ES");

  esClusterCollectionX_      = ps.getParameter<edm::InputTag>("ClusterCollectionX_ES");
  esClusterCollectionY_      = ps.getParameter<edm::InputTag>("ClusterCollectionY_ES");
  basicClusterCollection_ES_ = consumes<reco::BasicClusterCollection>(ps.getParameter<edm::InputTag>("basicClusterCollection_ES"));

  beamSpot_                  = ps.getParameter<edm::InputTag>("beamSpot");
  jets_                      = ps.getParameter<edm::InputTag>("jets");

  ethrEB_                    = ps.getParameter<double>("ethrEB");
  ethrEE_                    = ps.getParameter<double>("ethrEE");
  gainId_                    = ps.getParameter<double>("gainId");

  scEtThrEB_                 = ps.getParameter<double>("scEtThrEB");
  scEtThrEE_                 = ps.getParameter<double>("scEtThrEE");
  

  elePFlow_                  = consumes<reco::GsfElectronCollection>(ps.getParameter<edm::InputTag>("elePFlow"));
  elePFlowCand_              = consumes<reco::PFCandidateCollection>(ps.getParameter<edm::InputTag>("elePFlowCand"));
  esPFRec_                   = consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("esPFRec"));

  naiveId_ = 0;


  // histos 
  
  edm::Service<TFileService> fs;
  
  h_numberOfEvents = fs->make<TH1D>("h_numberOfEvents","h_numberOfEvents",10,0,10);
   
  h_PV_n = fs->make<TH1D>("h_PV_n","h_PV_n",50,0.,50.);
  h_PV_cut_n = fs->make<TH1D>("h_PV_cut_n","h_PV_cut_n",50,0.,50.);
  
  // ReducedRecHits ----------------------------------------------
  // ... barrel 
  h_redRecHits_EB_recoFlag = fs->make<TH1D>("h_redRecHits_EB_recoFlag","h_redRecHits_EB_recoFlag",16,-0.5,15.5);  
  // ... endcap 
  h_redRecHits_EE_recoFlag = fs->make<TH1D>("h_redRecHits_EE_recoFlag","h_redRecHits_EE_recoFlag",16,-0.5,15.5);  
  // ... all 
  h_redRecHits_recoFlag = fs->make<TH1D>("h_redRecHits_recoFlag","h_redRecHits_recoFlag",16,-0.5,15.5);  

  
  // RecHits ---------------------------------------------- 
  

  // Super Clusters ----------------------------------------------
  // ... endcap
  h_superClusters_EEP_energy = fs->make<TH1D>("h_superClusters_EEP_energy","h_superClusters_EEP_energy",2000,0.,400.);
  h_superClusters_EEP_rawEnergy = fs->make<TH1D>("h_superClusters_EEP_rawEnergy","h_superClusters_EEP_rawEnergy",2000,0.,400.);
  h_superClusters_EEP_rawEt = fs->make<TH1D>("h_superClusters_EEP_rawEt","h_superClusters_EEP_rawEt",2000,0.,400.);
  h_superClusters_EEP_occupancy     = fs->make<TH2D>("h_superClusters_EEP_occupancy","h_superClusters_EEP_occupancy",100,0.,100.,100,0.,100.);

  h_superClusters_EEM_energy = fs->make<TH1D>("h_superClusters_EEM_energy","h_superClusters_EEM_energy",2000,0.,400.);
  h_superClusters_EEM_rawEnergy = fs->make<TH1D>("h_superClusters_EEM_rawEnergy","h_superClusters_EEM_rawEnergy",2000,0.,400.);
  h_superClusters_EEM_rawEt = fs->make<TH1D>("h_superClusters_EEM_rawEt","h_superClusters_EEM_rawEt",2000,0.,400.);
  h_superClusters_EEM_occupancy     = fs->make<TH2D>("h_superClusters_EEM_occupancy","h_superClusters_EEM_occupancy",100,0.,100.,100,0.,100.);
  
  h_superClusters_occupancyPhiEta = fs->make<TH2D>("h_superClusters_occupancyPhiEta","h_superClusters_occupancyPhiEta",360,-3.1415927,3.1415927,150,-3.,3.);

  h_superClusters_eta        = fs->make<TH1D>("h_superClusters_eta","h_superClusters_eta",150,-3.,3.);
  
  // preshower

  h_gedEle_ESenergy = fs->make<TH1D>("h_gedEle_ESenergy", "", 1000,0.,50.);
  h_gedEle_ESenergy_plane1 = fs->make<TH1D>("h_gedEle_ESenergy_plane1", "", 1000,0.,50.);               
  h_gedEle_ESenergy_plane2 = fs->make<TH1D>("h_gedEle_ESenergy_plane2", "", 1000,0.,50.);
  h_gedEle_ESenergy_plane1_Allok = fs->make<TH1D>("h_gedEle_ESenergy_plane1_Allok", "", 1000,0.,50.);
  h_gedEle_ESenergy_plane1_Alldead = fs->make<TH1D>("h_gedEle_ESenergy_plane1_Alldead", "", 1000,0.,50.);
  h_gedEle_ESenergy_plane1_NOok = fs->make<TH1D>("h_gedEle_ESenergy_plane1_NOok", "", 1000,0.,50.);
  h_gedEle_ESenergy_plane2_Allok = fs->make<TH1D>("h_gedEle_ESenergy_plane2_Allok", "", 1000,0.,50.);
  h_gedEle_ESenergy_plane2_Alldead = fs->make<TH1D>("h_gedEle_ESenergy_plane2_Alldead", "", 1000,0.,50.);
  h_gedEle_ESenergy_plane2_NOok = fs->make<TH1D>("h_gedEle_ESenergy_plane2_NOok", "", 1000,0.,50.);

  h_gedEle_ESenergy_1corrected = fs->make<TH1D>("h_gedEle_ESenergy_1corrected", "", 1000,0.,50.);
  h_gedEle_ESenergy_2corrected = fs->make<TH1D>("h_gedEle_ESenergy_2corrected", "", 1000,0.,50.);
  h_gedEle_ESenergy_1overS = fs->make<TH1D>("h_gedEle_ESenergy_1overS", "", 1000,-10.,20.);
  h_gedEle_ESenergy_2overS = fs->make<TH1D>("h_gedEle_ESenergy_2overS", "", 1000,-10.,20.);

  h_ESenergy = fs->make<TH1D>("h_ESenergy", "", 1000,0.,50.);
  h_ESenergy_plane1 = fs->make<TH1D>("h_ESenergy_plane1", "", 1000,0.,50.);
  h_ESenergy_plane2 = fs->make<TH1D>("h_ESenergy_plane2", "", 1000,0.,50.);
  h_ESenergy_plane1_Allok = fs->make<TH1D>("h_ESenergy_plane1_Allok", "", 1000,0.,50.);
  h_ESenergy_plane1_Alldead = fs->make<TH1D>("h_ESenergy_plane1_Alldead", "", 1000,0.,50.);
  h_ESenergy_plane1_NOok = fs->make<TH1D>("h_ESenergy_plane1_NOok", "", 1000,0.,50.);
  h_ESenergy_plane2_Allok = fs->make<TH1D>("h_ESenergy_plane2_Allok", "", 1000,0.,50.);
  h_ESenergy_plane2_Alldead = fs->make<TH1D>("h_ESenergy_plane2_Alldead", "", 1000,0.,50.);
  h_ESenergy_plane2_NOok = fs->make<TH1D>("h_ESenergy_plane2_NOok", "", 1000,0.,50.);


  h_gedEle_ESenergy_2ok = fs->make<TH1D>("h_gedEle_ESenergy_2ok", "", 1000,0.,50.); 
  h_gedEle_ESenergy_2okKill1 = fs->make<TH1D>("h_gedEle_ESenergy_2okKill1", "", 1000,0.,50.); 
  h_gedEle_ESenergy_2okKill2 = fs->make<TH1D>("h_gedEle_ESenergy_2okKill2", "", 1000,0.,50.); 

  isZ = new std::vector<int>;
  eleEta = new std::vector<float>;
  elePhi = new std::vector<float>;
  elePt = new std::vector<float>;
  eleP = new std::vector<float>;
  eleEnergy = new std::vector<float>;
  elefBrem = new std::vector<float>;
  eletrackBrem = new std::vector<float>;
  eleSCBrem = new std::vector<float>;
  ele5x5 = new std::vector<float>;
  scEnergy = new std::vector<float>;
  scRawEnergy = new std::vector<float>;
  scEta = new std::vector<float>;
  scPhi = new std::vector<float>;
  deltaEtaEleClusterAtCalo = new std::vector<float>;
  deltaPhiEleClusterAtCalo = new std::vector<float>;
  deltaEtaEleClusterAtVtx = new std::vector<float>;
  deltaPhiEleClusterAtVtx = new std::vector<float>;

  isP1ok = new std::vector<int>;
  isP2ok = new std::vector<int>;

  eleES1 = new std::vector<float>;
  eleES2 = new std::vector<float>;

  ele_esRhFrac = new std::vector<float>;
  ele_esRhZ = new std::vector<int>;
  ele_esRhP = new std::vector<int>;
  ele_esRhX = new std::vector<int>;
  ele_esRhY = new std::vector<int>;
  ele_esRhS = new std::vector<int>;


  trkESz = new std::vector<float>;
  trkESp = new std::vector<float>;
  trkESx = new std::vector<float>;
  trkESy = new std::vector<float>;
  trkESs = new std::vector<float>;
  trkEEz = new std::vector<float>;
  trkEEx = new std::vector<float>;
  trkEEy = new std::vector<float>;
  trkEEEta = new std::vector<float>;
  trkEEPhi = new std::vector<float>;


  newT = fs->make<TTree>("newT","");
  newT->Branch("nEle", &nEle);
  newT->Branch("isZ", isZ);
  newT->Branch("eleEta", eleEta);
  newT->Branch("elePhi", elePhi);
  newT->Branch("elePt", elePt);
  newT->Branch("eleP", eleP);
  newT->Branch("eleEnergy", eleEnergy);
  newT->Branch("elefBrem", elefBrem);
  newT->Branch("eletrackBrem", eletrackBrem);
  newT->Branch("eleSCBrem", eleSCBrem);

  newT->Branch("ele5x5", ele5x5);
  newT->Branch("scEnergy", scEnergy);
  newT->Branch("scRawEnergy", scRawEnergy); //rawEnergy();
  newT->Branch("scEta", scEta);
  newT->Branch("scPhi", scPhi);
  newT->Branch("deltaEtaEleClusterAtCalo", deltaEtaEleClusterAtCalo);
  newT->Branch("deltaPhiEleClusterAtCalo", deltaPhiEleClusterAtCalo);
  newT->Branch("deltaEtaEleClusterAtVtx", deltaEtaEleClusterAtVtx);
  newT->Branch("deltaPhiEleClusterAtVtx", deltaPhiEleClusterAtVtx);

  //ES
  newT->Branch("isP1ok", isP1ok);
  newT->Branch("isP2ok", isP2ok);

  newT->Branch("eleES1", eleES1);
  newT->Branch("eleES2", eleES2);

  newT->Branch("ele_esRhFrac", ele_esRhFrac);
  newT->Branch("ele_esRhZ", ele_esRhZ);
  newT->Branch("ele_esRhP", ele_esRhP);
  newT->Branch("ele_esRhX", ele_esRhX);
  newT->Branch("ele_esRhY", ele_esRhY);
  newT->Branch("ele_esRhS", ele_esRhS);

  //track extrapolation
  newT->Branch("trkESz", trkESz);
  newT->Branch("trkESp", trkESp);
  newT->Branch("trkESx", trkESx);
  newT->Branch("trkESy", trkESy);
  newT->Branch("trkESs", trkESs);
  newT->Branch("trkEEz", trkEEz);
  newT->Branch("trkEEx", trkEEx);
  newT->Branch("trkEEy", trkEEy);
  newT->Branch("trkEEEta", trkEEEta);
  newT->Branch("trkEEPhi", trkEEPhi);
  
}



EcalValidationDigiless_ES::~EcalValidationDigiless_ES()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void EcalValidationDigiless_ES::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{

  //  std::cout << " in analyze event = " << std::endl;

//   float bx      = ev.bunchCrossing();
//   float ls      = ev.luminosityBlock();
//   float orbitNb = ev.orbitNumber();
 
  // Vertex Collection
  edm::Handle<reco::VertexCollection> vertexes;
  //  ev.getByLabel(PVTag_, vertexes); 
  ev.getByToken(PVTag_, vertexes); 
  if(vertexes->size() != 1) h_PV_n->Fill(vertexes->size());

  //  std::cout << " vertexes->size() = " << vertexes->size() << std::endl;
  //  if(vertexes->size() != 10.) return;
  h_PV_cut_n->Fill(vertexes->size());
  

  //Get the magnetic field
  edm::ESHandle<MagneticField> theMagField;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);

  //DataBase LaserCorrection
  edm::ESHandle<EcalLaserDbService> theLaser;
  iSetup.get<EcalLaserDbRecord>().get(theLaser);

  // ged electrons
  edm::Handle<reco::GsfElectronCollection> eleCand;
  ev.getByToken(elePFlow_, eleCand);
  const reco::GsfElectronCollection * electronCollection = eleCand.product();


  edm::Handle<reco::PFCandidateCollection> gedPFcand;
  ev.getByToken(elePFlowCand_, gedPFcand);
  const reco::PFCandidateCollection * pfCandidateCollection = gedPFcand.product();

  /*
  // Get PF RecHits
  edm::Handle<PFRecHitCollection> esPFRecHits;
  ev.getByToken(esPFRec_, esPFRecHits);
  //  const reco::PFRecHitCollection* esPFRecHitsCollection = esPFRecHits.product();
  */

  // EE supercluster
  edm::Handle<reco::SuperClusterCollection> superClusters_EE_h;
  ev.getByToken( superClusterCollection_EE_, superClusters_EE_h );
  const reco::SuperClusterCollection* theEndcapSuperClusters = superClusters_EE_h.product () ;
  if ( ! superClusters_EE_h.isValid() ) {
    std::cerr << "EcalValidation::analyze --> superClusters_EE_h not found" << std::endl; 
  }

  //ES channel status
  edm::ESHandle<ESChannelStatus> esChannelStatusHandle_;
  iSetup.get<ESChannelStatusRcd>().get(esChannelStatusHandle_);
  channelStatus_ = esChannelStatusHandle_.product();


  // ... ES basic cluster
  edm::Handle<reco::BasicClusterCollection> basicClusters_ES_h;
  ev.getByToken( basicClusterCollection_ES_, basicClusters_ES_h );
  //  const reco::BasicClusterCollection* theESBasicClusters = basicClusters_ES_h.product () ;
  if(!basicClusters_ES_h.isValid()){
    std::cout << "EcalValidation::analyze --> basicClusters_ES_h not found" << std::endl;
  }

  /*
  //tracks
  Handle<reco::TrackCollection> Tracks;
  if(ev.getByLabel(tracklabel_, Tracks)) {
    const reco::TrackCollection *track = Tracks.product();
  */


  //initialize
  isZ->clear();
  eleEta->clear();
  elePhi->clear();
  elePt->clear();
  eleP->clear();
  eleEnergy->clear();
  elefBrem->clear();
  eletrackBrem->clear();
  eleSCBrem->clear();
  ele5x5->clear();
  scEnergy->clear();
  scRawEnergy->clear();
  scEta->clear();
  scPhi->clear();
  deltaEtaEleClusterAtCalo->clear();
  deltaPhiEleClusterAtCalo->clear();
  deltaEtaEleClusterAtVtx->clear();
  deltaPhiEleClusterAtVtx->clear();

  isP1ok->clear();
  isP2ok->clear();

  eleES1->clear();
  eleES2->clear();
  ele_esRhFrac->clear();
  ele_esRhZ->clear();
  ele_esRhP->clear();
  ele_esRhX->clear();
  ele_esRhY->clear();
  ele_esRhS->clear();
  trkESz->clear();
  trkESp->clear();
  trkESx->clear();
  trkESy->clear();
  trkESs->clear();
  trkEEz->clear();
  trkEEx->clear();
  trkEEy->clear();
  trkEEEta->clear();
  trkEEPhi->clear();

  nEle = 0;
  //  if(int(isP1ok->size()) != nEle) std::cout << " >>> problemaAAAAAAAAAAAAAAAAAAAAAAA " << std::endl;

  ///////////// start analysis //////////////

  ++naiveId_;


  int countZ1 = 0;
  int countZ2 = 0;
  //selection on Zee
  //bool isZ = false;
  reco::GsfElectronCollection::const_iterator eleIt;
  reco::GsfElectronCollection::const_iterator iEle1 = electronCollection->begin();
  reco::GsfElectronCollection::const_iterator iEle2 = electronCollection->begin();
  TLorentzVector e1, e2;
  //  std::cout << " electronCollection->size() = " << electronCollection->size() << std::endl;
  for(eleIt=electronCollection->begin(); eleIt!=electronCollection->end(); eleIt++){

    if(eleIt->pt() > 10 && eleIt->pt() > iEle1->pt()) {
      e2.SetPtEtaPhiE(iEle1->pt(), iEle1->eta(), iEle1->phi(), iEle1->energy());
      e1.SetPtEtaPhiE(eleIt->pt(), eleIt->eta(), eleIt->phi(), eleIt->energy());
      iEle2 = iEle1;
      iEle1 = eleIt;
      countZ2 = countZ1;
      countZ1 = nEle;
    }


    eleEta->push_back(eleIt->eta());
    elePhi->push_back(eleIt->phi());
    elePt->push_back(eleIt->pt());
    eleEnergy->push_back(eleIt->energy());
    eleP->push_back(eleIt->p());
    elefBrem->push_back(eleIt->fbrem());
    eletrackBrem->push_back(eleIt->trackFbrem());
    eleSCBrem->push_back(eleIt->superClusterFbrem());
    ele5x5->push_back(eleIt->e5x5());
    scEnergy->push_back(eleIt->superCluster()->energy());
    scRawEnergy->push_back(eleIt->superCluster()->rawEnergy());
    scEta->push_back(eleIt->superCluster()->eta());
    scPhi->push_back(eleIt->superCluster()->phi());
    deltaEtaEleClusterAtCalo->push_back(eleIt->deltaEtaEleClusterTrackAtCalo());
    deltaPhiEleClusterAtCalo->push_back(eleIt->deltaPhiEleClusterTrackAtCalo());
    deltaEtaEleClusterAtVtx->push_back(eleIt->deltaEtaSuperClusterTrackAtVtx());
    deltaPhiEleClusterAtVtx->push_back(eleIt->deltaPhiSuperClusterTrackAtVtx());
    isZ->push_back(0);

    ////////////////
    //    std::cout << " entering ES part " << std::endl;

    bool okHaiFillato = false;
    //ES part
    if(fabs(eleIt->superCluster()->eta()) > 1.65 && fabs(eleIt->superCluster()->eta()) < 2.6 ){

      //gedEle
      h_gedEle_ESenergy->Fill(eleIt->superCluster()->preshowerEnergy());
      h_gedEle_ESenergy_plane1->Fill(eleIt->superCluster()->preshowerEnergyPlane1());
      h_gedEle_ESenergy_plane2->Fill(eleIt->superCluster()->preshowerEnergyPlane2());

      eleES1->push_back(eleIt->superCluster()->preshowerEnergyPlane1());
      eleES2->push_back(eleIt->superCluster()->preshowerEnergyPlane2());
	
      reco::GsfElectronRef eleRef(eleCand,nEle);
      reco::PFClusterRef clustops;
      
      for(unsigned int i=0; i<pfCandidateCollection->size(); i++){
	const reco::PFCandidate& iPFc = (*pfCandidateCollection)[i];
	if(iPFc.particleId()!=reco::PFCandidate::e) continue;
	
	float dEta = fabs(iPFc.eta() - eleIt->eta() );
	float dPhi = fabs(iPFc.phi() - eleIt->phi());
	if(dPhi > 3.14) dPhi = 3.1415926*2 - dPhi;
       
	//check if ES planes are working
	bool P1_Allok = true;
	bool P2_Allok = true;
	bool P1_Alldead = true;
	bool P2_Alldead = true;
	bool P1_ok = false;
	bool P2_ok = false;
	
	float fractionSum = 0.;
	float fractionOK = 0.;
	float highestFrac = 0.;
	int seedX = 0;
	int seedY = 0;
	int seedZ = 0;
	int seedP = 0;
	int seedS = 0;

	//std::cout << " P1_Allok = " << P1_Allok << " P2_Allok = " << P2_Allok
	// 	  << " P1_Alldead = " << P1_Alldead << " P2_Alldead = " << P2_Alldead
	// 	  << " P1_ok = " << P1_ok << " P2_ok = " << P2_ok
	// 	  << std::endl;
	

	//match gedEle wrt PFCandidate
	if(sqrt(pow(dEta,2) + pow(dPhi,2)) < 0.05){
	  //loop over PFClusters
	  for(unsigned b=0; b<iPFc.elementsInBlocks().size(); b++){
	    reco::PFBlockRef blockRef = iPFc.elementsInBlocks()[b].first;
	    unsigned elementIndex = iPFc.elementsInBlocks()[b].second;
	    if(blockRef.isNull()) continue;
	    const edm::OwnVector< reco::PFBlockElement >& elements = blockRef->elements();
	    const reco::PFBlockElement& pfbe(elements[elementIndex]); 
	    if(pfbe.type() == reco::PFBlockElement::PS1 || pfbe.type() == reco::PFBlockElement::PS2){
	      //	std::cout << " >>> PS cluster " << pfbe.type() << std::endl;
	      clustops = pfbe.clusterRef();		

	      // PFRecHits are not filled => match PFCluster with CaloCluster + extract RecHits from CaloCluster
	      std::vector< std::pair<DetId, float> > hits;
	      for(CaloCluster_iterator iES = eleIt->superCluster()->preshowerClustersBegin(); 
		  iES != eleIt->superCluster()->preshowerClustersEnd(); ++iES){
		if(clustops->energy() == (*iES)->energy() && clustops->eta() == (*iES)->eta()){
		  hits = (*iES)->hitsAndFractions();
		  break;
		}
	      }		


	      for(std::vector<std::pair<DetId,float> >::const_iterator rh = hits.begin(); rh!=hits.end(); ++rh){
		ESDetId strip1 = ESDetId((*rh).first);
		fractionSum += (*rh).second;
		if(strip1 != ESDetId(0)){
		  ESChannelStatusMap::const_iterator status = channelStatus_->getMap().find(strip1);

		  //getStatusCode() == 0 => active channel
		  // apply correction if all recHits are dead
		  if(status->getStatusCode() == 0) { //active     
		    if(pfbe.type() == reco::PFBlockElement::PS1) {
		      P1_Alldead = false;
		      P1_ok = true;
		    }
		    if(pfbe.type() == reco::PFBlockElement::PS2){
		      P2_Alldead = false;
		      P2_ok = true;
		    }
		    
		    fractionOK += (*rh).second;
		    if(highestFrac < (*rh).second){
		      highestFrac = (*rh).second;
		      
		      seedX = strip1.six();
		      seedY = strip1.siy();
		      seedZ = strip1.zside();
		      seedP = strip1.plane();
		      seedS = strip1.strip();
		    }		 
		  }
		  else{ // dead                                                                                                                                   
		    if(pfbe.type() == reco::PFBlockElement::PS1) {
		      P1_Allok = false;
		    }
		    if(pfbe.type() == reco::PFBlockElement::PS2){
		      P2_Allok = false;
		    }
		  }
		}
	      }//ES recHIts and fractions
	    } // pfClusters
	  } // blocks
	  
	  float ESenergyP1 = eleIt->superCluster()->preshowerEnergyPlane1();
	  float ESenergyP2 = eleIt->superCluster()->preshowerEnergyPlane2();
	  
	  if(P1_Allok) h_gedEle_ESenergy_plane1_Allok->Fill(ESenergyP1);
	  if(P2_Allok) h_gedEle_ESenergy_plane2_Allok->Fill(ESenergyP2);
	  if(P1_Alldead) h_gedEle_ESenergy_plane1_Alldead->Fill(ESenergyP1);
	  if(P2_Alldead) h_gedEle_ESenergy_plane2_Alldead->Fill(ESenergyP2);
	  if(P1_ok) h_gedEle_ESenergy_plane1_NOok->Fill(ESenergyP1);
	  if(P2_ok) h_gedEle_ESenergy_plane2_NOok->Fill(ESenergyP2);
	  
	  //	  std::cout << " P1_Allok = " << P1_Allok << " P2_Allok = " << P2_Allok
	  //  	    << " P1_Alldead = " << P1_Alldead << " P2_Alldead = " << P2_Alldead
	  //	  std::cout << " P1_ok = " << P1_ok << " P2_ok = " << P2_ok << std::endl;
	
	  if(P1_ok && P2_ok){
	    h_gedEle_ESenergy_2ok->Fill(ESenergyP1 + ESenergyP2);
	    h_gedEle_ESenergy_2okKill1->Fill(ESenergyP2);
	    h_gedEle_ESenergy_2okKill2->Fill(ESenergyP1);

	    h_gedEle_ESenergy_1corrected->Fill(1. / 0.0188 * 0.0434 * ESenergyP1);
	    h_gedEle_ESenergy_2corrected->Fill(1. / 0.0188 * 0.0287 * ESenergyP2);
	    if(ESenergyP1 != 0. && ESenergyP2 != 0.){
	      h_gedEle_ESenergy_1overS->Fill(ESenergyP1 / (ESenergyP1 + ESenergyP2) );	    
	      h_gedEle_ESenergy_2overS->Fill(ESenergyP2 / (ESenergyP1 + ESenergyP2) );	    
	    }
	  }

	  //	  std::cout << " fractionOK = " << fractionOK << std::endl;
	  //	  if(P1_ok == true && P1_ok != true) std::cout << << std::endl;	  

	  if(P1_ok == true) isP1ok->push_back(1);
	  else isP1ok->push_back(0);

	  if(P2_ok == true) isP2ok->push_back(1);
	  else isP2ok->push_back(0);

	  //	  if(P1_ok == true && P1_ok != true) std::cout << " >>>> SERIOOOOO " << std::endl; 

	  if(fractionSum != 0.) ele_esRhFrac->push_back(fractionOK/fractionSum);
	  else ele_esRhFrac->push_back(-100);
	  ele_esRhZ->push_back(seedZ);
	  ele_esRhP->push_back(seedP);
	  ele_esRhX->push_back(seedX);
	  ele_esRhY->push_back(seedY);
	  ele_esRhS->push_back(seedS);

	  okHaiFillato = true;
	  //	  if(int(isP1ok->size()) != nEle)  std::cout << " >>>> SERIOOOOO SIZEEEEEEEEEEEEEEEEE" << std::endl;

	}// == ele matched
      }//PFcand
    }//PS acceptance
    else{
      eleES1->push_back(-100); 
      eleES2->push_back(-100); 
    }
    
    if(!okHaiFillato){
      //      std::cout << " non ES part " << std::endl;  
      isP1ok->push_back(-100);
      isP2ok->push_back(-100);
      // eleES1->push_back(-100);  // already filled since in electron loop
      // eleES2->push_back(-100);
      ele_esRhFrac->push_back(-100);
      ele_esRhZ->push_back(-100);
      ele_esRhP->push_back(-100);
      ele_esRhX->push_back(-100);
      ele_esRhY->push_back(-100);
      ele_esRhS->push_back(-100);      
      //      if(int(isP1ok->size()) != nEle)  std::cout << " >>>> SERIOOOOO SIZEEEEEEEEEEEEEEEEE   fuori ES " << std::endl;
    } 
  
    /////////////////
  
    /*
    //////////////////  track extrapolation
    // nTrk = 0;
    // nGoodTrk = 0;
    //    reco::TrackCollection::const_iterator aTrk = eleIt->gsfTrack()->initialFreeState();
    //    FreeTrajectoryState aTrk = eleIt->gsfTrack()->initialFreeState();
    reco::GsfTrack aTrk = (*(eleIt->gsfTrack()));
    TrackDetectorAssociator::Direction direction = TrackDetectorAssociator::InsideOut;
      
    //cout<<"trk : "<<aTrk->pt()<<" "<<fabs(aTrk->vertex().rho())<<" "<<aTrk->found()<<endl;                                                                           
    TrackDetMatchInfo trk_info;
    //    trk_info = trackAssociator_.associate(ev, iSetup, *aTrk, parameters_, direction);

    trk_info = trackAssociator_.associate(ev, iSetup, aTrk, parameters_, direction);

    //position of track crossing the ES planes
    for(std::vector<DetId>::const_iterator trk_id = trk_info.crossedPreshowerIds.begin(); trk_id != trk_info.crossedPreshowerIds.end(); ++trk_id){
      //      GlobalPoint point = trk_info.getPosition(*trk_id);
      ESDetId trk_esid (trk_id->rawId()) ;
      trkESz->push_back(trk_esid.zside());
      trkESp->push_back(trk_esid.plane());
      trkESx->push_back(trk_esid.six());
      trkESy->push_back(trk_esid.siy());
      trkESs->push_back(trk_esid.strip());

      // trk_posx[nTrk][trk_esid.plane()-1]  = point.x();
      // trk_posy[nTrk][trk_esid.plane()-1]  = point.y();
      // trk_posz[nTrk][trk_esid.plane()-1]  = point.z();
    }

    for (std::vector<DetId>::const_iterator trk_eeid = trk_info.crossedEcalIds.begin(); trk_eeid != trk_info.crossedEcalIds.end(); ++trk_eeid) {
      GlobalPoint EE_point = trk_info.getPosition(*trk_eeid);
      EEDetId trk_EEid = EEDetId(trk_eeid->rawId()) ;
      trkEEz->push_back(trk_EEid.zside());
      trkEEx->push_back(trk_EEid.ix());
      trkEEy->push_back(trk_EEid.iy());
      trkEEEta->push_back(EE_point.eta());
      trkEEPhi->push_back(EE_point.phi());
    }

    //////////////////
    */

    ++nEle;
    // if(okHaiFillato == true && (int(isP1ok->size()) != nEle) ){
    //   std::cout << " >>> nEle = " << nEle << " isP1ok->size() = " << isP1ok->size() << std::endl;
    // }    
    if(int(eleES1->size()) != nEle)  std::cout << " >>> All nEle = " << nEle << " eleES1->size() = " << eleES1->size() << std::endl;
    if(isP1ok->size() != eleES1->size())  std::cout << " >>> All nEle = " << nEle << " eleES1->size() = " << eleES1->size() << " isP1ok->size() = " << isP1ok->size() << std::endl;
  }

  //  std::cout << "nEle = " << nEle << std::endl;
  // std::cout << "isZ[nEle] = " << isZ[nEle-1] << std::endl;


  //select Z candidates
  float Zmass = -1;
  if((iEle1 != electronCollection->begin() || iEle2 != electronCollection->begin()) && electronCollection->size() >= 2){
    TLorentzVector Z = e1+e2;
    Zmass = Z.M();
    //    std::cout << " @@@@@@@@@@@@>> Zmass = " << Zmass << std::endl;
    if(Zmass > 60. && Zmass < 120.) {

      isZ->at(countZ1) = 1;
      isZ->at(countZ2) = 1;
    }
  }
  
      /*
	if(electronCollection->size() >= 2 && isZ == false){
	TLorentzVector Z = e1+e2;    
	std::cout << " @@@@@@@@@@@@>> ele1 = " << e1.Pt() << " eta = " << e1.Eta() << " phi = " << e1.Phi() << " e = " << e1.Energy() << std::endl;
	std::cout << " @@@@@@@@@@@@>> ele2 = " << e2.Pt() << " eta = " << e2.Eta() << " phi = " << e2.Phi() << " e = " << e2.Energy() << std::endl;
	std::cout << " @@@@@@@@@@@@>> Zmass = " << Z.M() << std::endl;
	}
      */
      
      /*
	if(isZ == false || electronCollection->size() > 2){
	//    std::cout << " >>> isZ == false " << std::endl;
	return;
	}
      */
  
      // all clusters
      // ... endcap
  
      for(reco::SuperClusterCollection::const_iterator itSC = theEndcapSuperClusters->begin(); 
	  itSC != theEndcapSuperClusters->end(); ++itSC ) {
	
	double scEt = itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
	double scRawEt = itSC -> rawEnergy() * sin(2.*atan( exp(- itSC->position().eta() )));
	
	if (scEt < scEtThrEE_ ) continue;
	
	h_superClusters_eta       -> Fill( itSC -> eta() );
	h_superClusters_occupancyPhiEta -> Fill(itSC -> phi(),itSC -> eta());
	
	if  ( itSC -> z() > 0 ){
	  h_superClusters_EEP_energy -> Fill( itSC -> energy() );
	  h_superClusters_EEP_rawEnergy -> Fill( itSC -> rawEnergy() );
	  h_superClusters_EEP_rawEt -> Fill( scRawEt );
	}
	
	if  ( itSC -> z() < 0 ){
	  h_superClusters_EEM_energy -> Fill( itSC -> energy() );
	  h_superClusters_EEM_rawEnergy -> Fill( itSC -> rawEnergy() );
	  h_superClusters_EEM_rawEt -> Fill( scRawEt );
	}

	//--------------------------------------------------------
	//    while looping over EE supercluster	
	//Preshower RecHits
	
	if ( fabs(itSC->eta()) < 1.65 || fabs(itSC->eta()) > 2.6 ) continue;
	
	h_ESenergy->Fill(itSC->preshowerEnergy());
	h_ESenergy_plane1->Fill(itSC->preshowerEnergyPlane1());
	h_ESenergy_plane2->Fill(itSC->preshowerEnergyPlane2());
      	

	//check if ES planes are working   
	bool P1_Allok = true;
	bool P2_Allok = true;
	bool P1_Alldead = true;
	bool P2_Alldead = true;
	bool P1_ok = false;
	bool P2_ok = false;
	

	for(unsigned int icl = 0; icl < basicClusters_ES_h->size(); ++icl){
	  //Get the associated RecHits                       
	  const std::vector<std::pair<DetId,float> > & hits= (*basicClusters_ES_h)[icl].hitsAndFractions();

	  for(std::vector<std::pair<DetId,float> > ::const_iterator rh = hits.begin(); rh!=hits.end(); ++rh){
	    ESDetId rhid = ((*rh).first);
	    
	    if(rhid != ESDetId(0)){
	      ESChannelStatusMap::const_iterator status = channelStatus_->getMap().find(rhid);
	      // getStatusCode() == 1 => dead channel   
	      // apply correction if all recHits in dead region   
	      if(status->getStatusCode() == 0) { //active    
		if(rhid.plane() == 1) {
		  P1_Alldead = false;
		  P1_ok = true;
		}
		if(rhid.plane() == 2){
		  P2_Alldead = false;
		  P2_ok = true;
		}
	      }
	      else{ // dead   
		if(rhid.plane() == 1) {
		  P1_Allok = false;
		}
		if(rhid.plane() == 2){
		  P2_Allok = false;
		}
	      }
	    }
	  } // rh                  
	}
	
	if(P1_Allok) h_ESenergy_plane1_Allok->Fill(itSC->preshowerEnergyPlane1());
	if(P2_Allok) h_ESenergy_plane2_Allok->Fill(itSC->preshowerEnergyPlane2());
	if(P1_Alldead) h_ESenergy_plane1_Alldead->Fill(itSC->preshowerEnergyPlane1());
	if(P2_Alldead) h_ESenergy_plane2_Alldead->Fill(itSC->preshowerEnergyPlane2());
	if(P1_ok) h_ESenergy_plane1_NOok->Fill(itSC->preshowerEnergyPlane1());
	if(P2_ok) h_ESenergy_plane2_NOok->Fill(itSC->preshowerEnergyPlane2());
	
      }// end loop over EE superclusters


      // std::cout << " nEle = " << nEle                                                                           
      // 		<< " isZ = " << isZ->size()                                                                
      // 		<< " eleEta = " << eleEta->size() << " elePhi = " << elePhi->size()                        
      // 		<< " eleP = " << eleP->size() << " elePt = " << elePt->size()                              
      // 		<< " eleEnergy = " << eleEnergy->size() << " elefBrem = " << elefBrem->size()              
      // 		<< " eleSCBrem = " << eleSCBrem->size() << " ele5x5 = " << ele5x5->size()                  
      // 		<< " scEnergy = " << scEnergy->size() << " scRawEnergy = " << scRawEnergy->size()          
      // 		<< " scEta = " << scEta->size() << " scPhi = " << scPhi->size()                            
      // 		<< " deltaEtaEleClusterAtCalo = " << deltaEtaEleClusterAtCalo->size()                      
      // 		<< " deltaPhiEleClusterAtCalo = " << deltaPhiEleClusterAtCalo->size()                      
      // 		<< " deltaEtaEleClusterAtVtx = " << deltaEtaEleClusterAtVtx->size()                        
      // 		<< " deltaPhiEleClusterAtVtx = " << deltaPhiEleClusterAtVtx->size()                        
      // 		<< " eleES1 = " << eleES1->size() << " eleES2 = " << eleES2->size()                        
      // 		<< " ele_esRhFrac = " << ele_esRhFrac->size()                                              
      // 		<< " ele_esRhZ = " << ele_esRhZ->size() << " ele_esRhP = " << ele_esRhP->size()            
      // 		<< " ele_esRhX = " << ele_esRhX->size() << " ele_esRhY = " << ele_esRhY->size()            
      // 		<< " ele_esRhS = " << ele_esRhS->size()                                                    
      // 		<< " isP1ok = " << isP1ok->size() << " isP2ok = " << isP2ok->size()
      // 		<< std::endl;                                                                              


	// std::cout << " >>> isP1ok->size() = " << isP1ok->size() << std::endl;
	// std::cout << " >>> isP2ok->size() = " << isP2ok->size() << std::endl;
      
      // if(int(isP1ok->size()) != nEle) {
      // 	std::cout << " >>> problema " << std::endl;
      // 	std::cout << " >>> isP1ok->size() = " << isP1ok->size() << std::endl;
      // 	std::cout << " >>> ele_esRhFrac->size() = " << ele_esRhFrac->size() << std::endl;
      // 	std::cout << " >>> nEle = " << nEle << std::endl;
      // }

      newT->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
        void 
EcalValidationDigiless_ES::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EcalValidationDigiless_ES::endJob() {
  h_numberOfEvents ->Fill(0.,naiveId_);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EcalValidationDigiless_ES);
