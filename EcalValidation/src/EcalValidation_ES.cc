// -*- C++ -*-
//
// Package:    EcalValidation_ES
// Class:      EcalValidation_ES
// Original Author:  Martina Malberti
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

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
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
#include "CondFormats/ESObjects/interface/ESPedestals.h"
#include "CondFormats/DataRecord/interface/ESPedestalsRcd.h"

#include "Validation/EcalValidation/interface/EcalValidation_ES.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

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
EcalValidation_ES::EcalValidation_ES(const edm::ParameterSet& ps)
{
  //now do what ever initialization is needed
  recHitCollection_EB_       = ps.getParameter<edm::InputTag>("recHitCollection_EB");
  recHitCollection_EE_       = ps.getParameter<edm::InputTag>("recHitCollection_EE");
  recHitCollection_ES_       = ps.getParameter<edm::InputTag>("recHitCollection_ES");
  redRecHitCollection_ES_    = ps.getParameter<edm::InputTag>("redRecHitCollection_ES");
  esClusterCollectionX_      = ps.getParameter<edm::InputTag>("ClusterCollectionX_ES");
  esClusterCollectionY_      = ps.getParameter<edm::InputTag>("ClusterCollectionY_ES");
  //esDigiCollection_          = ps.getParameter<edm::InputTag>("esDigiCollection");

  EleTag_                     = ps.getParameter<edm::InputTag>("EleTag");
  PhoTag_                     = ps.getParameter<edm::InputTag>("PhoTag");
  saveDigis_                  = ps.getUntrackedParameter<bool>("saveDigis",   false);

  // histos
  edm::Service<TFileService> fs;
  h_numberOfEvents = fs->make<TH1D>("h_numberOfEvents","h_numberOfEvents",10,0,10);
  
  // PulseShape ----------------------------------------------
  h_PulseShape_EB               =  fs->make<TProfile>("h_PulseShape_EB","h_PulseShape_EB",20,-2,18);
  h_PulseShape_EE               =  fs->make<TProfile>("h_PulseShape_EE","h_PulseShape_EE",20,-2,18);
  h_PulseShape_ES               =  fs->make<TProfile>("h_PulseShape_ES","h_PulseShape_ES",20,-2,18);
  h_PulseShape_ESpedSub               =  fs->make<TProfile>("h_PulseShape_ESpedSub","h_PulseShape_ESpedSub",20,-2,18);
  h_PulseShape_ES_inTime        =  fs->make<TProfile>("h_PulseShape_ES_inTime","h_PulseShape_ES_inTime",20,-2,18);
  h_PulseShape_ES_eOOT          =  fs->make<TProfile>("h_PulseShape_ES_eOOT","h_PulseShape_ES_eOOT",20,-2,18);
  h_PulseShape_ES_lOOT          =  fs->make<TProfile>("h_PulseShape_ES_lOOT","h_PulseShape_ES_lOOT",20,-2,18);

  h_PulseShape_ES_kGood         =  fs->make<TProfile>("h_PulseShape_ES_kGood","h_PulseShape_ES_kGood",20,-2,18);
  h_PulseShape_ES_kGoodpedSub         =  fs->make<TProfile>("h_PulseShape_ES_kGoodpedSub","h_PulseShape_ES_kGoodpedSub",20,-2,18);
  h_PulseShape_ES_inTime_kGood        =  fs->make<TProfile>("h_PulseShape_ES_inTime_kGood","h_PulseShape_ES_inTime_kGood",20,-2,18);
  h_PulseShape_ES_eOOT_kGood          =  fs->make<TProfile>("h_PulseShape_ES_eOOT_kGood","h_PulseShape_ES_eOOT_kGood",20,-2,18);
  h_PulseShape_ES_lOOT_kGood          =  fs->make<TProfile>("h_PulseShape_ES_lOOT_kGood","h_PulseShape_ES_lOOT_kGood",20,-2,18);

  // preshower
  h_recHits_ES_size           = fs->make<TH1D>("h_recHits_ES_size","h_recHits_ES_size",1000,0.,10000);
  h_recHits_ES_size_F[0]      = fs->make<TH1D>("h_recHits_ES_size_F+","h_recHits_ES_size_F+",1000,0.,10000);
  h_recHits_ES_size_F[1]      = fs->make<TH1D>("h_recHits_ES_size_F-","h_recHits_ES_size_F-",1000,0.,10000);
  h_recHits_ES_size_R[0]      = fs->make<TH1D>("h_recHits_ES_size_R+","h_recHits_ES_size_R+",1000,0.,10000);
  h_recHits_ES_size_R[1]      = fs->make<TH1D>("h_recHits_ES_size_R-","h_recHits_ES_size_R-",1000,0.,10000);

  h_recHits_ES_energy         = fs->make<TH1D>("h_recHits_ES_energy","h_recHits_ES_energy",1000,0.,0.01);
  h_recHits_ES_energy_F[0]    = fs->make<TH1D>("h_recHits_ES_energy_F+","h_recHits_ES_energy_F+",1000,0.,0.01);
  h_recHits_ES_energy_F[1]    = fs->make<TH1D>("h_recHits_ES_energy_F-","h_recHits_ES_energy_F-",1000,0.,0.01);
  h_recHits_ES_energy_R[0]    = fs->make<TH1D>("h_recHits_ES_energy_R+","h_recHits_ES_energy_R+",1000,0.,0.01);
  h_recHits_ES_energy_R[1]    = fs->make<TH1D>("h_recHits_ES_energy_R-","h_recHits_ES_energy_R-",1000,0.,0.01);

  h_recHits_ES_energyMax      = fs->make<TH1D>("h_recHits_ES_energyMax","h_recHits_ES_energyMax",1000,0.,0.01);
  h_recHits_ES_energyMax_F[0] = fs->make<TH1D>("h_recHits_ES_energyMax_F+","h_recHits_ES_energyMax_F+",1000,0.,0.01);
  h_recHits_ES_energyMax_F[1] = fs->make<TH1D>("h_recHits_ES_energyMax_F-","h_recHits_ES_energyMax_F-",1000,0.,0.01);
  h_recHits_ES_energyMax_R[0] = fs->make<TH1D>("h_recHits_ES_energyMax_R+","h_recHits_ES_energyMax_R+",1000,0.,0.01);
  h_recHits_ES_energyMax_R[1] = fs->make<TH1D>("h_recHits_ES_energyMax_R-","h_recHits_ES_energyMax_R-",1000,0.,0.01);

  h_recHits_ES_time           = fs->make<TH1D>("h_recHits_ES_time","h_recHits_ES_time",400,-100.,100.);
  h_recHits_ES_time_F[0]      = fs->make<TH1D>("h_recHits_ES_time_F+","h_recHits_ES_time_F+",400,-100.,100.);
  h_recHits_ES_time_F[1]      = fs->make<TH1D>("h_recHits_ES_time_F-","h_recHits_ES_time_F-",400,-100.,100.);
  h_recHits_ES_time_R[0]      = fs->make<TH1D>("h_recHits_ES_time_R+","h_recHits_ES_time_R+",400,-100.,100.);
  h_recHits_ES_time_R[1]      = fs->make<TH1D>("h_recHits_ES_time_R-","h_recHits_ES_time_R-",400,-100.,100.);

  h_esClusters_energy_plane1 = fs->make<TH1D>("h_esClusters_energy_plane1","h_esClusters_energy_plane1",1000,0.,0.01);
  h_esClusters_energy_plane2 = fs->make<TH1D>("h_esClusters_energy_plane2","h_esClusters_energy_plane2",1000,0.,0.01);
  h_esClusters_energy_ratio  = fs->make<TH1D>("h_esClusters_energy_ratio","h_esClusters_energy_ratio",100,0.,10.);

  h_recHits_ES_recoFlag = fs->make<TH1D>("h_recHits_ES_recoFlag","h_recHits_ES_recoFlag",16,-0.5,15.5);  
  //  h_digiShape_ES_recoFlag = fs->make<TH1D>("h_digiShape_ES_recoFlag","h_digiShape_ES_recoFlag",20,-2,18);

  h_PulseShape_ES_kESGood = fs->make<TH1D>("h_PulseShape_ES_kESGood","h_PulseShape_ES_kESGood",20,-2,18);
  h_PulseShape_ES_kESDead = fs->make<TH1D>("h_PulseShape_ES_kESDead","h_PulseShape_ES_kESDead",20,-2,18);
  h_PulseShape_ES_kESHot = fs->make<TH1D>("h_PulseShape_ES_kESHot","h_PulseShape_ES_kESHot",20,-2,18);
  h_PulseShape_ES_kESPassBX = fs->make<TH1D>("h_PulseShape_ES_kESPassBX","h_PulseShape_ES_kESPassBX",20,-2,18);
  h_PulseShape_ES_kESTwoGoodRatios = fs->make<TH1D>("h_PulseShape_ES_kESTwoGoodRatios","h_PulseShape_ES_kESTwoGoodRatios",20,-2,18);
  h_PulseShape_ES_kESBadRatioFor12 = fs->make<TH1D>("h_PulseShape_ES_kESBadRatioFor12","h_PulseShape_ES_kESBadRatioFor12",20,-2,18);
  h_PulseShape_ES_kESBadRatioFor23Upper = fs->make<TH1D>("h_PulseShape_ES_kESBadRatioFor23Upper","h_PulseShape_ES_kESBadRatioFor23Upper",20,-2,18);
  h_PulseShape_ES_kESBadRatioFor23Lower = fs->make<TH1D>("h_PulseShape_ES_kESBadRatioFor23Lower","h_PulseShape_ES_kESBadRatioFor23Lower",20,-2,18);
  h_PulseShape_ES_kESTS1Largest = fs->make<TH1D>("h_PulseShape_ES_kESTS1Largest","h_PulseShape_ES_kESTS1Largest",20,-2,18);
  h_PulseShape_ES_kESTS3Largest = fs->make<TH1D>("h_PulseShape_ES_kESTS3Largest","h_PulseShape_ES_kESTS3Largest",20,-2,18);
  h_PulseShape_ES_kESTS3Negative = fs->make<TH1D>("h_PulseShape_ES_kESTS3Negative","h_PulseShape_ES_kESTS3Negative",20,-2,18);
  h_PulseShape_ES_kESSaturated = fs->make<TH1D>("h_PulseShape_ES_kESSaturated","h_PulseShape_ES_kESSaturated",20,-2,18);
  h_PulseShape_ES_kESTS2Saturated = fs->make<TH1D>("h_PulseShape_ES_kESTS2Saturated","h_PulseShape_ES_kESTS2Saturated",20,-2,18);
  h_PulseShape_ES_kESTS3Saturated = fs->make<TH1D>("h_PulseShape_ES_kESTS3Saturated","h_PulseShape_ES_kESTS3Saturated",20,-2,18);
  h_PulseShape_ES_kESTS13Sigmas = fs->make<TH1D>("h_PulseShape_ES_kESTS13Sigmas","h_PulseShape_ES_kESTS13Sigmas",20,-2,18);
  h_PulseShape_ES_kESTS15Sigmas = fs->make<TH1D>("h_PulseShape_ES_kESTS15Sigmas","h_PulseShape_ES_kESTS15Sigmas",20,-2,18);

  for(int i=0; i<100; ++i)  h_PulseShape_ES_pedSub[i] = fs->make<TH1D>(Form("h_PulseShape_ES_pedSub_%d",i),Form("h_PulseShape_ES_pedSub_%d",i),3,0,3);

  h_PulseShape_ES_flag = fs->make<TH1D>("h_PulseShape_ES_flag","h_PulseShape_ES_flag",20,-2,18);

  h_CentralStrip_P1_E = fs->make<TH1D>("h_CentralStrip_P1_E","h_CentralStrip_P1_E", 1000, 0., 0.05);
  h_CentralStrip_P2_E = fs->make<TH1D>("h_CentralStrip_P2_E","h_CentralStrip_P2_E", 1000, 0., 0.05);


  for(int i=0; i<31; ++i){
    plane1_stripsN_E[i] = fs->make<TH1D>(Form("plane1_stripsN_E_%d",i), "", 1000, 0., 0.05);
    plane2_stripsN_E[i] = fs->make<TH1D>(Form("plane2_stripsN_E_%d",i), "", 1000, 0., 0.05);
  }

}



EcalValidation_ES::~EcalValidation_ES()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void EcalValidation_ES::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  //*********** CALO GEOMETRY                                                                                                                                                                                           
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);
  const CaloGeometry *geometry = pGeometry.product();

  const CaloSubdetectorGeometry *geometryES = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  CaloSubdetectorTopology *topologyES = 0;
  if(geometryES) topologyES = new EcalPreshowerTopology(geometry);

  //Preshower RecHits and digis  
  edm::ESHandle<ESPedestals> ESPed;
  iSetup.get<ESPedestalsRcd>().get( ESPed );
  const ESPedestals *peds = ESPed.product();

  edm::Handle<ESRecHitCollection> recHitsES;
  ev.getByLabel (recHitCollection_ES_, recHitsES) ;
  const ESRecHitCollection* thePreShowerRecHits = recHitsES.product () ;

  if ( ! recHitsES.isValid() ) {
    std::cerr << "EcalValidation::analyze --> recHitsES not found" << std::endl;
  }

  h_recHits_ES_size -> Fill( recHitsES->size());


  //*********** ES REC HITS 
  ESrechits_map.clear();
  edm::Handle<EcalRecHitCollection> redrecHitsES;
  ev.getByLabel(redRecHitCollection_ES_, redrecHitsES);
  const EcalRecHitCollection* theESEcalRedRecHits = redrecHitsES.product () ;
  //  if(verbosity_)  std::cout << " theESEcalRedRecHits->size() =  " << theESEcalRedRecHits->size() << std::endl;
  if ( ! redrecHitsES.isValid() ) {
    std::cerr << "SimpleNtupleEoverP::analyze --> redrecHitsES not found" << std::endl;
  }
  else{
    // make the map of rechits
    EcalRecHitCollection::const_iterator it;
    for(it = theESEcalRedRecHits->begin(); it != theESEcalRedRecHits->end(); ++it) {
      // remove bad ES rechits  
      //      std::cout << " it->recoFlag() = " << it->recoFlag() << std::endl;
      if(it->recoFlag() == 1 || it->recoFlag() == 14 || (it->recoFlag() <= 10 && it->recoFlag() >= 5)) continue;
      //Make the map of DetID, EcalRecHit pairs 
      ESrechits_map.insert(std::make_pair(it->id(), *it));
      // std::cout << " it->energy() = " << it->energy() << std::endl;
    }
  }


  // if(saveDigis_){
  //   edm::Handle<ESDigiCollection> esDigis;
  //   ev.getByLabel (esDigiCollection_, esDigis) ;
  //   const ESDigiCollection* theEcalPreShowerDigis = esDigis.product () ;
  //   if (! (esDigis.isValid ()) ) {
  //     std::cerr << "EcalValidation_ES::analyze -->  esDigis not found" << std::endl;
  //   }
  // }
  

  float maxRecHitEnergyES = -999.;

  int nF[2]={0,0};
  int nR[2]={0,0};

  for (ESRecHitCollection::const_iterator esItr = thePreShowerRecHits->begin(); esItr != thePreShowerRecHits->end(); ++esItr)
    {

      h_recHits_ES_energy -> Fill(esItr->energy());
      h_recHits_ES_time   -> Fill(esItr->time());
      if (esItr -> energy() > maxRecHitEnergyES ) maxRecHitEnergyES = esItr -> energy() ;

      h_recHits_ES_recoFlag->Fill(esItr->recoFlag());

      ESDetId id = ESDetId(esItr->id());
      // front plane : id.plane()==1                                                                                                                       
      if ( id.plane()==1 && id.zside() > 0 ){
        h_recHits_ES_energy_F[0]-> Fill( esItr->energy() );
        h_recHits_ES_time_F[0]  -> Fill( esItr->time() );
        nF[0]++;
      }
      if ( id.plane()==1 && id.zside() < 0 ){
        h_recHits_ES_energy_F[1]-> Fill( esItr->energy() );
        h_recHits_ES_time_F[1]  -> Fill( esItr->time() );
        nF[1]++;
      }
      // rear plane : id.plane()==2 
      if ( id.plane()==2 && id.zside() > 0 ){
        h_recHits_ES_energy_R[0]-> Fill( esItr->energy() );
        h_recHits_ES_time_R[0]  -> Fill( esItr->time() );
        nR[0]++;
      }
      if ( id.plane()==2 && id.zside() < 0 ){
        h_recHits_ES_energy_R[1]->Fill( esItr->energy() );
        h_recHits_ES_time_R[1]  -> Fill( esItr->time() );
        nR[1]++;
      }
    
      if(saveDigis_){
    edm::Handle<ESDigiCollection> esDigis;
    ev.getByLabel (esDigiCollection_, esDigis) ;
    const ESDigiCollection* theEcalPreShowerDigis = esDigis.product () ;
    if (! (esDigis.isValid ()) ) {
      std::cerr << "EcalValidation_ES::analyze -->  esDigis not found" << std::endl;
    }
    

      //... digis ES from recHit  
      for(ESDigiCollection::const_iterator digiItr = theEcalPreShowerDigis->begin();
          digiItr != theEcalPreShowerDigis->end();
          ++digiItr)
	{
          if( digiItr->id() != id )continue;
          EcalDataFrame df = *digiItr;

	  ESPedestals::const_iterator it_ped = peds->find(digiItr->id());
	  float pedValue = it_ped->getMean();

	  //	  if(esItr->recoFlag() == EcalRecHit::kESGood) std::cout << " esItr->recoFlag() = " << esItr->recoFlag() << std::endl;

	  //	  if(esItr->recoFlag() == EcalRecHit::kESGood){
	  if(esItr->recoFlag() == 0){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kGoodpedSub->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_kGood->Fill(i, df.sample(i).adc());

	      h_PulseShape_ES_kESGood->Fill(i, df.sample(i).adc());
	      h_PulseShape_ES_flag->SetBinContent(1, esItr->recoFlag());

	      //inTime scenario  
	      if(df.sample(1).adc() > df.sample(0).adc() && df.sample(1).adc() > df.sample(2).adc()){
		for (int i=0; i < df.size(); i++ ){
		  h_PulseShape_ES_inTime_kGood->Fill(i, df.sample(i).adc() - pedValue);
		}
	      }
	      //early ootTime scenario
	      if(df.sample(0).adc() > df.sample(1).adc() && df.sample(1).adc() > df.sample(2).adc()){
		for (int i=0; i < df.size(); i++ ){
		  h_PulseShape_ES_eOOT_kGood->Fill(i, df.sample(i).adc() - pedValue);
		}
	      }

	      //late ootTime scenario
	      if(df.sample(2).adc() > df.sample(1).adc() && df.sample(1).adc() > df.sample(0).adc()){
		for (int i=0; i < df.size(); i++ ){
		  h_PulseShape_ES_lOOT_kGood->Fill(i, df.sample(i).adc() - pedValue);
		}
	      }

	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESDead){
	  if(esItr->recoFlag() == 1){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESDead->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(2, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESHot){
	  if(esItr->recoFlag() == 2){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESHot->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(3, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESPassBX){
	  if(esItr->recoFlag() == 3){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESPassBX->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(4, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESTwoGoodRatios){
	  if(esItr->recoFlag() == 4){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESTwoGoodRatios->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(5, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESBadRatioFor12){
	  if(esItr->recoFlag() == 5){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESBadRatioFor12->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(6, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESBadRatioFor23Upper){
	  if(esItr->recoFlag() == 6){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESBadRatioFor23Upper->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(7, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESBadRatioFor23Lower){
	  if(esItr->recoFlag() == 7){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESBadRatioFor23Lower->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(8, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESTS1Largest){
	  if(esItr->recoFlag() == 8){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESTS1Largest->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(9, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESTS3Largest){
	  if(esItr->recoFlag() == 9){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESTS3Largest->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(10, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESTS3Negative){
	  if(esItr->recoFlag() == 10){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESTS3Negative->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(11, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESSaturated){
	  if(esItr->recoFlag() == 11){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESSaturated->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(12, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESTS2Saturated){
	  if(esItr->recoFlag() == 12){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESTS2Saturated->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(13, esItr->recoFlag());
	    }
	  }

	  //	  if(esItr->recoFlag() == EcalRecHit::kESTS3Saturated){
	  if(esItr->recoFlag() == 13){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESTS3Saturated->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(14, esItr->recoFlag());
	    }
	  }


	  //	  if(esItr->recoFlag() == EcalRecHit::kESTS13Sigmas){
	  if(esItr->recoFlag() == 14){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESTS13Sigmas->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(15, esItr->recoFlag());
	    }
	  }


	  //	  if(esItr->recoFlag() == EcalRecHit::kESTS15Sigmas){
	  if(esItr->recoFlag() == 15){
	    for (int i=0; i < df.size(); i++ ){
	      h_PulseShape_ES_kESTS15Sigmas->Fill(i, df.sample(i).adc() - pedValue);
	      h_PulseShape_ES_flag->SetBinContent(16, esItr->recoFlag());
	    }
	  }
	}
      }//if saveDigis
    } // end loop over ES rec Hits                                                                                                                         

  h_recHits_ES_energyMax -> Fill(maxRecHitEnergyES );
  h_recHits_ES_size_F[0] -> Fill( nF[0] );
  h_recHits_ES_size_F[1] -> Fill( nF[1] );
  h_recHits_ES_size_R[0] -> Fill( nR[0] );
  h_recHits_ES_size_R[1] -> Fill( nR[1] );

  // ALL ES digis (not from recHits)
  int nEvt = 0;  
  if(saveDigis_){
    edm::Handle<ESDigiCollection> esDigis;
    ev.getByLabel (esDigiCollection_, esDigis) ;
    const ESDigiCollection* theEcalPreShowerDigis = esDigis.product () ;
    if (! (esDigis.isValid ()) ) {
      std::cerr << "EcalValidation_ES::analyze -->  esDigis not found" << std::endl;
    }

  for(ESDigiCollection::const_iterator digiItr = theEcalPreShowerDigis->begin();
      digiItr != theEcalPreShowerDigis->end(); ++digiItr){
    
    EcalDataFrame df = *digiItr;
    
    ESPedestals::const_iterator it_ped = peds->find(digiItr->id());
    float pedValue = it_ped->getMean();

    for (int i=0; i < df.size(); i++ ){
      h_PulseShape_ES->Fill(i, df.sample(i).adc());
      h_PulseShape_ESpedSub->Fill(i, df.sample(i).adc() - pedValue);
      if(nEvt < 100) h_PulseShape_ES_pedSub[nEvt]->SetBinContent(i+1, df.sample(i).adc() - pedValue);
    }
    ++nEvt;

    //inTime scenario
    if(df.sample(1).adc() > df.sample(0).adc() && df.sample(1).adc() > df.sample(2).adc()){
      for (int i=0; i < df.size(); i++ ){
	h_PulseShape_ES_inTime->Fill(i, df.sample(i).adc() - pedValue);
      }
    }
    //early ootTime scenario
    if(df.sample(0).adc() > df.sample(1).adc() && df.sample(1).adc() > df.sample(2).adc()){
      for (int i=0; i < df.size(); i++ ){
	h_PulseShape_ES_eOOT->Fill(i, df.sample(i).adc() - pedValue);
      }
    }

    //late ootTime scenario
    if(df.sample(2).adc() > df.sample(1).adc() && df.sample(1).adc() > df.sample(0).adc()){
      for (int i=0; i < df.size(); i++ ){
	h_PulseShape_ES_lOOT->Fill(i, df.sample(i).adc() - pedValue);
      }
    }

  }
  }//savedigis
 
  /////////// ES saturation 
  //************* ELECTRONS                                                                                                                                                               
  //Handle<View<reco::GsfElectron> > electronHandle;
  //ev.getByLabel(EleTag_,electronHandle);
  //View<reco::GsfElectron> electrons = *electronHandle;
  Handle<reco::GsfElectronCollection> electronHandle;
  ev.getByLabel(EleTag_, electronHandle);

  //************* CLUSTER LAZY TOOLS                                                                                                                             
  EcalClusterLazyTools lazyTools(ev,iSetup,recHitCollection_EB_,recHitCollection_EE_);

  const reco::GsfElectronCollection * electronCollection = electronHandle.product();

  //  for(reco::GsfElectronCollection::const_iterator eleIt=electrons.begin(); eleIt!=electrons.end(); ++eleIt){
  //  for(reco::GsfElectronCollection::const_iterator eleIt=electronHandle->begin(); eleIt!=electronHandle->end(); ++eleIt){
  for(reco::GsfElectronCollection::const_iterator eleIt=electronCollection->begin(); eleIt!=electronCollection->end(); eleIt++) {
    
    reco::SuperClusterRef sClRef = eleIt->superCluster();
    

    // ES variables   
    std::vector<float>    ele1_ESHits_plane1 = lazyTools.getESHits(sClRef->x(), sClRef->y(), sClRef->z(), ESrechits_map, geometry, topologyES, 0, 1);
    std::vector<float>    ele1_ESHits_plane2 = lazyTools.getESHits(sClRef->x(), sClRef->y(), sClRef->z(), ESrechits_map, geometry, topologyES, 0, 2);

    //    std::cout << " sClRef->y() = " << sClRef->y() << std::endl;
    //    std::cout << " sClRef->energy() = " << sClRef->energy() << std::endl;
    //    std::cout << " ele1_ESHits_plane1.at(0) = " << ele1_ESHits_plane1.at(0) << std::endl;
    if(ele1_ESHits_plane1.at(0) != 0.)
      h_CentralStrip_P1_E->Fill(ele1_ESHits_plane1.at(0));
    if(ele1_ESHits_plane2.at(0) != 0.)
    h_CentralStrip_P2_E->Fill(ele1_ESHits_plane2.at(0));

    for(unsigned int itr = 0; itr<ele1_ESHits_plane1.size(); ++itr){
      if(ele1_ESHits_plane1.at(itr) != 0.){
	plane1_stripsN_E[itr]->Fill(ele1_ESHits_plane1.at(itr));
      }
    }
    for(unsigned int itr = 0; itr<ele1_ESHits_plane2.size(); ++itr){
      if(ele1_ESHits_plane2.at(itr) != 0.){
	plane2_stripsN_E[itr]->Fill(ele1_ESHits_plane2.at(itr));
      }
    }
      
    if(sClRef->preshowerEnergy() != 0.) std::cout << " >>> energy >> 0 ######################## = " << sClRef->preshowerEnergy() << std::endl;
    //    h_ESenergy->Fill(sClRef->preshowerEnergy());
  }

}
   
  //--------------------------------------------------------
 


  // ---------- Do histos for pi0 peak

//   doPi0Barrel(geometry, topology, recHitsEB);
//   doPi0Endcap(geometry, topology, recHitsEE);

void 
EcalValidation_ES::beginJob()
{
  std::cout << " CIAO " << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EcalValidation_ES::endJob() {
  
  h_numberOfEvents->Fill(0.,naiveId_);
  
  // int nBins_noise = 30; 
  // TFile* noise_histos = new TFile("noise_histos.root","RECREATE");
  // noise_histos->cd();
  // noise_histos->Close();
   
  /// ---------- Compute and Fill RecHits occupancy deviations  





}



// ----------additional functions-------------------

void EcalValidation_ES::convxtalid(Int_t &nphi,Int_t &neta)
{
  // Barrel only
  // Output nphi 0...359; neta 0...84; nside=+1 (for eta>0), or 0 (for eta<0).
  // neta will be [-85,-1] , or [0,84], the minus sign indicates the z<0 side.
  
  if(neta > 0) neta -= 1;
  if(nphi > 359) nphi=nphi-360;
  
} //end of convxtalid

int EcalValidation_ES::diff_neta_s(Int_t neta1, Int_t neta2){
  Int_t mdiff;
  mdiff=(neta1-neta2);
  return mdiff;
}

// Calculate the distance in xtals taking into account the periodicity of the Barrel
int EcalValidation_ES::diff_nphi_s(Int_t nphi1,Int_t nphi2) {
  Int_t mdiff;
  if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) {
    mdiff=nphi1-nphi2;
  }
  else {
    mdiff=360-abs(nphi1-nphi2);
    if(nphi1>nphi2) mdiff=-mdiff;
  }
  return mdiff;
}

//define this as a plug-in
DEFINE_FWK_MODULE(EcalValidation_ES);
