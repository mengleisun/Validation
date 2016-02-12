#ifndef EcalValidationDigiless_ES_h
#define EcalValidationDigiless_ES_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "CondFormats/ESObjects/interface/ESChannelStatus.h"


// ROOT include
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <map>


// Less than operator for sorting EcalRecHits according to energy.
class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool>
{
public:
  bool operator()(EcalRecHit x, EcalRecHit y)
  {
    return (x.energy() > y.energy());
  }
};


class EcalValidationDigiless_ES : public edm::EDAnalyzer {
  
      public:
         explicit EcalValidationDigiless_ES(const edm::ParameterSet&);
	 ~EcalValidationDigiless_ES();
  
  
      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
	 virtual void endJob() ;

	 // ----------additional functions-------------------
	 void convxtalid(int & , int &);
	 int diff_neta_s(int,int);
	 int diff_nphi_s(int,int);

	 const ESChannelStatus* channelStatus_;

      protected:
	 void doPi0Barrel(const CaloGeometry *theCaloGeometry,
			  const CaloTopology *theCaloTopology,
			  edm::Handle<EcalRecHitCollection> recHitsEB);

	 void doPi0Endcap(const CaloGeometry *theCaloGeometry,
			  const CaloTopology *theCaloTopology,
			  edm::Handle<EcalRecHitCollection> recHitsEE);
  

	 // ----------member data ---------------------------
	 edm::EDGetTokenT<reco::VertexCollection> PVTag_;
	 //         edm::InputTag PVTag_;
	 edm::InputTag recHitCollection_EB_;
	 edm::InputTag recHitCollection_EE_;
         edm::InputTag redRecHitCollection_EB_;
         edm::InputTag redRecHitCollection_EE_;
         edm::InputTag basicClusterCollection_EB_;
	 edm::InputTag basicClusterCollection_EE_;
	 edm::InputTag superClusterCollection_EB_;
	 edm::EDGetTokenT<reco::SuperClusterCollection> superClusterCollection_EE_;
	 edm::InputTag esRecHitCollection_;
	 edm::InputTag esClusterCollectionX_ ;
	 edm::InputTag esClusterCollectionY_ ;
	 edm::EDGetTokenT<reco::BasicClusterCollection> basicClusterCollection_ES_;


	 edm::EDGetTokenT<reco::GsfElectronCollection> elePFlow_ ;
	 edm::EDGetTokenT<reco::PFCandidateCollection> elePFlowCand_ ;
	 edm::EDGetTokenT<reco::PFRecHitCollection> esPFRec_;

	 edm::InputTag beamSpot_ ;
         edm::InputTag jets_;
	 
	 double ethrEB_;
	 double ethrEE_;
         double gainId_;

	 double scEtThrEB_;
	 double scEtThrEE_;
         


	 // ------------- HISTOGRAMS ------------------------------------

	 int naiveId_;
	 
	 TH1D *h_numberOfEvents;
         TH1D* h_PV_n;
         TH1D* h_PV_cut_n;
         
	 // ReducedRecHits ----------------------------------------------
	 // ... barrel 
         TH1D *h_redRecHits_EB_recoFlag;
	 // ... endcap 
         TH1D *h_redRecHits_EE_recoFlag;
	 // ... all 
         TH1D *h_redRecHits_recoFlag;
	 
	 // RecHits ----------------------------------------------
         
	 // Super Clusters ----------------------------------------------
	 
	 // ... endcap
	 TH1D *h_superClusters_EEP_energy;
	 TH1D *h_superClusters_EEP_rawEnergy;
	 TH1D *h_superClusters_EEP_rawEt;
         TH2D *h_superClusters_EEP_occupancy;
	 
	 TH1D *h_superClusters_EEM_energy;
	 TH1D *h_superClusters_EEM_rawEnergy;
	 TH1D *h_superClusters_EEM_rawEt;
         TH2D *h_superClusters_EEM_occupancy;
	 
         TH2D *h_superClusters_occupancyPhiEta;
         
	 TH1D *h_superClusters_eta;
	 
	 // PRESHOWER ----------------------------------------------	 
	 TH1D* h_gedEle_ESenergy;
	 TH1D* h_gedEle_ESenergy_plane1;
	 TH1D* h_gedEle_ESenergy_plane2;

	 TH1D* h_gedEle_ESenergy_plane1_Allok;
	 TH1D* h_gedEle_ESenergy_plane1_Alldead;
	 TH1D* h_gedEle_ESenergy_plane1_NOok;
	 TH1D* h_gedEle_ESenergy_plane2_Allok;
	 TH1D* h_gedEle_ESenergy_plane2_Alldead;
	 TH1D* h_gedEle_ESenergy_plane2_NOok;

	 TH1D* h_gedEle_ESenergy_2ok;
	 TH1D* h_gedEle_ESenergy_2okKill1;
	 TH1D* h_gedEle_ESenergy_2okKill2;
	 TH1D* h_gedEle_ESenergy_1corrected;
	 TH1D* h_gedEle_ESenergy_2corrected;
	 TH1D* h_gedEle_ESenergy_1overS;
	 TH1D* h_gedEle_ESenergy_2overS;


	 TH1D* h_ESenergy;
	 TH1D* h_ESenergy_plane1;
	 TH1D* h_ESenergy_plane2;


	 TH1D* h_ESenergy_plane1_Allok;
	 TH1D* h_ESenergy_plane1_Alldead;
	 TH1D* h_ESenergy_plane1_NOok;
	 TH1D* h_ESenergy_plane2_Allok;
	 TH1D* h_ESenergy_plane2_Alldead;
	 TH1D* h_ESenergy_plane2_NOok;


	 TTree* newT;
	 int event;
	 int run;
	 int lumi;
	 int nEle;
	 std::vector<int>* isZ;
	 std::vector<float>* eleEta;
	 std::vector<float>* elePhi;
	 std::vector<float>* elePt;
	 std::vector<float>* eleP;
	 std::vector<float>* eleEnergy;
	 std::vector<float>* elefBrem;
	 std::vector<float>* eletrackBrem;
	 std::vector<float>* eleSCBrem;
	 std::vector<float>* ele5x5;
	 std::vector<float>* scEnergy;
	 std::vector<float>* scRawEnergy;
	 std::vector<float>* scEta;
	 std::vector<float>* scPhi;
	 std::vector<float>* deltaEtaEleClusterAtCalo;
	 std::vector<float>* deltaPhiEleClusterAtCalo;
	 std::vector<float>* deltaEtaEleClusterAtVtx;
	 std::vector<float>* deltaPhiEleClusterAtVtx;

	 //sc parte

	 //es part
	 std::vector<int>* isP1ok;
	 std::vector<int>* isP2ok;

	 std::vector<float>* eleES1;
	 std::vector<float>* eleES2;

	 std::vector<float>* ele_esRhFrac;
	 std::vector<int>* ele_esRhZ; 
	 std::vector<int>* ele_esRhP;
	 std::vector<int>* ele_esRhX; 
	 std::vector<int>* ele_esRhY;
	 std::vector<int>* ele_esRhS;

	 std::vector<float>* trkESz;
	 std::vector<float>* trkESp;
	 std::vector<float>* trkESx;
	 std::vector<float>* trkESy;
	 std::vector<float>* trkESs;
	 std::vector<float>* trkEEz;
	 std::vector<float>* trkEEx;
	 std::vector<float>* trkEEy;
	 std::vector<float>* trkEEEta;
	 std::vector<float>* trkEEPhi;


	 TrackDetectorAssociator trackAssociator_;
	 TrackAssociatorParameters parameters_;

};


#endif
