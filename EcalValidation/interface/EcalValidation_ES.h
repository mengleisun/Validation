#ifndef EcalValidation_ES_h
#define EcalValidation_ES_h

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


// ROOT include
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TProfile.h"
#include "TProfile2D.h"


// Less than operator for sorting EcalRecHits according to energy.
class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool>
{
public:
  bool operator()(EcalRecHit x, EcalRecHit y)
  {
    return (x.energy() > y.energy());
  }
};


class EcalValidation_ES : public edm::EDAnalyzer {
  
      public:
         explicit EcalValidation_ES(const edm::ParameterSet&);
	 ~EcalValidation_ES();
  
  
      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
	 virtual void endJob() ;

	 // ----------additional functions-------------------
	 void convxtalid(int & , int &);
	 int diff_neta_s(int,int);
	 int diff_nphi_s(int,int);

      protected:

	 // ----------member data ---------------------------
	 bool isMC_;
	 edm::InputTag recHitCollection_EB_;
	 edm::InputTag recHitCollection_EE_;
	 edm::InputTag recHitCollection_ES_;
         edm::InputTag redRecHitCollection_ES_;
	 edm::InputTag esClusterCollectionX_ ;
	 edm::InputTag esClusterCollectionY_ ;
	 edm::InputTag esDigiCollection_ ;

	 edm::InputTag EleTagMC_;
	 edm::InputTag EleTag_;
	 edm::InputTag PhoTag_;
	 bool saveDigis_;

	 std::map<DetId, EcalRecHit> ESrechits_map;

	 bool SaveSrFlag_;
         
	 double ethrEB_;
	 double ethrEE_;
         double gainId_;
	 

	 double scEtThrEB_;
	 double scEtThrEE_;


	 // for Pi0
	 PositionCalc posCalculator_ ;
	 
	 double clusSeedThr_;
	 double clusSeedThr_EE_;
	 int clusEtaSize_;
	 int clusPhiSize_;
	 
	 double seleXtalMinEnergy_;
	 double seleXtalMinEnergy_EE_;
	 
	 bool isMaskEB_;
	 bool isMaskEE_;
	 
	 /// Files with xtals masked
	 std::string maskEBFile_;
	 std::string maskEEFile_;
	 
	 /// Use Reco Flag from RH
	 bool useRecoFlag_;
	 
	 // ------------- HISTOGRAMS ------------------------------------

	 int naiveId_;
	 TH1D *h_numberOfEvents;
         
	//PulseShape
	TProfile* h_PulseShape_EB;
	TProfile* h_PulseShape_EE;
	TProfile* h_PulseShape_ES;
	TProfile* h_PulseShape_ESpedSub;
	TProfile* h_PulseShape_ES_inTime;
	TProfile* h_PulseShape_ES_eOOT;
	TProfile* h_PulseShape_ES_lOOT;
	TProfile* h_PulseShape_ES_kGood;
	TProfile* h_PulseShape_ES_kGoodpedSub;

	TProfile* h_PulseShape_ES_inTime_kGood;
	TProfile* h_PulseShape_ES_eOOT_kGood;
	TProfile* h_PulseShape_ES_lOOT_kGood;

	TH1D* h_PulseShape_ES_pedSub[100];

	 // PRESHOWER ----------------------------------------------
	 
	 TH1D *h_recHits_ES_size;
	 TH1D *h_recHits_ES_size_F[2];
	 TH1D *h_recHits_ES_size_R[2];

	 TH1D *h_recHits_ES_energy;
	 TH1D *h_recHits_ES_energy_F[2];
	 TH1D *h_recHits_ES_energy_R[2];

	 TH1D *h_recHits_ES_energyMax;
	 TH1D *h_recHits_ES_energyMax_F[2];
	 TH1D *h_recHits_ES_energyMax_R[2];
	
	 TH1D *h_recHits_ES_time;
	 TH1D *h_recHits_ES_time_F[2];
	 TH1D *h_recHits_ES_time_R[2];

	 TH1D *h_esClusters_energy_plane1;
	 TH1D *h_esClusters_energy_plane2;
	 TH1D *h_esClusters_energy_ratio;
	 
	 TH1D* h_recHits_ES_recoFlag;

	 TH1D *h_PulseShape_ES_kESGood;
	 TH1D *h_PulseShape_ES_kESDead;
	 TH1D *h_PulseShape_ES_kESHot;
	 TH1D *h_PulseShape_ES_kESPassBX;
	 TH1D *h_PulseShape_ES_kESTwoGoodRatios;
	 TH1D *h_PulseShape_ES_kESBadRatioFor12;
	 TH1D *h_PulseShape_ES_kESBadRatioFor23Upper;
	 TH1D *h_PulseShape_ES_kESBadRatioFor23Lower;
	 TH1D *h_PulseShape_ES_kESTS1Largest;
	 TH1D *h_PulseShape_ES_kESTS3Largest;
	 TH1D *h_PulseShape_ES_kESTS3Negative;
	 TH1D *h_PulseShape_ES_kESSaturated;
	 TH1D *h_PulseShape_ES_kESTS2Saturated;
	 TH1D *h_PulseShape_ES_kESTS3Saturated;
	 TH1D *h_PulseShape_ES_kESTS13Sigmas;
	 TH1D *h_PulseShape_ES_kESTS15Sigmas;
	 TH1D *h_PulseShape_ES_flag;

	 TH1D* h_CentralStrip_P1_E;
	 TH1D* h_CentralStrip_P2_E;

	 //	 TH1D* plane1_stripsN_E[31*500];
	 //	 TH1D* plane2_stripsN_E[31*500];


	 TH2D* plane2_vs_plane1;
	 TH1D* sclusterEnergy;
	 TH1D* ES_over_sclusterEnergy;
	 TH1D* nEvents;
	 TH1D* nElectronsPerEvents;

	 int totEvents;

	 TTree* tree_;
	 int event;
	 int run;
	 int lumi; 
	 int nRH;
	 int nt_rhES_z[200000];
	 int nt_rhES_p[200000];
	 int nt_rhES_x[200000];
	 int nt_rhES_y[200000];
	 int nt_rhES_s[200000];
	 float nt_rhES_t[200000];
	 float nt_rhES_e[200000];
	 float nt_rhES_posx[200000];
	 float nt_rhES_posy[200000];
	 float nt_rhES_posz[200000];
	 float nt_rhES_eta[200000];
	 float nt_rhES_phi[200000];
	 //int nt_rhES_ADC0;
	 //int nt_rhES_ADC1;
	 //iont nt_rhES_ADC2;
	 int nt_rhES_status[200000];
	 int nSC;
	 float nt_scEn[500];
	 float nt_scX[500];
	 float nt_scY[500];
	 float nt_scZ[500];
	 float nt_scEta[500];
	 float nt_scPhi[500];
	 float nt_scRawEn[500];
	 float nt_scESEn[500];
	 float nt_scEtaWidth[500];
	 float nt_scPhiWidth[500];

 	 float nt_p1_stripsN_E[15500];
 	 float nt_p2_stripsN_E[15500];

};


#endif
