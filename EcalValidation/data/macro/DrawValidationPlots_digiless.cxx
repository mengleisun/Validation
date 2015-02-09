//
// Macro to produce ECAL cosmic plots
//

// int Wait() {
//      cout << " Continue [<RET>|q]?  ";
//      char x;
//      x = getchar();
//      if ((x == 'q') || (x == 'Q')) return 1;
//      return 0;
// }

#include <algorithm>
#include <string>
#include <vector>
#include <map>

void DrawValidationPlots_digiless(Char_t* infile1 = 0, 
		     Char_t* infile2 = 0, 
		     Char_t* fileType = "png", 
		     Char_t* dirName = ".")
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  gStyle->SetOptStat(1110);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleBorderSize(0);

  if (!infile1 || !infile2) {
    cout << " No input file specified !" << endl;
    return;
  }
  
  cout << "Producing validation plots for: " << infile1 << " and " << infile2 << endl;

  TFile* _f[2]; 
  _f[0] = new TFile(infile2); 
  _f[1] = new TFile(infile1);

  TH1D *h_numberOfEvents[2];
  for (int i=0;i<2;i++)
    h_numberOfEvents[i] = (TH1D*)_f[i]->Get("ecalvalidation/h_numberOfEvents") ; 
  float nEvents0 = h_numberOfEvents[0]->GetBinContent(1);
  float nEvents1 = h_numberOfEvents[1]->GetBinContent(1);
  float s        = nEvents1/nEvents0;
  cout << "SCALE: " << nEvents1 << "/" << nEvents0 << endl;
  //s = 1;
  
  // Define list of object names
  const int nObj=106;
  const int nObj_2D=21;
  char *objName[nObj]={"ecalvalidation/h_recHits_EB_size",
		     "ecalvalidation/h_recHits_EEP_size",
		     "ecalvalidation/h_recHits_EEM_size",
		     "ecalvalidation/h_recHits_ES_size",
		     "ecalvalidation/h_recHits_EB_energy",
		     "ecalvalidation/h_recHits_EEP_energy",
		     "ecalvalidation/h_recHits_EEM_energy",
		     "ecalvalidation/h_recHits_ES_energy",
		     "ecalvalidation/h_recHits_EB_energyMax",
		     "ecalvalidation/h_recHits_EEP_energyMax",
		     "ecalvalidation/h_recHits_EEM_energyMax",
		     "ecalvalidation/h_recHits_ES_energyMax",
		     "ecalvalidation/h_recHits_EB_time",
		     "ecalvalidation/h_recHits_EEP_time",
		     "ecalvalidation/h_recHits_EEM_time",
		     "ecalvalidation/h_recHits_ES_time",
		     "ecalvalidation/h_recHits_eta",
		     "ecalvalidation/h_recHits_EB_phi",
		     "ecalvalidation/h_recHits_EE_phi",
		     "ecalvalidation/h_recHits_EB_Chi2",
		     "ecalvalidation/h_recHits_EEP_Chi2",
		     "ecalvalidation/h_recHits_EEM_Chi2",
		     "ecalvalidation/h_recHits_EB_E1oE4",
                     //Leo
                     "ecalvalidation/h_recHits_EB_recoFlag",
                     "ecalvalidation/h_recHits_EE_recoFlag",
                     "ecalvalidation/h_redRecHits_EB_recoFlag",
                     "ecalvalidation/h_redRecHits_EE_recoFlag",
                     "ecalvalidation/h_basicClusters_recHits_EB_recoFlag",
                     "ecalvalidation/h_basicClusters_recHits_EE_recoFlag",
                     //Leo
		     "ecalvalidation/h_basicClusters_EB_size",
		     "ecalvalidation/h_basicClusters_EEP_size",
		     "ecalvalidation/h_basicClusters_EEM_size",
		     "ecalvalidation/h_superClusters_EB_size",
		     "ecalvalidation/h_superClusters_EEP_size",
		     "ecalvalidation/h_superClusters_EEM_size",
		     "ecalvalidation/h_superClusters_EB_nXtals",
		     "ecalvalidation/h_superClusters_EEP_nXtals",
		     "ecalvalidation/h_superClusters_EEM_nXtals",
		     "ecalvalidation/h_superClusters_EB_nBC",
		     "ecalvalidation/h_superClusters_EEP_nBC",
		     "ecalvalidation/h_superClusters_EEM_nBC",
		     "ecalvalidation/h_superClusters_EB_energy",
		     "ecalvalidation/h_superClusters_EEP_energy",
		     "ecalvalidation/h_superClusters_EEM_energy",
		     "ecalvalidation/h_superClusters_eta",
		     "ecalvalidation/h_superClusters_EB_phi",
		     "ecalvalidation/h_superClusters_EE_phi",
		     "ecalvalidation/h_superClusters_EB_E1oE4",
		     "ecalvalidation/h_superClusters_EEP_E1oE4",
		     "ecalvalidation/h_superClusters_EEM_E1oE4",
		     "ecalvalidation/h_esClusters_energy_plane1",
		     "ecalvalidation/h_esClusters_energy_plane2",
		     "ecalvalidation/h_esClusters_energy_ratio",  
                     //Badder
                     "ecalvalidation/h_digis_EB_ped_mean",
                     "ecalvalidation/h_digisFromRechit_EB_ped_mean",
                     "ecalvalidation/h_digis_EEP_ped_mean",
                     "ecalvalidation/h_digis_EEM_ped_mean",
                     "ecalvalidation/h_digisFromRechit_EEP_ped_mean",
                     "ecalvalidation/h_digisFromRechit_EEM_ped_mean",
                     "ecalvalidation/h_digis_EB_ped_rms",
                     "ecalvalidation/h_digisFromRechit_EB_ped_rms",
                     "ecalvalidation/h_digis_EEP_ped_rms",
                     "ecalvalidation/h_digis_EEM_ped_rms",
                     "ecalvalidation/h_digisFromRechit_EEP_ped_rms",
                     "ecalvalidation/h_digisFromRechit_EEM_ped_rms",
                     "ecalvalidation/h_HF_noise_EB",
                     "ecalvalidation/h_LF_noise_EB",
                     "ecalvalidation/h_Total_noise_EB",
                     "ecalvalidation/h_HF_noise_FromRechit_EB",
                     "ecalvalidation/h_LF_noise_FromRechit_EB",          
                     "ecalvalidation/h_Total_noise_FromRechit_EB",
                     "ecalvalidation/h_HF_noise_EEP",
                     "ecalvalidation/h_LF_noise_EEP",
                     "ecalvalidation/h_Total_noise_EEP",
                     "ecalvalidation/h_HF_noise_EEM",
                     "ecalvalidation/h_LF_noise_EEM",
                     "ecalvalidation/h_Total_noise_EEM",
                     "ecalvalidation/h_HF_noise_FromRechit_EEP",
                     "ecalvalidation/h_LF_noise_FromRechit_EEP",
                     "ecalvalidation/h_Total_noise_FromRechit_EEP",
                     "ecalvalidation/h_HF_noise_FromRechit_EEM",
                     "ecalvalidation/h_LF_noise_FromRechit_EEM",
                     "ecalvalidation/h_Total_noise_FromRechit_EEM",
                     "ecalvalidation/h_HF_noise_FromRechit_vs_Eta",
                     "ecalvalidation/h_LF_noise_FromRechit_vs_Eta",
                     "ecalvalidation/h_Total_noise_FromRechit_vs_Eta", 
                     "ecalvalidation/h_HF_noise_vs_Eta",
                     "ecalvalidation/h_LF_noise_vs_Eta",
                     "ecalvalidation/h_Total_noise_vs_Eta",
                     "ecalvalidation/h_HF_noise_FromRechit_vs_Eta_ped",
                     "ecalvalidation/h_LF_noise_FromRechit_vs_Eta_ped",
                     "ecalvalidation/h_Total_noise_FromRechit_vs_Eta_ped", 
                     "ecalvalidation/h_recHits_EB_iPhiOccupancy",
                     "ecalvalidation/h_recHits_EB_iEtaOccupancy",
                     "ecalvalidation/h_Amplitude_vs_Eta",
                     "ecalvalidation/h_Amplitude_FromRechit_vs_Eta",
                     "ecalvalidation/h_recHits_EB_sumEt",  
                     "ecalvalidation/h_recHits_EEP_sumEt", 
                     "ecalvalidation/h_recHits_EEM_sumEt",  
                     "ecalvalidation/h_recHits_EEP_sumEtCut",  
                     "ecalvalidation/h_recHits_EEM_sumEtCut",
                     "ecalvalidation/h_HF_noise_vs_Eta_ped",
                     "ecalvalidation/h_LF_noise_vs_Eta_ped",
                     "ecalvalidation/h_Total_noise_vs_Eta_ped",
                     "ecalvalidation/h_recHits_EB_SRP",
                     "ecalvalidation/h_recHits_EE_SRP"};                       

 char *objName_2D[nObj_2D]={"ecalvalidation/h_recHits_EB_occupancy",
                          "ecalvalidation/h_recHits_EEP_occupancy",
                          "ecalvalidation/h_recHits_EEM_occupancy",
                          "ecalvalidation/h_digis_EB_occupancy",
                          "ecalvalidation/h_digis_EEP_occupancy",
                          "ecalvalidation/h_digis_EEM_occupancy",
                          "ecalvalidation/h_HF_noise_iphieta_EB",
                          "ecalvalidation/h_LF_noise_iphieta_EB",
                          "ecalvalidation/h_Total_noise_iphieta_EB",
                          "ecalvalidation/h_HF_noise_ixiy_EEP",
                          "ecalvalidation/h_LF_noise_ixiy_EEP",
                          "ecalvalidation/h_Total_noise_ixiy_EEP",
                          "ecalvalidation/h_HF_noise_ixiy_EEM",
                          "ecalvalidation/h_LF_noise_ixiy_EEM",
                          "ecalvalidation/h_Total_noise_ixiy_EEM",
                          "ecalvalidation/h_Amplitude_iphieta_EB",
                          "ecalvalidation/h_Amplitude_FromRechit_iphieta_EB",
                          "ecalvalidation/h_Amplitude_ixiy_EEP",
                          "ecalvalidation/h_Amplitude_ixiy_EEM",
                          "ecalvalidation/h_Amplitude_FromRechit_ixiy_EEP",
                          "ecalvalidation/h_Amplitude_FromRechit_ixiy_EEM"};
                          
               
 char *objTitle[nObj]={"Number of RecHits (EB)",
		     "Number of RecHits (EE+)",
		     "Number of RecHits (EE-)",
		     "Number of RecHits (ES)",
		     "RecHits Energy (EB)",
		     "RecHits Energy (EE+)",
		     "RecHits Energy (EE-)",
		     "RecHits Energy (ES)",
		     "RecHits Max Energy (EB)",
		     "RecHits Max Energy (EE+)",
		     "RecHits Max Energy (EE-)",
		     "RecHits Max Energy (ES)",
		     "RecHits Time (EB)",
		     "RecHits Time (EE+)",
		     "RecHits Time (EE-)",
		     "RecHits Time (ES)",
		     "RecHits Eta",
		     "RecHits Phi (EB)",
		     "RecHits Phi (EE)", 
		     "RecHits \ #chi^{2} (EB)",
		     "RecHits \ #chi^{2} (EE+)",
		     "RecHits \ #chi^{2} (EE-)",
		     "RecHits 1-E4/E1 (EB)",
                     //Leo
                     "Number of RecHits (EB)",
                     "Number of RecHits (EE)",
                     "Number of RedRecHits (EB)",
                     "Number of RedRecHits (EE)",
                     "Number of ClusRecHits (EB)",
                     "Number of ClusRecHits (EE)",
                     //Leo
                     "Number of Basic Clusters (EB)",
		     "Number of Basic Clusters (EE+)",
		     "Number of Basic Clusters (EE-)",
		     "Number of Superclusters (EB)",
		     "Number of Superclusters (EE+)",
		     "Number of Superclusters (EE-)",
		     "Number of Crystal per Supercluster (EB)",
		     "Number of Crystal per Supercluster (EE+)",
		     "Number of Crystal per Supercluster (EE-)",
		     "Number of Basic Clusters  per Supercluster (EB)",
		     "Number of Basic Clusters  per Supercluster (EE-)",
		     "Number of Basic Clusters  per Supercluster (EE+)",
		     "Supercluster Energy (EB)",
		     "Supercluster Energy (EE+)",
		     "Supercluster Energy (EE-)",
		     "Superclusters Eta",
		     "Superclusters Phi (EB)",
		     "Superclusters Phi (EE)",
		     "1-E4/E1",
		     "1-E4/E1",
		     "1-E4/E1",
		     "ES Clusters Energy - Plane 1",
		     "ES Clusters Energy - Plane 2",
		     "ES Clusters Energy - Ratio",
                     //Badder
                     "SimDigi Pedestal Mean (EB)",
                     "RecoDigi Pedestal Mean (EB)",
                     "SimDigi Pedestal Mean (EE+)",
                     "SimDigi Pedestal Mean (EE-)",
                     "RecoDigi Pedestal Mean (EE+)",
                     "RecoDigi Pedestal Mean (EE-)",
                     "SimDigi Pedestal RMS (EB)",
                     "RecoDigi Pedestal RMS (EB)",
                     "SimDigi Pedestal RMS (EE+)",
                     "SimDigi Pedestal RMS (EE-)",
                     "RecoDigi Pedestal RMS (EE+)",
                     "RecoDigi Pedestal RMS (EE-)",
                     "SimDigi HF-noise (EB)", 
                     "SimDigi LF-noise (EB)",
                     "SimDigi Total-noise (EB)",
                     "RecoDigi HF-noise (EB)", 
                     "RecoDigi LF-noise (EB)",
                     "RecoDigi Total-noise (EB)",
                     "SimDigi HF-noise (EE+)",
                     "SimDigi LF-noise (EE+)",
                     "SimDigi Total-noise (EE+)",
                     "SimDigi HF-noise (EE-)",
                     "SimDigi LF-noise (EE-)",
                     "SimDigi Total-noise (EE-)",
                     "RecoDigi HF-noise (EE+)", 
                     "RecoDigi LF-noise (EE+)",
                     "RecoDigi Total-noise (EE+)",
                     "RecoDigi HF-noise (EE-)", 
                     "RecoDigi LF-noise (EE-)",
                     "RecoDigi Total-noise (EE-)",
                     "RecoDigi HF-noise vs Eta",
                     "RecoDigi LF-noise vs Eta",
                     "RecoDigi Total noise vs Eta",
                     "SimDigi HF-noise vs Eta",
                     "SimDigi LF-noise vs Eta",
                     "SimDigi Total noise vs Eta",
                     "RecoDigi HF-noise vs Eta (pedestal)",
                     "RecoDigi LF-noise vs Eta (pedestal)",
                     "RecoDigi Total noise vs Eta (pedestal)",
                     "RecHits iPhiOccupancy",
                     "RecHits iEtaOccupancy",
                     "SimDigi Amplitude",
                     "RecoDigi Amplitude",
                     "RecHits SumEt (EB)",
                     "RecHits SumEt (EE+)",
                     "RecHits SumEt (EE-)", 
                     "RecHits SumEt (EE+, |eta| < 2.5)",
                     "RecHits SumEt (EE-, |eta| < 2.5)",
                     "SimDigi HF-noise vs Eta (pedestal)",
                     "SimDigi LF-noise vs Eta (pedestal)",
                     "SimDigi Total noise vs Eta (pedestal)",
                     "Rechit SRP (EB)",
                     "Rechit SRP (EE)"};
                     
  
  char *objTitle_2D[nObj_2D]={"RecHits Occupancy (EB)",
                            "RecHits Occupancy (EE+)",
                            "RecHits Occupancy (EE-)",
                            "Digis Occupancy (EB)",
                            "Digis Occupancy (EE+)",
                            "Digis Occupancy (EE-)",
                            "HF noise  vs iPhi-iEta (EB)",
                            "LF noise  vs iPhi-iEta (EB)",
                            "Total noise  vs iPphi-iEta (EB)",
                            "HF noise  vs iX-iY (EEP)",
                            "LF noise  vs iX-iY (EEP)",
                            "Total noise  vs iX-iY (EEP)",
                            "HF noise  vs iX-iY (EEM)",
                            "LF noise  vs iX-iY (EEM)",
                            "Total noise  vs iX-iY (EEM)",
                            "SimDigi Amplitude vs iPhi-iEta (EB)",
                            "RecoDigi Amplitude vs iPhi-iEta (EB)",
                            "SimDigi Amplitude  vs iX-iY (EEP)",
                            "SimDigi Amplitude  vs iX-iY (EEM)",
                            "RecoDigi Amplitude  vs iX-iY (EEP)",
                            "RecoDigi Amplitude  vs iX-iY (EEM)"};
                                         
  char *labelX[nObj]={"Number of RecHits/Event",
		    "Number of RecHits/Event",
		    "Number of RecHits/Event",
		    "Number of RecHits/Event",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "time(ns)",
		    "time(ns)",
		    "time(ns)",
		    "time(ns)",
		    "eta",
		    "phi",
		    "phi",
		    "#chi^{2}",
		    "#chi^{2}",
		    "#chi^{2}",
		    "1-E4/E1",
                    "recoFlag",
                    "recoFlag",
                    "recoFlag",
                    "recoFlag",
                    "recoFlag",
                    "recoFlag",
                    "BasicClusters/Event",
		    "BasicClusters/Event",
		    "BasicClusters/Event",
		    "Superclusters/Event",
		    "Superclusters/Event",
		    "Superclusters/Event",
		    "Crystals/Supercluster",
		    "Crystals/Supercluster",
		    "Crystals/Supercluster",
		    "BasicClusters/Supercluster",
		    "BasicClusters/Supercluster",
		    "BasicClusters/Supercluster",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "eta",
		    "phi",
		    "phi",
		    "1-E4/E1",
		    "1-E4/E1",
		    "1-E4/E1",
		    "Energy (GeV)",
		    "Energy (GeV)",
		    "EnergyPlane1/EnergyPlane1",
                    "ADC",
                    "ADC",
                    "ADC",
	            "ADC",
                    "ADC",
                    "ADC",
                    "ADC",
	            "ADC",
		    "ADC",
		    "ADC",
                    "ADC",
                    "ADC",
                    "ADC",
	            "ADC",
                    "ADC",
                    "ADC",
                    "ADC",
	            "ADC",
		    "ADC",
		    "ADC",
                    "ADC",
		    "ADC",
		    "ADC",
		    "ADC",
		    "ADC",
                    "ADC",
		    "ADC",
		    "ADC",
                    "ADC",
                    "ADC",
                    "Eta",
                    "Eta",
                    "Eta",
                    "Eta",
                    "Eta",
                    "Eta",
                    "Eta",
                    "Eta",
                    "Eta",
                    "iPhi",
                    "iEta",
                    "Eta",
                    "Eta",
                    "GeV",
                    "GeV",
                    "GeV",
                    "GeV",
                    "GeV",
                    "",
                    ""};
                    
                             
  char *labelY[1]={"Counts"};
  
  char *labelX_2D[nObj_2D]={"iPhi",
                          "iX",
                          "iX", 
                          "iPhi",
                          "iX",
                          "iX",
                          "iPhi",
                          "iPhi",
                          "iPhi",
                          "iX",
                          "iX",
                          "iX",
                          "iX",
                          "iX",
                          "iX",
                          "iPhi",
                          "iPhi",
                          "iX",
                          "iX",
                          "iX",
                          "iX"};

  char *labelY_2D[nObj_2D]={"iEta",
                          "iY",
                          "iY", 
                          "iEta",
                          "iY",
                          "iY",
                          "iEta",
                          "iEta",
                          "iEta",
                          "iY",
                          "iY",
                          "iY",
                          "iY",
                          "iY",
                          "iY",
                          "iEta",
                          "iEta",
                          "iY",
                          "iY",
                          "iY",
                          "iY"};
                
  char *recoFlagLabels[16]={"kGood",
                    "kPoorReco",
                    "kOutOfTime",
                    "kFaultyHardware",
                    "kNoisy",
                    "kPoorCalib",
                    "kSaturated",
                    "kLeadingEdgeRecovered",
                    "kNeighboursRecovered",
                    "kTowerRecovered",
                    "kDead",
                    "kKilled",
                    "kTPSaturated",
                    "kL1SpikeFlag",
                    "kWeird",
                    "kDiWeird"};
                    
//                     "kFake",
//                     "kFakeNeighbours",
//                     "kDead",
//                     "kKilled",
//                     "kTPSaturated",
//                     "kL1SpikeFlag"};

  double xMin[nObj]={0.,0.,0.,0., 
		   0,0,0,-10, 
		   0,1.,1.,1.,
		   -100.,-50.,-50.,-50.,
		   -3.,-3.2,-3.2,
		   0,0,0,
		   0,
                   0,0,
                   0,0,
                   0,0,
                   0,0,0,
		   0,0,0,
		   0,0,0,
		   0,0,0,
		   0,0,0,
		   -3.,-3.2,-3.2,
		   0,0,0,
		   0,0,0,
                   175,175,175,175,175,175,
		   -1,-1,-1,-1,-1,-1,
                   -12,-12,-12,-12,-12,
                   -12,-12,-12,-12,-12,
                   -12,-12,-12,-12,-12,
                   -12,-12,-12,
                   -3,-3,-3,
                   -3,-3,-3,
                   -3,-3,-3,
                   1.,-86,
                   -3,-3,
                   0,
                   0,0,
                   0,0,
                   -3,-3,-3,
                   0,0};
  
  double xMin_2D[nObj_2D]={1,0,0,
                           1,0,0, 
                           1,1,1, 
                           0,0,0, 
                           0,0,0,
                           1,1,
                           0,0,0,0};           
  
  double xMax[nObj]={10000.,1700.,1700.,8000.,  
		   300,300.,300.,0.06, 
		   300,300,300,0.05,
		   100.,50.,50.,50.,
		   3.,3.2,3.2,
		   70,70,70,
		   1.2,
                   32,32,
                   32,32,
                   32,32,
                   80,150,150,
		   60,60,60,
		   200,200,200,
		   20,20,20,
		   500,500,500,
		   3.,3.2,3.2,
		   1.2,1.2,1.2,
		   0.05,0.05,100,
                   220,220,220,220,220,220,
                   11,11,11,11,11,11,
                   12,12,12,12,12,
                   12,12,12,12,12,
                   12,12,12,12,12,12,12,12,
                   3,3,3,
                   3,3,3,
                   3,3,3,
                   361,86,
                   3,3,
                   1400,
                   150,150,
                   100,100,
                   3,3,3,
                   3,3};

  double xMax_2D[nObj_2D]={361,100,100,   
                           361,100,100,
                           361,361,361,
                           100,100,100,
                           100,100,100,
                           361,361,
                           100,100,100,100};
  
  double yMin_2D[nObj_2D]={-86,0,0,
                           -86,0,0,
                           -86,-86,-86,
                            0,0,0,
                            0,0,0,
                            -86,-86,
                            0,0,0,0};

  double yMax_2D[nObj_2D]={86,100,100,
                           86,100,100,
                           86,86,86,
                           100,100,100,
                           100,100,100,
                           86,86,
                           100,100,100,100};

  double zMin_2D[nObj_2D]={0,0,0,
                           1000,1000,1000,
                           0,0,0,
                           0,0,0,
                           0,0,0,
                           0,0,
                           0,0,0,0};
  
  double zMax_2D[nObj_2D]={5000,8000,8000,
                           1020,1020,1020,
                           6,6,6,
                           8,8,8,
                           8,8,8,
                           5,5,
                           8,8,8,8};

  int reBin[500]  = {20,2,2,2, 
		    20,20,20,8, 
		    4,4,4,4, 
		    1,1,1,1, 
		    5,5,5, 
		    5,5,5,
		    1,
                    1,1,
                    1,1,
                    1,1,
                    1,1,1,
		    1,1,1,
		    1,1,1,
		    1,1,1,
		    10,10,10,
		    5,5,5,
		    1,1,1,
		    10,10,1,
		    1,1,1,1,1,
		    1,1,1,1,1,
                    1,1,1,1,1,
                    1,1,1,1,1,
                    1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1};


  int optLogY[nObj] = {1,1,1,1,
		     1,1,1,1,
		     1,1,1,1,
		     1,0,0,0,
		     0,0,0,
		     1,1,1,
		     0,
                     1,1,
                     1,1,
                     1,1,
                     1,1,1,
		     1,1,1,
		     1,1,1,
		     1,1,1,
		     1,1,1,
		     0,0,0,
		     0,0,0,
		     1,1,1};


  TH1D* h[2][500]; 
  TCanvas *c[500];
  TPaveStats *st[500];
  
  bool isMissing = false;
  
  int iHisto = 0;
  while(iHisto<nObj){
    for (int ifile=0;ifile<2;ifile++){ 
    h[ifile][iHisto] = (TH1D*)_f[ifile]->Get(objName[iHisto]);
    if(h[ifile][iHisto] == 0)isMissing = true;
    if(h[ifile][iHisto] == 0) continue;
    //std::cout << "HistoName: " << h[ifile][iHisto]->GetName() << std::endl; 
    h[ifile][iHisto]->Rebin(reBin[iHisto]);
    if (ifile == 0) {
      // open a new canvas
      c[iHisto] = new TCanvas(objName[iHisto],objName[iHisto],50+iHisto*20,50+iHisto*5,500,400);
      c[iHisto]->cd();
      // customize and plot
      h[ifile][iHisto]->GetYaxis()->SetTitle(labelY[0]);
      std::string histo_name = std::string(h[ifile][iHisto]->GetName());
      if(histo_name.find("Eta") != std::string::npos) h[ifile][iHisto]->GetYaxis()->SetTitle("ADC");
      if ( iHisto > 22 && iHisto < 29 ) //Set reasonable labels for recoFlag plots
      {
        int nBin = h[ifile][iHisto] ->GetXaxis()-> GetNbins ();
        for ( int ibin = 1; ibin <= nBin; ibin++ ) 
          h[ifile][iHisto]->GetXaxis()->SetBinLabel( ibin, recoFlagLabels[ibin-1]);
      }
      else h[ifile][iHisto]->GetXaxis()->SetTitle(labelX[iHisto]);
      h[ifile][iHisto]->SetFillColor(kBlack+10);
      h[ifile][iHisto]->SetFillStyle(3002);
      h[ifile][iHisto]->SetTitle(objTitle[iHisto]);
      h[ifile][iHisto]->Scale(s);
      h[ifile][iHisto]->GetXaxis()->SetRangeUser(xMin[iHisto],xMax[iHisto]);
      h[ifile][iHisto]->Draw();

    }
    if (ifile == 1) {

      if(isMissing == true) continue;

      h[ifile][iHisto]->SetMarkerStyle(20);
      h[ifile][iHisto]->SetMarkerSize(0.7);
      h[ifile][iHisto]->SetMarkerColor(kRed);
      h[ifile][iHisto]->SetLineColor(kBlack);
      std::string histo_name = std::string(h[ifile][iHisto]->GetName());
      if(histo_name.find("Eta") != std::string::npos) h[ifile][iHisto]->Draw("psames");
      else h[ifile][iHisto]->Draw("epsames");
     
      // set the log-scale and range 
      float maxy = max (h[ifile][iHisto]->GetMaximum(),h[ifile-1][iHisto]->GetMaximum() );
      float miny = h[ifile][iHisto]->GetMinimum();
      if (optLogY[iHisto]) {
	c[iHisto]->SetLogy();
	h[0][iHisto]->SetMaximum(maxy*2);
	h[0][iHisto]->SetMinimum(0.1);
      }
      else  h[0][iHisto]->SetMaximum(maxy*1.1);
      
      c[iHisto]->Update();
      
      // stat. box
      st[iHisto]= (TPaveStats*)(h[ifile][iHisto]->GetListOfFunctions()->FindObject("stats"));
      st[iHisto]->SetY1NDC(0.72); //new y start position
      st[iHisto]->SetY2NDC(0.85); //new y end position
      st[iHisto]->SetTextColor(kRed);
      st[iHisto]->Draw();
     
      
    }
    
  }
  
  char *str = strtok (objName[iHisto],"/");
  str = strtok (NULL, "/");

  char myname[500];
  sprintf (myname,"%s/",dirName);
  strcat(myname,str);
  strcat(myname,".");
  strcat(myname,fileType);

  if(isMissing == false) c[iHisto]->Print(myname,fileType);
  isMissing = false;
  
  iHisto++;
  
  
  }

  /*std::vector<int> nEvents;
  nEvents.push_back(h[0][0]->GetEntries());
  nEvents.push_back(h[1][0]->GetEntries());*/
  
  TH2D* h2[2][100]; 
  TCanvas *c2[2][100];
  TPaveStats *st2[2][100];
  
  map<int, map<int,bool> > isMissing_2D;
  for(int ii = 0; ii < 2; ii++)
      for(int jj = 0; jj < nObj_2D; jj++)
          isMissing[ii][jj] = false;
  int iHisto = 0;
  while(iHisto<nObj_2D){
    for (int ifile=0;ifile<2;ifile++){ 
    h2[ifile][iHisto] = (TH2D*)_f[ifile]->Get(objName_2D[iHisto]);
    if(h2[ifile][iHisto] == 0) isMissing_2D[ifile][iHisto]  = true;
    if(h2[ifile][iHisto] == 0) continue;
    char c_name[500];
    sprintf (c_name,"%s_%d_%d",objName_2D[iHisto],iHisto,ifile);
    c2[ifile][iHisto] = new TCanvas(c_name,c_name,50+iHisto*20,50+iHisto*5,500,400);
    c2[ifile][iHisto]->cd();
    // customize and plot
    if(iHisto < 6) h2[ifile][iHisto]->Scale(h[0][0]->GetEntries()/h[ifile][0]->GetEntries());
    h2[ifile][iHisto]->GetXaxis()->SetTitle(labelX_2D[iHisto]);
    h2[ifile][iHisto]->GetYaxis()->SetTitle(labelY_2D[iHisto]);
    h2[ifile][iHisto]->GetXaxis()->SetRangeUser(xMin_2D[iHisto],xMax_2D[iHisto]);
    h2[ifile][iHisto]->GetYaxis()->SetRangeUser(yMin_2D[iHisto],yMax_2D[iHisto]);
    h2[ifile][iHisto]->GetZaxis()->SetRangeUser(zMin_2D[iHisto],zMax_2D[iHisto]);
    h2[ifile][iHisto]->Draw("colz");
    c2[ifile][iHisto]->Update();
    // stat. box
    st2[ifile][iHisto]= (TPaveStats*)(h2[ifile][iHisto]->GetListOfFunctions()->FindObject("stats"));
    /*st2[ifile][iHisto]->SetX1NDC(0.71); //new x start position
    st2[ifile][iHisto]->SetX2NDC(0.99); //new x end position
    st2[ifile][iHisto]->SetY1NDC(0.84); //new y start position
    st2[ifile][iHisto]->SetY2NDC(0.99); //new y end position*/
    st2[ifile][iHisto]->SetX1NDC(0.); //new x start position
    st2[ifile][iHisto]->SetX2NDC(0.); //new x end position
    st2[ifile][iHisto]->SetY1NDC(0.); //new y start position
    st2[ifile][iHisto]->SetY2NDC(0.); //new y end position
    //st[ifile][iHisto]->SetTextColor(kRed);
    st2[ifile][iHisto]->Draw();
    std::string Name = std::string(objName_2D[iHisto]);
    
    char *str = strtok (Name.c_str(),"/");
    str = strtok (NULL, "/");

    char myname[500];
    sprintf (myname,"%s/",dirName);
    strcat(myname,str);
    if(ifile == 0) strcat(myname,"_0");
    if(ifile == 1) strcat(myname,"_1");
    strcat(myname,".");
    strcat(myname,fileType);

    if(isMissing_2D[ifile][iHisto] == false) c2[ifile][iHisto]->Print(myname,fileType);
    isMissing_2D[ifile][iHisto] = true;

    }
  
  iHisto++;

  }
}

