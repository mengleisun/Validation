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
  float scaleEvents        = nEvents1/nEvents0;
  cout << "SCALE: " << nEvents1 << "/" << nEvents0 << endl;
  //s = 1;
  
  // Define list of object names
  const int nObj=130;
  const int nObj_2D=29;
  char *objName[nObj]={"ecalvalidation/h_recHits_EB_size",  //-1
		       "ecalvalidation/h_recHits_EEP_size",
		       "ecalvalidation/h_recHits_EEM_size",
		       //RA
		       "ecalvalidation/h_recHits_ES_size",
		       "ecalvalidation/h_recHits_ES_size_gr",
		       "ecalvalidation/h_recHits_ES_size_F+",
		       "ecalvalidation/h_recHits_ES_size_F-",
		       "ecalvalidation/h_recHits_ES_size_R+",
		       "ecalvalidation/h_recHits_ES_size_R-",
		       "ecalvalidation/h_recHits_EB_energy",
		       "ecalvalidation/h_recHits_EEP_energy",
		       "ecalvalidation/h_recHits_EEM_energy",//10
		       "ecalvalidation/h_recHits_ES_energy",
		       "ecalvalidation/h_recHits_ES_energy_gr",
		       "ecalvalidation/h_recHits_ES_energy_F+",
		       "ecalvalidation/h_recHits_ES_energy_F-",
		       "ecalvalidation/h_recHits_ES_energy_R+",
		       "ecalvalidation/h_recHits_ES_energy_R-",
		       "ecalvalidation/h_recHits_EB_energyMax",
		       "ecalvalidation/h_recHits_EEP_energyMax",
		       "ecalvalidation/h_recHits_EEM_energyMax",
		       "ecalvalidation/h_recHits_ES_energyMax",//20
		       "ecalvalidation/h_recHits_EB_time",
		       "ecalvalidation/h_recHits_EEP_time",
		       "ecalvalidation/h_recHits_EEM_time",
		       "ecalvalidation/h_recHits_ES_time", //24
		       "ecalvalidation/h_recHits_ES_time_gr",
		       "ecalvalidation/h_recHits_ES_time_F+",
		       "ecalvalidation/h_recHits_ES_time_F-",
		       "ecalvalidation/h_recHits_ES_time_R+",
		       "ecalvalidation/h_recHits_ES_time_R-", //29
		       "ecalvalidation/h_recHits_EB_size_cleaned",//30
		       "ecalvalidation/h_recHits_EB_energy_cleaned",
		       "ecalvalidation/h_recHits_EB_energyMax_cleaned",
		       "ecalvalidation/h_recHits_EB_time_cleaned",
		       "ecalvalidation/h_recHits_eta",
		       "ecalvalidation/h_recHits_EB_phi",
		       "ecalvalidation/h_recHits_EE_phi",
		       "ecalvalidation/h_recHits_EB_Chi2",
		       "ecalvalidation/h_recHits_EEP_Chi2",
		       "ecalvalidation/h_recHits_EEM_Chi2",
		       "ecalvalidation/h_recHits_EB_E1oE4",//40
		       //Leo
		       "ecalvalidation/h_recHits_EB_recoFlag",//41
		       "ecalvalidation/h_recHits_EE_recoFlag",
		       "ecalvalidation/h_redRecHits_EB_recoFlag",
		       "ecalvalidation/h_redRecHits_EE_recoFlag",
		       "ecalvalidation/h_basicClusters_recHits_EB_recoFlag",
		       "ecalvalidation/h_basicClusters_recHits_EE_recoFlag", //46
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
		       "ecalvalidation/h_basicClusters_EB_size_cleaned", //56
		       //"ecalvalidation/h_basicClusters_EB_nXtals_cleaned",
		       //"ecalvalidation/h_basicClusters_EB_energy_cleaned",		      
		       "ecalvalidation/h_superClusters_EB_size_cleaned",   
		       "ecalvalidation/h_superClusters_EB_nXtals_cleaned",
		       "ecalvalidation/h_superClusters_EB_nBC_cleaned",
		       "ecalvalidation/h_superClusters_EB_energy_cleaned",  //60
		       //"ecalvalidation/h_superClusters_EB_rawEnergy_cleaned",
		       //"ecalvalidation/h_superClusters_EB_rawEt_cleaned",    		       
		       "ecalvalidation/h_superClusters_EB_nBC",
		       "ecalvalidation/h_superClusters_EEP_nBC",
		       "ecalvalidation/h_superClusters_EEM_nBC",
		       "ecalvalidation/h_superClusters_EB_energy",
		       "ecalvalidation/h_superClusters_EEP_energy",
		       "ecalvalidation/h_superClusters_EEM_energy",
		       "ecalvalidation/h_superClusters_eta", //67
		       "ecalvalidation/h_superClusters_EB_phi",
		       "ecalvalidation/h_superClusters_EE_phi",
		       "ecalvalidation/h_superClusters_EB_E1oE4",//70
		       "ecalvalidation/h_superClusters_EEP_E1oE4",
		       "ecalvalidation/h_superClusters_EEM_E1oE4",
		       "ecalvalidation/h_esClusters_energy_plane1",
		       "ecalvalidation/h_esClusters_energy_plane2",
		       "ecalvalidation/h_esClusters_energy_ratio", //75   
		       //Badder
		       "ecalvalidation/h_digis_EB_ped_mean",
		       "ecalvalidation/h_digisFromRechit_EB_ped_mean",
		       "ecalvalidation/h_digis_EEP_ped_mean",
		       "ecalvalidation/h_digis_EEM_ped_mean",
		       "ecalvalidation/h_digisFromRechit_EEP_ped_mean", //80
		       "ecalvalidation/h_digisFromRechit_EEM_ped_mean",
		       "ecalvalidation/h_digis_EB_ped_rms",
		       "ecalvalidation/h_digisFromRechit_EB_ped_rms",
		       "ecalvalidation/h_digis_EEP_ped_rms",
		       "ecalvalidation/h_digis_EEM_ped_rms",
		       "ecalvalidation/h_digisFromRechit_EEP_ped_rms",
		       "ecalvalidation/h_digisFromRechit_EEM_ped_rms",
		       "ecalvalidation/h_HF_noise_EB",
		       "ecalvalidation/h_LF_noise_EB",
		       "ecalvalidation/h_Total_noise_EB",//90
		       "ecalvalidation/h_HF_noise_FromRechit_EB",
		       "ecalvalidation/h_LF_noise_FromRechit_EB",          
		       "ecalvalidation/h_Total_noise_FromRechit_EB",
		       "ecalvalidation/h_HF_noise_EEP",
		       "ecalvalidation/h_LF_noise_EEP",
		       "ecalvalidation/h_Total_noise_EEP",
		       "ecalvalidation/h_HF_noise_EEM",
		       "ecalvalidation/h_LF_noise_EEM",
		       "ecalvalidation/h_Total_noise_EEM",
		       "ecalvalidation/h_HF_noise_FromRechit_EEP", //100
		       "ecalvalidation/h_LF_noise_FromRechit_EEP",
		       "ecalvalidation/h_Total_noise_FromRechit_EEP",
		       "ecalvalidation/h_HF_noise_FromRechit_EEM",
		       "ecalvalidation/h_LF_noise_FromRechit_EEM",
		       "ecalvalidation/h_Total_noise_FromRechit_EEM",
		       "ecalvalidation/h_HF_noise_FromRechit_vs_Eta",
		       "ecalvalidation/h_LF_noise_FromRechit_vs_Eta",
		       "ecalvalidation/h_Total_noise_FromRechit_vs_Eta", 
		       "ecalvalidation/h_HF_noise_vs_Eta",
		       "ecalvalidation/h_LF_noise_vs_Eta",//110
		       "ecalvalidation/h_Total_noise_vs_Eta",
		       "ecalvalidation/h_HF_noise_FromRechit_vs_Eta_ped",
		       "ecalvalidation/h_LF_noise_FromRechit_vs_Eta_ped",
		       "ecalvalidation/h_Total_noise_FromRechit_vs_Eta_ped", 
		       "ecalvalidation/h_recHits_EB_iPhiOccupancy",
		       "ecalvalidation/h_recHits_EB_iEtaOccupancy",
		       "ecalvalidation/h_Amplitude_vs_Eta",
		       "ecalvalidation/h_Amplitude_FromRechit_vs_Eta",
		       "ecalvalidation/h_recHits_EB_sumEt",  
		       "ecalvalidation/h_recHits_EEP_sumEt", //120
		       "ecalvalidation/h_recHits_EEM_sumEt",  
		       "ecalvalidation/h_recHits_EEP_sumEtCut",  
		       "ecalvalidation/h_recHits_EEM_sumEtCut",
		       "ecalvalidation/h_HF_noise_vs_Eta_ped",
		       "ecalvalidation/h_LF_noise_vs_Eta_ped",
		       "ecalvalidation/h_Total_noise_vs_Eta_ped",
		       "ecalvalidation/h_recHits_EB_SRP",
		       "ecalvalidation/h_recHits_EE_SRP"};        //128               
  
  char *objName_2D[nObj_2D]={"ecalvalidation/h_recHits_EB_occupancy",
			     "ecalvalidation/h_recHits_EEP_occupancy",
			     "ecalvalidation/h_recHits_EEM_occupancy",
			     "ecalvalidation/h_recHits_ES_occupancy_F+",
			     "ecalvalidation/h_recHits_ES_occupancy_F-",
			     "ecalvalidation/h_recHits_ES_occupancy_R+",
			     "ecalvalidation/h_recHits_ES_occupancy_R-",
			     "ecalvalidation/h_recHits_ES_occupancy_F+_gr",
			     "ecalvalidation/h_recHits_ES_occupancy_F-_gr",
			     "ecalvalidation/h_recHits_ES_occupancy_R+_gr",
			     "ecalvalidation/h_recHits_ES_occupancy_R-_gr",
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
  
  
  char *objTitle[nObj]={"Number of RecHits (EB)",//-1
			"Number of RecHits (EE+)",
			"Number of RecHits (EE-)",
			"Number of RecHits (ES)",
			"Number of RecHits (ES goodRH)",
			"Number of RecHits (ES F+)",
			"Number of RecHits (ES F-)",
			"Number of RecHits (ES R+)",
			"Number of RecHits (ES R-)",
			"RecHits Energy (EB)",
			"RecHits Energy (EE+)",
			"RecHits Energy (EE-)",//10
			"RecHits Energy (ES)",
			"RecHits Energy (ES good RH)",
			"RecHits Energy (ES F+)",
			"RecHits Energy (ES F-)",
			"RecHits Energy (ES R+)",
			"RecHits Energy (ES R-)",
			"RecHits Max Energy (EB)",
			"RecHits Max Energy (EE+)",
			"RecHits Max Energy (EE-)",
			"RecHits Max Energy (ES)",//20
			"RecHits Time (EB)",
			"RecHits Time (EE+)",
			"RecHits Time (EE-)",
			"RecHits Time (ES)",
			"RecHits Time (ES good RH)",
			"RecHits Time (ES F+)",
			"RecHits Time (ES F-)",
			"RecHits Time (ES R+)",
			"RecHits Time (ES R-)",
			"Number of RecHits (EB spike cleaned)", //30
			"RecHits Energy (EB spike cleaned)",
			"RecHits Max Energy (EB spike cleaned)",
			"RecHits Time (EB spike cleaned)",
			"RecHits Eta",
			"RecHits Phi (EB)",
			"RecHits Phi (EE)", 
			"RecHits \#chi^{2} (EB)",
			"RecHits \#chi^{2} (EE+)",
			"RecHits \#chi^{2} (EE-)",
			"RecHits 1-E4/E1 (EB)",//40
			//Leo
			"Number of RecHits (EB)",//41
			"Number of RecHits (EE)",
			"Number of RedRecHits (EB)",
			"Number of RedRecHits (EE)",
			"Number of ClusRecHits (EB)",
			"Number of ClusRecHits (EE)",//46
			//Leo
			"Number of Basic Clusters (EB)",
			"Number of Basic Clusters (EE+)",
			"Number of Basic Clusters (EE-)",
			"Number of Superclusters (EB)",//50
			"Number of Superclusters (EE+)",
			"Number of Superclusters (EE-)",
			"Number of Crystal per Supercluster (EB)",
			"Number of Crystal per Supercluster (EE+)",
			"Number of Crystal per Supercluster (EE-)",
			"Number of Basic Clusters (EB spike cleaned)",//56
			"Number of Superclusters (EB spike cleaned)",
			"Number of Crystal per Supercluster (EB spike cleaned)",
			"Number of Basic Clusters  per Supercluster (EB spike cleaned)",
			"Supercluster Energy (EB spike cleaned)",//60
			"Number of Basic Clusters  per Supercluster (EB)",
			"Number of Basic Clusters  per Supercluster (EE-)",
			"Number of Basic Clusters  per Supercluster (EE+)",
			"Supercluster Energy (EB)",
			"Supercluster Energy (EE+)",
			"Supercluster Energy (EE-)",
			"Superclusters Eta",
			"Superclusters Phi (EB)",
			"Superclusters Phi (EE)",
			"1-E4/E1",//70
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
			"RecoDigi Pedestal Mean (EE+)",//80
			"RecoDigi Pedestal Mean (EE-)",
			"SimDigi Pedestal RMS (EB)",
			"RecoDigi Pedestal RMS (EB)",
			"SimDigi Pedestal RMS (EE+)",
			"SimDigi Pedestal RMS (EE-)",
			"RecoDigi Pedestal RMS (EE+)",
			"RecoDigi Pedestal RMS (EE-)",
			"SimDigi HF-noise (EB)", 
			"SimDigi LF-noise (EB)",
			"SimDigi Total-noise (EB)",//90
			"RecoDigi HF-noise (EB)", 
			"RecoDigi LF-noise (EB)",
			"RecoDigi Total-noise (EB)",
			"SimDigi HF-noise (EE+)",
			"SimDigi LF-noise (EE+)",
			"SimDigi Total-noise (EE+)",
			"SimDigi HF-noise (EE-)",
			"SimDigi LF-noise (EE-)",
			"SimDigi Total-noise (EE-)",
			"RecoDigi HF-noise (EE+)", //100
			"RecoDigi LF-noise (EE+)",
			"RecoDigi Total-noise (EE+)",
			"RecoDigi HF-noise (EE-)", 
			"RecoDigi LF-noise (EE-)",
			"RecoDigi Total-noise (EE-)",
			"RecoDigi HF-noise vs Eta",
			"RecoDigi LF-noise vs Eta",
			"RecoDigi Total noise vs Eta",
			"SimDigi HF-noise vs Eta",
			"SimDigi LF-noise vs Eta",//110
			"SimDigi Total noise vs Eta",
			"RecoDigi HF-noise vs Eta (pedestal)",
			"RecoDigi LF-noise vs Eta (pedestal)",
			"RecoDigi Total noise vs Eta (pedestal)",
			"RecHits iPhiOccupancy",
			"RecHits iEtaOccupancy",
			"SimDigi Amplitude",
			"RecoDigi Amplitude",
			"RecHits SumEt (EB)",
			"RecHits SumEt (EE+)",//120
			"RecHits SumEt (EE-)", 
			"RecHits SumEt (EE+, |eta| < 2.5)",
			"RecHits SumEt (EE-, |eta| < 2.5)",
			"SimDigi HF-noise vs Eta (pedestal)",
			"SimDigi LF-noise vs Eta (pedestal)",
			"SimDigi Total noise vs Eta (pedestal)",
			"Rechit SRP (EB)",
			"Rechit SRP (EE)"};//128
  
  
  char *objTitle_2D[nObj_2D]={"RecHits Occupancy (EB)", 
			      "RecHits Occupancy (EE+)",
			      "RecHits Occupancy (EE-)",
			      "RecHits Occupancy (ES F+)",
			      "RecHits Occupancy (ES F-)",
			      "RecHits Occupancy (ES R+)",
			      "RecHits Occupancy (ES R-)",
			      "RecHits Occupancy (ES F+ goodRH)",
			      "RecHits Occupancy (ES F- goodRH)",
			      "RecHits Occupancy (ES R+ goodRH)",
			      "RecHits Occupancy (ES R- goodRH)",
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
  
  char *labelX[nObj]={"Number of RecHits/Event",//-1
		      "Number of RecHits/Event",
		      "Number of RecHits/Event",
		      "Number of RecHits/Event",
		      "Number of RecHits/Event",
		      "Number of RecHits/Event",
		      "Number of RecHits/Event",
		      "Number of RecHits/Event",
		      "Number of RecHits/Event",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",//10
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",//20
		      "time(ns)",
		      "time(ns)",
		      "time(ns)",
		      "time(ns)",
		      "time(ns)",
		      "time(ns)",
		      "time(ns)",
		      "time(ns)",
		      "time(ns)",
		      "Number of RecHits/Event", //30
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "time(ns)",
		      "eta",
		      "phi",
		      "phi",
		      "#chi^{2}",
		      "#chi^{2}",
		      "#chi^{2}",
		      "1-E4/E1",//40
		      "recoFlag", //41
		      "recoFlag",
		      "recoFlag",
		      "recoFlag",
		      "recoFlag",
		      "recoFlag", //46
		      "BasicClusters/Event",
		      "BasicClusters/Event",
		      "BasicClusters/Event",
		      "Superclusters/Event", //50
		      "Superclusters/Event",
		      "Superclusters/Event",
		      "BasicClusters/Event",
		      "Superclusters/Event",
		      "Crystals/Supercluster",
		      "BasicClusters/Supercluster", //56
		      "Energy (GeV)",
		      "Crystals/Supercluster",
		      "Crystals/Supercluster",
		      "Crystals/Supercluster",  //60
		      "BasicClusters/Supercluster",
		      "BasicClusters/Supercluster",
		      "BasicClusters/Supercluster",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "eta",
		      "phi",
		      "phi",
		      "1-E4/E1", //70
		      "1-E4/E1",
		      "1-E4/E1",
		      "Energy (GeV)",
		      "Energy (GeV)",
		      "EnergyPlane1/EnergyPlane2", //75
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC", //80
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC", // inizio noise
		      "ADC",
		      "ADC",
		      "ADC",//90
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",//100
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "ADC",
		      "Eta",
		      "Eta",
		      "Eta",
		      "Eta",
		      "Eta",//110
		      "Eta",
		      "Eta",
		      "Eta",
		      "Eta",
		      "iPhi",
		      "iEta",
		      "Eta",
		      "Eta",
		      "GeV",
		      "GeV",//120
		      "GeV", 
		      "GeV",
		      "GeV",
		      "Eta",
		      "Eta",
		      "Eta",
		      "",
		      ""}; //128
  
  
  char *labelY[1]={"Counts"};
  
  char *labelX_2D[nObj_2D]={"iPhi",
			    "iX",
			    "iX",
			    "iX",
			    "iX",
			    "iX",
			    "iX",
			    "iX",
			    "iX",
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
			    "iY",
			    "iY",
			    "iY",
			    "iY",
			    "iY",
			    "iY",
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

  double xMin[nObj]={0.,0.,0.,0.,   //-1 
		     0.,0.,0.,0., 
		     0.,0.,
		     0.,0.,0.,0., 
		     0.,0.,0.,0., //13
		     0.,0.,0.,0.,
		     -100.,-50.,-50.,-50., //21
		     -50.,
		     -50.,-50.,-50.,-50.,
		     //cleaned
		     0., 0., 0., -100.,  //30
		     -3.,-3.2,-3.2,
		     0,0,0,
		     0, //40
		     0,0, //41
		     0,0,
		     0,0, //46
		     0,0,0,
		     0,0,0,//50
		     0,0,0,
		     0, 0, 0, 0, 0,
		     0,0,0,//61
		     0,0,0,
		     -3.,-3.2,-3.2,
		     0,0,0,//70
		     0,0,0, 
		     175,175,175,175,175,175,//76
		     -1,-1,-1,-1,-1,-1,//82
		     -12,-12,-12,-12,-12,
		     -12,-12,-12,-12,-12,//93
		     -12,-12,-12,-12,-12,
		     -12,-12,-12,//103
		     -3,-3,-3,
		     -3,-3,-3,//109
		     -3,-3,-3,
		     1.,-86,//115
		     -3,-3,
		     0,
		     0,0,//120
		     0,0,
		     -3,-3,-3,
		     0,0};//128
  
  double xMin_2D[nObj_2D]={1,0,0,
			   0,0,0,0,
			   0,0,0,0,
                           1,0,0, 
                           1,1,1, 
                           0,0,0, 
                           0,0,0,
                           1,1,
                           0,0,0,0};           
  
  double xMax[nObj]={10000.,1700.,1700.,8000.,   //-1
		     8000.,  
		     2000., 2000., 2000., 2000.,
		     300.,300.,300.,0.06, 
		     0.06, //12
		     0.06, 0.06, 0.06, 0.06, 
		     300.,300.,300.,0.06,
		     100.,50.,50.,50.,//21
		     50.,//25
		     50.,50.,50.,50.,
		     10000., 200., 200., 100.,//30
		     3.,3.2,3.2, 
		     70,70,70,
		     1.2,//40
		     32,32,//41
		     32,32,
		     32,32,
		     80,150,150, //47
		     60,60,60,
		     200,200,200,
		     80., 60., 200., 20., 400., 
		     20,20,20,//61
		     500,500,500,
		     3.,3.2,3.2,//67
		     1.2,1.2,1.2,
		     0.05,0.05,10,
		     220,220,220,220,220,220, //76
		     11,11,11,11,11,11,
		     12,12,12,12,12,//88
		     12,12,12,12,12,
		     12,12,12,12,12,12,12,12, 
		     3,3,3,//106
		     3,3,3,
		     3,3,3,//112
		     361,86,
		     3,3,
		     1400,//119
		     150,150,
		     100,100,
		     3,3,3,
		     3,3};//128
  
  double xMax_2D[nObj_2D]={361,100,100,   
			   41., 41.,41., 41.,
			   41., 41.,41., 41.,
                           361,100,100,
                           361,361,361,
                           100,100,100,
                           100,100,100,
                           361,361,
                           100,100,100,100};
  
  double yMin_2D[nObj_2D]={-86,0,0,
			   0,0,0,0,
			   0,0,0,0,
                           -86,0,0,
                           -86,-86,-86,
			   0,0,0,
			   0,0,0,
			   -86,-86,
			   0,0,0,0};
  
  double yMax_2D[nObj_2D]={86,100,100,
			   41., 41.,41., 41.,
			   41., 41.,41., 41.,
                           86,100,100,
                           86,86,86,
                           100,100,100,
                           100,100,100,
                           86,86,
                           100,100,100,100};

  double zMin_2D[nObj_2D]={0,0,0,
			   0,0,0,0,
			   0,0,0,0,
                           1000,1000,1000,
                           0,0,0,
                           0,0,0,
                           0,0,0,
                           0,0,
                           0,0,0,0};
  
  double zMax_2D[nObj_2D]={10000,10000,10000,
			   200000,200000,200000,200000,
			   200000,200000,200000,200000,
                           1020,1020,1020,
                           6,6,6,
                           8,8,8,
                           8,8,8,
                           5,5,
                           8,8,8,8};

  int reBin[500]  = {20,2,2,2,//-1
		     2,
		     2,2,2,2, 
		     20,20,20,8, 
		     8,//12
		     8, 8, 8, 8, 
		     4,4,4,4, 
		     4,//21
		     1,1,
		     4, //24 
		     4,
		     4,4,4, 
		     4, //29
		     2, 2, 8, //30
		     5,5,5, 
		     5,5,5,
		     1,//40
		     1,1,
		     1,1,
		     1,1,
		     1,1,1,//47
		     1,1,1,
		     1,1,1,
		     1,1,1,1,1,//56
		     1,1,1,//60
		     1,10,10,
		     5,5,5,
		     5,1,1,//70
		     1,10,10,
		     1,1,1,1,1,
		     1,1,1,1,1,//81
		     1,1,1,1,1,
		     1,1,1,1,1,//91
		     1,1,1,1,1,
		     1,1,1,1,1,1,1,1,1, //101
		     1,1,1,1,1,1,1,1,1,1,//110
		     1,1,1,1,1,4,1,1,1};//120
  

  int optLogY[nObj] = {1,1,1,1,//-1
		       1,1,1,1,
		       1,1,
		       1,1,1,1,
		       1,1,1,1,
		       1,1,1,1,//17
		       1,0,0,0,
		       0,//25
		       0,0,0,0,
		       1,1,1,1,
		       0,0,0,
		       1,1,1,
		       0,//40
		       1,1,
		       1,1,
		       1,1,
		       1,1,1,
		       1,1,1,1,1,//50
		       1,1,1,
		       1,1,1,
		       1,1,1,//61
		       1,1,1,
		       0,0,0,//67
		       0,0,0,
		       0,0,0,//73
		       0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0, //80
		       0,0,0,0,0,0,0,0,0,0, //90
 		       0,0,0,0,0,0,0,0,0,0, //100
		       0,0,0,0,0,0,0,0,0,0, //110
		       0,0,0,0,0,0,0,0,0}; //120
  

  TH1D* h[2][500]; 
  TCanvas *c[500];
  TPaveStats *st[500];
  
  bool isMissing = false;

  int iHisto = 0;
  while(iHisto<nObj){
    //      std::cout << " 1D = " << std::string(objName[iHisto]) << std::endl;
    for (int ifile=0;ifile<2;ifile++){ 

      std::cout << " ifile = " << ifile << " iHisto = " << iHisto << std::endl;
      h[ifile][iHisto] = (TH1D*)_f[ifile]->Get(objName[iHisto]);
      if(h[ifile][iHisto] == 0)isMissing = true;
      if(h[ifile][iHisto] == 0) continue;
          
      h[ifile][iHisto]->Rebin(reBin[iHisto]);
      
      if (ifile == 0) {
	// open a new canvas
	c[iHisto] = new TCanvas(objName[iHisto],objName[iHisto],50+iHisto*20,50+iHisto*5,500,400);
	c[iHisto]->cd();

	// customize and plot
	h[ifile][iHisto]->GetYaxis()->SetTitle(labelY[0]);
	std::string histo_name = std::string(h[ifile][iHisto]->GetName());
	if(histo_name.find("Eta") != std::string::npos) h[ifile][iHisto]->GetYaxis()->SetTitle("ADC");
	if ( iHisto > 41 && iHisto < 48 ) //Set reasonable labels for recoFlag plots
	  {
	    int nBin = h[ifile][iHisto]->GetXaxis()->GetNbins();
	    for ( int ibin = 1; ibin <= nBin; ibin++ ) 
	      h[ifile][iHisto]->GetXaxis()->SetBinLabel( ibin, recoFlagLabels[ibin-1]);
	  }
	else h[ifile][iHisto]->GetXaxis()->SetTitle(labelX[iHisto]);
	h[ifile][iHisto]->SetFillColor(kBlack+10);
	h[ifile][iHisto]->SetFillStyle(3002);
	h[ifile][iHisto]->SetTitle(objTitle[iHisto]);
	
	h[ifile][iHisto]->GetXaxis()->SetRangeUser(xMin[iHisto],xMax[iHisto]);
	h[ifile][iHisto]->Draw();
      }
      
      h[ifile][iHisto]->Scale(1. * h[0][iHisto]->GetEntries() / 1. / h[ifile][iHisto]->GetEntries() );
      std::cout << " ifile = " << ifile << " name = " << h[ifile][iHisto]->GetName() << " n " <<  h[ifile][iHisto]->GetEntries() << std::endl;
      
      if (ifile == 1){
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
  
  
  std::cout << " AND NOW 2D " << std::endl;

  TH2D* h2[2][100]; 
  TCanvas *c2[2][100];
  TPaveStats *st2[2][100];
  
  bool  isMissing_2D[2][100];
  for(int ii = 0; ii < 2; ii++)
    for(int jj = 0; jj < nObj_2D; jj++){
      isMissing_2D[ii][jj] = false;
    }
  
  int iHisto2D = 0;
  while(iHisto2D<nObj_2D){
    std::cout << " 2D = " << std::string(objName_2D[iHisto2D]) << std::endl;

    for (int ifile=0;ifile<2;ifile++){ 
      h2[ifile][iHisto2D] = (TH2D*)_f[ifile]->Get(objName_2D[iHisto2D]);
      if(h2[ifile][iHisto2D] == 0) isMissing_2D[ifile][iHisto2D]  = true;
      if(h2[ifile][iHisto2D] == 0) continue;

      std::cout << " >>> ciao1 " << std::endl;    

      std::cout << " >>> ciao7 " << std::endl;
      char c_name[500];
      sprintf (c_name,"%s_%d_%d",objName_2D[iHisto2D],iHisto2D,ifile);
      c2[ifile][iHisto2D] = new TCanvas(c_name,c_name,50+iHisto2D*20,50+iHisto2D*5,500,400);
      c2[ifile][iHisto2D]->cd();
      // customize and plot
      h2[ifile][iHisto2D]->Scale(h[0][0]->GetEntries()/h[ifile][0]->GetEntries());
      h2[ifile][iHisto2D]->GetXaxis()->SetTitle(labelX_2D[iHisto2D]);
      h2[ifile][iHisto2D]->GetYaxis()->SetTitle(labelY_2D[iHisto2D]);
      h2[ifile][iHisto2D]->GetXaxis()->SetRangeUser(xMin_2D[iHisto2D],xMax_2D[iHisto2D]);
      h2[ifile][iHisto2D]->GetYaxis()->SetRangeUser(yMin_2D[iHisto2D],yMax_2D[iHisto2D]);
      h2[ifile][iHisto2D]->GetZaxis()->SetRangeUser(zMin_2D[iHisto2D],zMax_2D[iHisto2D]);
      h2[ifile][iHisto2D]->Draw("colz");
      c2[ifile][iHisto2D]->Update();
      // stat. box
      st2[ifile][iHisto2D]= (TPaveStats*)(h2[ifile][iHisto2D]->GetListOfFunctions()->FindObject("stats"));
      /*
	st2[ifile][iHisto2D]->SetX1NDC(0.71); //new x start position
	st2[ifile][iHisto2D]->SetX2NDC(0.99); //new x end position
	st2[ifile][iHisto2D]->SetY1NDC(0.84); //new y start position
	st2[ifile][iHisto2D]->SetY2NDC(0.99); //new y end position
      */
      st2[ifile][iHisto2D]->SetX1NDC(0.); //new x start position
      st2[ifile][iHisto2D]->SetX2NDC(0.); //new x end position
      st2[ifile][iHisto2D]->SetY1NDC(0.); //new y start position
      st2[ifile][iHisto2D]->SetY2NDC(0.); //new y end position
      //st[ifile][iHisto2D]->SetTextColor(kRed);
      st2[ifile][iHisto2D]->Draw();

      std::cout << " >>> ciao8 " << std::endl;
      char dummyName_2D[500];
      strcat(dummyName_2D, objName_2D[iHisto2D]);
      char *str = strtok (dummyName_2D,"/");
      str = strtok (NULL, "/");

      char myname[500];
      sprintf (myname,"%s/",dirName);
      strcat(myname,str);
      
      if(ifile == 0) strcat(myname,"_0");
      if(ifile == 1) strcat(myname,"_1");
      
      strcat(myname,".");
      strcat(myname,fileType);
      
      if(isMissing_2D[ifile][iHisto2D] == false) c2[ifile][iHisto2D]->Print(myname,fileType);
      isMissing_2D[ifile][iHisto2D] = false;


      //ratio plots
      //if(ifile == 1 && (iHisto2D == 0 || iHisto2D == 1 || iHisto2D == 2 || iHisto2D == 7 || iHisto2D == 8 || iHisto2D == 9 || iHisto2D == 10)){
      if(ifile == 1 && (iHisto2D == 7 || iHisto2D == 8 || iHisto2D == 9 || iHisto2D == 10)){
	TH2D* histoRatio = (TH2D*)h2[0][iHisto2D]->Clone((std::string(objName_2D[iHisto2D])+"_DAoMC_ratio").c_str());
	histoRatio->Reset();
	histoRatio->Divide(h2[0][iHisto2D], h2[1][iHisto2D], 1, 1);

	std::cout << " 2Dratio histoRatio->GetName() = " << histoRatio->GetName() << std::endl;

	char c_nameRatio[500];
	sprintf (c_nameRatio,"%s_%d_%d_%s",objName_2D[iHisto2D],iHisto2D,ifile, "_DAoMC_ratio");
	c2[ifile][iHisto2D] = new TCanvas(c_nameRatio,c_nameRatio,50+iHisto2D*20,50+iHisto2D*5,500,400);
	c2[ifile][iHisto2D]->cd();
	std::cout << " >>> ciao2 " << std::endl;
	// customize and plot
	histoRatio->GetXaxis()->SetTitle(labelX_2D[iHisto2D]);
	histoRatio->GetYaxis()->SetTitle(labelY_2D[iHisto2D]);
	histoRatio->GetXaxis()->SetRangeUser(xMin_2D[iHisto2D],xMax_2D[iHisto2D]);
	histoRatio->GetYaxis()->SetRangeUser(yMin_2D[iHisto2D],yMax_2D[iHisto2D]);
	std::cout << " >>> ciao3 " << std::endl;
	double zMinR = 0.;
	double zMaxR = 10.;
	// if(iHisto2D == 0) histoRatio->GetZaxis()->SetRangeUser(zMin_2D[iHisto2D], 5.);
	// else histoRatio->GetZaxis()->SetRangeUser(zMin_2D[iHisto2D], 7.);

	//if(iHisto2D == 0 || iHisto2D == 1 || iHisto2D == 2)
	//histoRatio->GetZaxis()->SetRangeUser(0.5, 1.5);
	//else
	histoRatio->GetZaxis()->SetRangeUser(0., 10.);

	std::cout << " >>> ciao4 " << std::endl;
	
	histoRatio->Draw("colz");
	c2[ifile][iHisto2D]->Update();
	// stat. box
	st2[ifile][iHisto2D]= (TPaveStats*)(histoRatio->GetListOfFunctions()->FindObject("stats"));
     
	st2[ifile][iHisto2D]->SetX1NDC(0.); //new x start position
	st2[ifile][iHisto2D]->SetX2NDC(0.); //new x end position
	st2[ifile][iHisto2D]->SetY1NDC(0.); //new y start position
	st2[ifile][iHisto2D]->SetY2NDC(0.); //new y end position
	//st[ifile][iHisto2D]->SetTextColor(kRed);
	st2[ifile][iHisto2D]->Draw();
      	std::cout << " >>> ciao5 " << std::endl;

	
	char dummyName_2D_ratio[500];
	strcat(dummyName_2D_ratio, objName_2D[iHisto2D]);
	strcat(dummyName_2D_ratio, "_DAoMC_ratio");
	char *strR = strtok (dummyName_2D_ratio,"/");
	strR = strtok (NULL, "/");

	char mynameR[500];
	sprintf (mynameR,"%s/",dirName);
	strcat(mynameR,strR);
      
	strcat(mynameR,".");
	strcat(mynameR,fileType);
      	std::cout << " >>> ciao6 " << std::endl;
	c2[ifile][iHisto2D]->Print(mynameR,fileType);
	
      }


    }
  
    iHisto2D++;
  }

  std::cout << " in teoria finito " << std::endl;


}

