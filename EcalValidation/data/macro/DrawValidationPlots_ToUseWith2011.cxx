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

void DrawValidationPlots(Char_t* infile1 = 0, 
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


  // Define list of object names
  const int nObj=53;
  char *objName[nObj]={"ecalvalidation/h_recHits_EB_size",
		     "ecalvalidation/h_recHits_EEP_size",
		     "ecalvalidation/h_recHits_EEM_size",
		     "ecalvalidation/h_recHits_ES_size",
		     //"ecalvalidation/h_recHits_EB_energy_cleaned",
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
		     "ecalvalidation/h_esClusters_energy_ratio"  };

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
		     "ES Clusters Energy - Ratio"};

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
		    "EnergyPlane1/EnergyPlane1"};
                    
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

  char *labelY[1]={"Counts"};

  double xMin[nObj]={0.,0.,0.,0., 
		   0,0,0,-10, 
		   1.,1.,1.,1.,
		   -50.,-50.,-50.,-50.,
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
		   0,0,0};
  
  double xMax[nObj]={6000.,1700.,1700.,4000.,  
		   40,40.,40.,0.5, 
		   80,100,100,0.1,
		   50.,50.,50.,50.,
		   3.,3.2,3.2,
		   65,65,65,
		   1.2,
                   32,32,
                   32,32,
                   32,32,
                   80,80,80,
		   40,20,20,
		   100,100,100,
		   20,20,20,
		   300,300,300,
		   3.,3.2,3.2,
		   1.2,1.2,1.2,
		   0.02,0.02,30};

  int reBin[nObj]  = {2,2,2,2, 
		    1,1,1,8, 
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
		    10,10,1};


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


  TH1D* h[2][100]; 
  TCanvas *c[100];
  TPaveStats *st[100];
 
  int iHisto = 0;
  while(iHisto<nObj){
    for (int ifile=0;ifile<2;ifile++){ 
    h[ifile][iHisto] = (TH1D*)_f[ifile]->Get(objName[iHisto]);
    h[ifile][iHisto]->Rebin(reBin[iHisto]);
    if (ifile == 0) {
      // open a new canvas
      c[iHisto] = new TCanvas(objName[iHisto],objName[iHisto],50+iHisto*20,50+iHisto*5,500,400);
      c[iHisto]->cd();
      // customize and plot
      h[ifile][iHisto]->GetYaxis()->SetTitle(labelY[0]);
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

      h[ifile][iHisto]->SetMarkerStyle(20);
      h[ifile][iHisto]->SetMarkerSize(0.7);
      h[ifile][iHisto]->SetMarkerColor(kRed);
      h[ifile][iHisto]->Draw("epsames");
     
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
      st[iHisto]->SetY1NDC(0.72); //new x start position
      st[iHisto]->SetY2NDC(0.85); //new x end position
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

  c[iHisto]->Print(myname,fileType);
  
  
  iHisto++;
  
  
  }



 
}

