//
// Macro to produce ECAL Pi0 plots
//

// int Wait() {
//      cout << " Continue [<RET>|q]?  ";
//      char x;
//      x = getchar();
//      if ((x == 'q') || (x == 'Q')) return 1;
//      return 0;
// }

#include <algorithm>

void DrawValidationPlotsPi0(Char_t* infile1 = 0, 
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
  int nObj=2;
  char *objName[2]={"ecalvalidation/Pi0/h_Pi0_EB_mass",
		    "ecalvalidation/Pi0/h_Pi0_EE_mass"};

  char *objTitle[2]={"Pi0 peak (EB)",
		     "Pi0 peak (EE)"};

  char *labelX[2]={"mass (GeV)",
		   "mass (GeV)"};
		   
  char *labelY[1]={"Counts"};
  
  double xMin[2]={0.06,0.06};

  double xMax[2]={0.300,0.300};

  int reBin[2]  ={1,1};

  int optLogY[2] = {0,0};

  TH1D* h[2][10];
  TCanvas *c[10];
  TPaveStats *st[10];

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
      h[ifile][iHisto]->GetXaxis()->SetTitle(labelX[iHisto]);
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
      h[ifile][iHisto]->Draw("esames");
     
      // set the log-scale and range 
      float maxy = max (h[ifile][iHisto]->GetMaximum(),h[ifile-1][iHisto]->GetMaximum() );
      if (optLogY[iHisto]) {
	c[iHisto]->SetLogy();
	h[0][iHisto]->SetMaximum(maxy*2);
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

