#!/usr/bin/perl

# ----------------------------------------------------------------------------
#      MAIN PROGRAM
# ----------------------------------------------------------------------------

#PG lettura dei parametri da cfg file
#PG --------------------------------
print "reading ".$ARGV[0]."\n" ;

open (USERCONFIG,$ARGV[0]) ;

while (<USERCONFIG>)
  {
    chomp; 
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
#    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences{$var} = $value;
  }


$LOCALFolder      = $User_Preferences{"LOCALFolder"} ;
$CMSSWCfgTemplate = $User_Preferences{"CMSSWCfgTemplate"} ;
$RUNFILEList      = $User_Preferences{"RUNFILEList"} ;
$INPUTSAVEFolder  = $User_Preferences{"INPUTSAVEFolder"} ;
$OUTPUTSAVEPath   = $User_Preferences{"OUTPUTSAVEPath"} ;
$OUTPUTSAVEFolder = $User_Preferences{"OUTPUTSAVEFolder"} ;
$OUTPUTFileName   = $User_Preferences{"OUTPUTFileName"} ;
$JOBModulo        = $User_Preferences{"JOBModulo"} ;

print "LOCALFolder = "      .$LOCALFolder."\n" ;
print "CMSSWCfgTemplate = " .$CMSSWCfgTemplate."\n" ;
print "RUNFILEList = "      .$RUNFILEList."\n" ;
print "INPUTSAVEFolder = "  .$INPUTSAVEFolder."\n" ;
print "OUTPUTSAVEPath = "   .$OUTPUTSAVEPath."\n" ;
print "OUTPUTSAVEFolder = " .$OUTPUTSAVEFolder."\n" ;
print "OUTPUTFileName = "   .$OUTPUTFileName."\n" ;
print "JOBModulo = "        .$JOBModulo."\n" ;






##########################
# create folders in castor
##########################    

print ("\nrfmkdir ".$OUTPUTSAVEPath.$OUTPUTSAVEFolder."\n") ;
system ("rfmkdir ".$OUTPUTSAVEPath.$OUTPUTSAVEFolder."\n") ;
print ("rfchmod 777 ".$OUTPUTSAVEPath.$OUTPUTSAVEFolder."\n") ;
system ("rfchmod 777 ".$OUTPUTSAVEPath.$OUTPUTSAVEFolder."\n") ;





################
# loop over runs 
################

open (RUNFILEList,$RUNFILEList) ;
while (<RUNFILEList>)
{
  $run = $_;
  chomp ($run);
  print("\n----- RUN ".$run." -----\n");
  
  $currDir = `pwd` ;
  chomp ($currDir) ;
  $currDir = $currDir."/run".$run ;
  system ("mkdir ".$currDir) ;
  
  
   
  $tempo1 = "./tempo1" ;
  system ("cat ".$CMSSWCfgTemplate."   | sed -e s%LOCALFOLDER%".$LOCALFolder.
                                   "%g | sed -e s%RUN%".$run.
                                   "%g | sed -e s%OUTPUTFILENAME%".$OUTPUTFileName. 
                                   "%g > ".$tempo1) ;
  
  
  $LISTFile = $currDir."/list.txt";
  
  $command = "dbs search --query='find file where dataset=/SingleElectron/Run2011A-PromptReco-v4/RECO and run=".$run."' | grep root > ".$LISTFile;
  system($command);
  
  $listOfFiles;   
  $it = 0;
  open (LISTFile,$LISTFile) ;
  while (<LISTFile>)
  {
    chomp; 
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
    
    $file = $_ ;
    $listOfFiles = $listOfFiles."APICErfio:".$INPUTSAVEFolder."/".$file."APICE," ;
  }
  
  chop($listOfFiles) ;
  #print $listOfFiles."\n" ;
  
  $tempo2 = "./tempo2" ;    
  system ("cat ".$tempo1." | sed -e s%LISTOFFILES%".$listOfFiles."%g > ".$tempo2) ;
  $listOfFiles = "" ;
  
  $tempo3 = "./tempo3" ;
  system ("cat ".$tempo2." | sed -e s%APICE%\\'%g > ".$tempo3) ;
  
  
  
  
  
  
  $CMSSWCfgFile = $currDir."/".$run."_cfg.py" ;
  system ("mv ".$tempo3." ".$CMSSWCfgFile) ;
  system ("rm ./tempo*") ;
  
  
  
  $tempBjob = $currDir."/bjob_".$run.".sh" ;
  $command = "touch ".$tempBjob ;
  system ($command) ;
  $command = "chmod 777 ".$tempBjob ;
  system ($command) ;
  
  
  $command = "cd ".$currDir ;
  system ("echo ".$command." > ".$tempBjob) ;
  
  $command = "export SCRAM_ARCH=slc5_amd64_gcc434" ;
  system ("echo ".$command." >> ".$tempBjob) ;
  
  $command = "eval \\\`scramv1 runt -sh\\\`" ;
  system ("echo ".$command." >> ".$tempBjob) ;

  $command = "cd -" ;
  system ("echo ".$command." >> ".$tempBjob) ;
  
  $command = "cmsRun ".$CMSSWCfgFile ;
  system ("echo ".$command." >> ".$tempBjob) ;
  
  
  
  ######################
  # copy files to castor
  ######################    
  
  $command = "rfcp ./".$OUTPUTFileName."_run".$run.".root ".$OUTPUTSAVEPath.$OUTPUTSAVEFolder; 
  system ("echo ".$command." >> ".$tempBjob) ;

  
  
  ############
  # submit job
  ############
  
  system ("cd ".$currDir." \n") ;
  print ("bsub -q 1nd -cwd ".$currDir." ".$tempBjob."\n") ;
  system ("bsub -q 1nd -cwd ".$currDir." ".$tempBjob) ;
  
  print "\n" ;
}
  
