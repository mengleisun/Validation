#!/bin/bash

usage='Usage: -d <file1> -r <file2> -l <html dir> -o <http dir> -w <which dir>'



args=`getopt rd: -- "$@"`
if test $? != 0
     then
         echo $usage
         exit 1
fi



eval set -- "$args"
for i
  do
  case "$i" in
      -d) shift; file1=$2;shift;;
      -r) shift; file2=$2;shift;;
      -l) shift; htmldir=$2;shift;;
      -o) shift; httpdir=$2;shift;;
      -w) shift; whichdir=$2;shift;;
  esac
done

if [ "X"${file1} == "X" ]
    then
    echo "MISSING INPUT FILE NAME for the first dataset! Please give a valid one!"
    echo $usage
    exit
fi

if [ "X"${file2} == "X" ]
    then
    echo "MISSING INPUT FILE NAME for the second dataset! Please give a valid one!"
    echo $usage
    exit
fi

if [ "X"${htmldir} == "X" ]
    then
    htmldir=/afs/cern.ch/user/c/ccecal/scratch0/www/ValidationPlots/MC_vs_MC_noiseStudies
    echo "using default htmldir:" ${htmldir}
fi


if [ "X"${httpdir} == "X" ]
    then
    httpdir=http://test-ecal-cosmics.web.cern.ch/test-ecal-cosmics/ValidationPlots
    echo "using default httpdir:" ${httpdir}
fi

if [ "X"${whichdir} == "X" ]
    then
    whichdir=Data
    echo "using default httpdir:" ${httpdir}
fi


httpdir=${httpdir}/${whichdir}
htmldir=${htmldir}/${whichdir}


echo 'Preparing Validation Webpages' 

echo 'httpdir ' ${httpdir}
echo 'htmldir ' ${htmldir}

# specify directories here
cd ../../../
my_cmssw_base=`\pwd`
work_dir=${my_cmssw_base}/Validation/EcalValidation
echo $work_dir
cd -

out_dir=${file1}_vs_${file2};
echo ${out_dir}

plots_dir=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/ccecal/ValidationPlots/${whichdir}/${out_dir};

mkdir /afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/ccecal/ValidationPlots/${whichdir}
mkdir ${plots_dir}

cp ${work_dir}/test/crab/${file1}/EcalValidation_${file1}.root ${plots_dir}
cp ${work_dir}/test/crab/${file2}/EcalValidation_${file2}.root ${plots_dir}
cd ${my_cmssw_base}/Validation/EcalValidation/data/macro


echo
echo 'To make plots, run in ROOT:'

echo '.x '${my_cmssw_base}/'Validation/EcalValidation/data/macro/DrawValidationPlots_digiless.cxx("'${plots_dir}'/EcalValidation_'${file1}'.root","'${plots_dir}'/EcalValidation_'${file2}'.root" ,"png", "'${plots_dir}'")'
echo

echo '.x '${my_cmssw_base}/'Validation/EcalValidation/data/macro/DrawValidationPlotsPi0.cxx("'${plots_dir}'/EcalValidation_'${file1}'.root","'${plots_dir}'/EcalValidation_'${file2}'.root" ,"png", "'${plots_dir}'")'
echo

# run root command in batch

root -b <<!

.x ${my_cmssw_base}/Validation/EcalValidation/data/macro/DrawValidationPlots_digiless.cxx("${plots_dir}/EcalValidation_${file1}.root","${plots_dir}/EcalValidation_${file2}.root","png","${plots_dir}")

.x ${my_cmssw_base}/Validation/EcalValidation/data/macro/DrawValidationPlotsPi0.cxx("${plots_dir}/EcalValidation_${file1}.root","${plots_dir}/EcalValidation_${file2}.root","png","${plots_dir}")
.q
!



echo 'Making webpage for '${file1}' vs '${file2}''

data_set_1=`\grep datasetpath ${work_dir}/test/crab/${file1}/crab.cfg | tr "=" " " | awk '{print $2}'`;
echo ${data_set_1}

data_set_2=`\grep datasetpath ${work_dir}/test/crab/${file2}/crab.cfg | tr "=" " " | awk '{print $2}'`;
echo ${data_set_2}

cat > ${plots_dir}/index.html <<EOF


<HTML>

<HEAD><TITLE> ECAL VALIDATION PLOTS ${file1} vs ${file2} </TITLE></HEAD>
 
<BODY link="color: rgb(0, 0, 253);">
<FONT color="Black">

<Center>
<h1> ECAL RECO Validation </h1>
</Center>

<hr>

<Center>
<h3>  <FONT color="Red"> ${file1}  <FONT color="Black"> vs <FONT color="Grey"> ${file2}  </h3>
</center>


<FONT color="Black">


<h4> Datasets </h4>
<ul>
 <li> <a  href="https://cmsweb.cern.ch/dbs_discovery/getData?dbsInst=cms_dbs_prod_global&proc=${data_set_1}"> <FONT color="Red"> ${data_set_1} </a>  </FONT> <BR>
 <li> <a href="https://cmsweb.cern.ch/dbs_discovery/getData?dbsInst=cms_dbs_prod_global&proc=${data_set_2}">  <FONT color="Grey"> ${data_set_2} </a>    </FONT> <BR>
</ul> 


<h4> Validation Plots </h4>
<ul>
 <li> Rec Hits
 <ul>
 <li><A href="#RecHitsMultiplicity"> Rec Hits Multiplicity </A><BR>
 <li><A href="#RechitSRP"> Rechit SRP </A><BR>
 <li><A href="#RecHitsEnergy"> Rec Hits Energy</A><BR>
 <li><A href="#RecHitsEnergyMax"> Rec Hits Max Energy </A><BR>
 <li><A href="#RecHitsEtaPhi"> Rec Hits Eta/Phi </A><BR>
 <li><A href="#RecHitsTime"> Rec Hits Time </A><BR>
 <li><A href="#RecHitsChi2"> Rec Hits Chi<sup>2</sup></A><BR>
 <li><A href="#RecHitsE4"> Rec Hits 1-E4/E1</A><BR>
 <li><A href="#RecHitsRecoFlag"> Rec Hits Reco Flag</A><BR>
 <li><A href="#RedRecHitsRecoFlag"> Reduced Rec Hits Reco Flag</A><BR>
 <li><A href="#ClusRecHitsRecoFlag"> Cluster Rec Hits Reco Flag</A><BR>
 <li><A href="#RecHitsOccupancyiPhiiEta_0"> RecHits Occupancy iPhi-iEta (${file2}) </A><BR>
 <li><A href="#RecHitsOccupancyiPhiiEta_1"> RecHits Occupancy iPhi-iEta (${file1}) </A><BR>
 <li><A href="#RecHitsOccupancy"> RecHits Occupancy (EB) </A><BR>
 <li><A href="#RecHitsSumEt"> RecHits SumEt </A><BR>
 
 
</ul>
 <li> Clusters
 <ul>
 <li><A href="#NumberOfBasicClusters"> Number of BasicClusters </A><BR>
 <li><A href="#NumberOfSuperClusters"> Number of SuperClusters </A><BR>
 <li><A href="#EnergySC"> Super Clusters Energy </A><BR> 
 <li><A href="#SuperClustersEtaPhi"> Super Clusters Eta/Phi </A><BR>
 <li><A href="#NumberOfCrystalsInSC"> Number of Crystals per SuperCluster </A><BR>
 <li><A href="#NumberOfBCInSC"> Number of Basic Clusters per SuperCluster </A><BR>
 <li><A href="#1-E4/E1"> Supercluster Seed 1-E4/E1 </A><BR>
 <li><A href="#ESclusters"> ES clusters  </A><BR>
</ul>
<li><A href="#Pi0peak"> Pi0 peak  </A><BR>
</ul>
<h4> Root Files </h4> 
<ul>
 <li><A href="#RootFile1"> Root File for ${file1} </A><BR>
 <li><A href="#RootFile2"> Root File for ${file2}</A><BR>
</ul>


<hr>

<h3><A name="RecHitsMultiplicity"> Rec Hits Multiplicity  </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_size.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EB_size.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_size.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EEP_size.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_size.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EEM_size.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_ES_size.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_ES_size.png"> </A>

<hr>

<h3><A name="RechitSRP"> Rechit SRP </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_SRP.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EB_SRP.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EE_SRP.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EE_SRP.png"> </A>
</A>

<hr>

<h3><A name="RecHitsEnergy"> Rec Hits Energy </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_energy.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EB_energy.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_energy.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EEP_energy.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_energy.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EEM_energy.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_ES_energy.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_ES_energy.png"> </A>

<hr>



<h3><A name="RecHitsEnergyMax"> Rec Hits Max Energy </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_energyMax.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EB_energyMax.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_energyMax.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EEP_energyMax.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_energyMax.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EEM_energyMax.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_ES_energyMax.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_ES_energyMax.png"> </A>

<hr>


<h3><A name="RecHitsEtaPhi"> Rec Hits Eta/Phi </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_eta.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_eta.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_phi.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_phi.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EE_phi.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EE_phi.png"> </A>

<hr>



<h3><A name="RecHitsTime"> Rec Hits Time </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_time.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EB_time.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_time.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EEP_time.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_time.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_EEM_time.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_ES_time.png> <img height="350" src="${httpdir}/${out_dir}/h_recHits_ES_time.png"> </A>

<hr>



<h3><A name="RecHitsChi2"> Rec Hits Chi<sup>2</sup> </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_Chi2.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_Chi2.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_Chi2.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEP_Chi2.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_Chi2.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEM_Chi2.png"> </A>

<hr>



<h3><A name="RecHitsE4"> Rec Hits 1-E4/E1 </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_E1oE4.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_E1oE4.png"> </A>

<hr>



<h3><A name="RecHitsRecoFlag"> Rec Hits Reco Flag </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_recoFlag.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_recoFlag.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EE_recoFlag.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EE_recoFlag.png"> </A>
</A>

<hr>



<h3><A name="RedRecHitsRecoFlag"> Reduced Rec Hits Reco Flag </h3>

<A HREF=${httpdir}/${out_dir}/h_redRecHits_EB_recoFlag.png> <img height="300" src="${httpdir}/${out_dir}/h_redRecHits_EB_recoFlag.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_redRecHits_EE_recoFlag.png> <img height="300" src="${httpdir}/${out_dir}/h_redRecHits_EE_recoFlag.png"> </A>
</A>

<hr>



<h3><A name="ClusRecHitsRecoFlag"> Cluster Rec Hits Reco Flag </h3>

<A HREF=${httpdir}/${out_dir}/h_basicClusters_recHits_EB_recoFlag.png> <img height="300" src="${httpdir}/${out_dir}/h_basicClusters_recHits_EB_recoFlag.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_basicClusters_recHits_EE_recoFlag.png> <img height="300" src="${httpdir}/${out_dir}/h_basicClusters_recHits_EE_recoFlag.png"> </A>
</A>

<hr>



<h3><A name="RecHitsOccupancyiPhiiEta_0"> RecHits Occupancy iPhi-iEta (${file2}) </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_occupancy_0.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_occupancy_0.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_occupancy_0.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEP_occupancy_0.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_occupancy_0.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEM_occupancy_0.png"> </A>

<hr>


<h3><A name="RecHitsOccupancyiPhiiEta_1"> RecHits Occupancy iPhi-iEta (${file1}) </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_occupancy_1.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_occupancy_1.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_occupancy_1.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEP_occupancy_1.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_occupancy_1.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEM_occupancy_1.png"> </A>

<hr>


<h3><A name="RecHitsOccupancy"> RecHits Occupancy (EB) </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_iPhiOccupancy.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_iPhiOccupancy.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_h_recHits_EB_iEtaOccupancy.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_iEtaOccupancy.png"> </A>

<hr>

<h3><A name="RecHitsSumEt"> RecHits SumEt </h3>

<A HREF=${httpdir}/${out_dir}/h_recHits_EB_sumEt.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EB_sumEt.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_sumEt.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEP_sumEt.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEP_sumEtCut.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEP_sumEtCut.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_sumEt.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEM_sumEt.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_recHits_EEM_sumEtCut.png> <img height="300" src="${httpdir}/${out_dir}/h_recHits_EEM_sumEtCut.png"> </A>

<hr>


<h3><A name="NumberOfBasicClusters"> Number of BasicClusters </h3>

<A HREF=${httpdir}/${out_dir}/h_basicClusters_EB_size.png> <img height="300" src="${httpdir}/${out_dir}/h_basicClusters_EB_size.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_basicClusters_EEP_size.png> <img height="300" src="${httpdir}/${out_dir}/h_basicClusters_EEP_size.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_basicClusters_EEM_size.png> <img height="300" src="${httpdir}/${out_dir}/h_basicClusters_EEM_size.png"> </A>

<hr>



<h3><A name="NumberOfSuperClusters"> Number of SuperClusters </h3>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EB_size.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EB_size.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEP_size.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEP_size.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEM_size.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEM_size.png"> </A>

<hr>



<h3><A name="EnergySC">  Super Clusters Energy </h3>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EB_energy.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EB_energy.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEP_energy.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEP_energy.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEM_energy.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEM_energy.png"> </A>

<hr>


<h3><A name="SuperClustersEtaPhi">  Super Clusters Eta/Phi </h3>

<A HREF=${httpdir}/${out_dir}/h_superClusters_eta.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_eta.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EB_phi.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EB_phi.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EE_phi.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EE_phi.png"> </A>

<hr>



<h3><A name="NumberOfCrystalsInSC">  Number of Crystals per SuperCluster </h3>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EB_nXtals.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EB_nXtals.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEP_nXtals.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEP_nXtals.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEM_nXtals.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEM_nXtals.png"> </A>

<hr>



<h3><A name="NumberOfBCInSC">  Number of Basic Clusters per SuperCluster </h3>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EB_nBC.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EB_nBC.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEP_nBC.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEP_nBC.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEM_nBC.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEM_nBC.png"> </A>

<hr>



<h3><A name="1-E4/E1"> 1-E4/E1 </h3>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EB_E1oE4.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EB_E1oE4.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEP_E1oE4.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEP_E1oE4.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_superClusters_EEM_E1oE4.png> <img height="300" src="${httpdir}/${out_dir}/h_superClusters_EEM_E1oE4.png"> </A>

<hr>



<h3><A name="ESclusters"> ES clusters  </h3>

<A HREF=${httpdir}/${out_dir}/h_esClusters_energy_plane1.png> <img height="300" src="${httpdir}/${out_dir}/h_esClusters_energy_plane1.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_esClusters_energy_plane2.png> <img height="300" src="${httpdir}/${out_dir}/h_esClusters_energy_plane2.png"> </A>

<A HREF=${httpdir}/${out_dir}/h_esClusters_energy_ratio.png> <img height="300" src="${httpdir}/${out_dir}/h_esClusters_energy_ratio.png"> </A>

<hr>



<h3><A name="Pi0peak"> Pi0 peak  </h3>

<A HREF=${httpdir}/${out_dir}/h_Pi0_EB_mass.png> <img height="300" src="${httpdir}/${out_dir}/h_Pi0_EB_mass.png"> </A>
<A HREF=${httpdir}/${out_dir}/h_Pi0_EE_mass.png> <img height="300" src="${httpdir}/${out_dir}/h_Pi0_EE_mass.png"> </A>

<hr>

<br>
<h3> <A name="RootFile1"> ROOT File ${file1} </h3> <A HREF=${httpdir}/${out_dir}/EcalValidation_${file1}.root> EcalValidation_${file1}.root </A>

<br>
<h3> <A name="RootFile2"> ROOT File ${file2} </h3>
<A HREF=${httpdir}/${out_dir}/EcalValidation_${file2}.root> EcalValidation_${file2}.root </A>

</FONT>
</BODY>
</HTML>

EOF

echo "==============================================="
echo "==========Template for validation mail========="
echo "===============================================

"

echo "Dear all,

ECAL DATA validation results are available here:

http://test-ecal-cosmics.web.cern.ch/test-ecal-cosmics/ValidationPlots/${whichdir}/${file1}_vs_${file2}/

* samples used:

new -> ${data_set_1}
ref -> ${data_set_2} 

* status: SUCCESS

* comment: ---

Best regards,

Leonardo and Martina" 






