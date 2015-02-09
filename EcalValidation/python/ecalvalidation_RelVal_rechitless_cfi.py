import FWCore.ParameterSet.Config as cms

ecalvalidation = cms.EDAnalyzer("EcalValidationRelVal_rechitless",
    #superClusterCollection_EB = cms.InputTag("correctedHybridSuperClusters"),
    #superClusterCollection_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
    #basicClusterCollection_EE = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters"),
    #recHitCollection_EE       = cms.InputTag("hltEcalRecHitAll","EcalRecHitsEE"),
    #recHitCollection_EB       = cms.InputTag("hltEcalRecHitAll","EcalRecHitsEB"),
    #redRecHitCollection_EE    = cms.InputTag("reducedEcalRecHitsEE"),
    #redRecHitCollection_EB    = cms.InputTag("reducedEcalRecHitsEB"),
    #basicClusterCollection_EB = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    #recHitCollection_ES       = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
    #ClusterCollectionX_ES     = cms.InputTag("multi5x5SuperClustersWithPreshower","preshowerXClusters"),
    #ClusterCollectionY_ES     = cms.InputTag("multi5x5SuperClustersWithPreshower","preshowerYClusters"),
    ebDigiCollection          = cms.InputTag("simEcalDigis","ebDigis"),
    eeDigiCollection          = cms.InputTag("simEcalDigis","eeDigis"),
      
    #tracks                    = cms.InputTag("generalTracks"),
    beamSpot                  = cms.InputTag("hltOnlineBeamSpot"),
    #jets                      = cms.InputTag("hltAntiKT5CaloJets"),

                                    
    #ethrEB = cms.double(0.8),
    #ethrEE = cms.double(1.2),
    ethrEB = cms.double(0.0),
    ethrEE = cms.double(0.0),
    gainId = cms.double(1.0),

    scEtThrEB = cms.double(0.0),
    scEtThrEE = cms.double(0.0),

    #runId = cms.double(194533.0),
    #runId = cms.double(200519.0),
    #runId = cms.double(206859.0),
    runId = cms.double(1.0),
    
    # for pi0

    isMonEBpi0 = cms.untracked.bool(True),
    selePtGamma = cms.double(1 ),
    selePtPi0 = cms.double( 2. ),
    #seleMinvMaxPi0 = cms.double( 0.22 ),
    seleMinvMaxPi0 = cms.double( 0.35 ),
    seleMinvMinPi0 = cms.double( 0.06 ),
    seleS4S9Gamma = cms.double( 0.83 ),
    selePi0Iso = cms.double( 0.5 ),
    ptMinForIsolation = cms.double( 1 ),
    selePi0BeltDR = cms.double( 0.2 ),
    selePi0BeltDeta = cms.double( 0.05 ),

    isMonEEpi0 = cms.untracked.bool(True),
    selePtGamma_EE = cms.double( 0.8 ),
    selePtPi0_EE = cms.double( 3.0 ),
    seleS4S9Gamma_EE = cms.double( 0.9 ),
    seleMinvMaxPi0_EE = cms.double( 0.3 ),
    seleMinvMinPi0_EE = cms.double( 0.05 ),
    ptMinForIsolation_EE = cms.double( 0.5 ),
    selePi0Iso_EE = cms.double( 0.5 ),
    selePi0BeltDR_EE  = cms.double( 0.2 ),
    selePi0BeltDeta_EE  = cms.double( 0.05 ),

    region1_Pi0_EE = cms.double(2),
    selePtGammaPi0_EE_region1 = cms.double(0.7),
    selePtPi0_EE_region1 = cms.double(3),

    region2_Pi0_EE = cms.double(2.5),
    selePtGammaPi0_EE_region2 = cms.double(0.5),
    selePtPi0_EE_region2 = cms.double(2),

    selePtGammaPi0_EE_region3 = cms.double(0.3),
    selePtPi0_EE_region3 = cms.double(1.2),

    clusSeedThr = cms.double( 0.5 ),
    clusSeedThr_EE = cms.double( 1.0 ),
    clusEtaSize = cms.int32( 3 ),
    clusPhiSize = cms.int32( 3 ),
    seleXtalMinEnergy = cms.double( -0.15 ),
    seleXtalMinEnergy_EE = cms.double( -0.75 ),
                                    
    isMaskEB = cms.untracked.bool(True),
    isMaskEE = cms.untracked.bool(False),

    maskEBFile = cms.untracked.string('maskEB.txt'),
    maskEEFile = cms.untracked.string('maskEE.txt'),

    useRecoFlag = cms.untracked.bool(False),
      
    posCalcParameters = cms.PSet( 
      LogWeighted = cms.bool( True ),
      T0_barl = cms.double( 7.4 ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    )
                                 
)

