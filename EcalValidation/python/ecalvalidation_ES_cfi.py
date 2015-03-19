import FWCore.ParameterSet.Config as cms

ecalvalidationES = cms.EDAnalyzer("EcalValidation_ES",
    recHitCollection_EB       = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    recHitCollection_EE       = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    recHitCollection_ES       = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
    redRecHitCollection_ES    = cms.InputTag("reducedEcalRecHitsES"),

    ClusterCollectionX_ES     = cms.InputTag("multi5x5SuperClustersWithPreshower","preshowerXClusters"),
    ClusterCollectionY_ES     = cms.InputTag("multi5x5SuperClustersWithPreshower","preshowerYClusters"),
    #esDigiCollection          = cms.InputTag("simEcalUnsuppressedDigis",""),

    EleTag                    = cms.InputTag("gsfElectrons"),
    #EleTag                    = cms.InputTag("gedGsfElectrons"),
    PhoTag                    = cms.InputTag("photons"),
    saveDigis                 = cms.untracked.bool(False),
    SaveSrFlag                = cms.untracked.bool(True),
    tracks                    = cms.InputTag("generalTracks"),
    beamSpot                  = cms.InputTag("offlineBeamSpot"),
    jets                      = cms.InputTag("ak5CaloJets"),

    ethrEB = cms.double(0.0),
    ethrEE = cms.double(0.0),
    gainId = cms.double(1.0),

    scEtThrEB = cms.double(0.0),
    scEtThrEE = cms.double(0.0)    
    )

