import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

# Geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

# initialize magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")

# GT
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = "START71_V8::All"
#process.GlobalTag.globaltag = "MCRUN2_71_V1::All"
process.GlobalTag.globaltag = "FT_R_70_V1::All"


process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/work/a/amartell/Simulation/CMSSW_7_1_1/src/genPS/step4_RAW2DIGI_L1Reco_RECO.root'
#'file:/afs/cern.ch/work/a/amartell/Simulation/CMSSW_7_1_1/src/SimCalorimetry/EcalSimProducers/test/python/step3_RAW2DIGI_L1Reco_RECO.root'
## minBias runD
#'file:/afs/cern.ch/work/a/amartell/ES_Studies_GainSwitch/CMSSW_7_0_7/src/redigiDATA/step4DATA_RAW2DIGI_L1Reco_RECO.root'
'file:/afs/cern.ch/work/a/amartell/ES_Studies_GainSwitch/CMSSW_7_0_7/src/redigiDATA/step4DATA_RAW2DIGI_L1Reco_RECO_Zee.root' 
   )
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# filter on PhysDeclared bit
process.skimming = cms.EDFilter("PhysDecl",
    applyfilter = cms.untracked.bool(True)
)

# filter on bit 40 || 41 nad !(bit36 || bit37 || bit38 || bit39)
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

# FilterOutScraping
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

# Good Vertex Filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNumberOfTracks = cms.uint32(3) ,
                                           maxAbsZ = cms.double(15),	
                                           maxd0 = cms.double(2)	
                                           )


process.load("Validation.EcalValidation.ecalvalidation_Samples_cfi")

process.TFileService = cms.Service("TFileService",
         fileName = cms.string('dataRERECO_ES_Zee.root')
)

process.p = cms.Path(
     process.ecalvalidationSamples
    )

