import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.MessageLogger.cerr = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))

# Geometry
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

# initialize magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")

##noise    
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')  
                                                     
process.GlobalTag.globaltag = "START53_V14::All"
                                                                       
process.GlobalTag.toGet = cms.VPSet(
            cms.PSet(record = cms.string("EcalPedestalsRcd"),
            tag = cms.string("EcalPedestals_206766_200_mc"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
            ),
            cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
            tag = cms.string("EcalLaserAPDPNRatios_TL500_IL1E34_mc"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
            ),
            cms.PSet(record = cms.string('EcalLaserAlphasRcd'),
            tag = cms.string('EcalLaserAlphas_EB_sic1_btcp1_EE_sic1_btcp1'),
            connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
            )
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
    "root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_5_3_10_patch1/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG/START53_V7L_PhRD2_rundepMC_11May2013-v1/00000/FC060874-25BD-E211-9D6A-003048F1C58E.root",
    "root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_5_3_10_patch1/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG/START53_V7L_PhRD2_rundepMC_11May2013-v1/00000/DE9B3BEC-25BD-E211-87E0-003048F23FE8.root",
    "root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_5_3_10_patch1/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG/START53_V7L_PhRD2_rundepMC_11May2013-v1/00000/94DBD386-25BD-E211-9011-0025901D5DFA.root",
    "root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_5_3_10_patch1/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG/START53_V7L_PhRD2_rundepMC_11May2013-v1/00000/6AB76390-25BD-E211-ADEF-003048CFAD08.root",
    "root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_5_3_10_patch1/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG/START53_V7L_PhRD2_rundepMC_11May2013-v1/00000/54DB2494-25BD-E211-9173-003048FEC278.root"
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


#Good Vertex Filter (see GOODCOLL skim)
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 3 && abs(z) <= 24 && position.Rho <= 2"), 
  filter = cms.bool(True)   
)

# FilterOutScraping
process.noscraping = cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.25)
)

process.load("Validation.EcalValidation.ecalvalidation_RelVal_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EcalValidation_RelValZEE_CMSSW_5_3_10_patch1-START53_V7L_PhRD2_rundepMC_11May2013-v1_GEN-SIM-DIGI-RAW-HLTDEBUG_runD.root')
)

process.p = cms.Path(
    #process.skimming*
    #process.hltLevel1GTSeed*
    #process.noscraping*
    #process.primaryVertexFilter*
    process.ecalvalidation
    )

