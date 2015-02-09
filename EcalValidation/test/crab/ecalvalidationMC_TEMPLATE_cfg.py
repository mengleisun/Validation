import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

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
                                                     
process.GlobalTag.globaltag = "START62_V1::All"
                                                                       
process.GlobalTag.toGet = cms.VPSet(
            cms.PSet(record = cms.string("EcalPedestalsRcd"),
            tag = cms.string("EcalPedestals_206766_200_mc"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
            ),
            cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
            tag = cms.string("EcalLaserAPDPNRatios_20130130_447_p1_v2_run206859_mc"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
            ),
            cms.PSet(record = cms.string('EcalLaserAlphasRcd'),
            tag = cms.string('EcalLaserAlphas_EB_sic1_btcp152_EE_sic1_btcp116'),
            connect = cms.untracked.string('frontier://FrontierPrep/CMS_COND_ECAL')
            )
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_6_2_3/RelValDYJetsToLL/GEN-SIM-RECO/PU_START62_V1_rundepMC203002_dvmc-v2/00000/004136B2-9740-E311-849A-02163E008CCF.root'
    )
)



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#HLT selection
process.filter_1 = hlt.hltHighLevel.clone(
    HLTPaths = [ 'HLT_Ele*'],
    throw = False
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

process.load("Validation.EcalValidation.ecalvalidation_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EcalValidation.root')
)

process.p = cms.Path(
    #process.skimming*
    #process.hltLevel1GTSeed*
    #process.noscraping*
    #process.primaryVertexFilter*
    process.filter_1*
    process.ecalvalidation
    )

