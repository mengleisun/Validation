import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

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
#process.GlobalTag.globaltag = "FT_R_70_V1::All"
#process.GlobalTag.globaltag = "FT_R_53_V18::All" ## runD minBias 22Jan rereco
#process.GlobalTag.globaltag = "FT_R_70_V1::All" ## runD release compatibility tag

process.GlobalTag.globaltag = "START70_V7A::All" ## MC for 8TeV rereco

#readFiles = cms.untracked.vstring()
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
'/store/group/dpg_ecal/alca_ecalcalib/amartell/MinB_ESsat_MC/runD_ok/MinBias_TuneZ2star_8TeV-pythia6/crab_MinBias_MC_ESsat/150320_140818/0000/MinBias_MC_runD_8TeV_238.root'
#'/store/group/dpg_ecal/alca_ecalcalib/amartell/MinB_ESsat_MC/runD_ok/MinBias_TuneZ2star_8TeV-pythia6/crab_MinBias_MC_ESsat/150320_140818/0000/MinBias_MC_runD_8TeV_992.root'
    )
    #fileNames = readFiles
                            )
#readFiles.extend( [
#'/store/group/dpg_ecal/alca_ecalcalib/amartell/MinB_runD_ESsat/MinBias_runD_8TeV_1.root',
# ] );




process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# filter on PhysDeclared bit
process.skimming = cms.EDFilter("PhysDecl",
    applyfilter = cms.untracked.bool(True)
)


process.load("Validation.EcalValidation.ecalvalidation_ES_cfi")

process.TFileService = cms.Service("TFileService",
         fileName = cms.string('MCRERECO_ES_MinBias.root')
)

process.p = cms.Path(
    process.ecalvalidationES
    )

