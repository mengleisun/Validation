from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName   = 'RunIISpring15DR74_SingleNeutrino_Ecut_ok2'
config.General.transferLogs = False
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
#config.JobType.pluginName  = 'PrivateMC'
# Name of the CMSSW configuration file
config.JobType.psetName    = 'ecalvalidationDigiless_cfg.py'
#config.JobType.inputFiles = ['gbrv3ele_52x.root', 'gbrv3ph_52x.root']
config.JobType.outputFiles = ['EcalValidation_MinBias_MCRUN2_74_V8.root']

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
# For instance, this dataset will be named /CrabTestSingleMu/something/USER
config.Data.inputDataset = '/SingleNeutrino/RunIISpring15DR74-NhcalZSFlat10to30bx50_74X_mcRun2_startup_dataMC_comparison_v1-v1/GEN-SIM-RECO'
#config.Data.useParent = True
config.Data.inputDBS = 'global' #'phys03'
config.Data.splitting =  'FileBased'
config.Data.unitsPerJob = 5
config.Data.totalUnits = -1
config.Data.publication = False
# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
#config.Data.outLFN =  '/store/group/dpg_ecal/alca_ecalcalib/amartell/redigiHG/' # or '/store/group/<subdir>'   #'/store/group/dpg_ecal/alca_ecalcalib/amartell/'   #
config.Data.outLFNDirBase =  '/store/user/amartell/DATA_MC_RunB2015/MC_MinBias_MCRUN2_74_V8_Ecut_ok2/' # or '/store/group/<subdir>'   #'/store/group/dpg_ecal/alca_ecalcalib/amartell/'   #

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN' #'srm-eoscms.cern.ch'    #'T2_US_Nowhere'
#config.Site.whitelist = ['T2_*', 'T3_*']
