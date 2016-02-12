from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName   = 'DoubleEG_16Dec2015_v'
config.General.transferLogs = False
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
#config.JobType.pluginName  = 'PrivateMC'
# Name of the CMSSW configuration file
config.JobType.psetName    = 'ecalvalidationDigiless_ES_cfg.py'
#config.JobType.inputFiles = ['gbrv3ele_52x.root', 'gbrv3ph_52x.root']
config.JobType.outputFiles = ['DoubleEG_16Dec2015_DATA_runD.root']

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
# For instance, this dataset will be named /CrabTestSingleMu/something/USER
config.Data.inputDataset = '/DoubleEG/Run2015D-ZElectron-16Dec2015-v2/RAW-RECO'
#config.Data.useParent = True
config.Data.inputDBS = 'global' #'phys03'
config.Data.splitting =  'FileBased'
config.Data.unitsPerJob = 5
config.Data.totalUnits = -1
config.Data.publication = False
# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
#config.Data.outLFN =  '/store/group/dpg_ecal/alca_ecalcalib/amartell/redigiHG/' # or '/store/group/<subdir>'   #'/store/group/dpg_ecal/alca_ecalcalib/amartell/'   #
config.Data.outLFNDirBase =  '/store/group/dpg_ecal/alca_ecalcalib/amartell/eoy2015ES_v/' # or '/store/group/<subdir>'   #'/store/group/dpg_ecal/alca_ecalcalib/amartell/'   #

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN' #'srm-eoscms.cern.ch'    #'T2_US_Nowhere'
#config.Site.whitelist = ['T2_*', 'T3_*']
