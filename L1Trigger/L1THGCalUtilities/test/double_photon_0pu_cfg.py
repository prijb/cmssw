from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'double_photon_0pu'
config.General.workArea = 'crab_area'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'testHGCalL1T_RelValV11_cfg.py'
config.JobType.maxMemoryMB = 2500

config.Data.inputDataset = '/DoublePhoton_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-NoPU_111X_mcRun4_realistic_T15_v1-v1/FEVT'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 70
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/HGCalTPG_LLP_Test_v1/'
config.Data.publication = True
config.Data.outputDatasetTag = 'double_photon_0pu'

config.Site.storageSite = 'T2_UK_London_IC'