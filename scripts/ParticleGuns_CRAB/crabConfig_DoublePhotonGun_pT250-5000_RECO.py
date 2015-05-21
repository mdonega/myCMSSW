from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DoublePhotonGun_pT250-5000_17_RECO'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step_RECO.py' # 'pset_tutorial_MC_generation.py' 

config.Data.inputDataset = '/DoublePhotonGun/mdonega-DoublePhoton_pT250-5000-562b533eb6b95100b90218971f8215e6/USER'                            
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5


config.Data.outLFNDirBase = '/store/user/mdonega/' # or '/store/group/<subdir>'
config.Data.publication = True
config.Data.publishDataName = 'DoublePhotonGun_pT250-5000_RECO'



config.Site.storageSite = 'T2_CH_CSCS'#  <site where the user has permission to write>

#config.Site.whitelist = ["T2_US_Nebraska"]
config.Site.whitelist = ["T2_CH_CSCS"]

