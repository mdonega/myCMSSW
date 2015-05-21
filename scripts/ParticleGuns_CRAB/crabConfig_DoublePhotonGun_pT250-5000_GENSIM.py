from CRABClient.UserUtilities import config 
config = config()

config.General.requestName = 'DoublePhotonGun_pT250-5000_19'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'step_GENSIM.py' # 'pset_tutorial_MC_generation.py' 

config.Data.primaryDataset = 'DoublePhotonGun'

config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 200
NJOBS = 5000 # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/mdonega/' # or '/store/group/<subdir>'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'DoublePhoton_pT250-5000'

config.Data.ignoreLocality = True

config.Site.storageSite = 'T2_CH_CSCS'#  <site where the user has permission to write>

#config.Site.whitelist = ["T2_US_Nebraska"]
config.Site.whitelist = ["T2_CH_CSCS"]
