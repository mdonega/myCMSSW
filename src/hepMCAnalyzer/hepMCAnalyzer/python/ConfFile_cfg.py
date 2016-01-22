import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.MessageLogger = cms.Service("MessageLogger",
#                                    destinations   = cms.untracked.vstring('detailedInfo'),
#                                    categories      = cms.untracked.vstring('eventNumber'),
#                                    detailedInfo    = cms.untracked.PSet(eventNumber = cms.untracked.PSet(reportEvery = cms.untracked.int32(1)))
#                                    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/m/mdonega/work/CMSSW_7_4_7/src/eops/cms/store/mc/RunIIWinter15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v1/00000/B24DAAFB-90B0-E411-9FAA-0025905AF57E.root'
    )
)

process.demo = cms.EDAnalyzer('hepMCAnalyzer'
)


process.p = cms.Path(process.demo)
