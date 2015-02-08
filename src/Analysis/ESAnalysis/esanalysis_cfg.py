import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START61_V11::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

#process.load('Configuration/StandardSequences/Services_cff')
#process.load('Configuration/StandardSequences/MagneticField_38T_cff')
#process.load('Configuration/EventContent/EventContent_cff')

#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:root://eoscms//eos/cms/store/user/donega/ParticleGun610/Pi0/RECO/SinglePi0_cfi_py_RAW2DIGI_RECO_0.root',
#        'file:/afs/cern.ch/user/d/donega/work/CMSSW_6_1_0/test/trackerZeroMaterial/SingleGammaPt5-50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO.root'
        'file:/afs/cern.ch/user/d/donega/work/CMSSW_6_1_0/test/trackerZeroMaterial/SinglePi0Pt5-50_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO.root'
    )
)


process.demo = cms.EDAnalyzer('ESAnalysis',
                              verbose                 = cms.untracked.bool(False),
                              SCProducer              = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
                              ESrechitCollection      = cms.string('EcalRecHitsES'),
                              preshRecHitProducer     = cms.string('ecalPreshowerRecHit'),
                              #preshClusterProducer    = cms.string('multi5x5PreshowerClusterShape'),
                              outputFile              = cms.string('preshowerAnalyzer.root')
                              )


process.p = cms.Path(process.demo)

