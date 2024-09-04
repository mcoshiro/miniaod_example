import FWCore.ParameterSet.Config as cms

process = cms.Process("JetStudies")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

#one file of the QCD_HT100to200 2018 MiniAOD sample
process.source = cms.Source('PoolSource',
                             fileNames = cms.untracked.vstring(
                            'file:/net/cms26/cms26r0/oshiro/trees/002F46B9-4287-194A-BA9B-469CFB34D146.root'
                            )
                            )

process.jet_studies = cms.EDAnalyzer('jet_studies')
                              #tracks = cms.untracked.InputTag('generalTracks')
                              #)

process.p = cms.Path(process.jet_studies)
