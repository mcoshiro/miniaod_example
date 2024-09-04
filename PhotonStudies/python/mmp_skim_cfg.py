import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("PhotonSkim")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'ERROR'

files = FileUtils.loadListFromFile("data/datasetlist.txt")
process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(*files))

#process.source = cms.Source('PoolSource',
#                             fileNames = cms.untracked.vstring(
#                            #'file:/data2/oshiro/ntuples/DoubleMuon_sampleMiniAOD.root'
#                            'root://cmsxrootd.fnal.gov//store/data/Run2016D/DoubleMuon/MINIAOD/HIPM_UL2016_MiniAODv2-v1/130000/352D4168-8884-CF48-8990-B29F3A8CF4C2.root'
#                            #'root://cmsxrootd.fnal.gov//store/data/Run2016F/SingleMuon/MINIAOD/HIPM_UL2016_MiniAODv2-v2/140000/02A57141-01B6-4141-992A-A48786928471.root'
#                            )
#                            )

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("file:/homes/oshiro/test2.root")
#                               )

process.miniaod_mumuph_filter = cms.EDFilter('miniaod_mumuph_filter')

process.p = cms.Path(process.miniaod_mumuph_filter)

process.outp1=cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('file:/homes/oshiro/DoubleMuon_MiniAOD_mmpskim.root'),
        SelectEvents = cms.untracked.PSet(
                SelectEvents = cms.vstring('p')
                )
        )

process.ep = cms.EndPath(process.outp1)
