import FWCore.ParameterSet.Config as cms

process = cms.Process("PhotonStudies")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source('PoolSource',
                             fileNames = cms.untracked.vstring(
                            #'file:/data2/oshiro/ntuples/DoubleMuon_sampleMiniAOD.root'
                            'file:/homes/oshiro/DoubleMuon_MiniAOD_mmpskim.root'
                            #'root://cmsxrootd.fnal.gov//store/data/Run2016D/DoubleMuon/MINIAOD/HIPM_UL2016_MiniAODv2-v1/130000/352D4168-8884-CF48-8990-B29F3A8CF4C2.root'
                            #'root://cmsxrootd.fnal.gov//store/data/Run2016F/SingleMuon/MINIAOD/HIPM_UL2016_MiniAODv2-v2/140000/02A57141-01B6-4141-992A-A48786928471.root'
                            )
                            )

#run EGM post-reco tools to calculate photon MVA ID
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False, #no point in re-running them, they are already fine
                       runVID=True,
                       era='2016preVFP-UL',
                       eleIDModules=[],
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff']) 

process.miniaod_photon_studies = cms.EDAnalyzer('miniaod_photon_studies')

process.p = cms.Path(process.egammaPostRecoSeq * process.miniaod_photon_studies)
process.MessageLogger.cerr.threshold = 'ERROR'
