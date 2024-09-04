// Analysis script to skim MiniAOD files for analyzing 2016APV photon issue 
// in NanoAODs
//
// Requires 2 muons (H->Zg/H->ZZ signal criteria) and 1 loose photon (pT>20, 
// eta<1.566 or 1.4442<eta<2.5, pass CSEV)
//
// Original Author:  Michael Oshiro <oshiro@physics.ucsb.edu>
//         Created:  Sun, 25 Aug 2024


#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TSystem.h"
#include "TVector2.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

using std::cout;
using std::string;
using std::vector;

//main class and auxiliary methods
class miniaod_mumuph_filter : public edm::EDFilter  {
   public:
     explicit miniaod_mumuph_filter(const edm::ParameterSet&);
     ~miniaod_mumuph_filter();

     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
     virtual void beginJob();
     virtual bool filter(edm::Event&, const edm::EventSetup&);
     virtual void endJob();
     virtual bool beginRun(edm::Run&, edm::EventSetup const&);
     virtual bool endRun(edm::Run&, edm::EventSetup const&);
     virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
     virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

     edm::EDGetTokenT<pat::MuonCollection> token_muons_;
     edm::EDGetTokenT<pat::PhotonCollection> token_photons_;
};

miniaod_mumuph_filter::miniaod_mumuph_filter(const edm::ParameterSet& iConfig) :
  token_muons_(consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"))),
  token_photons_(consumes<pat::PhotonCollection>(edm::InputTag("slimmedPhotons"))) {}

miniaod_mumuph_filter::~miniaod_mumuph_filter() {}

void miniaod_mumuph_filter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void miniaod_mumuph_filter::beginJob() {}

void miniaod_mumuph_filter::endJob() {}

bool miniaod_mumuph_filter::beginRun(edm::Run&, edm::EventSetup const&) { return true; }

bool miniaod_mumuph_filter::endRun(edm::Run&, edm::EventSetup const&) { return true; }

bool miniaod_mumuph_filter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) { return true; }

bool miniaod_mumuph_filter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) { return true; }

// ------------ event loop  ------------
bool miniaod_mumuph_filter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(token_muons_, muons);
  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(token_photons_, photons);

  int nSignalMuon = 0;
  int nSignalPhoton = 0;

  //process muons
  for (const pat::Muon &mu : *muons) {
    //N.B. all muons with pT>5 are necessarily included in NanoAOD
    if (mu.pt() > 5.0 && fabs(mu.eta()) < 2.4 && 
      mu.passed(reco::Muon::CutBasedIdLoose)) {
      float mu_absdz = fabs(mu.dB(pat::Muon::PVDZ));
      float mu_absdxy = fabs(mu.dB(pat::Muon::PV2D));
      float mu_sip3d = fabs(mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D));
      if (mu_absdz < 1.0 && mu_absdxy < 0.5 && mu_sip3d < 4.0) {
        float mu_pfRelIso03_all = (mu.pfIsolationR03().sumChargedHadronPt + 
            fmax(mu.pfIsolationR03().sumNeutralHadronEt + 
            mu.pfIsolationR03().sumPhotonEt - 
            mu.pfIsolationR03().sumPUPt/2,0.0))/(mu.pt());
        if (mu_pfRelIso03_all < 0.35) {
          nSignalMuon += 1;
        }
      }
    }
  }

  //process photons
  for (const pat::Photon &ph : *photons) {
    float ph_origin_abseta = fabs(ph.superCluster()->eta());
    //N.B. all photons with pT>20 are necessarily included in NanoAOD
    if (ph.et() > 20.0 && ((ph_origin_abseta < 1.4442) || 
        (ph_origin_abseta > 1.566 && ph_origin_abseta < 2.5)) &&
        ph.passElectronVeto()) {
      nSignalPhoton++;
    }
  }

  return (nSignalMuon>=2 && nSignalPhoton>=1);
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniaod_mumuph_filter);
