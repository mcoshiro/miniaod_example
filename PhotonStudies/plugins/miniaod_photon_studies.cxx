// Analysis script to make plots from MiniAOD photons to check issue observed 
// in 2016preVFP NanoAOD
//
// Original Author:  Michael Oshiro <oshiro@physics.ucsb.edu>
//         Created:  Thu, 15 Aug 2024 13:35:42 GMT

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TVector2.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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

//helper functions
float deltaR(float eta1, float phi1, float eta2, float phi2){
  return hypot(TVector2::Phi_mpi_pi(phi2-phi1), eta2-eta1);
}

//main class and auxiliary methods
class miniaod_photon_studies : public edm::one::EDAnalyzer<>  {
   public:
     explicit miniaod_photon_studies(const edm::ParameterSet&);
     ~miniaod_photon_studies();

     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
     virtual void beginJob() override;
     virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
     virtual void endJob() override;

     edm::EDGetTokenT<pat::MuonCollection> token_muons_;
     edm::EDGetTokenT<pat::PhotonCollection> token_photons_;
     edm::EDGetTokenT<edm::TriggerResults> token_triggers_;
     TH1D* hist_photon_eta;
     TH1D* hist_photon_r9_etagood;
     TH1D* hist_photon_r9_etabad;
     //TH1D* hist_photon_mvaid_etagood;
     //TH1D* hist_photon_mvaid_etabad;
};

miniaod_photon_studies::miniaod_photon_studies(const edm::ParameterSet& iConfig) :
  token_muons_(consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"))),
  token_photons_(consumes<pat::PhotonCollection>(edm::InputTag("slimmedPhotons"))),
  token_triggers_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))) {}

miniaod_photon_studies::~miniaod_photon_studies() {}

void miniaod_photon_studies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ process setup  ------------
void miniaod_photon_studies::beginJob() { 
  std::cout << "Starting job.\n"; 
  hist_photon_eta = new TH1D("hist_photon_eta",
      "N_{#mu}=2, N_{#gamma}=1, 80<m_{#mu#mu#gamma}<100, HLT_IsoMu24, p_{T}^{lead #mu}>25 GeV; #eta; Events/bin",
      30,-2.5,2.5);
  hist_photon_r9_etagood = new TH1D("hist_photon_r9_etagood",
      "N_{#mu}=2, N_{#gamma}=1, 80<m_{#mu#mu#gamma}<100, HLT_IsoMu24, p_{T}^{lead #mu}>25 GeV, |#eta|<1.5 or |#eta|>2.0; R_{9}; Events/bin",
      30,0.0,1.0);
  hist_photon_r9_etabad = new TH1D("hist_photon_r9_etabad",
      "N_{#mu}=2, N_{#gamma}=1, 80<m_{#mu#mu#gamma}<100, HLT_IsoMu24, p_{T}^{lead #mu}>25 GeV, 1.5<|#eta|<2.0; R_{9}; Events/bin",
      30,0.0,1.0);
  //hist_photon_mvaid_etagood = new TH1D("hist_photon_mvaid_etagood",
  //    "N_{#mu}=2, N_{#gamma}=1, 80<m_{#mu#mu#gamma}<100, HLT_IsoMu24, p_{T}^{lead #mu}>25 GeV, |#eta|<1.5 or |#eta|>2.0; Fall17v2 MVA ID Score; Events/bin",
  //    30,0.0,1.0);
  //hist_photon_mvaid_etagood = new TH1D("hist_photon_mvaid_etabad",
  //    "N_{#mu}=2, N_{#gamma}=1, 80<m_{#mu#mu#gamma}<100, HLT_IsoMu24, p_{T}^{lead #mu}>25 GeV, 1.5<|#eta|<2.0; Fall17v2 MVA ID Score; Events/bin",
  //    30,0.0,1.0);
  std::cout << "Starting job.\n"; 
}

// ------------ process cleanip  ------------
void miniaod_photon_studies::endJob() { 
  TFile out_file("temp_output.root","RECREATE");
  hist_photon_eta->Write();
  hist_photon_r9_etagood->Write();
  hist_photon_r9_etabad->Write();
  //hist_photon_mvaid_etagood->Write();
  //hist_photon_mvaid_etabad->Write();
  out_file.Close();
  //delete hist_photon_eta;
  //delete hist_photon_r9_etagood;
  //delete hist_photon_r9_etabad;
  //delete hist_photon_mvaid_etagood;
  //delete hist_photon_mvaid_etabad;
  std::cout << "Ending job.\n"; 
}

// ------------ event loop  ------------
void miniaod_photon_studies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(token_muons_, muons);
  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(token_photons_, photons);
  edm::Handle<edm::TriggerResults> triggers;
  iEvent.getByToken(token_triggers_, triggers);

  int nSignalMuon = 0;
  int nSignalPhoton = 0;
  float Lead_SignalMuon_pt = -999.0;
  float Sublead_SignalMuon_pt = -999.0;
  float Lead_SignalPhoton_pt = -999.0;
  float Lead_SignalPhoton_eta = -999.0;
  //float Lead_SignalPhoton_mvaID = -999.0;
  float Lead_SignalPhoton_r9 = -999.0;
  TLorentzVector p_mu1, p_mu2, p_ph;
  vector<float> SignalMuon_eta;
  vector<float> SignalMuon_phi;
  bool HLT_IsoMu24 = false; 

  //process triggers
  const edm::TriggerNames &names = iEvent.triggerNames(*triggers);
  for (unsigned itrig = 0; itrig < triggers->size(); itrig++) {
    //if (names.triggerName(itrig).find("HLT_IsoMu24") != string::npos) {
    //note the below check includes DZ version
    if (names.triggerName(itrig).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL") != string::npos) { 
      if (triggers->accept(itrig)) {
        HLT_IsoMu24 = true;
      }
    }
  }

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
          if (nSignalMuon == 0)
            p_mu1.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(),0.106);
          else if (nSignalMuon == 1)
            p_mu2.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(),0.106);
          //is signal muon
          nSignalMuon += 1;
          if (mu.pt() > Lead_SignalMuon_pt) {
            Lead_SignalMuon_pt = mu.pt();
            Sublead_SignalMuon_pt = Lead_SignalMuon_pt;
          }
          else if (mu.pt() > Sublead_SignalMuon_pt) {
            Sublead_SignalMuon_pt = mu.pt();
          }
          SignalMuon_eta.push_back(mu.eta());
          SignalMuon_phi.push_back(mu.phi());
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
        ph.passElectronVeto() && ph.photonID("mvaPhoID-RunIIFall17-v2-wp80")) {
      float mindr = 999.0;
      for (unsigned imu = 0; imu < SignalMuon_eta.size(); imu++) {
        float dr = deltaR(SignalMuon_eta[imu], SignalMuon_phi[imu],
                          ph.eta(), ph.phi());
        if (dr < mindr)
          mindr = dr;
      }
      if (mindr > 0.3) {
        //is signal photon
        p_ph.SetPtEtaPhiM(ph.et(),ph.eta(),ph.phi(),0.0);
        nSignalPhoton++;
        if (ph.et() > Lead_SignalPhoton_pt) {
          Lead_SignalPhoton_pt = ph.et();
          Lead_SignalPhoton_eta = ph.eta();
          Lead_SignalPhoton_r9 = ph.full5x5_r9();
        }
      }
    }
  }

  float ZCandidate_mass = (p_mu1+p_mu2+p_ph).M();

  //fill plots
  if (nSignalMuon == 2 && nSignalPhoton == 1 && ZCandidate_mass > 80.0 && ZCandidate_mass < 100.0 && HLT_IsoMu24 && Lead_SignalMuon_pt > 20.0 && Sublead_SignalMuon_pt > 10.0) {
    hist_photon_eta->Fill(Lead_SignalPhoton_eta);
    if (fabs(Lead_SignalPhoton_eta) > 1.5 && fabs(Lead_SignalPhoton_eta) < 2.0) {
      hist_photon_r9_etabad->Fill(Lead_SignalPhoton_r9);
    }
    else {
      hist_photon_r9_etagood->Fill(Lead_SignalPhoton_r9);
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniaod_photon_studies);
