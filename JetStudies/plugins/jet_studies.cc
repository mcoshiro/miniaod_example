//
// Package:    Analysis/JetStudies
// Class:      jet_studies
//
/**\class jet_studies jet_studies.cc JetStudies/jet_studies/plugins/jet_studies.cc

 Description: Analyzer module for jet studies with miniAOD

 Implementation:
*/
//
// Original Author:  Michael Oshiro <oshiro@physics.ucsb.edu>
//         Created:  Sat, 27 Jan 2024 12:56:42 GMT
//
//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class jet_studies : public edm::one::EDAnalyzer<>  {
   public:
      explicit jet_studies(const edm::ParameterSet&);
      ~jet_studies();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::PackedGenParticleCollection> token_genpart_;
      edm::EDGetTokenT<pat::JetCollection> token_jet_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> token_pfcand_;
};

// helper functions, etc.
const float PI = 3.141593;

float delta_phi(float phi1, float phi2){
  float dphi = fmod(fabs(phi2-phi1), 2.0*PI);
  return dphi>PI ? 2.0*PI-dphi : dphi;
}

float delta_r(float eta1, float eta2, float phi1, float phi2) {
  float deta = eta1-eta2;
  float dphi = delta_phi(phi1, phi2);
  return sqrt(deta*deta+dphi*dphi);
}

struct SimpleParticle {
  double pt;
  double eta;
  double phi;
  int pdgId;
};

bool sort_SimpleParticle(SimpleParticle a, SimpleParticle b) {
  return (a.pt>b.pt);
}

// constants, enums and typedefs
const float jet_pt_cut = 30.0;

// static data member definitions

// constructors and destructor
jet_studies::jet_studies(const edm::ParameterSet& iConfig) :
  //initialize tokens
  token_genpart_(consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"))),
  token_jet_(consumes<pat::JetCollection>(edm::InputTag("slimmedJets"))),
  token_pfcand_(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates")))
{

}


jet_studies::~jet_studies() {

}


// member functions

// ------------ method called for each event  ------------
void
jet_studies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "New event.\n";

  //get collections from tokens
  edm::Handle<pat::PackedGenParticleCollection> genparts;
  iEvent.getByToken(token_genpart_, genparts);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(token_jet_, jets);
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(token_pfcand_, pfcands);

  std::vector<SimpleParticle> sorted_genparts;
  for (const pat::PackedGenParticle &genpart : *genparts) {
    sorted_genparts.push_back({genpart.pt(), genpart.eta(), genpart.phi(), genpart.pdgId()});
  }
  std::sort(sorted_genparts.begin(), sorted_genparts.end(), sort_SimpleParticle);

  for (const pat::Jet &jet : *jets) {
    if (jet.pt() > jet_pt_cut) {
      //std::cout << "-----------------------------------------------------------\n"
      //          << "AK4 CHS Jet pt: " << jet.pt() << ", eta: " << jet.eta()
      //          << ", phi: " << jet.phi() << ", mass: " << jet.mass()
      //          << "\n";

      //std::cout << "Linked PFCands: \n";
      std::vector<SimpleParticle> sorted_pfcands;
      for (unsigned ipf = 0; ipf < jet.numberOfDaughters(); ipf++) {
        sorted_pfcands.push_back({jet.daughter(ipf)->pt(), jet.daughter(ipf)->eta(), 
                                  jet.daughter(ipf)->phi(), jet.daughter(ipf)->pdgId()});
      }
      std::sort(sorted_pfcands.begin(), sorted_pfcands.end(), sort_SimpleParticle);
      std::cout << "-----------------------------------------------------------\n";
      std::cout << "var jet = {pt:" << jet.pt() << ",eta:" << jet.eta() << ",phi:" 
                << jet.phi() << "};\n";
      std::cout << "var parts = [\n";
      bool first = true;
      for (SimpleParticle &pfcand : sorted_pfcands) {
        //std::cout << "PFCand pt: " << pfcand.pt << ", eta: " << pfcand.eta
        //          << ", phi: " << pfcand.phi << ", id: " << pfcand.pdgId
        //          << "\n";
        if (first) {
          first = false;
        }
        else {
        std::cout << ",\n";
        }
        std::cout << "             {pt:" << pfcand.pt << ",eta:" << pfcand.eta
                  << ",phi:" << pfcand.phi << ",id:" << pfcand.pdgId << "}";
      }
      std::cout << "];\n";

      //std::cout << "Stable GenParts within dR of 0.4: \n";
      std::cout << "-----------------------------------------------------------\n";
      std::cout << "var jet = {pt:" << jet.pt() << ",eta:" << jet.eta() << ",phi:" 
                << jet.phi() << "};\n";
      std::cout << "var parts = [\n";
      first = true;
      for (SimpleParticle &genpart : sorted_genparts) {
        if (delta_r(jet.eta(), genpart.eta, jet.phi(), genpart.phi)<0.4) {
          //std::cout << "GenPart pt: " << genpart.pt << ", eta: " << genpart.eta
          //          << ", phi: " << genpart.phi << ", id: " << genpart.pdgId
          //          << "\n";
          if (first) {
            first = false;
          }
          else {
            std::cout << ",\n";
          }
          std::cout << "             {pt:" << genpart.pt << ",eta:" << genpart.eta
                    << ",phi:" << genpart.phi << ",id:" << genpart.pdgId << "}";
        }
      }
      std::cout << "];\n";
    }
  }

//   for(const auto& track : iEvent.get(tracksToken_) ) {
//      // do something with track parameters, e.g, plot the charge.
//      // int charge = track.charge();
//   }
//
//#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//   ESHandle<SetupData> pSetup;
//   iSetup.get<SetupRecord>().get(pSetup);
//#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
jet_studies::beginJob()
{
  std::cout << "Starting job.\n";
}

// ------------ method called once each job just after ending the event loop  ------------
void
jet_studies::endJob()
{
  std::cout << "Ending job.\n";
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
jet_studies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(jet_studies);
