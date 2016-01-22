// -*- C++ -*-
//
// Package:    hepMCAnalyzer/hepMCAnalyzer
// Class:      hepMCAnalyzer
// 
/**\class hepMCAnalyzer hepMCAnalyzer.cc hepMCAnalyzer/hepMCAnalyzer/plugins/hepMCAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mauro Donega
//         Created:  Thu, 21 Jan 2016 17:31:58 GMT
//
//


// system include files
#include <memory>
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoParticleFlow/Configuration/plugins/HepMCCopy.h"
#include "HepMC/GenEvent.h"

// ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TLorentzVector.h"

// miscellaneous  
#include <fstream>

using namespace std;

//
// class declaration
//

class hepMCAnalyzer : public edm::EDAnalyzer {
   public:
      explicit hepMCAnalyzer(const edm::ParameterSet&);
      ~hepMCAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  TFile   *f = 0;
  TNtuple *nt = 0;
  int     npho  = 0;
  int     npi0  = 0;
  int     nch   = 0;
  int     nother  = 0;
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
hepMCAnalyzer::hepMCAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  f = new TFile("nt.root","RECREATE");
  nt = new TNtuple("nt","data hepMC","pdgId:pt:eta:phi");

}


hepMCAnalyzer::~hepMCAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  f->Close();
  std::cout << "Npho = " << npho << " ; Npi0 = " << npi0 << " ; Nch = " << nch << " ; N other (decaying) = " << nother << std::endl;
}


//
// member functions
//

// ------------ method called for each event  ------------
void hepMCAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& iSetup)
{
  using namespace edm;
  edm::Handle<edm::HepMCProduct > genEvtHandle;
  e.getByLabel( "generator", genEvtHandle) ;
  const HepMC::GenEvent* Evt = genEvtHandle->GetEvent() ;
  
  
  //output ntuple
  //
  // this is an example loop over the hierarchy of vertices
  //
  for ( HepMC::GenEvent::vertex_const_iterator
	  itVtx=Evt->vertices_begin(); itVtx!=Evt->vertices_end(); ++itVtx )
    {
      //
      // this is an example loop over particles coming out of each vertex in the loop
      //
      for ( HepMC::GenVertex::particles_out_const_iterator itPartOut=(*itVtx)->particles_out_const_begin();
	    itPartOut!=(*itVtx)->particles_out_const_end(); ++itPartOut )
	{	  
	  
	  TLorentzVector v;
	  v.SetPxPyPzE((*itPartOut)->momentum().px(),(*itPartOut)->momentum().py(),(*itPartOut)->momentum().pz(),(*itPartOut)->momentum().e());
	  
	  if (v.Pt() > 0.0000001 ) { // ALWAYS remove stuff with no pT
	    //	    nt->Fill((*itPartOut)->pdg_id(), v.Pt(), v.Eta(), v.Phi());	   
	    if ((*itPartOut)->pdg_id() == 22  && v.Pt() > 0.1) npho++;
	    else if (    (*itPartOut)->pdg_id()  == 111  && v.Pt() > 0.2      ) npi0++;
	    else if (abs((*itPartOut)->pdg_id()) == 211  ) nch++; // pi+-
	    else if (abs((*itPartOut)->pdg_id()) == 2212 ) nch++; // p+-
	    else if (abs((*itPartOut)->pdg_id()) == 321  ) nch++; // k+-
	    else if (abs((*itPartOut)->pdg_id())  > 8    &&     //remove quarks/guons/W/Z and the previous particles
		     abs((*itPartOut)->pdg_id()) != 21   && 
		     abs((*itPartOut)->pdg_id()) != 22   && 
		     abs((*itPartOut)->pdg_id()) != 23   && 
		     abs((*itPartOut)->pdg_id()) != 24   &&
		     (*itPartOut)->pdg_id()      != 111  &&
		     abs((*itPartOut)->pdg_id()) != 211  && 
		     abs((*itPartOut)->pdg_id()) != 2212 &&
		     abs((*itPartOut)->pdg_id()) != 321)  {
	      nother++;
	      //	      std::cout << " pdgid/pt/eta/phi: " << (*itPartOut)->pdg_id() << "\t" << v.Pt() <<  "\t" << v.Eta() << "\t" << v.Phi() <<  std::endl;
	    }
	  }
	}
    }
  nt->Write();
  f ->Write();
  
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
hepMCAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
hepMCAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
hepMCAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
hepMCAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
hepMCAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
hepMCAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
hepMCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(hepMCAnalyzer);
