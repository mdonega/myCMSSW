// -*- C++ -*-
//
// Package:    ESAnalysis
// Class:      ESAnalysis
// 
/**\class ESAnalysis ESAnalysis.cc Analysis/ESAnalysis/src/ESAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mauro Donega,32 4-C16,+41227679448,
//         Created:  Tue Jun 18 15:20:35 CEST 2013
// $Id$
//
//

// http://cmslxr.fnal.gov/lxr/source/RecoEcal/EgammaClusterProducers/src/PreshowerPhiClusterProducer.cc?v=CMSSW_6_1_0#100
// http://cmslxr.fnal.gov/lxr/source/RecoEcal/EgammaClusterAlgos/src/PreshowerClusterAlgo.cc?v=CMSSW_6_1_0#013

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include <DataFormats/EcalDetId/interface/ESDetId.h>
#include <DataFormats/EcalDetId/interface/EEDetId.h>
#include <DataFormats/EcalDetId/interface/EBDetId.h>
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include "TH1.h"
#include "TH2.h"
#include "TTree.h"


#include "RecoEcal/EgammaClusterProducers/interface/PreshowerAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include "TFile.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "RecoEcal/EgammaClusterProducers/interface/PreshowerClusterProducer.h"

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "Math/Cartesian3D.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterShape.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "Rtypes.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/EgammaReco/interface/EgammaTrigger.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"

#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

using namespace std;


//
// class declaration
//

class ESAnalysis : public edm::EDAnalyzer 
{
public:
  explicit ESAnalysis(const edm::ParameterSet&);
  ~ESAnalysis();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  //virtual void endRun(edm::Run const&, edm::EventSetup const&);
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//edm::InputTag ESrechitCollection_;
edm::InputTag SCProducer_;
edm::InputTag preshHitProducer_;  


typedef std::map<DetId, EcalRecHit> RecHitsMap;

int nEvt_;         // internal counter of events

std::string preshRecHitProducer_;
std::string ESrechitCollection_;

std::string outputFile_; // output file
TFile*  rootFile_;

// three stripes of 101 strips
Float_t sxp[101];
Float_t sx0[101];
Float_t sxm[101];
Float_t syp[101];
Float_t sy0[101];
Float_t sym[101];
float ES_eseffsirir = 0 ;
float ES_eseffsixix = 0 ;
float ES_eseffsiyiy = 0 ;
float SC_E   = 0 ;
float SC_Et  = 0 ;
float SC_eta = 0 ;
float SC_phi = 0 ;
float SC_x = 0 ;
float SC_y = 0 ;
float SC_z = 0 ;

TTree *tree;

//
// static data member definitions
//

//
// constructors and destructor
//
ESAnalysis::ESAnalysis(const edm::ParameterSet& iConfig)

{

  SCProducer_          = iConfig.getParameter<edm::InputTag>("SCProducer");
  
  preshRecHitProducer_     = iConfig.getParameter<std::string>("preshRecHitProducer");
  ESrechitCollection_      = iConfig.getParameter<std::string>("ESrechitCollection");
  

  outputFile_   = iConfig.getParameter<std::string>("outputFile");
  rootFile_     = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms


  tree  = new TTree("tree","preshower");

  tree->Branch("sxp",sxp,"sxp[101]/F");
  tree->Branch("sx0",sx0,"sx0[101]/F");
  tree->Branch("sxm",sxm,"sxm[101]/F");

  tree->Branch("syp",syp,"syp[101]/F");
  tree->Branch("sy0",sy0,"sy0[101]/F");
  tree->Branch("sym",sym,"sym[101]/F");
    
  tree->Branch("ES_eseffsirir",&ES_eseffsirir,"ES_eseffsirir/F");
  tree->Branch("ES_eseffsixix",&ES_eseffsixix,"ES_eseffsixix/F");
  tree->Branch("ES_eseffsiyiy",&ES_eseffsiyiy,"ES_eseffsiyiy/F");

  tree->Branch("SC_E"  ,&SC_E  ,"SC_E/F");
  tree->Branch("SC_Et" ,&SC_Et ,"SC_Et/F");
  tree->Branch("SC_eta",&SC_eta,"SC_eta/F");
  tree->Branch("SC_phi",&SC_phi,"SC_phi/F");

  tree->Branch("SC_x",&SC_x,"SC_x/F");
  tree->Branch("SC_y",&SC_y,"SC_y/F");
  tree->Branch("SC_z",&SC_z,"SC_z/F");
  
  nEvt_ = 0; 
    
}


ESAnalysis::~ESAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  cout << "Close output file" << endl;
  rootFile_->Close();
  delete rootFile_;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ESAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm; // needed for all fwk related classe


  //   TH1F*hx  = new TH1F("hx","hx",41,-20,20);
  //   TH1F*hy  = new TH1F("hy","hy",41,-20,20);

  std::cout << " .......  Event # " << nEvt_+1 << " is analyzing ....... " << std::endl;

  
  
  // SC collection
  Handle<reco::SuperClusterCollection> SCHandle;
  iEvent.getByLabel(SCProducer_,SCHandle);
  //  int SCcoll_size = SCHandle->size(); 
  //  std::cout << " Number of SuperClusters = " << SCcoll_size << std::endl;



  //ES collection
  edm::Handle< EcalRecHitCollection >   pRecHits;
  iEvent.getByLabel(preshRecHitProducer_, ESrechitCollection_, pRecHits);
  // pointer to the object in the product
  const EcalRecHitCollection* rechits = pRecHits.product(); // EcalRecHitCollection hit_collection = *rhcHandle;  
  //  cout << "#preshower RecHits = "<< rechits->size() << endl;
  


  // Get the lazy tools
  EcalClusterLazyTools lazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"),edm::InputTag("reducedEcalRecHitsES"));  




  // make the map of rechits:
  RecHitsMap rechits_map;
  //  std::map<DetId, EcalRecHit> rechits_map;
  EcalRecHitCollection::const_iterator it;
  for (it = rechits->begin(); it != rechits->end(); it++) {
    // remove bad ES rechits
    //FIXME checkFlag();
    if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
    //Make the map of DetID, EcalRecHit pairs
    rechits_map.insert(std::make_pair(it->id(), *it));   
  }



  // get the ECAL geometry:
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  


  const CaloSubdetectorGeometry *geometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  const CaloSubdetectorGeometry *& geometry_p = geometry;
  


  CaloSubdetectorTopology * topology_p=0;
  if (geometry)
    topology_p  = new EcalPreshowerTopology(geoHandle);
  


  // Loop over SC 
  reco::SuperClusterCollection::const_iterator SCIter ;
  for( SCIter=SCHandle->begin() ; SCIter!=SCHandle->end() ; SCIter++ ){
        
    // SC variables
    SC_eta  = SCIter->eta();
    SC_E    = SCIter->energy();
    SC_Et   = SC_E/cosh(SCIter->eta());
    SC_phi  = SCIter->phi();

    std::cout << "E = " << SC_E << " ; Et = " << SC_Et << " ; eta = " << SC_eta << " ; phi = " << SC_phi << std::endl;
    

    ES_eseffsirir = lazyTools.eseffsirir(*SCIter);
    ES_eseffsixix = lazyTools.eseffsixix(*SCIter);
    ES_eseffsiyiy = lazyTools.eseffsiyiy(*SCIter);
    //    cout << ES_eseffsirir << " " << ES_eseffsixix << " " << ES_eseffsiyiy  << endl;

    // SC position in (x,y,z)
    const GlobalPoint pointSC(SCIter->x(),SCIter->y(),SCIter->z()); // get the centroid of the SC
    SC_x = pointSC.x();
    SC_y = pointSC.y();
    SC_z = pointSC.z();
    //    std::cout << "SC centroid = " << pointSC << std::endl;

    //
    // Get the cell in ES closest to the EE cluster
    //
    // like in: http://cmslxr.fnal.gov/lxr/source/RecoEcal/EgammaCoreTools/src/EcalClusterLazyTools.cc#665
    //
    DetId tmp1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(pointSC, 1);  //point, plane
    DetId tmp2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(pointSC, 2);
    ESDetId strip1 = (tmp1 == DetId(0)) ? ESDetId(0) : ESDetId(tmp1);
    ESDetId strip2 = (tmp2 == DetId(0)) ? ESDetId(0) : ESDetId(tmp2);     
    
    cout<<"plane 1 === "<<strip1.zside()<<" "<<strip1.plane()<<" "<<strip1.six()<<" "<<strip1.siy()<<" "<<strip1.strip()<<endl;
    cout<<"plane 2 === "<<strip2.zside()<<" "<<strip2.plane()<<" "<<strip2.six()<<" "<<strip2.siy()<<" "<<strip2.strip()<<endl;
    


    //Make a navigator, and set it to the strip cell.
    EcalPreshowerNavigator navigator(strip1, topology_p);

    // central strip in X = strip1
    navigator.setHome(strip1);    
    RecHitsMap::iterator strip_it = rechits_map.find(strip1);   


    int halfWindow =50;

    //    std::cout << " Central Strip plane1 = " << strip_it->second.energy() << std::endl;
    //  hx->SetBinContent(21,strip_it->second.energy());
    if(strip_it !=rechits_map.end()) sx0[halfWindow]=strip_it->second.energy(); 
    else sx0[halfWindow]=0; 


    // CENTRAL STRIPE IN X

    int iNav      = 1;
    navigator.setHome(strip1);  
    while ( iNav < halfWindow+1 ){
      ESDetId strip_east = navigator.east();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_east);   
      //      std::cout << " E east i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hx->SetBinContent(halfWindow-iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sx0[halfWindow-iNav] = strip_i->second.energy();
      else sx0[halfWindow-iNav]=0; 
      ++iNav;
    }
    iNav      = 1;
    navigator.setHome(strip1);    
    while ( iNav < halfWindow+1 ){
      ESDetId strip_west = navigator.west();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_west);   
      //      std::cout << " E west i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hx->SetBinContent(halfWindow+1+iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sx0[halfWindow+iNav] = strip_i->second.energy();
      else sx0[halfWindow+iNav] =0;
      ++iNav;
    }    

    // TOP STRIPE in X

    navigator.setHome(strip1);    
    ESDetId strip_top = navigator.north();
    RecHitsMap::iterator strip_it_top = rechits_map.find(strip_top);   
    if(strip_it_top !=rechits_map.end()) sxp[halfWindow]=strip_it_top->second.energy();

    iNav      = 1;
    while ( iNav < halfWindow+1 ){
      ESDetId strip_east = navigator.east();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_east);   
      //      std::cout << " E east i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hx->SetBinContent(halfWindow-iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sxp[halfWindow-iNav] = strip_i->second.energy();
      else sxp[halfWindow-iNav] =0;
      ++iNav;
    }
    iNav      = 1;
    navigator.setHome(strip_top);    
    while ( iNav < halfWindow+1 ){
      ESDetId strip_west = navigator.west();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_west);   
      //      std::cout << " E west i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hx->SetBinContent(halfWindow+1+iNav,strip_i->second.energy());
    if(strip_i !=rechits_map.end()) sxp[halfWindow+iNav] = strip_i->second.energy();
    else sxp[halfWindow+iNav] =0;
      ++iNav;
    }    

    // BOTTOM STRIPE in X

    navigator.setHome(strip1);    
    ESDetId strip_bottom = navigator.south();
    RecHitsMap::iterator strip_it_bottom = rechits_map.find(strip_bottom);   
    if(strip_it_bottom !=rechits_map.end()) sxm[halfWindow]=strip_it_bottom->second.energy();
    else sxm[halfWindow] =0;

    iNav      = 1;
    while ( iNav < halfWindow+1 ){
      ESDetId strip_east = navigator.east();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_east);   
      //      std::cout << " E east i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hx->SetBinContent(halfWindow-iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sxm[halfWindow-iNav] = strip_i->second.energy();
      else sxm[halfWindow-iNav] =0;
      ++iNav;
    }
    iNav      = 1;
    navigator.setHome(strip_bottom);    
    while ( iNav < halfWindow+1 ){
      ESDetId strip_west = navigator.west();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_west);   
      //      std::cout << " E west i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hx->SetBinContent(halfWindow+1+iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sxm[halfWindow+iNav] = strip_i->second.energy();
      else sxm[halfWindow+iNav] =0;
      ++iNav;
    }    


    // check for plane 2


    //Make a navigator, and set it to the strip cell.
    navigator.setHome(strip2);    
    strip_it = rechits_map.find(strip2);       
    //    std::cout << " Central Strip plane2 = " << strip_it->second.energy() << std::endl;
    //    hy->SetBinContent(21,strip_it->second.energy());
    if(strip_it !=rechits_map.end()) sy0[halfWindow]=strip_it->second.energy();
    else sy0[halfWindow]=0;


    // CENTRAL STRIPE IN Y

    iNav      = 1;
    navigator.setHome(strip2);  
    while ( iNav < halfWindow+1 ){
      ESDetId strip_north = navigator.north();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_north);   
      //      std::cout << " E north i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hy->SetBinContent(halfWindow-iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sy0[halfWindow-iNav] = strip_i->second.energy();
      else sy0[halfWindow-iNav] =0;
      ++iNav;
    }
    iNav      = 1;
    navigator.setHome(strip2);    
    while ( iNav < halfWindow+1 ){
      ESDetId strip_south = navigator.south();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_south);   
      //      std::cout << " E south i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hy->SetBinContent(halfWindow+1+iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sy0[halfWindow+iNav] = strip_i->second.energy();
      else sy0[halfWindow+iNav] =0;
      ++iNav;
    }    

    // LEFT STRIPE in Y

    navigator.setHome(strip2);    
    ESDetId strip_left = navigator.east();
    RecHitsMap::iterator strip_it_left = rechits_map.find(strip_left);   
    if(strip_it_left !=rechits_map.end()) sym[halfWindow]=strip_it_left->second.energy();
    else sym[halfWindow] =0;
    //cout << halfWindow<<"/"<<sym[halfWindow] << " ";
    iNav      = 1;
    while ( iNav < halfWindow+1 ){
      ESDetId strip_north = navigator.north();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_north);   
      //      std::cout << " E north i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hy->SetBinContent(halfWindow-iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sym[halfWindow-iNav] = strip_i->second.energy();
      else sym[halfWindow-iNav] =0;
      //cout << halfWindow-iNav<<"/"<<sym[halfWindow-iNav] << " ";
      ++iNav;
    }
    iNav      = 1;
    navigator.setHome(strip_left);    
    while ( iNav < halfWindow+1 ){
      ESDetId strip_south = navigator.south();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_south);   
      //      std::cout << " E south i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hy->SetBinContent(halfWindow+1+iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) sym[halfWindow+iNav] = strip_i->second.energy();
      else sym[halfWindow+iNav] =0;
      //cout << halfWindow+iNav<<"/"<<sym[halfWindow+iNav] << " ";
      ++iNav;
    }    
    //    cout << endl;

    //     for (int i = 0; i<102; ++i) cout << i<<"/"<<sym[i] << " " ;
    //     cout << endl;
    
    // RIGHT STRIPE in Y

    navigator.setHome(strip2);    
    ESDetId strip_right = navigator.west();
    RecHitsMap::iterator strip_it_right = rechits_map.find(strip_right);   
    if(strip_it_right !=rechits_map.end()) syp[halfWindow]=strip_it_right->second.energy();
    else syp[halfWindow] =0;

    iNav      = 1;
    while ( iNav < halfWindow+1 ){
      ESDetId strip_north = navigator.north();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_north);   
      //      std::cout << " E north i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hy->SetBinContent(halfWindow-iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) syp[halfWindow-iNav] = strip_i->second.energy();
      else syp[halfWindow-iNav] =0;
      ++iNav;
    }
    iNav      = 1;
    navigator.setHome(strip_right);    
    while ( iNav < halfWindow+1 ){
      ESDetId strip_south = navigator.south();
      RecHitsMap::iterator strip_i = rechits_map.find(strip_south);   
      //      std::cout << " E south i = " << iNav << " " << strip_i->second.energy() << std::endl;    
      //      hy->SetBinContent(halfWindow+1+iNav,strip_i->second.energy());
      if(strip_i !=rechits_map.end()) syp[halfWindow+iNav] = strip_i->second.energy();
      else syp[halfWindow+iNav] =0;
      ++iNav;
    }    

    tree->Fill();
    //    cout << "CHECK" << endl;
  }
  
  //   hx->Write();
  //   hy->Write();
  
  ++nEvt_;
  

}


// ------------ method called once each job just before starting event loop  ------------
void 
ESAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ESAnalysis::endJob() 
{
  cout << "Write the tree" << endl;
  
  tree->Write();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
ESAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ESAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ESAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ESAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ESAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ESAnalysis);
