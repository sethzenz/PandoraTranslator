#include "runPandora.h"
#include "HGCal/PandoraTranslator/interface/CMSBFieldPlugin.h"
#include "HGCal/PandoraTranslator/interface/CMSPseudoLayerPlugin.h"
#include "LCContent.h"
#include "HGCal/PandoraTranslator/interface/CMSTemplateAlgorithm.h"
#include "PandoraMonitoringApi.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"

#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

// Addition for HGC geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/FlatTrd.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//We need the speed of light
#include "CLHEP/Units/PhysicalConstants.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/RootAutoLibraryLoader/interface/RootAutoLibraryLoader.h"

#include "TGClient.h"
#include "TVirtualX.h"
#include "TROOT.h"
#include "TRint.h"
 
#include <fstream>
#include <iostream>
#include <cstdlib>

#include <algorithm>
#include <utility>
#include <iterator> 
#include <map> 

using namespace edm;
using namespace reco;
using namespace pandora;  
using namespace lc_content;
using namespace cms_content;

enum { D_HGCEE, D_HGCHEF, D_HGCHEB };
enum { EM, HAD };

namespace cms_content {
  pandora::StatusCode RegisterBasicPlugins(const pandora::Pandora &pandora)
  {
    LC_ENERGY_CORRECTION_LIST(PANDORA_REGISTER_ENERGY_CORRECTION);
    LC_PARTICLE_ID_LIST(PANDORA_REGISTER_PARTICLE_ID);

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(pandora, new cms_content::CMSPseudoLayerPlugin));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetShowerProfilePlugin(pandora, new lc_content::LCShowerProfilePlugin));

    return pandora::STATUS_CODE_SUCCESS;
  }
}

pandora::Pandora * runPandora::m_pPandora = NULL;
//
// constructors and destructor
//
runPandora::runPandora(const edm::ParameterSet& iConfig) 
{
  //now do what ever initialization is needed
  m_pPandora = new pandora::Pandora();

  //stm = new steerManager("energyWeight.txt");

  //mFileNames  = iConfig.getParameter<std::vector<std::string> > ("filenames"); 
  inputTagEcalRecHitsEB_ = iConfig.getParameter<InputTag>("ecalRecHitsEB");
  inputTagHcalRecHitsHBHE_ = iConfig.getParameter<InputTag>("hcalRecHitsHBHE");
  inputTagHGCEErechit_ = iConfig.getParameter<InputTag>("HGCEErechitCollection");
  inputTagHGCHEFrechit_ = iConfig.getParameter<InputTag>("HGCHEFrechitCollection");
  inputTagHGCHEBrechit_ = iConfig.getParameter<InputTag>("HGCHEBrechitCollection");
  inputTagGeneralTracks_ = iConfig.getParameter< std::vector < InputTag > >("generaltracks");
  //  inputTagtPRecoTrackAsssociation_ = iConfig.getParameter<InputTag>("tPRecoTrackAsssociation");
  //  iEvent.getByLabel(inputTagtPRecoTrackAsssociation_, rectosimCollection);
  //  const reco::RecoToSimCollection pRecoToSim = *(rectosimCollection.product());

  inputTagGenParticles_ = iConfig.getParameter<InputTag>("genParticles");
  m_pandoraSettingsXmlFile = iConfig.getParameter<edm::FileInPath>("inputconfigfile");

  m_calibrationParameterFile = iConfig.getParameter<edm::FileInPath>("calibrParFile");
  m_energyCorrMethod = iConfig.getParameter<std::string>("energyCorrMethod");
  m_energyWeightingFilename  = iConfig.getParameter<edm::FileInPath>("energyWeightFile");
  _outputFileName = iConfig.getParameter<std::string>("outputFile");

  stm = new steerManager(m_energyWeightingFilename.fullPath().c_str());
  
// NS // SHOWER PROFILE CALCULATOR
  
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterAlgorithms(*m_pPandora));

  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, cms_content::RegisterBasicPlugins(*m_pPandora));

  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetBFieldPlugin(*m_pPandora, new CMSBFieldPlugin()));    
  
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*m_pPandora, "Template", new CMSTemplateAlgorithm::Factory));

  // prepareGeometry(iSetup);
}

runPandora::~runPandora()
{
 

   delete stm;

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void runPandora::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   resetVariables();

  std::cout << "Analyzing events" << std::endl ; 
  if ( firstEvent_ ) { 
    firstEvent_ = false ; 
    std::cout << "At the first event...preparing geometry" << std::endl ; 
    prepareGeometry(iSetup) ; 
    std::cout << "Done with Geometry setup...moving along" << std::endl ; 
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_pandoraSettingsXmlFile.fullPath()));

  }

  // std::cout << "Analyzing events 1 " << std::endl ;

  // Get the primary vertex
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
  reco::Vertex pv = pvHandle->at(0) ; 
  
  // std::cout << "Analyzing events 2 " << std::endl ;
  
  // get the Calorimeter rechit collections
  edm::Handle<EcalRecHitCollection> ecalRecHitHandleEB;
  edm::Handle<HBHERecHitCollection> hcalRecHitHandleHBHE;
  edm::Handle<HGCRecHitCollection> HGCeeRecHitHandle;
  edm::Handle<HGCRecHitCollection> HGChefRecHitHandle;
  edm::Handle<HGCRecHitCollection> HGChebRecHitHandle;
  
  // std::cout << iEvent.getByLabel(inputTagHGCEErechit_, HGCeeRecHitHandle) << " " 
  //         << iEvent.getByLabel(inputTagHGCHEFrechit_, HGChefRecHitHandle) << " " 
  //         << iEvent.getByLabel(inputTagHGCHEBrechit_, HGChebRecHitHandle) << std::endl ; 

  bool found = iEvent.getByLabel(inputTagEcalRecHitsEB_, ecalRecHitHandleEB) && 
    iEvent.getByLabel(inputTagHcalRecHitsHBHE_, hcalRecHitHandleHBHE) && 
    iEvent.getByLabel(inputTagHGCEErechit_, HGCeeRecHitHandle) && 
    iEvent.getByLabel(inputTagHGCHEFrechit_, HGChefRecHitHandle) && 
    iEvent.getByLabel(inputTagHGCHEBrechit_, HGChebRecHitHandle);

  edm::Handle<reco::RecoToSimCollection > rectosimCollection;
  iEvent.getByLabel(inputTagtPRecoTrackAsssociation_, rectosimCollection);
  const reco::RecoToSimCollection pRecoToSim = *(rectosimCollection.product());
    
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  //Why PF uses global point (0,0,0) for all events?
  math::XYZVector B_(math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0))));

  edm::Handle<std::vector<reco::GenParticle> > genpart;
  iEvent.getByLabel(inputTagGenParticles_,genpart);
  
  if(!found ) {
    std::ostringstream err;
    err<<"cannot find rechits: "<< HGCeeRecHitHandle.isValid() << "," << HGChefRecHitHandle.isValid() << "," << HGChebRecHitHandle.isValid() ;
    LogError("runPandora")<<err.str()<<std::endl;
    throw cms::Exception( "MissingProduct", err.str());
  } 

  //  prepareTrack(B_,pRecoToSim,iEvent,iSetup);
  prepareTrack(B_,iEvent,iSetup);
  preparemcParticle(genpart);
  prepareHits( ecalRecHitHandleEB,hcalRecHitHandleHBHE,HGCeeRecHitHandle,HGChefRecHitHandle,HGChebRecHitHandle,pv,iEvent,iSetup );
  //preparemcParticle(genpart); //put before prepareHits() to have mc info, for mip calib check
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,PandoraApi::ProcessEvent(*m_pPandora));
  //  preparePFO(iEvent,iSetup);
  preparePFO(iEvent,iSetup,HGCeeRecHitHandle,HGChefRecHitHandle,HGChebRecHitHandle,pv);
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,PandoraApi::Reset(*m_pPandora));

}

void runPandora::initPandoraCalibrParameters()
{
   m_firstMCpartEta = -99999.;
   m_firstMCpartPhi = -99999.;
   m_secondMCpartEta = -99999.;
   m_secondMCpartPhi = -99999.;


   m_Calibr_ADC2GeV_EE  = 0.000012 ; //w/o absorber thickness correction
   m_Calibr_ADC2GeV_HEF = 0.0000176; //w/o absorber thickness correction
   m_Calibr_ADC2GeV_HEB = 0.0003108; //w/o absorber thickness correction

   m_EM_addCalibrEE = 10.;
   m_EM_addCalibrHEF= 10.;
   m_EM_addCalibrHEB= 10.;
   m_HAD_addCalibrEE = 10.;
   m_HAD_addCalibrHEF= 10.;
   m_HAD_addCalibrHEB= 10.;


   m_eCalThresBarrel    = 0.;
   m_eCalThresEndCap    = 27.55e-6; //EE
   m_hCalThresEndCapHEF = 42.50e-6;
   m_hCalThresEndCapHEB = 742.2e-6;
   m_hCalThresBarrel    = 0.;


   m_eCalMipThresBarrel    = 0.5;
   m_eCalMipThresEndCap    = 0.5;
   m_hCalMipThresEndCapHEF = 0.5;
   m_hCalMipThresEndCapHEB = 0.5;
   m_hCalMipThresBarrel    = 0.5;

   m_eCalToMipBarrel       = 3.3333333;
   m_eCalToMipEndCap       = 18149.;
   m_hCalToMipEndCapHEF    = 11765.;
   m_hCalToMipEndCapHEB    = 667.4;
   m_hCalToMipBarrel       = 3.3333333;

   m_eCalToEMGeVEndCap     = 1.;
   m_eCalToEMGeVBarrel     = 1.;
   m_hCalToEMGeVEndCapHEF  = 1.;
   m_hCalToEMGeVEndCapHEB  = 1.;
   m_hCalToEMGeVBarrel     = 1.;

   m_eCalToHadGeVEndCap    = 1.;
   m_eCalToHadGeVBarrel    = 1.;
   m_hCalToHadGeVEndCapHEF = 1.;
   m_hCalToHadGeVEndCapHEB = 1.;
   m_hCalToHadGeVBarrel    = 1.;
   m_muonToMip             = 1.;

   return;
}

void runPandora::readCalibrParameterFile()
{

  std::ifstream calibrParFile(m_calibrationParameterFile.fullPath().c_str() , std::ifstream::in );

   if (!calibrParFile.is_open()) {
      std::cout << "runPandora::readCalibrParameterFile: calibrParFile does not exist ("
         << m_calibrationParameterFile << ")" << std::endl;
      return;
   }

   while ( !calibrParFile.eof() ) {
      std::string linebuf;
      getline( calibrParFile, linebuf );
      if (linebuf.substr(0,1) == "#") continue;
      if (linebuf.substr(0,2) == "//") continue;

      if (linebuf.empty()) continue;

      std::string paraName;
      double paraValue;
      std::stringstream ss(linebuf);
      ss >> paraName >> paraValue;
      std::cout << "reading calibr parameter " << paraName << " ";

      if (paraName=="Calibr_ADC2GeV_EE"     ) {m_Calibr_ADC2GeV_EE       = paraValue; std::cout <<  m_Calibr_ADC2GeV_EE     << std::endl;}
      if (paraName=="Calibr_ADC2GeV_HEF"    ) {m_Calibr_ADC2GeV_HEF      = paraValue; std::cout <<  m_Calibr_ADC2GeV_HEF    << std::endl;}
      if (paraName=="Calibr_ADC2GeV_HEB"    ) {m_Calibr_ADC2GeV_HEB      = paraValue; std::cout <<  m_Calibr_ADC2GeV_HEB    << std::endl;}

      if (paraName=="EMaddCalibrEE"         ) {m_EM_addCalibrEE          = paraValue; std::cout <<  m_EM_addCalibrEE        << std::endl;}
      if (paraName=="EMaddCalibrHEF"        ) {m_EM_addCalibrHEF         = paraValue; std::cout <<  m_EM_addCalibrHEF       << std::endl;}
      if (paraName=="EMaddCalibrHEB"        ) {m_EM_addCalibrHEB         = paraValue; std::cout <<  m_EM_addCalibrHEB       << std::endl;}
      if (paraName=="HADaddCalibrEE"        ) {m_HAD_addCalibrEE         = paraValue; std::cout <<  m_HAD_addCalibrEE       << std::endl;}
      if (paraName=="HADaddCalibrHEF"       ) {m_HAD_addCalibrHEF        = paraValue; std::cout <<  m_HAD_addCalibrHEF      << std::endl;}
      if (paraName=="HADaddCalibrHEB"       ) {m_HAD_addCalibrHEB        = paraValue; std::cout <<  m_HAD_addCalibrHEB      << std::endl;}


      if (paraName=="ECalThresBarrel"       ) {m_eCalThresBarrel         = paraValue; std::cout <<  m_eCalThresBarrel       << std::endl;}
      if (paraName=="ECalThresEndCap"       ) {m_eCalThresEndCap         = paraValue; std::cout <<  m_eCalThresEndCap       << std::endl;}
      if (paraName=="HCalThresEndCapHEF"    ) {m_hCalThresEndCapHEF      = paraValue; std::cout <<  m_hCalThresEndCapHEF    << std::endl;}
      if (paraName=="HCalThresEndCapHEB"    ) {m_hCalThresEndCapHEB      = paraValue; std::cout <<  m_hCalThresEndCapHEB    << std::endl;}
      if (paraName=="HCalThresBarrel"       ) {m_hCalThresBarrel         = paraValue; std::cout <<  m_hCalThresBarrel       << std::endl;}


      if (paraName=="ECalMipThresEndCap"    ) {m_eCalMipThresEndCap      = paraValue; std::cout <<  m_eCalMipThresEndCap    << std::endl;}
      if (paraName=="ECalMipThresBarrel"    ) {m_eCalMipThresBarrel      = paraValue; std::cout <<  m_eCalMipThresBarrel    << std::endl;}
      if (paraName=="HCalMipThresEndCapHEF" ) {m_hCalMipThresEndCapHEF   = paraValue; std::cout <<  m_hCalMipThresEndCapHEF << std::endl;}
      if (paraName=="HCalMipThresEndCapHEB" ) {m_hCalMipThresEndCapHEB   = paraValue; std::cout <<  m_hCalMipThresEndCapHEB << std::endl;}
      if (paraName=="HCalMipThresBarrel"    ) {m_hCalMipThresBarrel      = paraValue; std::cout <<  m_hCalMipThresBarrel    << std::endl;}

      if (paraName=="ECalToMipEndCap"       ) {m_eCalToMipEndCap         = paraValue; std::cout <<  m_eCalToMipEndCap       << std::endl;}
      if (paraName=="ECalToMipBarrel"       ) {m_eCalToMipBarrel         = paraValue; std::cout <<  m_eCalToMipBarrel       << std::endl;}
      if (paraName=="HCalToMipEndCapHEF"    ) {m_hCalToMipEndCapHEF      = paraValue; std::cout <<  m_hCalToMipEndCapHEF    << std::endl;}
      if (paraName=="HCalToMipEndCapHEB"    ) {m_hCalToMipEndCapHEB      = paraValue; std::cout <<  m_hCalToMipEndCapHEB    << std::endl;}
      if (paraName=="HCalToMipBarrel"       ) {m_hCalToMipBarrel         = paraValue; std::cout <<  m_hCalToMipBarrel       << std::endl;}

      if (paraName=="ECalToEMGeVEndCap"     ) {m_eCalToEMGeVEndCap       = paraValue; std::cout <<  m_eCalToEMGeVEndCap     << std::endl;}
      if (paraName=="ECalToEMGeVBarrel"     ) {m_eCalToEMGeVBarrel       = paraValue; std::cout <<  m_eCalToEMGeVBarrel     << std::endl;}
      if (paraName=="HCalToEMGeVEndCapHEF"  ) {m_hCalToEMGeVEndCapHEF    = paraValue; std::cout <<  m_hCalToEMGeVEndCapHEF  << std::endl;}
      if (paraName=="HCalToEMGeVEndCapHEB"  ) {m_hCalToEMGeVEndCapHEB    = paraValue; std::cout <<  m_hCalToEMGeVEndCapHEB  << std::endl;}
      if (paraName=="HCalToEMGeVBarrel"     ) {m_hCalToEMGeVBarrel       = paraValue; std::cout <<  m_hCalToEMGeVBarrel     << std::endl;}

      if (paraName=="ECalToHadGeVEndCap"    ) {m_eCalToHadGeVEndCap      = paraValue; std::cout <<  m_eCalToHadGeVEndCap    << std::endl;}
      if (paraName=="ECalToHadGeVBarrel"    ) {m_eCalToHadGeVBarrel      = paraValue; std::cout <<  m_eCalToHadGeVBarrel    << std::endl;}
      if (paraName=="HCalToHadGeVEndCapHEF" ) {m_hCalToHadGeVEndCapHEF   = paraValue; std::cout <<  m_hCalToHadGeVEndCapHEF << std::endl;}
      if (paraName=="HCalToHadGeVEndCapHEB" ) {m_hCalToHadGeVEndCapHEB   = paraValue; std::cout <<  m_hCalToHadGeVEndCapHEB << std::endl;}
      if (paraName=="HCalToHadGeVBarrel"    ) {m_hCalToHadGeVBarrel      = paraValue; std::cout <<  m_hCalToHadGeVBarrel    << std::endl;}

      if (paraName=="MuonToMip"             ) {m_muonToMip               = paraValue; std::cout <<  m_muonToMip             << std::endl;}
   }

   calibrParFile.close();

   return;
}


void runPandora::readEnergyWeight()
{
   std::cout << "runPandora::readEnergyWeight" << std::endl;

   //FIXME : for the moment, everything is given in unit of MIP

   stm->addArrayParameter("layerSet_EE");
   stm->addArrayParameter("layerSet_HEF");
   stm->addArrayParameter("layerSet_HEB");
   stm->addArrayParameter("energyWeight_EM_EE");
   stm->addArrayParameter("energyWeight_EM_HEF");
   stm->addArrayParameter("energyWeight_EM_HEB");
   stm->addArrayParameter("energyWeight_Had_EE");
   stm->addArrayParameter("energyWeight_Had_HEF");
   stm->addArrayParameter("energyWeight_Had_HEB");
   
   stm->addSingleParameter("offSet_EM");
   stm->addSingleParameter("offSet_Had");

   stm->read();
   stm->printPars();

   offSet_EM  = stm->getSinglePara("offSet_EM");
   offSet_Had = stm->getSinglePara("offSet_Had");

   //notneeded layerSet_EE          = stm->getArrayPara("layerSet_EE");
   //notneeded layerSet_HEF         = stm->getArrayPara("layerSet_HEF");
   //notneeded layerSet_HEB         = stm->getArrayPara("layerSet_HEB");
   //notneeded energyWeight_EM_EE   = stm->getArrayPara("energyWeight_EM_EE");
   //notneeded energyWeight_EM_HEF  = stm->getArrayPara("energyWeight_EM_HEF");
   //notneeded energyWeight_EM_HEB  = stm->getArrayPara("energyWeight_EM_HEB");
   //notneeded energyWeight_Had_EE  = stm->getArrayPara("energyWeight_Had_EE");
   //notneeded energyWeight_Had_HEF = stm->getArrayPara("energyWeight_Had_HEF");
   //notneeded energyWeight_Had_HEB = stm->getArrayPara("energyWeight_Had_HEB");


   return;
}



void runPandora::prepareGeometry(const edm::EventSetup& iSetup){ // function to setup a geometry for pandora

  // std::cout << "I am preparing my geometry!!!" << std::endl ; 

  // 
  // Add ECAL/HCAL parameters to geometry
  //
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
  // Get the ecal/hcal barrel, endcap geometry
  const CaloSubdetectorGeometry *ebtmp = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  const CaloSubdetectorGeometry *hbtmp = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);

  // Additions
  edm::ESHandle<HGCalGeometry> hgceeGeoHandle ; 
  edm::ESHandle<HGCalGeometry> hgchefGeoHandle ; 
  edm::ESHandle<HGCalGeometry> hgchebGeoHandle ; 
  
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeGeoHandle) ; 
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hgchefGeoHandle) ; 
  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",hgchebGeoHandle) ; 

  const HGCalGeometry &hgceetmp = *hgceeGeoHandle ; 
  const HGCalGeometry &hgcheftmp = *hgchefGeoHandle ; 
  const HGCalGeometry &hgchebtmp = *hgchebGeoHandle ; 

  // std::cout << "Basic check: " << hgceetmp.producerTag() << " " << hgcheftmp.producerTag() << " " << hgchebtmp.producerTag() << std::endl ; 
  
  std::vector<DetId> ecalBarrelCells = geoHandle->getValidDetIds(DetId::Ecal, EcalBarrel);
  std::vector<DetId> ecalEndcapCells = hgceeGeoHandle->getValidDetIds(DetId::Forward, HGCEE);
  std::vector<DetId> hcalBarrelCells = geoHandle->getValidDetIds(DetId::Hcal, HcalBarrel);
  std::vector<DetId> hcalEndcapCellsFront = hgchefGeoHandle->getValidDetIds(DetId::Forward, HGCHEF);
  std::vector<DetId> hcalEndcapCellsBack  = hgchebGeoHandle->getValidDetIds(DetId::Forward, HGCHEB);
  
  const EcalBarrelGeometry* ecalBarrelGeometry = dynamic_cast< const EcalBarrelGeometry* > (ebtmp);
  const HcalGeometry* hcalBarrelGeometry = dynamic_cast< const HcalGeometry* > (hbtmp);
  
  assert( ecalBarrelGeometry );
  assert( hcalBarrelGeometry );
  
  const HGCalGeometry &HGCEEGeometry  = dynamic_cast< const HGCalGeometry& > (hgceetmp);
  const HGCalGeometry &HGCHEFGeometry = dynamic_cast< const HGCalGeometry& > (hgcheftmp);
  const HGCalGeometry &HGCHEBGeometry = dynamic_cast< const HGCalGeometry& > (hgchebtmp);
  assert( &HGCEEGeometry );
  assert( &HGCHEFGeometry );
  assert( &HGCHEBGeometry );

  PandoraApi::Geometry::SubDetector::Parameters *ebParameters = new PandoraApi::Geometry::SubDetector::Parameters();
  PandoraApi::Geometry::SubDetector::Parameters *eeParameters = new PandoraApi::Geometry::SubDetector::Parameters();
  PandoraApi::Geometry::SubDetector::Parameters *hbParameters = new PandoraApi::Geometry::SubDetector::Parameters();
  PandoraApi::Geometry::SubDetector::Parameters *heParameters = new PandoraApi::Geometry::SubDetector::Parameters();

  PandoraApi::Geometry::LayerParameters *ebLayerParameters    = new PandoraApi::Geometry::LayerParameters();
  PandoraApi::Geometry::LayerParameters *hbLayerParameters    = new PandoraApi::Geometry::LayerParameters();

  std::vector<PandoraApi::Geometry::LayerParameters*> hgcEELayerParameters;
  std::vector<PandoraApi::Geometry::LayerParameters*> hgcHEFLayerParameters;
  std::vector<PandoraApi::Geometry::LayerParameters*> hgcHEBLayerParameters;

  //todo: un-hardcode these values
  unsigned int nHGCeeLayers = 32, nHGChefLayers = 32, nHGChebLayers = 21 ; 
  std::vector<double> min_innerR_depth_ee, min_innerZ_depth_ee ; 
  std::vector<double> min_innerR_depth_hef, min_innerZ_depth_hef ; 
  std::vector<double> min_innerR_depth_heb, min_innerZ_depth_heb ;
  for (unsigned int i=0; i<nHGCeeLayers; i++) { 
    PandoraApi::Geometry::LayerParameters *eeLayerParameters;
    eeLayerParameters = new PandoraApi::Geometry::LayerParameters();
    hgcEELayerParameters.push_back( eeLayerParameters ) ; 
    min_innerR_depth_ee.push_back( 99999.0 ) ; 
    min_innerZ_depth_ee.push_back( 99999.0 ) ; 

    PandoraApi::Geometry::LayerParameters *hefLayerParameters;
    hefLayerParameters = new PandoraApi::Geometry::LayerParameters();
    hgcHEFLayerParameters.push_back( hefLayerParameters ) ; 
    min_innerR_depth_hef.push_back( 99999.0 ) ; 
    min_innerZ_depth_hef.push_back( 99999.0 ) ; 

    if ( i < nHGChebLayers ) { 
      PandoraApi::Geometry::LayerParameters *hebLayerParameters;
      hebLayerParameters = new PandoraApi::Geometry::LayerParameters();
      hgcHEBLayerParameters.push_back( hebLayerParameters ) ; 
      min_innerR_depth_heb.push_back( 99999.0 ) ; 
      min_innerZ_depth_heb.push_back( 99999.0 ) ; 
    }
  }

  // dummy vectors for CalculateCornerSubDetectorParameters function
  std::vector<double> min_innerR_depth_eb, min_innerZ_depth_eb ; 
  std::vector<double> min_innerR_depth_hb, min_innerZ_depth_hb ; 
  
  // To be enabled
  // PandoraApi::Geometry::LayerParameters *hgceeLayerParameters  = new PandoraApi::Geometry::LayerParameters();
  // PandoraApi::Geometry::LayerParameters *hgchefLayerParameters = new PandoraApi::Geometry::LayerParameters();
  // PandoraApi::Geometry::LayerParameters *hgchebLayerParameters = new PandoraApi::Geometry::LayerParameters();
  
  SetDefaultSubDetectorParameters("EcalBarrel", pandora::ECAL_BARREL, *ebParameters);
  SetDefaultSubDetectorParameters("EcalEndcap", pandora::ECAL_ENDCAP, *eeParameters);
  SetDefaultSubDetectorParameters("HcalBarrel", pandora::HCAL_BARREL, *hbParameters);
  SetDefaultSubDetectorParameters("HcalEndcap", pandora::HCAL_ENDCAP, *heParameters);

  //corner parameters for EB
  double min_innerRadius = 99999.0 ; double max_outerRadius = 0.0 ;
  double min_innerZ = 99999.0 ; double max_outerZ = 0.0 ;
  CalculateCornerSubDetectorParameters(*ecalBarrelGeometry, ecalBarrelCells, pandora::ECAL_BARREL, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       false, min_innerR_depth_eb, min_innerZ_depth_eb);
  SetCornerSubDetectorParameters(*ebParameters, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ); // One ECAL layer
  SetSingleLayerParameters(*ebParameters,*ebLayerParameters);
  ebParameters->m_nLayers = 1 ; // One ECAL layer
  
  //corner parameters for HB
  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;
  CalculateCornerSubDetectorParameters(*hcalBarrelGeometry, hcalBarrelCells, pandora::HCAL_BARREL, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       false, min_innerR_depth_hb, min_innerZ_depth_hb);
  SetCornerSubDetectorParameters(*hbParameters, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ);
  SetSingleLayerParameters(*hbParameters, *hbLayerParameters);
  hbParameters->m_nLayers = 1 ; //todo: include HCAL layers
  
  //corner & layer parameters for EE
  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;
  int nLayers = 0;
  CalculateCornerSubDetectorParameters(HGCEEGeometry, ecalEndcapCells, pandora::ECAL_ENDCAP, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       true, min_innerR_depth_ee, min_innerZ_depth_ee);
  SetCornerSubDetectorParameters(*eeParameters, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ);
  SetMultiLayerParameters(*eeParameters, hgcEELayerParameters, min_innerR_depth_ee, min_innerZ_depth_ee, nHGCeeLayers, nLayers);
  eeParameters->m_nLayers = nLayers ; // HACK(?) to account for the invalid layer 0(???)

  //corner & layer parameters for HE
  //consider both HEF and HEB together
  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;
  nLayers = 0;
  CalculateCornerSubDetectorParameters(HGCHEFGeometry, hcalEndcapCellsFront, pandora::HCAL_ENDCAP, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       true, min_innerR_depth_hef, min_innerZ_depth_hef);
  CalculateCornerSubDetectorParameters(HGCHEBGeometry, hcalEndcapCellsBack, pandora::HCAL_ENDCAP, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       true, min_innerR_depth_heb, min_innerZ_depth_heb);
  SetCornerSubDetectorParameters(*heParameters, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ);
  SetMultiLayerParameters(*heParameters, hgcHEFLayerParameters, min_innerR_depth_hef, min_innerZ_depth_hef, nHGChefLayers, nLayers);
  SetMultiLayerParameters(*heParameters, hgcHEBLayerParameters, min_innerR_depth_heb, min_innerZ_depth_heb, nHGChebLayers, nLayers);
  heParameters->m_nLayers = nLayers ; // HACK

  // std::cout << "before set GEO" << std::endl;
  // std::cout << "Idle check: " << geometryParameters.m_InnerRCoordinate << std::endl ; 
  PandoraApi::Geometry::SubDetector::Create(*m_pPandora, *ebParameters);
  PandoraApi::Geometry::SubDetector::Create(*m_pPandora, *hbParameters);
  PandoraApi::Geometry::SubDetector::Create(*m_pPandora, *eeParameters);
  PandoraApi::Geometry::SubDetector::Create(*m_pPandora, *heParameters);
  // std::cout << "after set GEO" << std::endl;

}

void runPandora::SetDefaultSubDetectorParameters(const std::string &subDetectorName, const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const {
  //identification
  parameters.m_subDetectorName = subDetectorName;
  parameters.m_subDetectorType = subDetectorType;
  
  // Phi Coordinate is when start drawing the detector, wrt x-axis.  
  // Assuming this is 0 since CMS ranges from -pi to pi
  parameters.m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  parameters.m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  
  // Symmetry order is how you draw the "polygon" detector.  
  // Circle approximation for now (0), but can be configured to match N(cells)
  parameters.m_innerSymmetryOrder = 0 ; 
  parameters.m_outerSymmetryOrder = 0 ; 
  
  parameters.m_isMirroredInZ = true ; // Duplicate detector +/- z
}

void runPandora::CalculateCornerSubDetectorParameters(const CaloSubdetectorGeometry& geom,  const std::vector<DetId>& cells, const pandora::SubDetectorType subDetectorType, 
                                                      double& min_innerRadius, double& max_outerRadius, double& min_innerZ, double& max_outerZ,
                                                      bool doLayers, std::vector<double>& min_innerR_depth, std::vector<double>& min_innerZ_depth) const
{
  //barrel:
  // Inner radius taken as average magnitude of (x,y) for corners 0-3
  // Outer radius taken as average magnitude of (x,y) for corners 4-7  
  double ci_barrel[4] = {0,1,2,3};
  double co_barrel[4] = {4,5,6,7};
  
  //endcap:
  // Inner radius taken as average magnitude of (x,y) for corners 0,3,4,7 
  // Outer radius taken as average magnitude of (x,y) for corners 1,2,5,6  
  double ci_endcap[4] = {0,3,4,7};
  double co_endcap[4] = {1,2,5,6};
  
  double *ci, *co;
  if(subDetectorType==pandora::ECAL_BARREL || subDetectorType==pandora::HCAL_BARREL){
    ci = ci_barrel;
    co = co_barrel;
  }
  else { //if(subDetectorType==pandora::ECAL_ENDCAP || subDetectorType==pandora::HCAL_ENDCAP){
    ci = ci_endcap;
    co = co_endcap;
  }
  
  // Determine: inner/outer detector radius
  for (std::vector<DetId>::const_iterator ib=cells.begin(); ib!=cells.end(); ib++) {
    const CaloCellGeometry *thisCell = geom.getGeometry(*ib);
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
    
    //kind of hacky
    unsigned int layer = 0;
    if(doLayers && subDetectorType==pandora::ECAL_ENDCAP) layer = (unsigned int) ((HGCEEDetId)(*ib)).layer() ;
    else if(doLayers && subDetectorType==pandora::HCAL_ENDCAP) layer = (unsigned int) ((HGCHEDetId)(*ib)).layer() ;    
    
    //inner radius calculation
    double avgX_inner = 0.25 * (corners[ci[0]].x() + corners[ci[1]].x() + corners[ci[2]].x() + corners[ci[3]].x()) ;
    double avgY_inner = 0.25 * (corners[ci[0]].y() + corners[ci[1]].y() + corners[ci[2]].y() + corners[ci[3]].y()) ;
    double innerRadius = sqrt( avgX_inner * avgX_inner + avgY_inner * avgY_inner ) ;
    if ( innerRadius < min_innerRadius ) min_innerRadius = innerRadius ;
    if ( doLayers && innerRadius < min_innerR_depth.at(layer) ) min_innerR_depth.at(layer) = innerRadius ; 
    
    //outer radius calculation
    double avgX_outer = 0.25 * (corners[co[0]].x() + corners[co[1]].x() + corners[co[2]].x() + corners[co[3]].x()) ;
    double avgY_outer = 0.25 * (corners[co[0]].y() + corners[co[1]].y() + corners[co[2]].y() + corners[co[3]].y()) ;
    double outerRadius = sqrt( avgX_outer * avgX_outer + avgY_outer * avgY_outer ) ;
    if ( outerRadius > max_outerRadius ) max_outerRadius = outerRadius ;
    
    //z calculation
    for( unsigned int isubcell = 0; isubcell<8; isubcell++){
        if( fabs(corners[isubcell].z()) < min_innerZ ) min_innerZ = fabs(corners[isubcell].z()) ;
        if ( corners[isubcell].z() > max_outerZ ) max_outerZ = corners[isubcell].z();
        if ( doLayers && fabs(corners[isubcell].z()) < min_innerZ_depth.at(layer) ) min_innerZ_depth.at(layer) = fabs(corners[isubcell].z()) ;
    }
  }

}

void runPandora::SetCornerSubDetectorParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, 
                                                const double& min_innerRadius, const double& max_outerRadius, const double& min_innerZ, const double& max_outerZ) const 
{
  parameters.m_innerRCoordinate = min_innerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  parameters.m_outerRCoordinate = max_outerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  parameters.m_innerZCoordinate = min_innerZ * 10.0 ; // CMS units cm, Pandora expects mm
  parameters.m_outerZCoordinate = max_outerZ * 10.0 ; // CMS units cm, Pandora expects mm
}

void runPandora::SetSingleLayerParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, PandoraApi::Geometry::LayerParameters &layerParameters) const {
  layerParameters.m_closestDistanceToIp = parameters.m_innerRCoordinate; 
  layerParameters.m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  layerParameters.m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
  parameters.m_layerParametersList.push_back(layerParameters) ; 
}

void runPandora::SetMultiLayerParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, std::vector<PandoraApi::Geometry::LayerParameters*> &layerParameters,
                                         std::vector<double>& min_innerR_depth, std::vector<double>& min_innerZ_depth, const unsigned int& nTotalLayers, int& nLayers) const 
{
  for (unsigned int i=0; i<nTotalLayers; i++) { 
    double distToIP = 10.0 * sqrt(min_innerR_depth.at(i)*min_innerR_depth.at(i) + min_innerZ_depth.at(i)*min_innerZ_depth.at(i)) ; 
    layerParameters.at(i)->m_closestDistanceToIp = distToIP ; 
    layerParameters.at(i)->m_nInteractionLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
    layerParameters.at(i)->m_nRadiationLengths = 0.0 ; // No idea what to say here.  Include the tracker material?
    if ( distToIP < 10000.0 ) { 
      nLayers++ ; 
      parameters.m_layerParametersList.push_back(*(layerParameters.at(i))) ; 
    }
  }
}

//void runPandora::prepareTrack( math::XYZVector B_, const reco::RecoToSimCollection pRecoToSim, const edm::Event& iEvent, const edm::EventSetup& iSetup){ // function to setup tracks in an event for pandora
void runPandora::prepareTrack( math::XYZVector B_, const edm::Event& iEvent, const edm::EventSetup& iSetup){ // function to setup tracks in an event for pandora
  PandoraApi::Track::Parameters trackParameters;
  //We need the speed of light
  double speedoflight = (CLHEP::c_light*CLHEP::mm)/CLHEP::ns;
  std::cout<< speedoflight << " mm/ns" << std::endl;

 std::cout << "prepareTrack 1 " << inputTagGeneralTracks_.size() << std::endl ;

  for (unsigned int istr=0; istr<inputTagGeneralTracks_.size();istr++){
    
 std::cout << "prepareTrack 2 " << std::endl;

    //Track collection
    // edm::Handle<reco::TrackCollection> tkRefCollection;
    edm::Handle<edm::View<reco::Track> > tkRefCollection;
    bool found1 = iEvent.getByLabel(inputTagGeneralTracks_[istr], tkRefCollection);
    if(!found1 ) {
      std::ostringstream err;
      err<<"cannot find generalTracks: "<< inputTagGeneralTracks_[istr];
      LogError("runPandora")<<err.str()<<std::endl;
      throw cms::Exception( "MissingProduct", err.str());
    } 
        

 std::cout << "prepareTrack 3 " << tkRefCollection->size() << std::endl;

    for(reco::TrackCollection::size_type i=0; i<tkRefCollection->size(); i++) {
      
      std::cout << "prepareTrack 5 " << std::endl;

      const reco::Track * track = &(*tkRefCollection)[i];
    
      //For the d0 = -dxy
      trackParameters.m_d0 = track->d0() * 10. ; //in mm
      //For the z0
      trackParameters.m_z0 = track->dz() * 10. ; //in mm
      //For the Track momentum at the 2D distance of closest approach
      //For tracks reconstructed in the CMS Tracker, the reference position is the point of closest approach to the centre of CMS. (math::XYZPoint posClosest = track->referencePoint();)
      // According to TrackBase.h the momentum() method returns momentum vector at innermost (reference) point on track
      const pandora::CartesianVector momentumAtClosestApproach(track->momentum().x(),track->momentum().y(),track->momentum().z()); //in GeV
      trackParameters.m_momentumAtDca = momentumAtClosestApproach;
 
      //For the track of the state at the start in mm and GeV
      const pandora::CartesianVector positionAtStart(track->innerPosition().x()* 10.,track->innerPosition().y()* 10., track->innerPosition().z() * 10. );
      const pandora::CartesianVector momentumAtStart(track->innerMomentum().x(),track->innerMomentum().y(), track->innerMomentum().z() );
      trackParameters.m_trackStateAtStart = pandora::TrackState(positionAtStart,momentumAtStart);
      //For the track of the state at the end in mm and GeV
      const pandora::CartesianVector positionAtEnd(track->outerPosition().x() * 10.,track->outerPosition().y() * 10., track->outerPosition().z() * 10.);
      const pandora::CartesianVector momentumAtEnd(track->outerMomentum().x(),track->outerMomentum().y(), track->outerMomentum().z() );
      trackParameters.m_trackStateAtEnd = pandora::TrackState(positionAtEnd,momentumAtEnd);
      //For the charge
      double charge = track->charge();
      // std::cout << "charge " << charge << std::endl;
      trackParameters.m_charge = charge;
      //Associate the reconstructed Track (in the Tracker) with the corresponding MC true simulated particle
      trackParameters.m_particleId = 211; // INITIALIZATION // NS

      edm::RefToBase<reco::Track> tr(tkRefCollection, i);
      std::vector<std::pair<TrackingParticleRef, double> > tp;
      TrackingParticleRef tpr; 

      if(pRecoToSim.find(tr) != pRecoToSim.end()){
    tp = pRecoToSim[tr];
    std::cout << "Reco Track pT: "  << track->pt() <<  " matched to " << tp.size() << " MC Tracks" << " associated with quality: " << tp.begin()->second << std::endl;
    tpr = tp.begin()->first;
      }
      


      //So the pdg code of the track is 
      if(pRecoToSim.find(tr) != pRecoToSim.end()) {

       std::cout << " EDW  KAI PALI tpr->pdgId() " << tpr->pdgId() << std::endl;

       trackParameters.m_particleId = (tpr->pdgId());
       std::cout << "the pdg id of this track is " << (tpr->pdgId()) << std::endl;
       //The parent vertex (from which this track was produced) has daughter particles.
       //These are the desire siblings of this track which we need. 
       TrackingParticleRefVector simSiblings = getTpSiblings(tpr);

       const TrackingParticle * sib; 
       int numofsibs = 0;
       std::vector<int> pdgidofsibs; pdgidofsibs.clear();

        
      if (simSiblings.isNonnull()) {
    for(TrackingParticleRefVector::iterator si = simSiblings.begin(); si != simSiblings.end(); si++){
      //Check if the current sibling is the track under study
      if ( (*si) ==  tpr  ) {continue;}
      sib = &(**si);
      pdgidofsibs.push_back(sib->pdgId());
      ++numofsibs;
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora, track, sib)); 
    }
    std::cout<< "This track has " << numofsibs << " sibling tracks with pdgids:" << std::endl; 
    for (std::vector<int>::iterator sib_pdg_it = pdgidofsibs.begin(); sib_pdg_it != pdgidofsibs.end(); sib_pdg_it++){
      std::cout << (*sib_pdg_it) << std::endl;
    }
      } else {
    std::cout << "Particle pdgId = "<< (tpr->pdgId()) << " produced at rho = " << (tpr->vertex().Rho()) << ", z = " << (tpr->vertex().Z()) << ", has NO siblings!"<< std::endl;
      }
     
      //Now the track under study has daughter particles. To find them we study the decay vertices of the track
      TrackingParticleRefVector simDaughters = getTpDaughters(tpr);
      const TrackingParticle * dau; 
      int numofdaus = 0;
      std::vector<int> pdgidofdaus; pdgidofdaus.clear();

      if (simDaughters.isNonnull()) {
    for(TrackingParticleRefVector::iterator di = simDaughters.begin(); di != simDaughters.end(); di++){
      //We have already checked that simDaughters don't contain the track under study
      dau = &(**di);
      pdgidofdaus.push_back(dau->pdgId());
      ++numofdaus;
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(*m_pPandora, track, dau)); 
    }
    std::cout<< "This track has " << numofdaus << " daughter tracks with pdgids:" << std::endl; 
    for (std::vector<int>::iterator dau_pdg_it = pdgidofdaus.begin(); dau_pdg_it != pdgidofdaus.end(); dau_pdg_it++){
      std::cout << (*dau_pdg_it) << std::endl;
    }
      } else {
    std::cout << "Particle pdgId = "<< (tpr->pdgId()) << " produced at rho = " << (tpr->vertex().Rho()) << ", z = " << (tpr->vertex().Z()) << ", has NO daughters!"<< std::endl;
      }

} // KLEINW TO IF pRecoToSim.find(tr) != pRecoToSim.end() NS FIX


      //The mass 
      trackParameters.m_mass = pandora::PdgTable::GetParticleMass(trackParameters.m_particleId.Get());
 

      //For the ECAL entrance
      // Starting from outermost hit position of the track and propagating to ECAL Entrance
      double pfoutenergy = track->outerMomentum().Mag2();

      //the input BaseParticlePropagator needs to be cm and GeV
      BaseParticlePropagator theOutParticle = BaseParticlePropagator( RawParticle(XYZTLorentzVector(track->outerMomentum().x(),
                                                    track->outerMomentum().y(),
                                                    track->outerMomentum().z(),
                                                    pfoutenergy),
                                          XYZTLorentzVector(track->outerPosition().x(),
                                                    track->outerPosition().y(),
                                                    track->outerPosition().z(),
                                                    0.)),
                                      0.,0.,B_.z());

      theOutParticle.setCharge(track->charge());
      math::XYZPoint theOutParticle_position = math::XYZPoint(theOutParticle.vertex());
      math::XYZTLorentzVector theOutParticle_momentum = theOutParticle.momentum();
      std::cout << "magnetic field z " << B_.z() << std::endl;
      std::cout << "theOutParticle x position before propagation in cm "<< theOutParticle_position.x()<< std::endl;
      std::cout << "theOutParticle x momentum before propagation in cm "<< theOutParticle_momentum.x()<< std::endl;
      theOutParticle.propagateToEcalEntrance(false);
      bool reachesCalorimeter = false;
      bool isonendcap = false;

      //Set position and momentum after propagation to ECAL
      theOutParticle_position = math::XYZPoint(theOutParticle.vertex());
      theOutParticle_momentum = theOutParticle.momentum();
 
      std::cout << "theOutParticle x position after propagation to ECAL "<< theOutParticle_position.x()<< std::endl;
      std::cout << "theOutParticle x momentum after propagation to ECAL "<< theOutParticle_momentum.x()<< std::endl;
     
      if(theOutParticle.getSuccess()!=0){
    // std::cout<< "!!!Reached ECAL!!! "<< std::endl;
    reachesCalorimeter = true;
      }
      if( abs(theOutParticle.getSuccess()) == 2){
    // std::cout<< "It is on the endcaps "<< std::endl;
    isonendcap = true;
      }

      trackParameters.m_reachesCalorimeter = reachesCalorimeter;
      if (reachesCalorimeter){
    const pandora::CartesianVector positionAtCalorimeter(theOutParticle_position.x() * 10.,theOutParticle_position.y() * 10.,theOutParticle_position.z() * 10.);//in mm
    const pandora::CartesianVector momentumAtCalorimeter(theOutParticle_momentum.x(),theOutParticle_momentum.y(),theOutParticle_momentum.z());
    trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(positionAtCalorimeter, momentumAtCalorimeter);
    // For the time at calorimeter we need the speed of light
    //This is in BaseParticlePropagator c_light() method in mm/ns but is protected (299.792458 mm/ns)
    //So we take it from CLHEP
    trackParameters.m_timeAtCalorimeter = positionAtCalorimeter.GetMagnitude() / speedoflight; // in ns
      
      } else { 
    trackParameters.m_trackStateAtCalorimeter = trackParameters.m_trackStateAtEnd.Get();
    //trackParameters.m_timeAtCalorimeter = std::numeric_limits<double>::max();
    trackParameters.m_timeAtCalorimeter = 999999999.;
      }

/*
      trackParameters.m_isProjectedToEndCap = isonendcap; 

      bool canFormPfo = false; 
      bool canFormClusterlessPfo = false;
 
      //Add more criteria here
      if (trackParameters.m_reachesCalorimeter.Get()){
    canFormPfo = true;
    canFormClusterlessPfo = true;
    std::cout<< "Yes, this track can form pfo" << std::endl;
      }
      trackParameters.m_canFormPfo = canFormPfo;
      trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;
*/


      trackParameters.m_isProjectedToEndCap = isonendcap;

      bool canFormPfo = false;
      bool canFormClusterlessPfo = false;

      double mom = sqrt((theOutParticle_momentum.x()*theOutParticle_momentum.x()+theOutParticle_momentum.y()*theOutParticle_momentum.y()+theOutParticle_momentum.z()*theOutParticle_momentum.z()));

      if(trackParameters.m_reachesCalorimeter.Get()>0 && mom>2. && track->quality(reco::TrackBase::tight) ){ // SCZ change in hard coded pt cut from diff_runpandora
        canFormPfo = true;
        canFormClusterlessPfo = true;
      }
      trackParameters.m_canFormPfo = canFormPfo;
      trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;


     //The parent address
      trackParameters.m_pParentAddress =  (void *) track;

      //Some cout
      // std::cout <<  track->innerDetId() << std::endl;
      // std::cout <<  track->outerPx() << std::endl;
      // std::cout <<  track->d0() << std::endl;
      // std::cout <<  track->dz() << std::endl;
      // std::cout <<  pfrt->pdgCode() << std::endl;


 
       PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));
    }

    // PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraMonitoringApi::VisualizeTracks(  &(*tkRefCollection)  , "currentTrackList", AUTO, false, true  ) );
    // PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraMonitoringApi::VisualizeTracks(  &(*tkRefCollection)  , "currentTrackList",  true  ) );
    // PANDORA_MONITORING_API(VisualizeTracks(  &(*tkRefCollection)  , "currentTrackList", AUTO, false, true  ) );
    // PANDORA_MONITORING_API(ViewEvent() );


  }

}

void runPandora::prepareHits( edm::Handle<EcalRecHitCollection> ecalRecHitHandleEB,
                  edm::Handle<HBHERecHitCollection> hcalRecHitHandleHBHE, 
                  edm::Handle<HGCRecHitCollection> HGCeeRecHitHandle,
                  edm::Handle<HGCRecHitCollection> HGChefRecHitHandle,
                  edm::Handle<HGCRecHitCollection> HGChebRecHitHandle,
                  reco::Vertex& pv, 
                  const edm::Event& iEvent, const edm::EventSetup& iSetup){

  PandoraApi::RectangularCaloHitParameters caloHitParameters;

  double speedoflight = (CLHEP::c_light/CLHEP::cm)/CLHEP::ns;
  std::cout<< speedoflight << " cm/ns" << std::endl;

  sumallhitsEMenergyB4Clustering_[0] = 0;   // for - side
  sumallhitsEMenergyB4Clustering_[1] = 0;   // for + side

  sumallhitsHGCalenergyB4Clustering_[0] = 0;   // for - side
  sumallhitsHGCalenergyB4Clustering_[1] = 0;   // for + side

  double sumCaloEnergy = 0.;
  double sumCaloEnergyEM = 0.;
  double sumCaloEnergyHAD = 0.;

  double simDir_sumCaloEnergyEM = 0.;
  double simDir_sumCaloEnergyHAD = 0.;


  double sumCaloECALEnergyEM = 0.;
  double sumCaloHCALEnergyEM  = 0.;
  double sumCaloECALEnergyHAD= 0.;
  double sumCaloECALEnergyHAD_unc= 0.;
  double sumCaloHCALEnergyHAD = 0.;

  // 
  // Process ECAL barrel rechits 
  // 
  for(unsigned i=0; i<ecalRecHitHandleEB->size(); i++) {
    continue;      
    const EcalRecHit * erh = &(*ecalRecHitHandleEB)[i];
    const DetId& detid = erh->detid();
    double energy = erh->energy();

    if (energy < m_eCalThresBarrel) continue;

    double time = erh->time();
    // std::cout << "energy " << energy <<  " time " << time <<std::endl;
    EcalSubdetector esd=(EcalSubdetector)detid.subdetId();
    if (esd != 1) continue;
      
    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
    // get the ecalBarrel geometry
    const CaloSubdetectorGeometry *ebtmp = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  
    const EcalBarrelGeometry* ecalBarrelGeometry = dynamic_cast< const EcalBarrelGeometry* > (ebtmp);
    assert( ecalBarrelGeometry );

    // get the ecalBarrel topology
    EcalBarrelTopology ecalBarrelTopology(geoHandle);

    const CaloSubdetectorGeometry* geom = ecalBarrelGeometry;
    
    const CaloCellGeometry *thisCell = geom->getGeometry(detid);
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraECAL") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      continue;
    }
        
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
 //   for (unsigned int i=0; i<8; i++) { 
 //     std::cout << "Corners " << i << ": x " << corners[i].x() << " y " << corners[i].y() << " z " << corners[i].z() << std::endl ; 
 //   }
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    // dist = cm, c = cm/nsec, rechit t in psec
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; 

    caloHitParameters.m_hitType = pandora::ECAL;
    caloHitParameters.m_hitRegion = pandora::BARREL;
    caloHitParameters.m_inputEnergy = energy;
    caloHitParameters.m_electromagneticEnergy = m_eCalToEMGeVBarrel * energy; 
    caloHitParameters.m_mipEquivalentEnergy = m_eCalToMipBarrel * energy; // = energy; // HTAN 0 NS AP


    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_eCalMipThresBarrel) {
       //std::cout << "EcalBarrel MIP threshold rejected" << std::endl;
       continue;
    }

    sumCaloEnergy += energy;
    sumCaloEnergyEM += energy  ;//FIXME * absorberCorrectionEM;
    sumCaloEnergyHAD += energy ;//FIXME * absorberCorrectionHAD;
    sumCaloECALEnergyEM  += energy  ;//* absorberCorrectionEM;
    sumCaloECALEnergyHAD += energy  ;//* absorberCorrectionHAD;
 
    caloHitParameters.m_hadronicEnergy = m_eCalToHadGeVBarrel * energy; // = energy; 
    caloHitParameters.m_layer = 1.;//PFLayer::ECAL_BARREL;
    caloHitParameters.m_nCellRadiationLengths = 0.0; // 6.;
    caloHitParameters.m_nCellInteractionLengths = 0.0; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) erh;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 
    
  }    

  //
  // process HCAL Barrel Hits
  //
  for(unsigned i=0; i<hcalRecHitHandleHBHE->size(); i++) {
        continue;       
    const HBHERecHit * hrh = &(*hcalRecHitHandleHBHE)[i];
    const DetId& detid = hrh->detid();
    double energy = hrh->energy();
    if (energy < m_hCalThresBarrel) continue;
    double time = hrh->time();

    HcalSubdetector hsd=(HcalSubdetector)detid.subdetId();
    if (hsd != 1) continue; // HCAL BARREL ONLY

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
    // get the hcal geometry
    const CaloSubdetectorGeometry *hbtmp = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
    const HcalGeometry* hcalBarrelGeometry = dynamic_cast< const HcalGeometry* > (hbtmp);
    assert( hcalBarrelGeometry );

    const CaloSubdetectorGeometry* hbgeom = hcalBarrelGeometry;
    const CaloCellGeometry *thisCell = hbgeom->getGeometry(detid) ;
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraHCAL") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      continue;
    }
        
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    // dist = cm, c = cm/nsec, rechit t in psec
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; 

    caloHitParameters.m_hitType = pandora::HCAL;
    caloHitParameters.m_hitRegion = pandora::BARREL ;
    caloHitParameters.m_inputEnergy = energy;
    caloHitParameters.m_electromagneticEnergy = m_hCalToEMGeVBarrel * energy; 
    caloHitParameters.m_mipEquivalentEnergy =  m_hCalToMipBarrel * energy; // = energy; // HTAN 0 NS AP  

    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_hCalMipThresBarrel) {
       //std::cout << "HcalBarrel MIP threshold rejected" << std::endl;
       continue;
    }
    sumCaloEnergy += energy;
    sumCaloEnergyEM  += energy;// * absorberCorrectionEM ;
    sumCaloEnergyHAD += energy;// * absorberCorrectionHAD;
    sumCaloHCALEnergyEM  += energy  ;//* absorberCorrectionEM;
    sumCaloHCALEnergyHAD += energy  ;//* absorberCorrectionHAD;


    caloHitParameters.m_hadronicEnergy = m_hCalToHadGeVBarrel * energy; // = energy; 

    caloHitParameters.m_layer = 1.;
    caloHitParameters.m_nCellRadiationLengths = 0.0; // 6.;
    caloHitParameters.m_nCellInteractionLengths = 0.0; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) hrh;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
  }

  //
  // Process HGC EE rec hits 
  // 

  int nNotFoundEE = 0, nCaloHitsEE = 0 ; 
  for(unsigned i=0; i<HGCeeRecHitHandle->size(); i++) {
      
    const HGCRecHit *eerh = &(*HGCeeRecHitHandle)[i];
    const HGCEEDetId& detid = eerh->id();
    //double energy = eerh->energy();
    double energy = eerh->energy() * m_Calibr_ADC2GeV_EE;
    if (energy<m_eCalThresEndCap) continue;
    //energy *= m_EM_addCalibrEE; //this must be after threshold rejection
    double time = eerh->time();

    double nCellInteractionLengths = 0.;
    double nCellRadiationLengths = 0.;
    double absorberCorrectionEM = 1.;
    double absorberCorrectionHAD = 1.;
    getLayerPropertiesEE(eerh, detid.layer(),
          nCellInteractionLengths, nCellRadiationLengths,
          absorberCorrectionEM, absorberCorrectionHAD);


     //std::cout << "HGC EE rechit cell " << detid.cell() << ", sector " << detid.sector() 
     //          << ", subsector " << detid.subsector() << ", layer " << detid.layer() << ", z " << detid.zside() 
     //          << " energy " << energy <<  " and time " << time << std::endl;

    ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
    if (thesubdet != 3) continue; // HGC EE
      
    math::XYZVector axis;

    edm::ESHandle<HGCalGeometry> hgceeGeoHandle ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeGeoHandle) ; 
    const HGCalGeometry &hgceetmp = *hgceeGeoHandle ; 

    const HGCalGeometry &hgceeGeometry = dynamic_cast< const HGCalGeometry& > (hgceetmp);
    assert( &hgceeGeometry );

    // Get the HGC topology
    // edm::ESHandle<HGCalTopology> hgceeTopoHandle ; 
    // iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeTopoHandle) ; 
    // const HGCalTopology &hgceeTopology = *hgceeTopoHandle ; 
    
    // const HGCalGeometry geom = hgceeGeometry;
    const FlatTrd *thisCell = static_cast<const FlatTrd*>(hgceeGeometry.getGeometry(detid));    
    // const CaloCellGeometry *thisCell = hgceeGeometry.getGeometry(detid);
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraHGCEE") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      nNotFoundEE++ ; 
      continue;
    }
        
    const HGCalGeometry::CornersVec corners = ( std::move( hgceeGeometry.getCorners( detid ) ) );
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    //Hieu
    //caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    //caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; // dist = cm, c = cm/nsec, rechit t in psec
    
    caloHitParameters.m_hitType = pandora::ECAL;
    caloHitParameters.m_hitRegion = pandora::ENDCAP;
    caloHitParameters.m_inputEnergy = energy;
    if (m_energyCorrMethod == "ABSCORR")
       caloHitParameters.m_electromagneticEnergy = m_eCalToEMGeVEndCap * energy * absorberCorrectionEM * m_EM_addCalibrEE; 
    else if (m_energyCorrMethod == "WEIGHTING") {
       //std::cout << "ENERGY CORR: WEIGHTING" << std::endl;
       //int x_layer = caloHitParameters.m_layer.Get()+1;

       //std::cout << "11111111111111111111111111111111" << std::endl;
       //double erjcorr =  getEnergyWeight(D_HGCEE, x_layer, EM);
       //std::cout << "22222222222222222222222222222222" << std::endl;
       //caloHitParameters.m_electromagneticEnergy = m_eCalToEMGeVEndCap * energy * erjcorr; 
       caloHitParameters.m_electromagneticEnergy = m_eCalToEMGeVEndCap * energy * getEnergyWeight(D_HGCEE, detid.layer(), EM) * m_eCalToMipEndCap; 
    }
    caloHitParameters.m_mipEquivalentEnergy =   m_eCalToMipEndCap   * energy; // = energy; // HTAN 0 NS AP  
    
    double angleCorrectionMIP(1.); 
    double hitR = distToFrontFace; 
    double hitZ = zf;
    angleCorrectionMIP = hitR/hitZ; 

    //---choose hits with eta-phi close to init. mc particle
    TVector3 hit3v(xf,yf,zf);
    double hitEta = hit3v.PseudoRapidity();
    double hitPhi = hit3v.Phi();
    h_hit_Eta -> Fill(hitEta);
    h_hit_Phi -> Fill(hitPhi);
    //std::cout << "TEST: m_firstMCpartEta = " << m_firstMCpartEta
    //   << ", hitEta = " << hitEta << std::endl;
    if (std::fabs(hitEta-m_firstMCpartEta) < 0.05
          && std::fabs(hitPhi-m_firstMCpartPhi) < 0.05) {
       h_MIP_EE -> Fill(caloHitParameters.m_mipEquivalentEnergy.Get());
       h_MIP_Corr_EE -> Fill(caloHitParameters.m_mipEquivalentEnergy.Get()*angleCorrectionMIP);
    }
    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_eCalMipThresEndCap) {
       //std::cout << "EE MIP threshold rejected" << std::endl;
       continue;
    }
    sumCaloEnergy    += energy;
    sumCaloEnergyEM  += energy * absorberCorrectionEM * m_EM_addCalibrEE;
    sumCaloEnergyHAD += energy * absorberCorrectionHAD * m_HAD_addCalibrEE;
//    if ( ( std::fabs(hitEta-m_firstMCpartEta) < 0.5
//             || std::fabs(hitEta+m_firstMCpartEta) < 0.5 )
//          && std::fabs(hitPhi-m_firstMCpartPhi)< 0.5) 
    if ( (std::fabs(hitEta-m_firstMCpartEta) < 0.2
          && std::fabs(hitPhi-m_firstMCpartPhi) < 0.2)
          ||
          (std::fabs(hitEta-m_secondMCpartEta) < 0.2
           && std::fabs(hitPhi-m_secondMCpartPhi) < 0.2)
       )
    {
       simDir_sumCaloEnergyEM  += energy * absorberCorrectionEM * m_EM_addCalibrEE;
       simDir_sumCaloEnergyHAD += energy * absorberCorrectionHAD * m_HAD_addCalibrEE;
    }

    if(z<0) sumallhitsEMenergyB4Clustering_[0] += energy* absorberCorrectionEM * m_EM_addCalibrEE;
    if(z>=0) sumallhitsEMenergyB4Clustering_[1] += energy* absorberCorrectionEM * m_EM_addCalibrEE;

    if(z<0) sumallhitsHGCalenergyB4Clustering_[0] += energy* absorberCorrectionEM * m_EM_addCalibrEE;
    if(z>=0) sumallhitsHGCalenergyB4Clustering_[1] += energy* absorberCorrectionEM * m_EM_addCalibrEE;

    sumCaloECALEnergyEM  += energy  * absorberCorrectionEM * m_EM_addCalibrEE;
    sumCaloECALEnergyHAD += energy  * absorberCorrectionHAD * m_HAD_addCalibrEE;
    sumCaloECALEnergyHAD_unc += energy ;


    if (m_energyCorrMethod == "ABSCORR")
    caloHitParameters.m_hadronicEnergy = m_eCalToHadGeVEndCap * energy * absorberCorrectionHAD * m_HAD_addCalibrEE;
    else if (m_energyCorrMethod == "WEIGHTING")
       caloHitParameters.m_hadronicEnergy = m_eCalToHadGeVEndCap * energy * getEnergyWeight(D_HGCEE, detid.layer(), HAD) * m_eCalToMipEndCap; 

    caloHitParameters.m_layer = detid.layer() ;//PFLayer::ECAL_BARREL;
    caloHitParameters.m_nCellRadiationLengths = nCellRadiationLengths; // 6.;
    caloHitParameters.m_nCellInteractionLengths = nCellInteractionLengths; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) eerh;

    //std::cout << "TEST input caloHitParameters: "
    //   << ", hadronicEnergy = " << caloHitParameters.m_hadronicEnergy.Get()
    //   << ", hitType = " << caloHitParameters.m_hitType.Get()
    //   << ", m_hitRegion = " << caloHitParameters.m_hitRegion.Get()
    //   << std::endl;

    // std::cout << "Parameters for input: " << std::endl ; 
    // std::cout << "position vector X,Y,Z: " << caloHitParameters.m_positionVector.GetX() << " " 
    //           << caloHitParameters.m_positionVector.GetY() << " " << caloHitParameters.m_positionVector.GetZ() << std::endl ; 
    // std::cout << "expected direction X,Y,Z: " << caloHitParameters.m_expectedDirection.GetX() << " " 
    //           << caloHitParameters.m_expectedDirection.GetY() << " " << caloHitParameters.m_expectedDirection.GetZ() << std::endl ; 
    

    // std::cout << "Processing HGC EE rec hit with position (in cm) " << x << " " << y << " " << z << std::endl ; 

    m_hitEperLayer_EM [caloHitParameters.m_layer.Get()] += caloHitParameters.m_electromagneticEnergy.Get();
    m_hitEperLayer_HAD[caloHitParameters.m_layer.Get()] += caloHitParameters.m_hadronicEnergy.Get();
    //std::cout << "HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << std::endl;
    //std::cout << "EM ENERGY: " << caloHitParameters.m_electromagneticEnergy.Get()
    //   << std::endl;
    
    m_hitEperLayer_EM_EE [caloHitParameters.m_layer.Get()] += caloHitParameters.m_electromagneticEnergy.Get();
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 
    nCaloHitsEE++ ; 
     //std::cout << "BMD: This calo hit has been created" << std::endl ;   
  }


  //
  // Process HGC HEF rec hits 
  // 
  int nNotFoundHEF = 0, nCaloHitsHEF = 0 ; 
  for(unsigned i=0; i<HGChefRecHitHandle->size(); i++) {
      
    const HGCRecHit *hefrh = &(*HGChefRecHitHandle)[i];
    const HGCHEDetId& detid = hefrh->id();
    double energy = hefrh->energy() * m_Calibr_ADC2GeV_HEF;
    if (energy<m_hCalThresEndCapHEF)
       continue;
    //energy *= m_addCalibrHEF;
    double time = hefrh->time();

    double nCellInteractionLengths = 0.;
    double nCellRadiationLengths = 0.;
    double absorberCorrectionEM = 1.;
    double absorberCorrectionHAD = 1.;
    getLayerPropertiesHEF(hefrh, detid.layer(),
      nCellInteractionLengths, nCellRadiationLengths,
      absorberCorrectionEM, absorberCorrectionHAD);


    // std::cout << "HGC HEF rechit cell " << detid.cell() << ", sector " << detid.sector() 
    //           << ", subsector " << detid.subsector() << ", layer " << detid.layer() << ", z " << detid.zside() 
    //           << " energy " << energy <<  " and time " << time << std::endl;

    ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
    if (thesubdet != 4) continue; // HGC HEF
      
    math::XYZVector axis;

    edm::ESHandle<HGCalGeometry> hgchefGeoHandle ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hgchefGeoHandle) ; 
    const HGCalGeometry &hgcheftmp = *hgchefGeoHandle ; 

    const HGCalGeometry &hgchefGeometry = dynamic_cast< const HGCalGeometry& > (hgcheftmp);
    assert( &hgchefGeometry );

    // Get the HGC topology
    // edm::ESHandle<HGCalTopology> hgceeTopoHandle ; 
    // iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeTopoHandle) ; 
    // const HGCalTopology &hgceeTopology = *hgceeTopoHandle ; 
    
    // const HGCalGeometry geom = hgceeGeometry;
    const FlatTrd *thisCell = static_cast<const FlatTrd*>(hgchefGeometry.getGeometry(detid));    
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraHGCHEF") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      nNotFoundHEF++ ; 
      continue;
    }
        
    const HGCalGeometry::CornersVec corners = ( std::move( hgchefGeometry.getCorners( detid ) ) );
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; // dist = cm, c = cm/nsec, rechit t in psec
    
    caloHitParameters.m_hitType = pandora::HCAL;
    caloHitParameters.m_hitRegion = pandora::ENDCAP;
    caloHitParameters.m_inputEnergy = energy;
    if (m_energyCorrMethod == "ABSCORR")
       caloHitParameters.m_electromagneticEnergy = m_hCalToEMGeVEndCapHEF * energy * m_EM_addCalibrHEF; 
    else if (m_energyCorrMethod == "WEIGHTING")
       caloHitParameters.m_electromagneticEnergy = m_hCalToEMGeVEndCapHEF * energy * getEnergyWeight(D_HGCHEF, detid.layer(), EM); 

    caloHitParameters.m_mipEquivalentEnergy =  m_hCalToMipEndCapHEF * energy; // = energy; // HTAN 0 NS AP  

    //---choose hits with eta-phi close to init. mc particle
    TVector3 hit3v(xf,yf,zf);
    double hitEta = hit3v.PseudoRapidity();
    double hitPhi = hit3v.Phi();
    if (std::fabs(hitEta-m_firstMCpartEta) < 0.05
          && std::fabs(hitPhi-m_firstMCpartPhi)< 0.05) {

       h_MIP_HEF->Fill(caloHitParameters.m_mipEquivalentEnergy.Get());
       double angleCorrectionMIP(1.);
       double hitR = distToFrontFace; 
       double hitZ = zf;
       angleCorrectionMIP = hitR/hitZ;
       h_MIP_Corr_HEF->Fill(caloHitParameters.m_mipEquivalentEnergy.Get()*angleCorrectionMIP);
    }

    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_hCalMipThresEndCapHEF) {
       //std::cout << "HEF MIP threshold rejected" << std::endl;
       continue;
    }
    sumCaloEnergy += energy;
    sumCaloEnergyEM  += energy * absorberCorrectionEM * m_EM_addCalibrHEF;
    sumCaloEnergyHAD += energy * absorberCorrectionHAD * m_HAD_addCalibrHEF;

    if(z<0) sumallhitsEMenergyB4Clustering_[0] += energy * absorberCorrectionEM * m_EM_addCalibrHEF;
    if(z>=0) sumallhitsEMenergyB4Clustering_[1] += energy * absorberCorrectionEM * m_EM_addCalibrHEF;

    if(z<0) sumallhitsHGCalenergyB4Clustering_[0] += energy * absorberCorrectionEM * m_EM_addCalibrHEF;
    if(z>=0) sumallhitsHGCalenergyB4Clustering_[1] += energy * absorberCorrectionEM * m_EM_addCalibrHEF;

//    if ( (std::fabs(hitEta-m_firstMCpartEta) < 0.5
//             || std::fabs(hitEta+m_firstMCpartEta) < 0.5 )
//          && std::fabs(hitPhi-m_firstMCpartPhi)< 0.5) 
    if ( (std::fabs(hitEta-m_firstMCpartEta) < 0.2
          && std::fabs(hitPhi-m_firstMCpartPhi) < 0.2)
          ||
          (std::fabs(hitEta-m_secondMCpartEta) < 0.2
           && std::fabs(hitPhi-m_secondMCpartPhi) < 0.2)
       )
    {
       simDir_sumCaloEnergyEM  += energy * absorberCorrectionEM * m_EM_addCalibrHEF;
       simDir_sumCaloEnergyHAD += energy * absorberCorrectionHAD * m_HAD_addCalibrHEF;
    }
    sumCaloHCALEnergyEM  += energy  * absorberCorrectionEM * m_EM_addCalibrHEF;
    sumCaloHCALEnergyHAD += energy  * absorberCorrectionHAD * m_HAD_addCalibrHEF;

    if (m_energyCorrMethod == "ABSCORR")
    caloHitParameters.m_hadronicEnergy = m_hCalToHadGeVEndCapHEF * energy * m_HAD_addCalibrHEF; // = energy; 
    else if (m_energyCorrMethod == "WEIGHTING")
    caloHitParameters.m_hadronicEnergy = m_hCalToHadGeVEndCapHEF * energy * getEnergyWeight(D_HGCHEF, detid.layer(), HAD) * m_hCalToMipEndCapHEF;

    caloHitParameters.m_layer = detid.layer() ;
    caloHitParameters.m_nCellRadiationLengths = nCellRadiationLengths; // 6.;
    caloHitParameters.m_nCellInteractionLengths = nCellInteractionLengths; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) hefrh;

    // std::cout << "Parameters for input: " << std::endl ; 
    // std::cout << "position vector X,Y,Z: " << caloHitParameters.m_positionVector.GetX() << " " 
    //           << caloHitParameters.m_positionVector.GetY() << " " << caloHitParameters.m_positionVector.GetZ() << std::endl ; 
    // std::cout << "expected direction X,Y,Z: " << caloHitParameters.m_expectedDirection.GetX() << " " 
    //           << caloHitParameters.m_expectedDirection.GetY() << " " << caloHitParameters.m_expectedDirection.GetZ() << std::endl ; 
    

    // std::cout << "Processing HGC HEF rec hit with position (in cm) " << x << " " << y << " " << z << std::endl ; 
    m_hitEperLayer_EM [caloHitParameters.m_layer.Get()+31] += caloHitParameters.m_electromagneticEnergy.Get();
    m_hitEperLayer_HAD[caloHitParameters.m_layer.Get()+31] += caloHitParameters.m_hadronicEnergy.Get();
    m_hitEperLayer_EM_HEF[caloHitParameters.m_layer.Get()] += caloHitParameters.m_electromagneticEnergy.Get();
  
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 

    // std::cout << "BMD: This calo hit has been created" << std::endl ;   
    nCaloHitsHEF++ ; 
  }


  //
  // Process HGC HEB rec hits 
  // 
  int nNotFoundHEB = 0, nCaloHitsHEB = 0 ; 
  for(unsigned i=0; i<HGChebRecHitHandle->size(); i++) {
      
    const HGCRecHit *hebrh = &(*HGChebRecHitHandle)[i];
    const HGCHEDetId& detid = hebrh->id();
    double energy = hebrh->energy() * m_Calibr_ADC2GeV_HEB;
    if (energy<m_hCalThresEndCapHEB) continue;
    //energy *= m_addCalibrHEB;
    double time = hebrh->time();
    double nCellInteractionLengths = 0.;
    double nCellRadiationLengths = 0.;
    double absorberCorrectionEM = 1.;
    double absorberCorrectionHAD = 1.;
    getLayerPropertiesHEB(hebrh, detid.layer(),
      nCellInteractionLengths, nCellRadiationLengths,
      absorberCorrectionEM, absorberCorrectionHAD);


     //std::cout << "HGC HEB rechit cell " << detid.cell() << ", sector " << detid.sector() 
     //          << ", subsector " << detid.subsector() << ", layer " << detid.layer() << ", z " << detid.zside() 
     //          << " energy " << energy <<  " and time " << time << std::endl;

    ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
    if (thesubdet != 5) continue; // HGC HEB
      
    math::XYZVector axis;

    edm::ESHandle<HGCalGeometry> hgchebGeoHandle ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",hgchebGeoHandle) ; 
    const HGCalGeometry &hgchebtmp = *hgchebGeoHandle ; 

    const HGCalGeometry &hgchebGeometry = dynamic_cast< const HGCalGeometry& > (hgchebtmp);
    assert( &hgchebGeometry );

    // Get the HGC topology
    // edm::ESHandle<HGCalTopology> hgceeTopoHandle ; 
    // iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeTopoHandle) ; 
    // const HGCalTopology &hgceeTopology = *hgceeTopoHandle ; 
    
    // const HGCalGeometry geom = hgceeGeometry;
    const FlatTrd *thisCell = static_cast<const FlatTrd*>(hgchebGeometry.getGeometry(detid));    
    // const CaloCellGeometry *thisCell = hgchebGeometry.getGeometry(detid);
  
    // find rechit geometry
    if(!thisCell) {
      LogError("runPandoraHGCHEB") << "warning detid "<<detid.rawId() <<" not found in geometry"<<std::endl;
      nNotFoundHEB++ ; 
      continue;
    }
        
    const HGCalGeometry::CornersVec corners = ( std::move( hgchebGeometry.getCorners( detid ) ) );
    assert( corners.size() == 8 );
    const pandora::CartesianVector NECorner( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector SECorner( corners[1].x(), corners[1].y(), corners[1].z() );

    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;

    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();

    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; // dist = cm, c = cm/nsec, rechit t in psec
    
    caloHitParameters.m_hitType = pandora::HCAL;
    caloHitParameters.m_hitRegion = pandora::ENDCAP;
    caloHitParameters.m_inputEnergy = energy;
    if (m_energyCorrMethod == "ABSCORR")
    caloHitParameters.m_electromagneticEnergy = m_hCalToEMGeVEndCapHEB * energy * m_EM_addCalibrHEB; 
    else if (m_energyCorrMethod == "WEIGHTING")
    caloHitParameters.m_electromagneticEnergy = m_hCalToEMGeVEndCapHEB * energy * getEnergyWeight(D_HGCHEB, detid.layer(), EM) * m_hCalToMipEndCapHEB; 
    caloHitParameters.m_mipEquivalentEnergy = m_hCalToMipEndCapHEB * energy ; /// m_addCalibrHEB; // = energy; // HTAN 0 NS AP  

    //---choose hits with eta-phi close to init. mc particle
    TVector3 hit3v(xf,yf,zf);
    double hitEta = hit3v.PseudoRapidity();
    double hitPhi = hit3v.Phi();
    if (std::fabs(hitEta-m_firstMCpartEta) < 0.05
          && std::fabs(hitPhi-m_firstMCpartPhi)< 0.05) {

     h_MIP_HEB->Fill(caloHitParameters.m_mipEquivalentEnergy.Get());
    }
    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_hCalMipThresEndCapHEB) {
       //std::cout << "HEB MIP threshold rejected" << std::endl;
       continue;
    }
    sumCaloEnergy += energy;
    sumCaloEnergyEM  += energy * absorberCorrectionEM * m_EM_addCalibrHEB;
    sumCaloEnergyHAD += energy * absorberCorrectionHAD * m_HAD_addCalibrHEB;

    if(z<0) sumallhitsEMenergyB4Clustering_[0] += energy * absorberCorrectionEM * m_EM_addCalibrHEB;
    if(z>=0) sumallhitsEMenergyB4Clustering_[1] += energy * absorberCorrectionEM * m_EM_addCalibrHEB;

    if(z<0) sumallhitsHGCalenergyB4Clustering_[0] += energy * absorberCorrectionEM * m_EM_addCalibrHEB;
    if(z>=0) sumallhitsHGCalenergyB4Clustering_[1] += energy * absorberCorrectionEM * m_EM_addCalibrHEB;

//    if ( (std::fabs(hitEta-m_firstMCpartEta) < 0.5
//          || std::fabs(hitEta-m_firstMCpartEta) < 0.5 )
//          && std::fabs(hitPhi-m_firstMCpartPhi)< 0.5) 
    if ( (std::fabs(hitEta-m_firstMCpartEta) < 0.2
             && std::fabs(hitPhi-m_firstMCpartPhi) < 0.2)
          ||
          (std::fabs(hitEta-m_secondMCpartEta) < 0.2
           && std::fabs(hitPhi-m_secondMCpartPhi) < 0.2)
       )   
    {
       simDir_sumCaloEnergyEM  += energy * absorberCorrectionEM * m_EM_addCalibrHEB;
       simDir_sumCaloEnergyHAD += energy * absorberCorrectionHAD * m_HAD_addCalibrHEB;
    }

    sumCaloHCALEnergyEM  += energy  * absorberCorrectionEM * m_EM_addCalibrHEB;
    sumCaloHCALEnergyHAD += energy  * absorberCorrectionHAD * m_HAD_addCalibrHEB;

    if (m_energyCorrMethod == "ABSCORR")
    caloHitParameters.m_hadronicEnergy = m_hCalToHadGeVEndCapHEB * energy * m_HAD_addCalibrHEB; // = energy; 
    else if (m_energyCorrMethod == "WEIGHTING")
    caloHitParameters.m_hadronicEnergy = m_hCalToHadGeVEndCapHEB * energy * getEnergyWeight(D_HGCHEB, detid.layer(), HAD) * m_hCalToMipEndCapHEB;
    //caloHitParameters.m_layer = 31 + detid.layer();//PFLayer::ECAL_BARREL;
    caloHitParameters.m_layer = detid.layer();//PFLayer::ECAL_BARREL;
    caloHitParameters.m_nCellRadiationLengths = nCellRadiationLengths; // 6.;
    caloHitParameters.m_nCellInteractionLengths = nCellInteractionLengths; // 6.;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *) hebrh;

    // std::cout << "Parameters for input: " << std::endl ; 
    // std::cout << "position vector X,Y,Z: " << caloHitParameters.m_positionVector.GetX() << " " 
    //           << caloHitParameters.m_positionVector.GetY() << " " << caloHitParameters.m_positionVector.GetZ() << std::endl ; 
    // std::cout << "expected direction X,Y,Z: " << caloHitParameters.m_expectedDirection.GetX() << " " 
    //           << caloHitParameters.m_expectedDirection.GetY() << " " << caloHitParameters.m_expectedDirection.GetZ() << std::endl ; 
    

    // std::cout << "Processing HGC HEB rec hit with position (in cm) " << x << " " << y << " " << z << std::endl ; 
    m_hitEperLayer_EM [caloHitParameters.m_layer.Get()+31+12] += caloHitParameters.m_electromagneticEnergy.Get();
    m_hitEperLayer_HAD[caloHitParameters.m_layer.Get()+31+12] += caloHitParameters.m_hadronicEnergy.Get();
    m_hitEperLayer_EM_HEB[caloHitParameters.m_layer.Get()] += caloHitParameters.m_electromagneticEnergy.Get();
  
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 

    nCaloHitsHEB++ ; 
    // std::cout << "BMD: This calo hit has been created" << std::endl ;   
  }

  h_sumCaloE->Fill(sumCaloEnergy);
  h_sumCaloEM->Fill(sumCaloEnergyEM);
  h_sumCaloHad->Fill(sumCaloEnergyHAD);
  h_simDir_sumCaloEM ->Fill(simDir_sumCaloEnergyEM );
  h_simDir_sumCaloHad->Fill(simDir_sumCaloEnergyHAD);

  h2_Calo_EM_hcalEecalE->Fill(sumCaloECALEnergyEM, sumCaloHCALEnergyEM);
  h2_Calo_Had_hcalEecalE->Fill(sumCaloECALEnergyHAD, sumCaloHCALEnergyHAD);
  h_sumEcalEEM->Fill(sumCaloECALEnergyEM);
  h_sumHcalEEM->Fill(sumCaloHCALEnergyEM);

  h_sumEcalEHad->Fill(sumCaloECALEnergyHAD);
  h_sumEcalEHad_unc->Fill(sumCaloECALEnergyHAD_unc);
  if(sumCaloECALEnergyHAD>=0.5) h_sumEcalEHadc->Fill(sumCaloECALEnergyHAD);
  if(sumCaloECALEnergyHAD<0.5) h_sumHcalEHadc->Fill(sumCaloHCALEnergyHAD);
  h_sumHcalEHad->Fill(sumCaloHCALEnergyHAD);
  h_sumEHad->Fill(sumCaloHCALEnergyHAD+sumCaloECALEnergyHAD);

  for (int ilay=0; ilay<100; ilay++) {
     int ibin = ilay+1;
     h_hitEperLayer_EM->SetBinContent(ibin,m_hitEperLayer_EM[ilay]+ h_hitEperLayer_EM->GetBinContent(ibin));
     h_hitEperLayer_HAD->SetBinContent(ibin,m_hitEperLayer_HAD[ilay]+ h_hitEperLayer_HAD->GetBinContent(ibin));

     h_hitEperLayer_EM_EE ->SetBinContent(ibin,m_hitEperLayer_EM_EE [ilay]+ h_hitEperLayer_EM_EE ->GetBinContent(ibin));
     h_hitEperLayer_EM_HEF->SetBinContent(ibin,m_hitEperLayer_EM_HEF[ilay]+ h_hitEperLayer_EM_HEF->GetBinContent(ibin));
     h_hitEperLayer_EM_HEB->SetBinContent(ibin,m_hitEperLayer_EM_HEB[ilay]+ h_hitEperLayer_EM_HEB->GetBinContent(ibin));
  }



  std::cout << "sumCaloEnergy = " << sumCaloEnergy << std::endl;
  std::cout << "sumCaloEnergyEM  = " << sumCaloEnergyEM  << std::endl;
  std::cout << "sumCaloEnergyHAD = " << sumCaloEnergyHAD << std::endl;

  std::cout << "prepareHits HGC summary: " << std::endl ; 
  std::cout << "HGC Calo Hits               : " << nCaloHitsEE << " (HGC EE) " 
        << nCaloHitsHEF << " (HGC HEF) " << nCaloHitsHEB << " (HGC HEB) " << std::endl ;
  std::cout << "DetIDs not found in geometry: " << nNotFoundEE << " (HGC EE) " 
        << nNotFoundHEF << " (HGC HEF) " << nNotFoundHEB << " (HGC HEB) " << std::endl ;

}


void runPandora::preparemcParticle(edm::Handle<std::vector<reco::GenParticle> > genpart){ // function to setup a mcParticle for pandora
  // PandoraPFANew/v00-09/include/Api/PandoraApi.h
  //class MCParticleParameters
  //{
  //public:
  //    pandora::InputFloat             m_energy;                   ///< The energy of the MC particle, units GeV
  //    pandora::InputCartesianVector   m_momentum;                 ///< The momentum of the MC particle, units GeV
  //    pandora::InputCartesianVector   m_vertex;                   ///< The production vertex of the MC particle, units mm
  //    pandora::InputCartesianVector   m_endpoint;                 ///< The endpoint of the MC particle, units mm
  //    pandora::InputInt               m_particleId;               ///< The MC particle's ID (PDG code)
  //    pandora::InputAddress           m_pParentAddress;           ///< Address of the parent MC particle in the user framework
  //};
  // for(std::vector<reco::GenParticle>::const_iterator cP = genpart->begin();  cP != genpart->end(); cP++ ) {
      
   const GenParticle * firstMCp = &(*genpart)[0];
   if (firstMCp) {
      m_firstMCpartEta = firstMCp->eta();
      m_firstMCpartPhi = firstMCp->phi();
   }
   if (genpart->size()>=2) {
      const GenParticle * secondMCp = &(*genpart)[1];
      if (secondMCp) {
         m_secondMCpartEta = secondMCp->eta();
         m_secondMCpartPhi = secondMCp->phi();
      }
   } 
  
   //h_MCp_Eta->Fill(m_firstMCpartEta);
   //h_MCp_Phi->Fill(m_firstMCpartPhi);
      
  
  RminVtxDaughter[0] = 999999.; //initialise for each event
  RminVtxDaughter[1] = 999999.; //initialise for each event
                                //FIXME Attention will crash for one particle sample
  ZminVtxDaughter[0] = 999999.; //initialise for each event
  ZminVtxDaughter[1] = 999999.; //initialise for each event
                                //FIXME Attention will crash for one particle sample
  
     
  isDecayedBeforeCalo[0] = 0;
  isDecayedBeforeCalo[1] = 0;

  for(size_t i = 0; i < genpart->size(); ++ i) {
    const GenParticle * pa = &(*genpart)[i];
    PandoraApi::MCParticle::Parameters parameters;
    parameters.m_energy = pa->energy(); 
    parameters.m_momentum = pandora::CartesianVector(pa->px() , pa->py(),  pa->pz() );
    parameters.m_vertex = pandora::CartesianVector(pa->vx() * 10. , pa->vy() * 10., pa->vz() * 10. ); //in mm
        
        
    // parameters.m_endpoint = pandora::CartesianVector(position.x(), position.y(), position.z());
    // Definition of the enpoint depends on the application that created the particle, e.g. the start point of the shower in a calorimeter.
    // If the particle was not created as a result of a continuous process where the parent particle continues, i.e.
    // hard ionization, Bremsstrahlung, elastic interactions, etc. then the vertex of the daughter particle is the endpoint.
    parameters.m_endpoint = pandora::CartesianVector(pa->vx() * 10. , pa->vy() * 10., pa->vz() * 10. ); //IS THIS CORRECT?! //NO, should be where it starts to decay
    parameters.m_particleId = pa->pdgId();
    parameters.m_mcParticleType = pandora::MCParticleType::MC_3D;
    parameters.m_pParentAddress = (void*) pa;
    if(i==0) std::cout << "The mc particle pdg id " << pa->pdgId() << " with energy " << pa->energy() << std::endl;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, parameters));  
        
    // Create parent-daughter relationships
    // HepMC::GenParticle * theParticle = new HepMC::GenParticle(HepMC::FourVector(0.,0.,0.,0.),12,1);
     size_t n = pa->numberOfDaughters();
    //std::cout << "The mc particle pdg id " << pa->pdgId() << " with energy " << pa->energy() << " and " << n << " daughters " <<  std::endl;
        
    //Bool_t MCEndpointFound = kFALSE;
    //if (n
        
    for(size_t j = 0; j < n; ++ j) {
      const Candidate * d = pa->daughter( j );
      //if we want to keep it also in GenParticle uncomment here
      const GenParticle * da = NULL;
      //We need to check if this daughter has an integer charge
      bool integercharge = ( ( (int) d->charge() ) - (d->charge()) ) == 0 ? true : false;
      da = new GenParticle( d->charge(), d->p4() , d->vertex() , d->pdgId() , d->status() , integercharge);
        
        
    double RaVeDa = 10 * std::sqrt(da->vx()*da->vx()+da->vy()*da->vy()+da->vz()*da->vz());
         
    if (i<2) {
       if (RminVtxDaughter[i]>RaVeDa)
          RminVtxDaughter[i] = RaVeDa;
       if (ZminVtxDaughter[i]>da->vz())
          ZminVtxDaughter[i] = da->vz();
    }
  
  
  
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*m_pPandora, pa , da));
  
    }
  
   
  }
    
  if (ZminVtxDaughter[0] < 3170) isDecayedBeforeCalo[0] = 1;
  if (ZminVtxDaughter[1] < 3170) isDecayedBeforeCalo[1] = 1;
    
}

// void runPandora::preparePFO(const edm::Event& iEvent, const edm::EventSetup& iSetup){
void runPandora::preparePFO(const edm::Event& iEvent, const edm::EventSetup& iSetup,edm::Handle<HGCRecHitCollection> 
			    HGCeeRecHitHandle,
			    edm::Handle<HGCRecHitCollection> HGChefRecHitHandle,
			    edm::Handle<HGCRecHitCollection> HGChebRecHitHandle,
			    reco::Vertex& pv){

    // PandoraPFANew/v00-09/include/Pandora/PandoraInternal.h
    // typedef std::set<ParticleFlowObject *> PfoList;  
    //     PandoraPFANew/v00-09/include/Api/PandoraContentApi.h
        //    class ParticleFlowObject
        //    {
        //    public:
        //        /**
        //         *  @brief  Parameters class
        //         */
        //        class Parameters
        //        {
        //        public:
        //            pandora::InputInt               m_particleId;       ///< The particle flow object id (PDG code)
        //            pandora::InputInt               m_charge;           ///< The particle flow object charge
        //            pandora::InputFloat             m_mass;             ///< The particle flow object mass
        //            pandora::InputFloat             m_energy;           ///< The particle flow object energy
        //            pandora::InputCartesianVector   m_momentum;         ///< The particle flow object momentum
        //            pandora::ClusterList            m_clusterList;      ///< The clusters in the particle flow object
        //            pandora::TrackList              m_trackList;        ///< The tracks in the particle flow object
        //        };
        //        /**
        //         *  @brief  Create a particle flow object
        //         * 
        //         *  @param  algorithm the algorithm creating the particle flow object
        //         *  @param  particleFlowObjectParameters the particle flow object parameters
        //         */
        //        static pandora::StatusCode Create(const pandora::Algorithm &algorithm, const Parameters &parameters);
        //    };
        //    typedef ParticleFlowObject::Parameters ParticleFlowObjectParameters;

    // const pandora::CartesianVector momentum(1., 2., 3.);
    // for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO){
    //   (*itPFO)->SetParticleId();
    //   (*itPFO)->SetCharge();
    //   (*itPFO)->SetMass();
    //   (*itPFO)->SetEnergy();
    //   (*itPFO)->SetMomentum();
    // }
  
  memset(&EvtInfo,0x00,sizeof(EvtInfo));
  memset(&GenInfo,0x00,sizeof(GenInfo));
  memset(&GenJetInfo,0x00,sizeof(GenJetInfo));
  memset(&PFOInfo,0x00,sizeof(PFOInfo));
  memset(&FJOInfo,0x00,sizeof(FJOInfo));
  memset(&GenFJOInfo,0x00,sizeof(GenFJOInfo));
  memset(&HitInfo,0x00,sizeof(HitInfo));

  //  const pandora::PfoList *pPfoList = NULL;
  const pandora::PfoList *pPfoListB4Sorting = NULL;
  pandora::PfoList pPfoList_tmp;

  // PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));

  //  const pandora::StatusCode statusCode(PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));  
  const pandora::StatusCode statusCode(PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoListB4Sorting));  

  //Sort PFO by descending ET-- start --
  int Npfob4sorting = 0;
  std::vector< std::pair<double, int> > ValueByET;
  std::vector<ParticleFlowObject *> pPfoV;
  std::vector<ParticleFlowObject *> pPfoVsorted;
  for (pandora::PfoList::const_iterator itPFO = pPfoListB4Sorting->begin(), itPFOEnd = pPfoListB4Sorting->end(); itPFO != itPFOEnd; ++itPFO){ 
        double pT = 
	  sqrt(((*itPFO)->GetMomentum()).GetX()*((*itPFO)->GetMomentum()).GetX()+((*itPFO)->GetMomentum()).GetY()*((*itPFO)->GetMomentum()).GetY()) ; 

	ValueByET.push_back( std::pair<double, int>(pT, Npfob4sorting) );
	pPfoV.push_back((ParticleFlowObject*)(*itPFO));
	Npfob4sorting++;
  }
  sort(ValueByET.begin(), ValueByET.end());
  reverse(ValueByET.begin(), ValueByET.end());

  for(int inpfo = 0; inpfo < Npfob4sorting; inpfo++)
    pPfoVsorted.push_back((ParticleFlowObject*)(pPfoV.at(ValueByET.at(inpfo).second)));
  const std::vector<ParticleFlowObject *> *pPfoList = new std::vector<ParticleFlowObject *>(pPfoVsorted);
  //Sort PFO by descending ET-- end --

  if (pandora::STATUS_CODE_SUCCESS != statusCode){throw pandora::StatusCodeException(statusCode);}

  edm::Handle<std::vector<reco::GenParticle> > genpart;
  iEvent.getByLabel(inputTagGenParticles_,genpart); //NS ADD
  std::cout << " GENPART SIZE IS " << genpart->size() << std::endl;

  //int found_recoIndx = -1;
  Double_t ene_all_true = 0;

  // Fast jet for genetic particles
  std::vector<fastjet::PseudoJet> genfjInputs;
  genfjInputs.clear();

  for(size_t i = 0; i < genpart->size(); ++ i) { // 1 (gen-particle loop)

    const GenParticle * pa = &(*genpart)[i];
    GenInfo.ene[GenInfo.Size]    = pa->energy();

    if (pa->numberOfDaughters()==0)  ene_all_true = ene_all_true+GenInfo.ene[GenInfo.Size];

    GenInfo.charge[GenInfo.Size] = 0;
    GenInfo.pT[GenInfo.Size]     = sqrt(pa->px()*pa->px()+pa->py()*pa->py());
    GenInfo.pZ[GenInfo.Size]     = pa->pz();
    GenInfo.eta[GenInfo.Size]    = 0.5 * log((sqrt(GenInfo.pT[GenInfo.Size]*GenInfo.pT[GenInfo.Size]+GenInfo.pZ[GenInfo.Size]*GenInfo.pZ[GenInfo.Size])+GenInfo.pZ[GenInfo.Size])/(sqrt(GenInfo.pT[GenInfo.Size]*GenInfo.pT[GenInfo.Size]+GenInfo.pZ[GenInfo.Size]*GenInfo.pZ[GenInfo.Size])-GenInfo.pZ[GenInfo.Size]));
    TLorentzVector GenVector;
    GenVector.SetPxPyPzE(pa->px(), pa->py(), pa->pz(), GenInfo.ene[GenInfo.Size]);
    GenInfo.phi[GenInfo.Size]  = GenVector.Phi();       
    GenInfo.eta_check[GenInfo.Size]   = pa->eta();
    GenInfo.phi_check[GenInfo.Size]   = pa->phi();

    GenInfo.pid[GenInfo.Size]    = pa->pdgId();
    GenInfo.mass[GenInfo.Size]   = pa->mass();
    GenInfo.status[GenInfo.Size]   = pa->status();


    if(pa->pdgId()>0) GenInfo.charge[GenInfo.Size] = 1;
    if(pa->pdgId()<0) GenInfo.charge[GenInfo.Size] = -1;
    if(pa->pdgId()==22) GenInfo.charge[GenInfo.Size] = 0;
    GenInfo.charge[GenInfo.Size] = pa->charge();

    if(pa->status()==1)// stable particle (handed over to Geant4 for detector simulation)
      genfjInputs.push_back(fastjet::PseudoJet(
					       pa->px(), 
					       pa->py(),
					       pa->pz(),
					       GenInfo.ene[GenInfo.Size] 
					       ));
    GenInfo.Size++;
  } // 1  (gen-particle loop)
  std::cout << " ENERGY ALL TRUE " << ene_all_true << std::endl;

  edm::Handle<std::vector<reco::GenJet> > genjet;
  iEvent.getByLabel("ak5GenJets",genjet); //Jacky ADD
  std::cout << " GENJET SIZE IS " << genjet->size() << std::endl;
  for(size_t i = 0; i < genjet->size(); ++ i) { // 1 (gen-particle loop)

    const reco::GenJet * pa = &(*genjet)[i];

    TLorentzVector GenVector;
    GenVector.SetPtEtaPhiM(pa->pt(), pa->eta(), pa->phi(), pa->mass());

    //GenJetInfo.ene[GenJetInfo.Size]    = GenVector.Energy();
    GenJetInfo.ene[GenJetInfo.Size]    = pa->energy();
    GenJetInfo.pT[GenJetInfo.Size]     = pa->pt();
    GenJetInfo.pZ[GenJetInfo.Size]     = GenVector.Pz();
    GenJetInfo.eta[GenJetInfo.Size]    = pa->eta();
    GenJetInfo.phi[GenJetInfo.Size]    = pa->phi();       

    GenJetInfo.pid[GenJetInfo.Size]    = pa->pdgId();
    GenJetInfo.mass[GenJetInfo.Size]   = pa->mass();

    // get all constituents
    std::vector<const GenParticle*> theJetConstituents = pa->getGenConstituents();
    GenJetInfo.Nconstituents[GenJetInfo.Size]   = theJetConstituents.size();

    double EMene = 0.;
    double ChargedHadene = 0.;
    double NeutralHadene = 0.;
    for( std::vector<const GenParticle*>::const_iterator aCandidate = theJetConstituents.begin();
	 aCandidate != theJetConstituents.end();
	 ++aCandidate)
      {
	int pdgId = std::abs((*aCandidate)->pdgId());
	const Candidate* theCandidate = static_cast<const Candidate*>(*aCandidate);
	//filter nus
	if ((pdgId>=11 && pdgId<=15) || (pdgId==22))// EM particle
	  {
	    EMene += theCandidate->energy(); 
	  } else 
	  {
	    if(theCandidate->charge()==0){
	      NeutralHadene += theCandidate->energy(); 
	    }else{
	      ChargedHadene += theCandidate->energy(); 
	    }
	  }
      }
    GenJetInfo.EMF[GenJetInfo.Size] = (EMene)/pa->energy();
    GenJetInfo.CHF[GenJetInfo.Size] = (ChargedHadene)/pa->energy();
    GenJetInfo.NHF[GenJetInfo.Size] = (NeutralHadene)/pa->energy();
    //GenJetInfo.charge[GenJetInfo.Size] = pa->charge();
    GenJetInfo.Size++;
  } // 1  (gen-particle loop)

  edm::Handle<std::vector<reco::GenJet> > genjet;
  iEvent.getByLabel("ak5GenJets",genjet); //Jacky ADD
  std::cout << " GENJET SIZE IS " << genjet->size() << std::endl;
  for(size_t i = 0; i < genjet->size(); ++ i) { // 1 (gen-particle loop)

    const reco::GenJet * pa = &(*genjet)[i];

    TLorentzVector GenVector;
    GenVector.SetPtEtaPhiM(pa->pt(), pa->eta(), pa->phi(), pa->mass());

    //GenJetInfo.ene[GenJetInfo.Size]    = GenVector.Energy();
    GenJetInfo.ene[GenJetInfo.Size]    = pa->energy();
    GenJetInfo.pT[GenJetInfo.Size]     = pa->pt();
    GenJetInfo.pZ[GenJetInfo.Size]     = GenVector.Pz();
    GenJetInfo.eta[GenJetInfo.Size]    = pa->eta();
    GenJetInfo.phi[GenJetInfo.Size]    = pa->phi();       

    GenJetInfo.pid[GenJetInfo.Size]    = pa->pdgId();
    GenJetInfo.mass[GenJetInfo.Size]   = pa->mass();

    // get all constituents
    std::vector<const GenParticle*> theJetConstituents = pa->getGenConstituents();
    GenJetInfo.Nconstituents[GenJetInfo.Size]   = theJetConstituents.size();

    double EMene = 0.;
    double ChargedHadene = 0.;
    double NeutralHadene = 0.;
    for( std::vector<const GenParticle*>::const_iterator aCandidate = theJetConstituents.begin();
	 aCandidate != theJetConstituents.end();
	 ++aCandidate)
      {
	int pdgId = std::abs((*aCandidate)->pdgId());
	const Candidate* theCandidate = static_cast<const Candidate*>(*aCandidate);
	//filter nus
	if ((pdgId>=11 && pdgId<=15) || (pdgId==22))// EM particle
	  {
	    EMene += theCandidate->energy(); 
	  } else 
	  {
	    if(theCandidate->charge()==0){
	      NeutralHadene += theCandidate->energy(); 
	    }else{
	      ChargedHadene += theCandidate->energy(); 
	    }
	  }
      }
    GenJetInfo.EMF[GenJetInfo.Size] = (EMene)/pa->energy();
    GenJetInfo.CHF[GenJetInfo.Size] = (ChargedHadene)/pa->energy();
    GenJetInfo.NHF[GenJetInfo.Size] = (NeutralHadene)/pa->energy();
    //GenJetInfo.charge[GenJetInfo.Size] = pa->charge();
    GenJetInfo.Size++;
  } // 1  (gen-particle loop)


  edm::Handle<std::vector<reco::GenParticle> > genpart;
  iEvent.getByLabel(inputTagGenParticles_,genpart); //NS ADD
  std::cout << " GENPART SIZE IS " << genpart->size() << std::endl;

  Double_t found_energy = -1;
  Double_t ene_all_true = 0;

  for(size_t i = 0; i < genpart->size(); ++ i) { // 1 
          
    const GenParticle * pa = &(*genpart)[i];
    ene_true    = pa->energy();
     
    if (pa->numberOfDaughters()==0)  ene_all_true = ene_all_true+ene_true;

    charge_true = 0;
    pT_true     = sqrt(pa->px()*pa->px()+pa->py()*pa->py());
    pid_true    = pa->pdgId();
    mass_true   = pa->mass();

    isDecBefCal = isDecayedBeforeCalo[i];

    if(pa->pdgId()>0) charge_true = 1;
    if(pa->pdgId()<0) charge_true = -1;
    if(pa->pdgId()==22) charge_true = 0;
     //std::cout << " IN GENPART LOOP charge min " << std::endl;
      double diff_min  = 1e9;
      ene_match       = 0;
      mass_match      = 0;
      pid_match       = 0;
      pT_match        = 0;
      charge_match    = 0;
      ene_match_em    = 0;
      ene_match_had   = 0;
      ene_match_track = 0;
 
      nbPFOs = 0;
      for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO){ // 4
        nbPFOs=nbPFOs+1;
        double charge = (*itPFO)->GetCharge() ;
        double energy = (*itPFO)->GetEnergy();
        double pid    = (*itPFO)->GetParticleId();
 
        //FIXME
        //NOTYET TLorentzVector Pfo4vec(0.,0.,0.,0.);
        //NOTYET if (m_energyCorrMethod == "WEIGHTING") {
        //NOTYET    if (pid == 22)
        //NOTYET       energy += 2.095/18148.82;
        //NOTYET    Pfo4vec.SetPxPyPzE(((*itPFO)->GetMomentum()).GetX(),
        //NOTYET          ((*itPFO)->GetMomentum()).GetY(),
        //NOTYET          ((*itPFO)->GetMomentum()).GetZ(),
        //NOTYET          energy);
        //NOTYET }

        double mass   = (*itPFO)->GetMass() ;
        double pT     = sqrt(((*itPFO)->GetMomentum()).GetX()*((*itPFO)->GetMomentum()).GetX()+((*itPFO)->GetMomentum()).GetY()*((*itPFO)->GetMomentum()).GetY()) ; // EDW GetX

        const ClusterList &clusterList((*itPFO)->GetClusterList());
        //std::cout << " size of cluster list " << clusterList.size() << std::endl;
        ClusterVector clusterVector(clusterList.begin(), clusterList.end());
        ene_em  =0;
        ene_had =0;
 
        unsigned int firstLayer = std::numeric_limits<unsigned int>::max();
        unsigned int lastLayer = 0;

    for (ClusterVector::const_iterator clusterIter = clusterVector.begin(), clusterIterEnd = clusterVector.end();
        clusterIter != clusterIterEnd; ++clusterIter)
    {
        Cluster *pCluster = (*clusterIter);
//        ene_em  = pCluster->GetElectromagneticEnergy();
//        ene_had = pCluster->GetHadronicEnergy();
//     std::cout << " clusterEMenergy INI " << pCluster->GetElectromagneticEnergy()  << std::endl;  
// hits
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList pCaloHitList;
        orderedCaloHitList.GetCaloHitList(pCaloHitList);
      
       for (CaloHitList::const_iterator hitIter = pCaloHitList.begin(), hitIterEnd = pCaloHitList.end(); hitIter != hitIterEnd; ++hitIter)
        {
        const CaloHit *pCaloHit = (*hitIter);
        // Determing extremal pseudolayers
        const unsigned int pseudoLayer(pCaloHit->GetPseudoLayer());
        if (pseudoLayer > lastLayer)
            lastLayer = pseudoLayer;
        if (pseudoLayer < firstLayer)
            firstLayer = pseudoLayer;
        }
        
     }
      

    
        //Hieu, 22.08.14
        //FIXME: for the moment this is a fast fix, should be merged with similar part below
        //       or simply remove that part (just for a few control histograms)
        double clusterEMenergy  = 0.;
        double clusterHADenergy = 0.;
                 
        const pandora::ClusterAddressList clusterAddressList((*itPFO)->GetClusterAddressList());
        for (pandora::ClusterAddressList::const_iterator itCluster = clusterAddressList.begin(), itClusterEnd = clusterAddressList.end();
            itCluster != itClusterEnd; ++itCluster)
        {
           const unsigned int nHitsInCluster((*itCluster).size());
              
           for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
           {
              const HGCRecHit *hgcHit = (HGCRecHit*)((*itCluster)[iHit]);
              //const HGCEEDetId& detidEE = hgcHit->id();
              double nCellInteractionLengths = 0.;
              double nCellRadiationLengths = 0.;
              double absorberCorrectionEM = 1.;
              double absorberCorrectionHAD = 1.;
                  
              const DetId& detid = hgcHit->id();
              if (!detid)
                 continue;
        
              //const HGCHEDetId& detidHE = hgcHit->id();
           
              ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
              //ForwardSubdetector thesubdetHE = (ForwardSubdetector)detidHE.subdetId();
              if (thesubdet == 3) {
                 int layer = (int) ((HGCEEDetId)(detid)).layer() ;
                 if (m_energyCorrMethod=="ABSCORR") {
                    getLayerPropertiesEE(hgcHit, layer,
                          nCellInteractionLengths, nCellRadiationLengths,
                          absorberCorrectionEM, absorberCorrectionHAD);
                    clusterEMenergy += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * absorberCorrectionEM * m_EM_addCalibrEE;
                    clusterHADenergy += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * absorberCorrectionHAD * m_HAD_addCalibrEE;
                 }
                 else if (m_energyCorrMethod == "WEIGHTING") {
                    clusterEMenergy += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * getEnergyWeight(D_HGCEE, layer, EM) * m_eCalToMipEndCap;
                    clusterHADenergy += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * getEnergyWeight(D_HGCEE, layer,HAD) * m_eCalToMipEndCap;
                 }
              }
              else if (thesubdet == 4) {
                 int layer = (int) ((HGCHEDetId)(detid)).layer() ;
                 if (m_energyCorrMethod=="ABSCORR") {
                    clusterEMenergy += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToEMGeVEndCapHEF  * absorberCorrectionEM  * m_EM_addCalibrHEF ;
                    clusterHADenergy += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToHadGeVEndCapHEF * absorberCorrectionHAD * m_HAD_addCalibrHEF;
                 }
                 else if (m_energyCorrMethod=="WEIGHTING") {
                    clusterEMenergy += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToEMGeVEndCapHEF  * getEnergyWeight(D_HGCHEF, layer, EM) * m_hCalToMipEndCapHEF;
                    clusterHADenergy += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToHadGeVEndCapHEF * getEnergyWeight(D_HGCHEF, layer,HAD) * m_hCalToMipEndCapHEF;
                 }
           
              }
              else if (thesubdet == 5) {
                 int layer = (int) ((HGCHEDetId)(detid)).layer() ;
                 if (m_energyCorrMethod=="ABSCORR") {
                    clusterEMenergy  += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * absorberCorrectionEM  * m_EM_addCalibrHEB ;
                    clusterHADenergy += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * absorberCorrectionHAD * m_HAD_addCalibrHEB;
                 }
                 else if (m_energyCorrMethod=="WEIGHTING") {
                    clusterEMenergy  += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * getEnergyWeight(D_HGCHEB, layer, EM) * m_hCalToMipEndCapHEB;
                    clusterHADenergy += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * getEnergyWeight(D_HGCHEB, layer,HAD) * m_hCalToMipEndCapHEB;
                 }  
              }
              else {
              }
       
           }
        }
        
        ene_em = clusterEMenergy;
        ene_had = clusterHADenergy;

    const TrackList &trackList((*itPFO)->GetTrackList());
        TrackVector trackVector(trackList.begin(), trackList.end());

    ene_track=0; 

     for (TrackVector::const_iterator trackIter = trackVector.begin(), trackIterEnd = trackVector.end();
        trackIter != trackIterEnd; ++trackIter)
    {
         pandora::Track *pPandoraTrack = (*trackIter);
        // Extract pandora track states
        const TrackState &trackState(pPandoraTrack->GetTrackStateAtStart());
        const CartesianVector &momentum(trackState.GetMomentum());
        ene_track = momentum.GetMagnitude();
     }

    
        double diff = 1e10;  
        if(charge_true == charge && energy != found_energy){
          diff = abs(energy-ene_true);
          if(diff<diff_min){ // 3
           diff_min      = diff;
           ene_match     = energy;
           mass_match    = mass;
           charge_match  = charge;
           pT_match      = pT;
           pid_match     = pid;
           ene_match_em  = ene_em;
           ene_match_had = ene_had;
           ene_match_track = ene_track;
           first_layer_match = double(firstLayer*1.);
           last_layer_match  = double(lastLayer*1.);
           if(pid_match==11  && last_layer_match > 30) pid_match = -211;
           if(pid_match==-11 && last_layer_match > 30) pid_match =  211;
          } // 3
        }

      } // 4
        

      if(ene_match>0){ // 2
       runno = iEvent.id().run();
       eventno = iEvent.id().event();
       lumi = iEvent.luminosityBlock();
       found_energy = ene_match;
       mytree->Fill();
       //std::cout << " FOUND MATCH " << ene_match << " to " << ene_true << std::endl;
       //std::cout << " FOUND MATCH " << pT_match << " to " << pT_true << std::endl;
       //std::cout << " FOUND MATCH " << pid_match << " to " << pid_true << std::endl;
       //std::cout << " FOUND MATCH " << mass_match << " to " << mass_true << std::endl;
      } // 2
 } // 1

  std::cout << " ENERGY ALL TRUE " << ene_all_true << std::endl;

  double _sumPFOEnergy(0.);
  double sumClustEMEcalE(0.); //PFO cluster energy in Ecal
  double sumClustEMHcalE(0.); //PFO cluster energy in Hcal
  double sumClustHADEcalE(0.); //PFO cluster energy in Ecal
  double sumClustHADHcalE(0.); //PFO cluster energy in Hcal

  double ene_all = 0;
  int nbPFOs(0);
  for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO){
    nbPFOs++;
    std::cout << "Particle Id: " << (*itPFO)->GetParticleId() << std::endl;
    //std::cout << "Charge: " << (*itPFO)->GetCharge() << std::endl;
    //std::cout << "Mass: " << (*itPFO)->GetMass() << std::endl;
    std::cout << "Energy: " << (*itPFO)->GetEnergy() << std::endl;
    _sumPFOEnergy += (*itPFO)->GetEnergy();

    ene_all = ene_all + (*itPFO)->GetEnergy() ;

    //For the cluster we will deal with it after the finishing of the calo hit
    const pandora::ClusterAddressList clusterAddressList((*itPFO)->GetClusterAddressList());
    const pandora::TrackAddressList trackAddressList((*itPFO)->GetTrackAddressList());
    const pandora::TrackList trackList((*itPFO)->GetTrackList());
   
    //TGClient *gclient  = NULL;
    //PandoraMonitoringApi::VisualizeTracks(  &trackList  , "currentTrackList", AUTO);
    //PANDORA_MONITORING_API(ViewEvent() );

    PANDORA_MONITORING_API(VisualizeTracks(  trackAddressList  , "currentTrackList", AUTO, false, true  ) );    
    PANDORA_MONITORING_API(ViewEvent() );

    for (pandora::ClusterAddressList::const_iterator itCluster = clusterAddressList.begin(), itClusterEnd = clusterAddressList.end(); itCluster != itClusterEnd; ++itCluster){
      const unsigned int nHitsInCluster((*itCluster).size());
      
      //const pandora::CaloHitAddressList &caloHitAddressList(*itCluster);
      //EcalRecHit * hgcrh = NULL;
      //HBHERecHit * hrh = NULL;

      
      //int nbNonEHcalHit = 0;
      for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
      {

         const HGCRecHit *hgcHit = (HGCRecHit*)((*itCluster)[iHit]);
         //const HGCEEDetId& detidEE = hgcHit->id();
         double nCellInteractionLengths = 0.;
         double nCellRadiationLengths = 0.;
         double absorberCorrectionEM = 1.;
         double absorberCorrectionHAD = 1.;

         const DetId& detid = hgcHit->id();
         if (!detid)
            continue;

         //const HGCHEDetId& detidHE = hgcHit->id();

         ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
         //ForwardSubdetector thesubdetHE = (ForwardSubdetector)detidHE.subdetId();
         if (thesubdet == 3) {
            int layer = (int) ((HGCEEDetId)(detid)).layer() ;
            if (m_energyCorrMethod=="ABSCORR") {
               getLayerPropertiesEE(hgcHit, layer,
                     nCellInteractionLengths, nCellRadiationLengths,
                     absorberCorrectionEM, absorberCorrectionHAD);
               sumClustEMEcalE  += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * absorberCorrectionEM * m_EM_addCalibrEE;
               sumClustHADEcalE += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * absorberCorrectionHAD * m_HAD_addCalibrEE;
            }
            else if (m_energyCorrMethod == "WEIGHTING") {
               sumClustEMEcalE  += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * getEnergyWeight(D_HGCEE, layer, EM) * m_eCalToMipEndCap;
               sumClustHADEcalE += hgcHit->energy() * m_Calibr_ADC2GeV_EE * m_eCalToEMGeVEndCap * getEnergyWeight(D_HGCEE, layer,HAD) * m_eCalToMipEndCap;
            }
         }
         else if (thesubdet == 4) {
            int layer = (int) ((HGCHEDetId)(detid)).layer() ;
            if (m_energyCorrMethod=="ABSCORR") {
               sumClustEMHcalE  += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToEMGeVEndCapHEF  * absorberCorrectionEM  * m_EM_addCalibrHEF ;
               sumClustHADHcalE += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToHadGeVEndCapHEF * absorberCorrectionHAD * m_HAD_addCalibrHEF;
            }
            else if (m_energyCorrMethod=="WEIGHTING") {
               sumClustEMHcalE  += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToEMGeVEndCapHEF  * getEnergyWeight(D_HGCHEF, layer, EM) * m_hCalToMipEndCapHEF;
               sumClustHADHcalE += hgcHit->energy() * m_Calibr_ADC2GeV_HEF * m_hCalToHadGeVEndCapHEF * getEnergyWeight(D_HGCHEF, layer,HAD) * m_hCalToMipEndCapHEF;
            }

         }
         else if (thesubdet == 5) {
           int layer = (int) ((HGCHEDetId)(detid)).layer() ;
           if (m_energyCorrMethod=="ABSCORR") {
            sumClustEMHcalE  += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * absorberCorrectionEM  * m_EM_addCalibrHEB ;
            sumClustHADHcalE += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * absorberCorrectionHAD * m_HAD_addCalibrHEB;
           }
            else if (m_energyCorrMethod=="WEIGHTING") {
            sumClustEMHcalE  += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * getEnergyWeight(D_HGCHEB, layer, EM) * m_hCalToMipEndCapHEB;
            sumClustHADHcalE += hgcHit->energy() * m_Calibr_ADC2GeV_HEB * m_hCalToEMGeVEndCapHEB * getEnergyWeight(D_HGCHEB, layer,HAD) * m_hCalToMipEndCapHEB;
            }
         }
         else {
         }

      }

//      for (pandora::CaloHitAddressList::const_iterator hIter = caloHitAddressList.begin(), hIterEnd = caloHitAddressList.end(); hIter != hIterEnd; ++hIter){    
//    pandora::CaloHit * ch = (pandora::CaloHit *) (*hIter);
//    EcalRecHit * hgcrh = NULL;
//    HBHERecHit * hrh = NULL;
//
//    if (ch->GetHitType() ==  pandora::ECAL) { 
//      hgcrh = (EcalRecHit *) (*hIter);
//      std::cout << "EcalRecHit energy " << hgcrh->energy() <<  std::endl;
//     sumClustEcalE += hgcrh->energy();
//    } else if (ch->GetHitType() ==  pandora::HCAL) {  
//      hrh = (HBHERecHit *) (*hIter); 
//      std::cout << "HcalRecHit energy " << hrh->energy() <<  std::endl;          
//     sumClustHcalE += hrh->energy();
//    }
//    else {
//      std::cout << " No ECAL or HCAL??? What is this? " << ch->GetHitType() << std::endl;
//      nbNonEHcalHit++;
//    }
    
//      }
//      std::cout << "nbNonEHcalHit: " << nbNonEHcalHit << std::endl;
 

      // for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit){
      //     EVENT::CalorimeterHit *pCalorimeterHit = (CalorimeterHit*)((*itCluster)[iHit]);
      // }
 
    }


    for (pandora::TrackAddressList::const_iterator itTrack = trackAddressList.begin(), itTrackEnd = trackAddressList.end();itTrack != itTrackEnd; ++itTrack){
      reco::Track * track =  (reco::Track *) (*itTrack);
      std::cout<< "Track from pfo charge " << track->charge() << std::endl;
      std::cout<< "Track from pfo transverse momentum " << track->pt() << std::endl;

    }

  std::cout << " ENERGY  is  " << ene_all << std::endl;

  }

  h2_EM_hcalEecalE->Fill(sumClustEMEcalE,sumClustEMHcalE);
  h2_Had_hcalEecalE->Fill(sumClustHADEcalE,sumClustHADHcalE);
  h_sumPfoE->Fill(_sumPFOEnergy);
  h_nbPFOs->Fill(nbPFOs);

    
}

//Get the track siblings
TrackingParticleRefVector runPandora::getTpSiblings(TrackingParticleRef tp){

  if (tp.isNonnull() && tp->parentVertex().isNonnull() && !tp->parentVertex()->daughterTracks().empty()) {
    return tp->parentVertex()->daughterTracks();
  } else {
    return TrackingParticleRefVector();
  }

}
//Get the track daughters
TrackingParticleRefVector runPandora::getTpDaughters(TrackingParticleRef tp){

  TrackingVertexRefVector trvertexes;
  TrackingParticleRefVector trdaughter;

  if (tp.isNonnull() && tp->decayVertices().isNonnull() ) {
    trvertexes = tp->decayVertices();
    //Loop on vector of TrackingVertex objects where the TrackingParticle decays. 
    for(TrackingVertexRefVector::iterator vi = trvertexes.begin(); vi != trvertexes.end(); vi++){
      //loop on all daughter tracks 
      for(TrackingParticleRefVector::iterator di = (**vi).daughterTracks().begin(); di != (**vi).daughterTracks().end(); di++){
    //Check if the daughter is the same as our mother tp particle
    if ( (*di) == tp  ) {continue;}
    trdaughter.push_back( (*di) );
      }//end on loop over daughter
    }//end on loop over vertices
    return trdaughter;
  } else {
    return TrackingParticleRefVector();
  }

}



// ------------ method called once each job just before starting event loop  ------------
void runPandora::beginJob()
{   
  std::cout << "I am beginning my job...LOOK!!!" << std::endl ; 
  firstEvent_ = true ;


// AP
  const char *pDisplay(::getenv("DISPLAY"));
  if (NULL == pDisplay) {
    std::cout << "DISPLAY environment not set" << std::endl;
  }  else {
    std::cout << "DISPLAY environment set to " << pDisplay << std::endl;
  }
  int argc = 0;
  char* argv = (char *)"";
  TApplication *m_pApplication;
  m_pApplication = gROOT->GetApplication();
  std::cout << "In runPandora::beginJob gVirtualX->GetDisplay()" << gVirtualX->GetDisplay() << std::endl;
  if(!m_pApplication){
    std::cout << "In if of m_pApplication in runPandora::beginJob " << std::endl;
    m_pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
  } 
// END AP

  file = new TFile(_outputFileName.c_str(),"recreate");

  const bool oldAddDir = TH1::AddDirectoryStatus(); 
  TH1::AddDirectory(true); 

  h_sumCaloE = new TH1F("sumCaloE","sum hit E in Calos",1000,0,400);
  h_sumCaloEM = new TH1F("sumCaloEM","sum hit E in Calos",1000,0,400);
  h_sumCaloHad = new TH1F("sumCaloHad","sum hit E in Calos",1000,0,400);
  h_simDir_sumCaloEM  = new TH1F("sumCaloEMsimDir","sum hit E in Calos",1000,0,400);
  h_simDir_sumCaloHad = new TH1F("sumCaloHadsimDir","sum hit E in Calos",1000,0,400);


  h_sumEcalEEM   = new TH1F("sumEcalEEM","sum hit EM E in Ecal",1000,0,350);
  h_sumHcalEEM   = new TH1F("sumHcalEEM","sum hit EM E in Hcal",1000,0,350);
  h_sumEcalEHad  = new TH1F("sumEcalEHad","sum hit Had E in Ecal",1000,0,400);
  h_sumEcalEHad_unc = new TH1F("sumEcalEHad_unc","sum hit Had E in Ecal",1000,0,400);
  h_sumHcalEHad  = new TH1F("sumHcalEHad","sum hit Had E in Hcal",1000,0,400);
  h_sumHcalEHadc = new TH1F("sumHcalEHadc","sum hit Had E in Hcal",1000,0,400);
  h_sumEcalEHadc = new TH1F("sumEcalEHadc","sum hit Had E in Hcal",1000,0,400);
  h_sumEHad      = new TH1F("sumEHad","sum hit Had E ",1000,0,400);

  h_sumPfoE = new TH1F("hsumPfoE","sumPFOenergy",1000,0.,1000.);
  h_nbPFOs = new TH1F("hnbPfos","nb of rec PFOs",30,0.,30.);

  h2_Calo_EM_hcalEecalE = new TH2F("CalohcalEecalEem","",1000,0,400,1000,0,400);
  h2_Calo_Had_hcalEecalE = new TH2F("CalohcalEecalEhad","",1000,0,400,1000,0,400);


  h2_EM_hcalEecalE = new TH2F("hcalEecalEem","",1000,0,400,1000,0,400);
  h2_Had_hcalEecalE = new TH2F("hcalEecalEhad","",1000,0,400,400,0,400);

  h_MCp_Eta = new TH1F("MCp_Eta","MCp_Eta",300,-3.5,3.5);
  h_MCp_Phi = new TH1F("MCp_Phi","MCp_Phi",300,-3.5,3.5);
  h_hit_Eta = new TH1F("hit_Eta","hit_Eta",300,-3.5,3.5);
  h_hit_Phi = new TH1F("hit_Phi","hit_Phi",300,-3.5,3.5);

  h_hitEperLayer_EM  = new TH1F("hitEperLayer_EM" ,"sum hit E per layer",100,-0.5,99.5);
  h_hitEperLayer_HAD = new TH1F("hitEperLayer_HAD","sum hit E per layer",100,-0.5,99.5);
  h_hitEperLayer_EM_EE  = new TH1F("hitEperLayer_EM_EE"  ,"sum hit E per layer EE",100,-0.5,99.5);
  h_hitEperLayer_EM_HEF = new TH1F("hitEperLayer_EM_HEF" ,"sum hit E per layer HEF",100,-0.5,99.5);
  h_hitEperLayer_EM_HEB = new TH1F("hitEperLayer_EM_HEB" ,"sum hit E per layer HEB",100,-0.5,99.5);

  h_MIP_EE  = new TH1F("MIPEE" ,"Mip in EE ",1000,0,10);
  h_MIP_HEF = new TH1F("MIPHEF","Mip in HEF",1000,0,10);
  h_MIP_HEB = new TH1F("MIPHEB","Mip in HEB",1000,0,10);

  h_MIP_Corr_EE  = new TH1F("MIPCorrEE" ,"Mip corrected in EE ",1000,0,10);
  h_MIP_Corr_HEF = new TH1F("MIPCorrHEF","Mip corrected in HEF",1000,0,10);
  h_MIP_Corr_HEB = new TH1F("MIPCorrHEB","Mip corrected in HEB",1000,0,10);

  mytree = new TTree("mytree","mytree");
  mytree->Branch("ene_true",&ene_true);
  mytree->Branch("mass_true",&mass_true);
  mytree->Branch("pT_true",&pT_true);
  mytree->Branch("charge_true",&charge_true);
  mytree->Branch("pid_true",&pid_true);
  mytree->Branch("ene_match",&ene_match);
  mytree->Branch("ene_match_em",&ene_match_em);
  mytree->Branch("ene_match_had",&ene_match_had);
  mytree->Branch("ene_match_track",&ene_match_track);
  mytree->Branch("mass_match",&mass_match);
  mytree->Branch("pT_match",&pT_match);
  mytree->Branch("charge_match",&charge_match);
  mytree->Branch("pid_match",&pid_match);
  mytree->Branch("isDecBefCal",&isDecBefCal);
  mytree->Branch("first_layer_match",&first_layer_match);
  mytree->Branch("last_layer_match",&last_layer_match);
  mytree->Branch("runno",&runno);
  mytree->Branch("eventno",&eventno);
  mytree->Branch("lumi",&lumi);
  mytree->Branch("nbPFOs",&nbPFOs);

  TH1::AddDirectory(oldAddDir); 

  // read in calibration parameters
  initPandoraCalibrParameters();
  readCalibrParameterFile();
  if (m_energyCorrMethod == "WEIGHTING")
     readEnergyWeight();

  m_hitEperLayer_EM  = new double[100];
  m_hitEperLayer_HAD = new double[100];

  m_hitEperLayer_EM_EE  = new double[100];
  m_hitEperLayer_EM_HEF = new double[100];
  m_hitEperLayer_EM_HEB = new double[100];

}


void runPandora::getLayerPropertiesEE (const HGCRecHit *eerh, int layer,
      double & nCellInteractionLengths, double & nCellRadiationLengths,
      double & absorberCorrectionEM, double & absorberCorrectionHAD
      )
{
   double absorberThickness_Pb[3] = { 1.6, 3.3, 5.6}; //mm
   double absorberThickness_Cu[3] = { 3. , 3. , 3. }; //mm

   double absorberRadLength_Pb   = 0.561 * 10; //mm
   double absorberRadLength_Cu   = 1.436 * 10; //mm

   double absorberIntLength_Pb = 17.59 * 10; //mm
   double absorberIntLength_Cu = 15.32 * 10; //mm

   double firstLayerCu = 3.0; //mm
   double firstLayerX0   = firstLayerCu / absorberRadLength_Cu;
   double firstLayerLdaI = firstLayerCu / absorberIntLength_Cu;

   int layerSet = 0;
   if ( 2 <= layer && layer <= 11 ) layerSet = 0;
   else if ( 12<= layer && layer <= 21 ) layerSet = 1;
   else if ( 22<= layer && layer <= 31 ) layerSet = 2;
   else {
      //std::cout << "OOUPSSSSSSSSSSSSSS, unknown layer!!!!" << std::endl;
      return;
   } //if layer == 1 then take predefined values (should be 1 & 0)

   double thisLayerX0 = absorberThickness_Pb[layerSet] / absorberRadLength_Pb
      + absorberThickness_Cu[layerSet] / absorberRadLength_Cu;
   double thisLayerLdaI = absorberThickness_Pb[layerSet] / absorberIntLength_Pb
      + absorberThickness_Cu[layerSet] / absorberIntLength_Cu;

   absorberCorrectionEM  = thisLayerX0   / firstLayerX0  ;
   absorberCorrectionHAD = thisLayerLdaI / firstLayerLdaI;

   nCellInteractionLengths = thisLayerLdaI;
   nCellRadiationLengths   = thisLayerX0  ;

   return;
}
void runPandora::getLayerPropertiesHEF (const HGCRecHit *eerh, int layer,
      double & nCellInteractionLengths, double & nCellRadiationLengths,
      double & absorberCorrectionEM, double & absorberCorrectionHAD
      )
{
   double absorberRadLength_Cu = 1.436 * 10; //mm
   double absorberIntLength_Cu = 15.32 * 10; //mm

   double absorberRadLength_Br = 0.561 * 10; //lead instead, Valeri's mail, 16/7/2014
   double absorberIntLength_Br = 17.59 * 10; //FIXME v5; v4 of geometry: may change to Brass

   double absorberThickness_Cu = 3.  ; //mm
   double absorberThickness_Br = 52. ; //mm

   double layerX0 = absorberThickness_Cu / absorberRadLength_Cu
      + absorberThickness_Br / absorberRadLength_Br;
   double layerLdaI = absorberThickness_Cu / absorberIntLength_Cu
      + absorberThickness_Br / absorberIntLength_Br;


   nCellInteractionLengths = layerLdaI;
   nCellRadiationLengths   = layerX0  ;

   return;
}
void runPandora::getLayerPropertiesHEB (const HGCRecHit *eerh, int layer,
      double & nCellInteractionLengths, double & nCellRadiationLengths,
      double & absorberCorrectionEM, double & absorberCorrectionHAD
      )
{
   double absorberThickness_Br = 34.5 + 9; //mm

   double absorberRadLength_Br = 0.561 * 10;
   double absorberIntLength_Br = 17.59 * 10;

   double layerX0   = absorberThickness_Br / absorberRadLength_Br;
   double layerLdaI = absorberThickness_Br / absorberIntLength_Br;

   nCellInteractionLengths = layerLdaI;
   nCellRadiationLengths   = layerX0  ;

   return;
}

/*****************************************************************************/ 
double runPandora::getEnergyWeight(int subDet, int layer, int showerType)
{
   double weight = 1.;

   //std::cout << "runPandora::getEnergyWeight: "
   //   << ", subDet: " << subDet << ", layer: " << layer
   //   << "showerType: " << showerType
   //   << std::endl;

   if (showerType == EM)
      switch (subDet) {
         case D_HGCEE :
            weight = stm->getCorrectionAtPoint(layer+1,"layerSet_EE","energyWeight_EM_EE");
            break;
         case D_HGCHEF:
            weight = stm->getCorrectionAtPoint(layer+1,"layerSet_HEF","energyWeight_EM_HEF");
            break;
         case D_HGCHEB:
            weight = stm->getCorrectionAtPoint(layer+1,"layerSet_HEB","energyWeight_EM_HEB");
            break;
         default :
            break;
      }
   else 
      switch (subDet) {
         case D_HGCEE :
            weight = stm->getCorrectionAtPoint(layer+1,"layerSet_EE","energyWeight_Had_EE");
            break;
         case D_HGCHEF:
            weight = stm->getCorrectionAtPoint(layer+1,"layerSet_HEF","energyWeight_Had_HEF");
            break;
         case D_HGCHEB:
            weight = stm->getCorrectionAtPoint(layer+1,"layerSet_HEB","energyWeight_Had_HEB");
            break;
         default :
            break;
      }

   //std::cout << "weight: " << weight << std::endl;
  
   return weight;
}


// ------------ method called once each job just after ending the event loop  ------------
void 
runPandora::endJob() 
{
  file->Write();
  file->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
runPandora::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
runPandora::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
runPandora::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
runPandora::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
runPandora::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(runPandora);

