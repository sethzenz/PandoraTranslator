// -*- C++ -*-
//
// Package:    PFCal/runPandora
// Class:      runPandora
// 
/**\class runPandora runPandora.cc PFCal/runPandora/plugins/runPandora.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Andreas Psallidas
//         Created:  Mon, 11 Nov 2013 15:11:14 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <map>
#include <unordered_map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Pandora/Pandora.h"
#include "Pandora/StatusCodes.h"
#include "Api/PandoraApi.h"
#include "TLorentzVector.h"
#include "Objects/ParticleFlowObject.h"
#include "Pandora/PandoraInputTypes.h"
#include "Pandora/PandoraInternal.h"
#include "Objects/Cluster.h"
#include "Objects/Track.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"

#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/Vector3D.h"
//DQM services for histogram
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "HGCal/PandoraTranslator/interface/steerManager.h"

#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include "TH1F.h"

//
// class declaration
//
class CalibCalo;
class CaloGeometryRecord;
class IdealGeometryRecord;
class CaloGeometry;
class HGCalGeometry;
class MagneticField;

// namespace pandora {class Pandora;}

class PandoraCMSPFCandProducer : public edm::EDProducer {
public:
  //enum for subdet types
  enum subdet { EB = 1, HB = 2, EE = 3, HEF = 4, HEB = 5 };

  explicit PandoraCMSPFCandProducer(const edm::ParameterSet&);
  ~PandoraCMSPFCandProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static pandora::Pandora        *m_pPandora;

  void prepareTrack(edm::Event& iEvent);
  void prepareHits(edm::Event& iEvent);
  void preparemcParticle(edm::Event& iEvent);
  void ProcessRecHits(edm::Handle<reco::PFRecHitCollection> PFRecHitHandle, int subdet, const CaloSubdetectorGeometry* geom, CalibCalo* calib, int& nCaloHits, int& nNotFound, reco::Vertex& pv,
                     const pandora::HitType hitType, const pandora::HitRegion hitRegion, PandoraApi::RectangularCaloHitParameters& caloHitParameters);

  
  void preparePFO(edm::Event& iEvent);
  void prepareGeometry();
  void SetDefaultSubDetectorParameters(const std::string &subDetectorName, const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const;
  void CalculateCornerSubDetectorParameters(const CaloSubdetectorGeometry* geom,  const std::vector<DetId>& cells, const pandora::SubDetectorType subDetectorType, 
                                            double& min_innerRadius, double& max_outerRadius, double& min_innerZ, double& max_outerZ,
                                            bool doLayers, std::vector<double>& min_innerR_depth, std::vector<double>& min_innerZ_depth) const;
  void SetCornerSubDetectorParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, 
                                      const double& min_innerRadius, const double& max_outerRadius, const double& min_innerZ, const double& max_outerZ) const;
  void SetSingleLayerParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, PandoraApi::Geometry::LayerParameters &layerParameters) const;
  void SetMultiLayerParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, std::vector<PandoraApi::Geometry::LayerParameters*> &layerParameters,
                               std::vector<double>& min_innerR_depth, std::vector<double>& min_innerZ_depth, const unsigned int& nTotalLayers, int& nLayers) const;
                                         
  TrackingParticleRefVector getTpSiblings(TrackingParticleRef tp);
  TrackingParticleRefVector getTpDaughters(TrackingParticleRef tp);

  std::string _outputFileName;
  edm::FileInPath m_pandoraSettingsXmlFile;

  edm::FileInPath m_calibrationParameterFile;
  edm::FileInPath m_energyWeightingFilename;

  std::string m_energyCorrMethod; //energy correction method


  void initPandoraCalibrParameters();
  void readCalibrParameterFile();
  void readEnergyWeight(); //FIXME part of calibration, to be merged to readCalibrParameterFile??

private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;

  virtual void resetVariables() {
     for (int il=0; il<100; il++) {
        m_hitEperLayer_EM[subdet::EE][il] = 0.;
        m_hitEperLayer_EM[subdet::HEF][il] = 0.;
        m_hitEperLayer_EM[subdet::HEB][il] = 0.;
        m_hitEperLayer_HAD[subdet::EE][il] = 0.;
        m_hitEperLayer_HAD[subdet::HEF][il] = 0.;
        m_hitEperLayer_HAD[subdet::HEB][il] = 0.;    
     }
  };

  steerManager * stm;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  // ----------access to event data
  edm::InputTag    inputTagEcalRecHitsEB_ ; 
  edm::InputTag    inputTagHcalRecHitsHBHE_ ; 
  edm::InputTag    inputTagHGCEErechit_;
  edm::InputTag    inputTagHGCHEFrechit_;
  edm::InputTag    inputTagHGCHEBrechit_;
  edm::InputTag    inputTagtPRecoTrackAsssociation_;
  edm::InputTag    inputTagGenParticles_;
  std::vector<edm::InputTag>  inputTagGeneralTracks_;
  //std::vector<std::string> mFileNames;
  
  // hash tables to translate back to CMSSW collection index from Pandora
  std::unordered_map<const void*,unsigned int> recHitMap;
  
  //geometry handles
  edm::ESHandle<CaloGeometry> geoHandle;
  edm::ESHandle<HGCalGeometry> hgceeGeoHandle ; 
  edm::ESHandle<HGCalGeometry> hgchefGeoHandle ; 
  edm::ESHandle<HGCalGeometry> hgchebGeoHandle ; 
  edm::ESHandle<MagneticField> magneticField;

  TFile * file;
  TTree *mytree;
  double ene_track, ene_match_track,ene_match_em,ene_match_had, ene_had,ene_em,ene_match,mass_match,pid_match,pT_match,charge_match;
  double ene_true,mass_true,pid_true,pT_true,charge_true;
  double first_layer,last_layer,first_layer_match,last_layer_match;
  int runno, eventno, lumi , nbPFOs;

  int isDecBefCal;
  double RminVtxDaughter[2];
  double ZminVtxDaughter[2];
  int isDecayedBeforeCalo[2];

  TH1F * Epfos;
  TH1F * Egenpart;
  TH1F * Energy_res;

  TH1F * h_sumPfoE;
  TH1F * h_sumPfoEEM;
  TH1F * h_sumPfoEHad;
  TH1F * h_nbPFOs;

  TH1F * h_sumCaloE;
  TH1F * h_sumEcalEEM;
  TH1F * h_sumHcalEEM;
  TH1F * h_sumEcalEHad;
  TH1F * h_sumEcalEHad_unc;
  TH1F * h_sumHcalEHad;
  TH1F * h_sumEcalEHadc;
  TH1F * h_sumHcalEHadc;
  TH1F * h_sumEHad;

  TH1F * h_sumCaloEM;
  TH1F * h_sumCaloHad;

  TH1F * h_simDir_sumCaloEM ;//take only hits in sim part. direction
  TH1F * h_simDir_sumCaloHad;//take only hits in sim part. direction

  std::map<subdet,TH1F*> h_MIP ;
  std::map<subdet,TH1F*> h_MIP_Corr ;
  TH1F * h_MCp_Eta;
  TH1F * h_MCp_Phi;
  TH1F * h_hit_Eta;
  TH1F * h_hit_Phi;

  std::map<subdet,TH1F*> h_hitEperLayer_EM;
  std::map<subdet,double*> m_hitEperLayer_EM;
  std::map<subdet,TH1F*> h_hitEperLayer_HAD;
  std::map<subdet,double*> m_hitEperLayer_HAD;

  TH2F * h2_Calo_EM_hcalEecalE;
  TH2F * h2_Calo_Had_hcalEecalE;

  TH2F * h2_EM_hcalEecalE;
  TH2F * h2_Had_hcalEecalE;

  double sumCaloEnergy;
  double sumCaloEnergyEM;
  double sumCaloEnergyHAD;

  double simDir_sumCaloEnergyEM;
  double simDir_sumCaloEnergyHAD;

  double sumCaloECALEnergyEM;
  double sumCaloHCALEnergyEM;
  double sumCaloECALEnergyHAD;
  double sumCaloECALEnergyHAD_unc;
  double sumCaloHCALEnergyHAD;
  
  //-------------- energy weighting ------------------
  unsigned int nHGCeeLayers, nHGChefLayers, nHGChebLayers; 
  CalibCalo *m_calibEB, *m_calibHB, *m_calibEE, *m_calibHEF, *m_calibHEB;
  
  double  offSet_EM;
  double  offSet_Had;

  double m_firstMCpartEta;
  double m_firstMCpartPhi;
  double m_secondMCpartEta;
  double m_secondMCpartPhi;

  double m_muonToMip;

  short _debugLevel;
  double speedoflight;

};

//helper class to store calibration info
class CalibCalo {
  public:
    //constructor
    CalibCalo(PandoraCMSPFCandProducer::subdet id) : m_id(id) {}
    //destructor
    virtual ~CalibCalo() {}
    
    //helper functions
    virtual void initialize() {};
    
    //accessors
    virtual double GetADC2GeV() { return 1.; }
    virtual double GetEMCalib(unsigned int layer) { return m_CalToEMGeV; }
    virtual double GetHADCalib(unsigned int layer) { return m_CalToHADGeV; }
    
    //member variables
    PandoraCMSPFCandProducer::subdet m_id;
    double m_CalThresh;
    double m_CalMipThresh;
    double m_CalToMip;
    double m_CalToEMGeV;
    double m_CalToHADGeV;
    double m_Calibr_ADC2GeV;
    double m_EM_addCalibr;
    double m_HAD_addCalibr;
    std::vector<double> nCellInteractionLengths;
    std::vector<double> nCellRadiationLengths;
};

//extension for HGC
class CalibHGC : public CalibCalo {
  public:
    //constructor
    CalibHGC(PandoraCMSPFCandProducer::subdet id, int totalLayers, std::string energyCorrMethod, steerManager *stm)
            : CalibCalo(id), m_TotalLayers(totalLayers), m_energyCorrMethod(energyCorrMethod), m_stm(stm) {}
    //destructor
    virtual ~CalibHGC() {}
    
    //accessors
    virtual double GetADC2GeV() { return m_Calibr_ADC2GeV; }
    virtual double GetEMCalib(unsigned int layer) { 
      if(layer>m_TotalLayers) return 1.;
      
      if(m_energyCorrMethod == "ABSCORR"){
        return m_CalToEMGeV * m_absorberCorrectionEM[layer] * m_EM_addCalibr;
      }
      else if(m_energyCorrMethod == "WEIGHTING"){
        return m_CalToEMGeV * m_energyWeightEM[layer] * m_CalToMip;
      }
      else return 1.;
    }
    virtual double GetHADCalib(unsigned int layer) {
      if(layer>m_TotalLayers) return 1.;
      
      if(m_energyCorrMethod == "ABSCORR"){
        return m_CalToHADGeV * m_absorberCorrectionHAD[layer] * m_HAD_addCalibr;
      }
      else if(m_energyCorrMethod == "WEIGHTING"){
        return m_CalToHADGeV * m_energyWeightHAD[layer] * m_CalToMip;
      }
      else return 1.;
    }
    
    //helper functions
    virtual void initialize() {
      nCellInteractionLengths.reserve(m_TotalLayers);
      nCellRadiationLengths.reserve(m_TotalLayers);
      m_absorberCorrectionEM.reserve(m_TotalLayers);
      m_absorberCorrectionHAD.reserve(m_TotalLayers);
      m_energyWeightEM.reserve(m_TotalLayers);
      m_energyWeightHAD.reserve(m_TotalLayers);
      getLayerProperties();
    }
    virtual void getLayerProperties() {}
    
    //member variables
    unsigned int m_TotalLayers;
    std::string m_energyCorrMethod;
    steerManager * m_stm;
    std::vector<double> m_absorberCorrectionEM;
    std::vector<double> m_absorberCorrectionHAD;
    std::vector<double> m_energyWeightEM;
    std::vector<double> m_energyWeightHAD;

};

//extension for HGCEE
class CalibHGCEE : public CalibHGC {
  public:
    //constructor
    CalibHGCEE(int totalLayers, std::string energyCorrMethod, steerManager *stm)
              : CalibHGC(PandoraCMSPFCandProducer::subdet::EE, totalLayers, energyCorrMethod, stm) {}
    //destructor
    virtual ~CalibHGCEE() {}
              
    //helper functions
    virtual void getLayerProperties() {
      double absorberThickness_Pb[3] = { 1.6, 3.3, 5.6}; //mm
      double absorberThickness_Cu[3] = { 3. , 3. , 3. }; //mm
      
      double absorberRadLength_Pb   = 0.561 * 10; //mm
      double absorberRadLength_Cu   = 1.436 * 10; //mm
      
      double absorberIntLength_Pb = 17.59 * 10; //mm
      double absorberIntLength_Cu = 15.32 * 10; //mm
      
      double firstLayerCu = 3.0; //mm
      double firstLayerX0   = firstLayerCu / absorberRadLength_Cu;
      double firstLayerLdaI = firstLayerCu / absorberIntLength_Cu;
      
      for(unsigned layer = 0; layer < m_TotalLayers; layer++){
        int layerSet = 0;
        if( layer <= 1){ //if layer == 1 then take predefined values (should be 1 & 0)
          m_absorberCorrectionEM.push_back(1.);
          m_absorberCorrectionHAD.push_back(1.);
        
          nCellInteractionLengths.push_back(0.);
          nCellRadiationLengths.push_back(0.);
        }
        if ( 2 <= layer && layer <= 11 ) layerSet = 0;
        else if ( 12<= layer && layer <= 21 ) layerSet = 1;
        else if ( 22<= layer && layer <= 31 ) layerSet = 2;
        else {
           //std::cout << "OOUPSSSSSSSSSSSSSS, unknown layer!!!!" << std::endl;
           continue;
        } 
        
        double thisLayerX0 = absorberThickness_Pb[layerSet] / absorberRadLength_Pb + absorberThickness_Cu[layerSet] / absorberRadLength_Cu;
        double thisLayerLdaI = absorberThickness_Pb[layerSet] / absorberIntLength_Pb + absorberThickness_Cu[layerSet] / absorberIntLength_Cu;
        
        m_absorberCorrectionEM.push_back(thisLayerX0   / firstLayerX0);
        m_absorberCorrectionHAD.push_back(thisLayerLdaI / firstLayerLdaI);
        
        nCellInteractionLengths.push_back(thisLayerLdaI);
        nCellRadiationLengths.push_back(thisLayerX0);
        
        m_energyWeightEM.push_back(m_stm->getCorrectionAtPoint(layer+1,"layerSet_EE","energyWeight_EM_EE"));
        m_energyWeightHAD.push_back(m_stm->getCorrectionAtPoint(layer+1,"layerSet_EE","energyWeight_HAD_EE"));
      }
    }
};

//extension for HGCHEF
class CalibHGCHEF : public CalibHGC {
  public:
    //constructor
    CalibHGCHEF(int totalLayers, std::string energyCorrMethod, steerManager *stm)
              : CalibHGC(PandoraCMSPFCandProducer::subdet::HEF, totalLayers, energyCorrMethod, stm) {}
    //destructor
    virtual ~CalibHGCHEF() {}
              
    //helper functions
    virtual void getLayerProperties() {   
      double absorberRadLength_Cu = 1.436 * 10; //mm
      double absorberIntLength_Cu = 15.32 * 10; //mm
      
      double absorberRadLength_Br = 0.561 * 10; //lead instead, Valeri's mail, 16/7/2014
      double absorberIntLength_Br = 17.59 * 10; //FIXME v5; v4 of geometry: may change to Brass
      
      double absorberThickness_Cu = 3.  ; //mm
      double absorberThickness_Br = 52. ; //mm
      
      double layerX0 = absorberThickness_Cu / absorberRadLength_Cu + absorberThickness_Br / absorberRadLength_Br;
      double layerLdaI = absorberThickness_Cu / absorberIntLength_Cu + absorberThickness_Br / absorberIntLength_Br;
      
      for(unsigned layer = 0; layer < m_TotalLayers; layer++){     
        m_absorberCorrectionEM.push_back(1.);
        m_absorberCorrectionHAD.push_back(1.);
        
        nCellInteractionLengths.push_back(layerLdaI);
        nCellRadiationLengths.push_back(layerX0);
        
        m_energyWeightEM.push_back(m_stm->getCorrectionAtPoint(layer+1,"layerSet_HEF","energyWeight_EM_HEF"));
        m_energyWeightHAD.push_back(m_stm->getCorrectionAtPoint(layer+1,"layerSet_HEF","energyWeight_HAD_HEF"));
      }
    }
};

//extension for HGCHEB
class CalibHGCHEB : public CalibHGC {
  public:
    //constructor
    CalibHGCHEB(int totalLayers, std::string energyCorrMethod, steerManager *stm)
              : CalibHGC(PandoraCMSPFCandProducer::subdet::HEB, totalLayers, energyCorrMethod, stm) {}
    //destructor
    virtual ~CalibHGCHEB() {}
              
    //helper functions
    virtual void getLayerProperties() {   
      double absorberThickness_Br = 34.5 + 9; //mm
      
      double absorberRadLength_Br = 0.561 * 10;
      double absorberIntLength_Br = 17.59 * 10;
      
      double layerX0   = absorberThickness_Br / absorberRadLength_Br;
      double layerLdaI = absorberThickness_Br / absorberIntLength_Br;
      
      for(unsigned layer = 0; layer < m_TotalLayers; layer++){     
        m_absorberCorrectionEM.push_back(1.);
        m_absorberCorrectionHAD.push_back(1.);
        
        nCellInteractionLengths.push_back(layerLdaI);
        nCellRadiationLengths.push_back(layerX0);
        
        m_energyWeightEM.push_back(m_stm->getCorrectionAtPoint(layer+1,"layerSet_HEB","energyWeight_EM_HEB"));
        m_energyWeightHAD.push_back(m_stm->getCorrectionAtPoint(layer+1,"layerSet_HEB","energyWeight_HAD_HEB"));
      }
    }
};
