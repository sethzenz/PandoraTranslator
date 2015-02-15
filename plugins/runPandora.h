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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"

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
namespace reco {class PFRecHit;}

// namespace pandora {class Pandora;}

class runPandora : public edm::EDAnalyzer {
public:
  explicit runPandora(const edm::ParameterSet&);
  ~runPandora();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static pandora::Pandora        *m_pPandora;

  //  void prepareTrack(math::XYZVector B_,const reco::RecoToSimCollection pRecoToSim,const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void prepareTrack(math::XYZVector B_,const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void prepareHits(edm::Handle<EcalRecHitCollection> ecalRecHitHandleEB,edm::Handle<HBHERecHitCollection> hcalRecHitHandleHBHE,edm::Handle<HGCeeRecHitCollection> HGCeeRecHitHandle,edm::Handle<HGChefRecHitCollection> HGChefRecHitHandle,edm::Handle<HGChebRecHitCollection> HGChebRecHitHandle, reco::Vertex& pv, const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void preparemcParticle(edm::Handle<std::vector<reco::GenParticle> > genpart);

  void preparePFO(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void prepareGeometry(const edm::EventSetup& iSetup);
  void SetDefaultSubDetectorParameters(const std::string &subDetectorName, const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const;
  void CalculateCornerSubDetectorParameters(const CaloSubdetectorGeometry& geom,  const std::vector<DetId>& cells, const pandora::SubDetectorType subDetectorType, 
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
  void getLayerPropertiesEE (const HGCRecHit *eerh, int layer,
        double & nCellInteractionLengths, double & nCellRadiationLengths,
        double & absorberCorrectionEM, double & absorberCorrectionHAD
        );
  void getLayerPropertiesHEF (const HGCRecHit *eerh, int layer,
        double & nCellInteractionLengths, double & nCellRadiationLengths,
        double & absorberCorrectionEM, double & absorberCorrectionHAD
        );
  void getLayerPropertiesHEB (const HGCRecHit *eerh, int layer,
        double & nCellInteractionLengths, double & nCellRadiationLengths,
        double & absorberCorrectionEM, double & absorberCorrectionHAD
        );


  virtual double getEnergyWeight(int subDet, int layer, int showerType);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void resetVariables() {
     for (int il=0; il<100; il++) {
        m_hitEperLayer_EM [il] = 0.;
        m_hitEperLayer_HAD[il] = 0.;

        m_hitEperLayer_EM_EE [il] = 0.;
        m_hitEperLayer_EM_HEF[il] = 0.;
        m_hitEperLayer_EM_HEB[il] = 0.;
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

  TFile * file;
  TTree *mytree;
  double ene_track, ene_match_track,ene_match_em,ene_match_had, ene_had,ene_em,ene_match,mass_match,pid_match,pT_match,charge_match;
  double ene_true,mass_true,pid_true,pT_true,charge_true;
  double first_layer,last_layer,first_layer_match,last_layer_match;
  int runno, eventno, lumi , nbPFOs;


  GenInfoBranches GenInfo;
  GenJetInfoBranches GenJetInfo;
  //double ene_true,mass_true,pid_true,pT_true,charge_true, pZ_true, eta_true, phi_true, eta_check, phi_check;
  double first_layer,last_layer,first_layer_match,last_layer_match;

  EvtInfoBranches EvtInfo;
  //int runno, eventno, lumi , nbPFOs, nbPFOsInSameSide;
  
  HitInfoBranches HitInfo;
  //int hit_size; // hit_size_before9MIP
  //int hit_size_Total;
  //int hit_size_after9MIP;
  //int hit_layer[MAX_NHIT];
  //double hit_EMenergy[MAX_NHIT];
  //double hit_Hadenergy[MAX_NHIT];
  //double hit_Rawenergy[MAX_NHIT];
  //double hit_Eta[MAX_NHIT];
  //double hit_Phi[MAX_NHIT];
  //double hit_X[MAX_NHIT];
  //double hit_Y[MAX_NHIT];
  //double hit_Z[MAX_NHIT];
  //double sumallhitsEMenergyB4Clustering; 
  //double sumallhitsHGCalenergyB4Clustering; 

  //int hit_size_dummy;
  //int hit_layer_dummy[MAX_NHIT];
  //double hit_EMenergy_dummy[MAX_NHIT];
  //double hit_Hadenergy_dummy[MAX_NHIT];
  //double hit_Rawenergy_dummy[MAX_NHIT];

  double sumallhitsEMenergyB4Clustering_[2];   // 0 : -Z; 1 : +Z 
  double sumallhitsHGCalenergyB4Clustering_[2];   // 0 : -Z; 1 : +Z 

  // Hit information for recovery of hits
  double DRmaxMips[MAX_LAYER];
  double EmaxMips[MAX_LAYER];
  double EtamaxMips[MAX_LAYER];
  double PhimaxMips[MAX_LAYER];
  double XmaxMips[MAX_LAYER];
  double YmaxMips[MAX_LAYER];
  double PUforDRmaxMips1[MAX_LAYER];
  double PUforDRmaxMips2[MAX_LAYER];
  double PUforDRmaxMips12[MAX_LAYER];
  double PUforDRmaxMips3[MAX_LAYER];

  TH1F *hDRCut1p0MIP;
  TH1F *hDRCut0p5MIP;

  PFOInfoBranches PFOInfo;
  // PFOs for each gen-particle
  //int PFO_size; 
  //double PFOtrack_dx[MAX_NPFO];
  //double PFOtrack_dy[MAX_NPFO];
  //double PFOtrack_dz[MAX_NPFO];
  //double PFOcharge[MAX_NPFO];
  //int PFOstart[MAX_NPFO];
  //int PFOend[MAX_NPFO];
  //double PFOpid[MAX_NPFO]; 
  //double PFOeta[MAX_NPFO]; 
  //double PFOphi[MAX_NPFO]; 
  //double PFOene[MAX_NPFO]; 
  //double PFOene_track[MAX_NPFO];
  //double PFOene_em_cluster[MAX_NPFO];
  //double PFOene_em_clusterPUsubtraction1[MAX_NPFO];
  //double PFOene_em_clusterPUsubtraction2[MAX_NPFO];
  //double PFOene_em_clusterPUsubtraction12[MAX_NPFO];
  //double PFOene_em_clusterPUsubtraction3[MAX_NPFO];
  //double PFO_PUsubtracted_ene_em_cluster[MAX_NPFO];
  //double PFO_PUsubtracted_ene_em_first3layers[MAX_NPFO];
  //double PFOene_had_cluster[MAX_NPFO];
  //double PFOene_em_allhits[MAX_NPFO]; 
  //double PFOene_em_cluster_woRecovery[MAX_NPFO]; 
  //double PFOene_em_first3layers[MAX_NPFO];
  //int PFOnclusters[MAX_NPFO]; 
  //int PFOnhits[MAX_NPFO]; 
  //int PFOnhitsRecovered[MAX_NPFO]; 
  //// Information for FJO
  //double PFOene_ForAllType[MAX_NPFO];// including energy from EM particle, neutral hadrons, and charged hadrons
  //int PFOene_type[MAX_NPFO];            // -1 : no; 0 : EM particle; 1 : neutral hadrons; 2 : charged hadrons

  FJOInfoBranches FJOInfo;
  GenFJOInfoBranches GenFJOInfo;
  //// Fast Jet object using Pandora particles
  //int FJO_size; 
  //double FJO_ene[MAX_NFJO];
  //double FJO_pT[MAX_NFJO];
  //double FJO_eta[MAX_NFJO];
  //double FJO_phi[MAX_NFJO];
  //int    FJO_Nconstituents[MAX_NFJO];
  //double FJO_ene_true[MAX_NFJO];
  //double FJO_CHF[MAX_NFJO];// Fraction of energy in charged hadrons
  //double FJO_NHF[MAX_NFJO];// Fraction of energy in neutral hadrons
  //double FJO_EMF[MAX_NFJO];// Fraction of energy in EM

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


  TH1F * h_MIP_EE ;
  TH1F * h_MIP_HEF;
  TH1F * h_MIP_HEB;
  TH1F * h_MIP_Corr_EE ;
  TH1F * h_MIP_Corr_HEF;
  TH1F * h_MIP_Corr_HEB;
  TH1F * h_MCp_Eta;
  TH1F * h_MCp_Phi;
  TH1F * h_hit_Eta;
  TH1F * h_hit_Phi;

  TH1F * h_hitEperLayer_EM ;
  TH1F * h_hitEperLayer_HAD;
  double * m_hitEperLayer_EM ;
  double * m_hitEperLayer_HAD;

  TH1F * h_hitEperLayer_EM_EE ;
  TH1F * h_hitEperLayer_EM_HEF;
  TH1F * h_hitEperLayer_EM_HEB;
  double * m_hitEperLayer_EM_EE ;
  double * m_hitEperLayer_EM_HEF;
  double * m_hitEperLayer_EM_HEB;


  TH2F * h2_Calo_EM_hcalEecalE;
  TH2F * h2_Calo_Had_hcalEecalE;

  TH2F * h2_EM_hcalEecalE;
  TH2F * h2_Had_hcalEecalE;



  //-------------- energy weighting ------------------
  double * layerSet_EE ;
  double * layerSet_HEF;
  double * layerSet_HEB;
  double * energyWeight_EM_EE ;
  double * energyWeight_EM_HEF;
  double * energyWeight_EM_HEB;
  double * energyWeight_Had_EE ;
  double * energyWeight_Had_HEF;
  double * energyWeight_Had_HEB;
  double  offSet_EM;
  double  offSet_Had;



  double m_firstMCpartEta;
  double m_firstMCpartPhi;
  double m_secondMCpartEta;
  double m_secondMCpartPhi;


  double m_Calibr_ADC2GeV_EE     ;
  double m_Calibr_ADC2GeV_HEF    ;
  double m_Calibr_ADC2GeV_HEB    ;

  double m_EM_addCalibrEE ;
  double m_EM_addCalibrHEF;
  double m_EM_addCalibrHEB;
  double m_HAD_addCalibrEE ;
  double m_HAD_addCalibrHEF;
  double m_HAD_addCalibrHEB;

  double m_hCalThresBarrel       ;
  double m_hCalThresEndCapHEF    ;
  double m_hCalThresEndCapHEB    ;
  double m_eCalThresBarrel       ;
  double m_eCalThresEndCap       ;

  double m_hCalMipThresBarrel    ;
  double m_hCalMipThresEndCapHEF ;
  double m_hCalMipThresEndCapHEB ;
  double m_eCalMipThresBarrel    ;
  double m_eCalMipThresEndCap    ;

  double m_eCalToMipEndCap       ;
  double m_eCalToMipBarrel       ;
  double m_hCalToMipEndCapHEF    ;
  double m_hCalToMipEndCapHEB    ;
  double m_hCalToMipBarrel       ;

  double m_eCalToEMGeVEndCap     ;
  double m_eCalToEMGeVBarrel     ;
  double m_hCalToEMGeVEndCapHEF  ;
  double m_hCalToEMGeVEndCapHEB  ;
  double m_hCalToEMGeVBarrel     ;

  double m_eCalToHadGeVEndCap    ;
  double m_eCalToHadGeVBarrel    ;
  double m_hCalToHadGeVEndCapHEF ;
  double m_hCalToHadGeVEndCapHEB ;
  double m_hCalToHadGeVBarrel    ;

  double m_muonToMip             ;

  bool firstEvent_ ; 
  short _debugLevel;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

