import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi")

# No data of type "HGCalGeometry" with label "HGCalEESensitive" in record "IdealGeometryRecord"
process.load('Configuration.Geometry.GeometryExtended2023HGCalV4MuonReco_cff')

#The below three lines were added to solve an error Exception Message:
#No "CaloGeometryRecord" record found in the EventSetup.
### process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
### process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
### process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

#Add the next two lines to solve an error Exception Message:
#No "IdealMagneticFieldRecord" record found in the EventSetup.
### process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#Add the next three lines to solve an error Exception Message:
#No "TransientTrackRecord" record found in the EventSetup.
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
# process.GlobalTag.globaltag = 'START70_V1::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# this is for the event display 
#process.EveService = cms.Service("EveService")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/work/v/vandreev/public/RecHits/SLHC14/step3_elecPt35.root'
#        'file:/afs/cern.ch/work/a/apsallid/public/step3_SinglePiPt100_EtaEE.root'
#	 'file:/afs/cern.ch/work/a/apsallid/public/step3_SingleElectronPt35_EtaEE.root'
#        'file:/afs/cern.ch/work/a/apsallid/public/step3_SingleMuPt10_EtaEE.root'
#	 'file:/afs/cern.ch/work/a/apsallid/public/step3_SinglePiPt10_EtaEE.root'

#	 'file:/afs/cern.ch/work/a/apsallid/public/step3_SinglePiPt10.root'
#        'file:/afs/cern.ch/work/a/apsallid/public/step3_SinglePiPt20.root'
#        'file:/afs/cern.ch/work/a/apsallid/public/step3_SingleElectronPt35.root'
#        'file:/afs/cern.ch/work/a/apsallid/public/step3_SingleGammaPt35.root'
        'file:INPUT'
#         'root://eoscms//eos/cms/store/group/phys_b2g/apsallid/hg/SinglePiPt10/Step3Files/step3_2.root' 
#        'file:/afs/cern.ch/work/a/apsallid/public/step3_SingleElectronPt50.root'
#        'file:/afs/cern.ch/work/a/apsallid/public/step3_SinglePi0E20.root'
#	 'file:/afs/cern.ch/user/l/lgray/work/public/CMSSW_6_2_X_SLHC_2014-07-17-0200/src/matrix_tests/140_pu/step3.root'
    )
)

process.source.skipEvents = cms.untracked.uint32(0)

process.pandorapfanew = cms.EDAnalyzer('runPandora',
    ecalRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    hcalRecHitsHBHE = cms.InputTag("reducedHcalRecHits","hbheUpgradeReco"),
    HGCEErechitCollection  = cms.InputTag('HGCalRecHit','HGCEERecHits'), 
    HGCHEFrechitCollection = cms.InputTag('HGCalRecHit','HGCHEFRecHits'), 
    HGCHEBrechitCollection = cms.InputTag('HGCalRecHit','HGCHEBRecHits'), 
    generaltracks = cms.VInputTag(cms.InputTag("generalTracks")),
    tPRecoTrackAsssociation= cms.InputTag("trackingParticleRecoTrackAsssociation"),
    genParticles= cms.InputTag("genParticles"),
#    inputconfigfile = cms.string('PandoraSettingsDefault_WithoutMonitoring.xml'),
#    inputconfigfile = cms.string('PandoraSettingsDefault.xml'),
#    inputconfigfile = cms.string('PandoraSettingsBasic.xml'),
    inputconfigfile = cms.string('PandoraSettingsBasic_WithoutMonitoring.xml'),
#    inputconfigfile = cms.string('PandoraSettingsMuon.xml')

    energyCorrMethod = cms.string('ABSCORR'),
#   absorber thickness correction
#   energyCorrMethod = cms.string('WEIGHTING'),
    energyWeightFile = cms.string('energyWeight.txt'),

    calibrParFile = cms.string('pandoraCalibrPars.txt'),
    outputFile = cms.string('pandoraoutput.root')
)


process.p = cms.Path(process.pandorapfanew)
