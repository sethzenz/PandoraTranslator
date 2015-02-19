import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi")

# No data of type "HGCalGeometry" with label "HGCalEESensitive" in record "IdealGeometryRecord"
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')

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

process.TrackAssociatorRecord = cms.ESSource("EmptyESSource",
        recordName = cms.string('TrackAssociatorRecord'),
        iovIsRunNotTime = cms.bool(True),
        firstValid = cms.vuint32(1)
)
process.load('SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# this is for the event display 
#process.EveService = cms.Service("EveService")

infiles="""
Events_130_100_19.root
Events_130_100_2.root
Events_130_100_20.root
Events_130_100_21.root
Events_130_100_22.root
Events_130_100_24.root
Events_130_100_25.root
Events_130_100_26.root
Events_130_100_27.root
Events_130_100_28.root
Events_130_100_29.root
Events_130_100_3.root
Events_130_100_32.root
Events_130_100_34.root
Events_130_100_35.root
Events_130_100_37.root
Events_130_100_38.root
Events_130_100_39.root
Events_130_100_4.root
Events_130_100_40.root
Events_130_100_41.root
Events_130_100_42.root
Events_130_100_44.root
Events_130_100_45.root
Events_130_100_47.root
Events_130_100_48.root
Events_130_100_49.root
Events_130_100_5.root
Events_130_100_50.root
Events_130_100_6.root
Events_130_100_7.root
Events_130_100_9.root
""".split()

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:simple_jets.root'
    )
)

#["/store/cmst3/group/hgcal/CMSSW/Single130-FixE_CMSSW_6_2_0_SLHC23_patch2/%s"%a for a in infiles]

process.source.skipEvents = cms.untracked.uint32(0)

process.load('SimGeneral.MixingModule.mixNoPU_cfi')                                                                                                                               
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("RecoParticleFlow/PFClusterProducer/particleFlowRecHitHGCEE_cfi")
process.load("RecoParticleFlow/PFClusterProducer/particleFlowRecHitECAL_cfi")
process.load("RecoParticleFlow/PFClusterProducer/particleFlowRecHitHBHE_cfi")

process.trackingParticleRecoTrackAsssociation = cms.EDProducer(
    "TrackAssociatorEDProducer",
    label_tr = cms.InputTag("generalTracks"),
    associator = cms.string('quickTrackAssociatorByHits'),
    label_tp = cms.InputTag("mix","MergedTrackTruth"),
    ignoremissingtrackcollection = cms.untracked.bool(False)
    )

process.pandorapfanew = cms.EDProducer('PandoraCMSPFCandProducer',
    ecalRecHitsEB = cms.InputTag("particleFlowRecHitECAL",""),
    hcalRecHitsHBHE = cms.InputTag("particleFlowRecHitHBHE",""),
    HGCrechitCollection  = cms.InputTag("particleFlowRecHitHGCEE",""), 
    generaltracks = cms.VInputTag(cms.InputTag("generalTracks")),
    tPRecoTrackAsssociation= cms.InputTag("trackingParticleRecoTrackAsssociation"),
    genParticles= cms.InputTag("genParticles"),
#    inputconfigfile = cms.string('PandoraSettingsDefault_WithoutMonitoring.xml'),
#    inputconfigfile = cms.string('PandoraSettingsDefault.xml'),
#    inputconfigfile = cms.string('PandoraSettingsBasic.xml'),
    inputconfigfile = cms.FileInPath('HGCal/PandoraTranslator/data/PandoraSettingsBasic_cms.xml'),
#    inputconfigfile = cms.string('PandoraSettingsMuon.xml')

    energyCorrMethod = cms.string('ABSCORR'),
#   absorber thickness correction
#   energyCorrMethod = cms.string('WEIGHTING'),
    energyWeightFile = cms.FileInPath('HGCal/PandoraTranslator/data/energyWeight.txt'),

    calibrParFile = cms.FileInPath('HGCal/PandoraTranslator/data/pandoraCalibrPars.txt'),
    layerDepthFile = cms.FileInPath('HGCal/PandoraTranslator/data/HGCmaterial_v5.root'),
    overburdenDepthFile = cms.FileInPath('RecoParticleFlow/PFClusterProducer/data/HGCMaterialOverburden.root'),
    useOverburdenCorrection = cms.bool(False), #disabled until the overburden values make sense
    outputFile = cms.string('pandoraoutput.root')
)

process.load('UserCode.HGCanalysis.hgcTrackerInteractionsFilter_cfi')

process.reconstruction_step = cms.Path(#process.trackerIntFilter*
                                       process.particleFlowRecHitHGCEE*
                                       process.particleFlowRecHitHBHE*
                                       process.particleFlowRecHitECAL*
                                       process.trackingParticleRecoTrackAsssociation*
                                       process.pandorapfanew)
process.schedule = cms.Schedule(process.reconstruction_step)
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCalMuon
process = cust_2023HGCalMuon(process)

infiles = """
Events_130_20_1.root
Events_130_20_10.root
Events_130_20_11.root
Events_130_20_14.root
Events_130_20_15.root
Events_130_20_16.root
Events_130_20_17.root
Events_130_20_19.root
Events_130_20_2.root
Events_130_20_21.root
Events_130_20_22.root
Events_130_20_24.root
Events_130_20_25.root
Events_130_20_27.root
Events_130_20_28.root
Events_130_20_29.root
Events_130_20_30.root
Events_130_20_31.root
Events_130_20_33.root
Events_130_20_34.root
Events_130_20_35.root
Events_130_20_37.root
Events_130_20_38.root
Events_130_20_39.root
Events_130_20_4.root
Events_130_20_40.root
Events_130_20_41.root
Events_130_20_42.root
Events_130_20_43.root
Events_130_20_44.root
Events_130_20_45.root
Events_130_20_46.root
Events_130_20_47.root
Events_130_20_49.root
Events_130_20_50.root
Events_130_20_6.root
Events_130_20_7.root
Events_130_20_8.root
Events_130_20_9.root
""".split()
