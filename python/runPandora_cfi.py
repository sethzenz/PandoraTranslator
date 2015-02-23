import FWCore.ParameterSet.Config as cms

pandorapfanew = cms.EDProducer('PandoraCMSPFCandProducer',
    debugPrint = cms.bool(False), #for cout statements
    debugHisto = cms.bool(False), #for diagnostic/calibration histograms
    useRecoTrackAsssociation = cms.bool(False), #needed to turn off for 140PU                                                               
    HGCrechitCollection  = cms.InputTag("particleFlowRecHitHGCEE",""), 
    generaltracks = cms.InputTag("HGCalTrackCollection","TracksInHGCal"),
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
