// S. Zenz, 12 February 2015
//
// Splits a track collection into two, based on whether they propagate to the HGCal or not
// Tracks with bad pt resolution (suspected fakes) are dropped and not in either collection

#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

//#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"


// for track propagation through HGC  
// N.B. we are only propogating to first layer, so check these later
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

//geometry records
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include<unordered_map>

class HGCalTrackCollectionProducer : public edm::EDProducer {

public:
  HGCalTrackCollectionProducer( const edm::ParameterSet & );
private:
  bool goodPtResolution( const reco::TrackRef& trackref ) const;
  void produce( edm::Event &, const edm::EventSetup & ) override;

  edm::EDGetTokenT<edm::View<reco::PFRecTrack> > _src;
  
  // variables needed for copied function
  // need to go back and figure out sensible values
  bool _debug;
  const std::vector<double> _DPtovPtCut;
  const std::vector<unsigned> _NHitCut;
  const bool _useIterTracking;
};

HGCalTrackCollectionProducer::HGCalTrackCollectionProducer(const edm::ParameterSet & iConfig) :
  _src(consumes<edm::View<reco::PFRecTrack> >(iConfig.getParameter<edm::InputTag> ("src"))),
  _DPtovPtCut(iConfig.getParameter<std::vector<double> >("DPtOverPtCuts_byTrackAlgo")),
  _NHitCut(iConfig.getParameter<std::vector<unsigned> >("NHitCuts_byTrackAlgo")),
  _useIterTracking(iConfig.getParameter<bool>("useIterativeTracking"))
{
  _debug = true; // That's right, I hard-coded debug-mode.

  produces<reco::PFRecTrackCollection>("TracksInHGCal");
  produces<reco::PFRecTrackCollection>("TracksNotInHGCal");

}

void HGCalTrackCollectionProducer::produce(edm::Event & evt, const edm::EventSetup & iSetup) {

  edm::Handle<edm::View<reco::PFRecTrack> > trackHandle;
  evt.getByToken(_src,trackHandle);
  //  const reco::PFRecTrackRefVector& tracks = trackHandle->refVector();
  //  const edm::RefToBaseVector<reco::PFRecTrack>& tracks = trackHandle->refVector();
  const edm::PtrVector<reco::PFRecTrack>& tracks = trackHandle->ptrVector();

  std::auto_ptr<reco::PFRecTrackCollection> outputInHGCal(new reco::PFRecTrackCollection);
  std::auto_ptr<reco::PFRecTrackCollection> outputNotInHGCal(new reco::PFRecTrackCollection);

  for ( unsigned int i = 0 ; i < tracks.size() ; i++) {
    bool isGood = goodPtResolution(tracks[i]->trackRef());
    if (_debug) std::cout << "Track number " << i << " has a goodPtResolution result of" << isGood << std::endl;

  }

  evt.put(outputInHGCal,"TracksInHGCal");
  evt.put(outputNotInHGCal,"TracksNotInHGCal");
}

// Copied from https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/RecoParticleFlow/PFProducer/plugins/importers/GeneralTracksImporter.cc#L140
bool HGCalTrackCollectionProducer::
goodPtResolution( const reco::TrackRef& trackref) const {

  const double P = trackref->p();
  const double Pt = trackref->pt();
  const double DPt = trackref->ptError();
  const unsigned int NHit = 
    trackref->hitPattern().trackerLayersWithMeasurement();
  const unsigned int NLostHit = 
    trackref->hitPattern().trackerLayersWithoutMeasurement();
  const unsigned int LostHits = trackref->numberOfLostHits();
  const double sigmaHad = sqrt(1.20*1.20/P+0.06*0.06) / (1.+LostHits);

  // iteration 1,2,3,4,5 correspond to algo = 1/4,5,6,7,8,9
  unsigned int Algo = 0; 
  switch (trackref->algo()) {
  case reco::TrackBase::ctf:
  case reco::TrackBase::iter0:
  case reco::TrackBase::iter1:
  case reco::TrackBase::iter2:
    Algo = 0;
    break;
  case reco::TrackBase::iter3:
    Algo = 1;
    break;
  case reco::TrackBase::iter4:
    Algo = 2;
    break;
  case reco::TrackBase::iter5:
    Algo = 3;
    break;
  case reco::TrackBase::iter6:
    Algo = 4;
    break;
  default:
    Algo = _useIterTracking ? 5 : 0;
    break;
  }

  // Protection against 0 momentum tracks
  if ( P < 0.05 ) return false;

  // Temporary : Reject all tracking iteration beyond 5th step. 
  if ( Algo > 4 ) return false;
 
  if (_debug) std::cout << " PFBlockAlgo: PFrecTrack->Track Pt= "
			<< Pt << " DPt = " << DPt << std::endl;
  if ( ( _DPtovPtCut[Algo] > 0. && 
	 DPt/Pt > _DPtovPtCut[Algo]*sigmaHad ) || 
       NHit < _NHitCut[Algo] ) { 
    // (Algo >= 3 && LostHits != 0) ) {
    if (_debug) std::cout << " PFBlockAlgo: skip badly measured track"
			  << ", P = " << P 
			  << ", Pt = " << Pt 
			  << " DPt = " << DPt 
			  << ", N(hits) = " << NHit << " (Lost : " << LostHits << "/" << NLostHit << ")"
			  << ", Algo = " << Algo
			  << std::endl;
    if (_debug) std::cout << " cut is DPt/Pt < " << _DPtovPtCut[Algo] * sigmaHad << std::endl;
    if (_debug) std::cout << " cut is NHit >= " << _NHitCut[Algo] << std::endl;
    /*
      std::cout << "Track REJECTED : ";
      std::cout << ", P = " << P 
      << ", Pt = " << Pt 
      << " DPt = " << DPt 
      << ", N(hits) = " << NHit << " (Lost : " << LostHits << "/" << NLostHit << ")"
      << ", Algo = " << Algo
      << std::std::endl;
    */
    return false;
  }
  /*
    std::cout << "Track Accepted : ";
    std::cout << ", P = " << P 
    << ", Pt = " << Pt 
    << " DPt = " << DPt 
    << ", N(hits) = " << NHit << " (Lost : " << LostHits << "/" << NLostHit << ")"
    << ", Algo = " << Algo
    << std::std::endl;
  */
  return true;
}

