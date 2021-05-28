#include "RecoMuon/TrackerSeedGenerator/interface/SeedMvaEstimator.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CommonTools/MVAUtils/interface/GBRForestTools.h"

using namespace std;

SeedMvaEstimator::SeedMvaEstimator(const edm::FileInPath& weightsfile, std::vector<double> scale_mean, std::vector<double> scale_std) {
  gbrForest_  = createGBRForest(weightsfile);
  scale_mean_ = scale_mean;
  scale_std_  = scale_std;
}

SeedMvaEstimator::~SeedMvaEstimator() {}

namespace {
  enum inputIndexes {
    kTsosErr0,         // 0
    kTsosErr2,         // 1
    kTsosErr5,         // 2
    kTsosDxdz,         // 3
    kTsosDydz,         // 4
    kTsosQbp,          // 5
    kDRdRL1SeedP,      // 6
    kDPhidRL1SeedP,    // 7
    kLastL1,           // 8

    kDRdRL2SeedP = 8,  // 8
    kDPhidRL2SeedP,    // 9
    kLastL2,           // 10
  };
}  // namespace

void SeedMvaEstimator::getL1MuonVariables(
  const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::MuonBxCollection> h_L1Muon,
  float& dRdRL1SeedP,
  float& dPhidRL1SeedP,
  float& dRdPhiL1SeedX,
  float& dPhidPhiL1SeedX ) const {

  for(int ibx = h_L1Muon->getFirstBX(); ibx<=h_L1Muon->getLastBX(); ++ibx)
  {
    if(ibx != 0) continue; // -- only take when ibx == 0 -- //
    for(auto it=h_L1Muon->begin(ibx); it!=h_L1Muon->end(ibx); it++)
    {
      l1t::MuonRef ref_L1Mu(h_L1Muon, distance(h_L1Muon->begin(h_L1Muon->getFirstBX()), it) );

      // FIXME: 7 should be configurable
      if(ref_L1Mu->hwQual() < 7)
        continue;

      float dR_L1SeedP_AtVtx   = reco::deltaR( ref_L1Mu->etaAtVtx(), ref_L1Mu->phiAtVtx(), global_p.eta(), global_p.phi());
      float dPhi_L1SeedP_AtVtx = reco::deltaPhi( ref_L1Mu->phiAtVtx(), global_p.phi());
      float dR_L1SeedX_AtVtx   = reco::deltaR( ref_L1Mu->etaAtVtx(), ref_L1Mu->phiAtVtx(), global_x.eta(), global_x.phi());
      float dPhi_L1SeedX_AtVtx = reco::deltaPhi( ref_L1Mu->phiAtVtx(), global_x.phi());

      if( dR_L1SeedP_AtVtx < dRdRL1SeedP ) {
        dRdRL1SeedP = dR_L1SeedP_AtVtx;
        dPhidRL1SeedP = dPhi_L1SeedP_AtVtx;
      }
      if( fabs(dPhi_L1SeedX_AtVtx) < fabs(dPhidPhiL1SeedX) ) {
        dRdPhiL1SeedX = dR_L1SeedX_AtVtx;
        dPhidPhiL1SeedX = dPhi_L1SeedX_AtVtx;
      }
    }
  }
}

void SeedMvaEstimator::getL2MuonVariables(
  const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon,
  float& dRdRL2SeedP,
  float& dPhidRL2SeedP,
  float& dRdPhiL2SeedX,
  float& dPhidPhiL2SeedX ) const {

  for( unsigned int i_L2=0; i_L2<h_L2Muon->size(); i_L2++)
  {
    reco::RecoChargedCandidateRef ref_L2Mu(h_L2Muon, i_L2);

    float dR_L2SeedP   = reco::deltaR( *ref_L2Mu, global_p);
    float dPhi_L2SeedP = reco::deltaPhi( ref_L2Mu->phi(), global_p.phi());
    float dR_L2SeedX   = reco::deltaR( *ref_L2Mu, global_x);
    float dPhi_L2SeedX = reco::deltaPhi( ref_L2Mu->phi(), global_x.phi());

    if( dR_L2SeedP < dRdRL2SeedP ) {
      dRdRL2SeedP = dR_L2SeedP;
      dPhidRL2SeedP = dPhi_L2SeedP;
    }
    if( fabs(dPhi_L2SeedX) < fabs(dPhidPhiL2SeedX) ) {
      dRdPhiL2SeedX = dR_L2SeedX;
      dPhidPhiL2SeedX = dPhi_L2SeedX;
    }
  }
}

float SeedMvaEstimator::computeMva( const TrajectorySeed& seed,
  GlobalVector global_p,
  GlobalPoint  global_x,
  edm::Handle<l1t::MuonBxCollection> h_L1Muon,
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon,
  bool IsFromL1
) const {

  if(IsFromL1){
    
    float var[kLastL1]{};

    var[kTsosErr0]   = seed.startingState().error(0);
    // var[kTsosErr1]   = seed.startingState().error(1);
    var[kTsosErr2]   = seed.startingState().error(2);
    // var[kTsosErr3]   = seed.startingState().error(3);
    // var[kTsosErr4]   = seed.startingState().error(4);
    var[kTsosErr5]   = seed.startingState().error(5);
    // var[kTsosErr6]   = seed.startingState().error(6);
    // var[kTsosErr7]   = seed.startingState().error(7);
    // var[kTsosErr8]   = seed.startingState().error(8);
    // var[kTsosErr9]   = seed.startingState().error(9);
    // var[kTsosErr10]  = seed.startingState().error(10);
    // var[kTsosErr11]  = seed.startingState().error(11);
    // var[kTsosErr12]  = seed.startingState().error(12);
    // var[kTsosErr13]  = seed.startingState().error(13);
    // var[kTsosErr14]  = seed.startingState().error(14);
    var[kTsosDxdz]   = seed.startingState().parameters().dxdz();
    var[kTsosDydz]   = seed.startingState().parameters().dydz();
    var[kTsosQbp]    = seed.startingState().parameters().qbp();
    // var[kTsosCharge] = seed.startingState().parameters().charge();

    // FIXME: should be configurable
    float initDRdPhi = 99999.;

    float dRdRL1SeedP = initDRdPhi;
    float dPhidRL1SeedP = initDRdPhi;
    float dRdPhiL1SeedX = initDRdPhi;
    float dPhidPhiL1SeedX = initDRdPhi;
    getL1MuonVariables( seed, global_p, global_x, h_L1Muon, dRdRL1SeedP, dPhidRL1SeedP, dRdPhiL1SeedX, dPhidPhiL1SeedX );

    float dRdRL2SeedP = initDRdPhi;
    float dPhidRL2SeedP = initDRdPhi;
    float dRdPhiL2SeedX = initDRdPhi;
    float dPhidPhiL2SeedX = initDRdPhi;
    getL2MuonVariables( seed, global_p, global_x, h_L2Muon, dRdRL2SeedP, dPhidRL2SeedP, dRdPhiL2SeedX, dPhidPhiL2SeedX );

    var[kDRdRL1SeedP]     = dRdRL1SeedP;
    var[kDPhidRL1SeedP]   = dPhidRL1SeedP;
    // var[kDRdPhiL1SeedX]   = dRdPhiL1SeedX;
    // var[kDPhidPhiL1SeedX] = dPhidPhiL1SeedX;
    // var[kDRdRL2SeedP]     = dRdRL2SeedP;
    // var[kDPhidRL2SeedP]   = dPhidRL2SeedP;
    // var[kDRdPhiL2SeedX]   = dRdPhiL2SeedX;
    // var[kDPhidPhiL2SeedX] = dPhidPhiL2SeedX;

    for(int iv=0; iv<kLastL1; ++iv) {
      var[iv] = (var[iv] - scale_mean_.at(iv)) / scale_std_.at(iv);
    }

    return gbrForest_->GetResponse( var );
  }
  else{

    float var[kLastL2]{};

    var[kTsosErr0]   = seed.startingState().error(0);
    var[kTsosErr2]   = seed.startingState().error(2);
    var[kTsosErr5]   = seed.startingState().error(5);
    var[kTsosDxdz]   = seed.startingState().parameters().dxdz();
    var[kTsosDydz]   = seed.startingState().parameters().dydz();
    var[kTsosQbp]    = seed.startingState().parameters().qbp();

    float initDRdPhi = 99999.;

    float dRdRL1SeedP = initDRdPhi;
    float dPhidRL1SeedP = initDRdPhi;
    float dRdPhiL1SeedX = initDRdPhi;
    float dPhidPhiL1SeedX = initDRdPhi;
    getL1MuonVariables( seed, global_p, global_x, h_L1Muon, dRdRL1SeedP, dPhidRL1SeedP, dRdPhiL1SeedX, dPhidPhiL1SeedX );

    float dRdRL2SeedP = initDRdPhi;
    float dPhidRL2SeedP = initDRdPhi;
    float dRdPhiL2SeedX = initDRdPhi;
    float dPhidPhiL2SeedX = initDRdPhi;
    getL2MuonVariables( seed, global_p, global_x, h_L2Muon, dRdRL2SeedP, dPhidRL2SeedP, dRdPhiL2SeedX, dPhidPhiL2SeedX );

    var[kDRdRL1SeedP]     = dRdRL1SeedP;
    var[kDPhidRL1SeedP]   = dPhidRL1SeedP;
    var[kDRdRL2SeedP]     = dRdRL2SeedP;
    var[kDPhidRL2SeedP]   = dPhidRL2SeedP;

    for(int iv=0; iv<kLastL2; ++iv) {
      var[iv] = (var[iv] - scale_mean_.at(iv)) / scale_std_.at(iv);
    }

    return gbrForest_->GetResponse( var );
  }
}
