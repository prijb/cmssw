#include <iostream>
#include <memory>
#include <cmath>
#include <numeric>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //for pat
#include "DataFormats/Candidate/interface/Candidate.h" // for reco

#include "PhysicsTools/NearestNeighbors/interface/NearestNeighbors.h"

#include "ABCNetMakeInputs.h"

using namespace abcnet;
using PointCloud = cms::nanoflann::PointCloud<float, 2>; //default value for the DIM parameter is 3, but want to gather using 2 dims only (eta, phi). See template definition L11 of PhysicsTools/NearestNeighbors/interface/NearestNeighbors.h

template <typename A, typename B> void zip(const std::vector<A> &a, const std::vector<B> &b, std::vector<std::pair<A,B>> &zipped) {
  for(size_t i=0; i<a.size(); ++i) {
    zipped.push_back(std::make_pair(a[i], b[i]));
  }
}

template <typename A, typename B> void unzip( const std::vector<std::pair<A, B>> &zipped, std::vector<A> &a, std::vector<B> &b) {
  for(size_t i=0; i<a.size(); i++) {
    a[i] = zipped[i].first;
    b[i] = zipped[i].second;
  }
}

template <typename A> std::vector<size_t> get_indices(std::vector<A> & v) {
  //initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  //first, get indices mapping sorted to unsorted elements
  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
  //then, invert the mapping to get unsorted ---> sorted mapping
  std::vector<size_t> idx_inverse(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) idx_inverse[idx[i]] = i;
  
  return idx_inverse;
}

std::tuple< std::unordered_map<std::string, std::vector<float>>, std::vector<float> > ABCNetMakeInputs::makeFeatureMap (const reco::CandidateView * pfCol, std::vector<size_t> & indices, int n_pf_cands, int n_knns, bool debug) {
  
  //store PF candidates and their pts into vectors
  std::vector<pat::PackedCandidate> PFCands;
  std::vector<float> Pts;
  
  for (auto const & aPF : *pfCol ) {
    const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&aPF);    
    PFCands.push_back(*lPack);
    Pts.push_back(lPack->pt());
  } // end loop over PF candidates
  
  //get indices mapping unsorted and sorted PF candidates
  indices = get_indices(Pts);
  
  //sort PF candidates based on their Pt
  //zip the vectors
  std::vector<std::pair<pat::PackedCandidate,float>> zipped_vec;
  zip(PFCands, Pts, zipped_vec);
  // Sort the vector of pairs
  std::stable_sort(std::begin(zipped_vec), std::end(zipped_vec), [&](const auto& a, const auto& b) { return a.second > b.second; });
  // Write the sorted pairs back to the original vectors
  unzip(zipped_vec, PFCands, Pts);

  //fill feature map
  std::unordered_map<std::string, std::vector<float>> fts;
  PointCloud::Points points;
  
  for (auto const & aPF : PFCands) {
    
    points.push_back({{float(aPF.eta()), float(aPF.phi())}});
    
    fts["PFCandEta"].push_back(aPF.eta()); //f0
    fts["PFCandPhi"].push_back(aPF.phi()); //f1
    fts["PFCandLogPt"].push_back(log(aPF.pt())); //f2
    fts["PFCandLogE"].push_back(log(aPF.energy())); //f3c
    fts["PFCandCharge"].push_back(aPF.charge()); //f4
    if (abs(aPF.pdgId()) == 11) fts["PFCandIsEle"].push_back(1.0); else fts["PFCandIsEle"].push_back(0.0); //f5
    if (abs(aPF.pdgId()) == 13) fts["PFCandIsMu"].push_back(1.0); else fts["PFCandIsMu"].push_back(0.0); //f6
    if (abs(aPF.pdgId()) == 211 || abs(aPF.pdgId()) == 2 || abs(aPF.pdgId()) == 1) fts["PFCandIsHad"].push_back(1.0); else fts["PFCandIsHad"].push_back(0.0); //f7
    if (abs(aPF.pdgId()) == 130) fts["PFCandIsNeutralHad"].push_back(1.0); else fts["PFCandIsNeutralHad"].push_back(0.0); //f8
    if (abs(aPF.pdgId()) == 22) fts["PFCandIsPhoton"].push_back(1.0); else fts["PFCandIsPhoton"].push_back(0.0); //f9
    fts["PFCandNumHits"].push_back(aPF.numberOfHits()); //f14
    fts["PFCandNumLayersHit"].push_back(aPF.trackerLayersWithMeasurement()); //f15
    fts["PFCandFromPV"].push_back(aPF.fromPV()); //f16
    fts["PFCandTrackHighPurity"].push_back(aPF.trackHighPurity()); //f17
    if (aPF.bestTrack()) {
      if ( isinf(aPF.dz()/aPF.dzError()) ) fts["PFCandDZSig"].push_back(999.0); else fts["PFCandDZSig"].push_back(aPF.dz()/aPF.dzError()); //f10
      if ( isinf(aPF.dxy()/aPF.dxyError()) ) fts["PFCandDXYSig"].push_back(999.0); else fts["PFCandDXYSig"].push_back(aPF.dxy()/aPF.dxyError()); //f11
      const auto *trk = aPF.bestTrack();
      fts["PFCandNormChi2"].push_back(trk->normalizedChi2()); //f13
      fts["PFCandQuality"].push_back(trk->qualityMask()); //f18
    }
    else {
      fts["PFCandDZSig"].push_back(0.0);
      fts["PFCandDXYSig"].push_back(0.0);
      fts["PFCandNormChi2"].push_back(999.0);
      fts["PFCandQuality"].push_back(0.0);
    }
    if (aPF.pdgId() == 130 || aPF.pdgId() == 1) fts["PFCandHCalFrac"].push_back(aPF.hcalFraction()); else if (aPF.isIsolatedChargedHadron()) fts["PFCandHCalFrac"].push_back(aPF.rawHcalFraction()); else fts["PFCandHCalFrac"].push_back(0.0); //f12

    //make sure you don't store more than n_pf_cands (4000 in Run2), as the network is expecting exactly n_pf_cands. If you have less, 0-pad later, if you have more, cut them away now
    if ((int)(fts["PFCandEta"].size()) == n_pf_cands) break;
  }
  if (debug) {
    std::cout << "*** NOW CHECKING THE SORTING OF PF CANDIDATES ***" << std::endl;
    if ( !std::is_sorted(fts["PFCandLogPt"].begin(), fts["PFCandLogPt"].end(), std::greater<float>() ) ) std::cout << "*** WARNING *** PFCandidates are not sorted in Pt!!!" << std::endl;
    else std::cout << "No issues with the sorting of PF candidates detected" << std::endl;
  }
  
  // method to gather kNN indices
  auto knn = [&](size_t num_neighbors, size_t max_support_size, size_t max_query_size = 0) {
    if (max_query_size == 0)
      max_query_size = max_support_size;

    PointCloud::Points supports(points.begin(), points.begin() + std::min(max_support_size, points.size()));
    PointCloud::Points queries(points.begin(), points.begin() + std::min(max_query_size, points.size()));
    auto result = PointCloud::knn<float>(supports, queries, num_neighbors);  // queries.size() * num_neighbors
    std::vector<float> output(max_query_size * num_neighbors, max_support_size - 1);
    std::copy(result.begin(), result.begin() + std::min(result.size(), output.size()), output.begin());

    //std::cout << "[knn] k=" << num_neighbors << ", num_points=" << points.size()
    //<< ", max_support_size=" << max_support_size << ", max_query_size=" << max_query_size
    //<< ", knn_result_size=" << result.size() << ", final_output_size=" << output.size() << std::endl;
    return output;
  };

  std::vector<float> KNNs = knn(n_knns+1, points.size(), points.size()); //gather n_knns+1 k-nearest neighbors. This vector has size (n_particles * (n_knns+1))
  //std::cout << "size of KNNs is " << KNNs.size() << std::endl;
  //for each particle, the first KNN is the particle itself. The training was done excluding this first KNN. Exclude it during evaluation as well
  for (int i = points.size()-1; i >= 0; i--) KNNs.erase(KNNs.begin()+i*(n_knns+1)); //remove elements starting from the end
  return {fts,KNNs};

};
