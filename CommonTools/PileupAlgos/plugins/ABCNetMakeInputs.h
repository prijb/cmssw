#ifndef ABCNetMakeInputs_h
#define ABCNetMakeInputs_h

#include <string>
#include <vector>
#include <unordered_map>

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"

namespace abcnet {

  class ABCNetMakeInputs {
  public:
    ABCNetMakeInputs() {}
    ~ABCNetMakeInputs() {}

    static std::tuple< std::unordered_map<std::string, std::vector<float>>, std::vector<float> >  makeFeatureMap( const reco::CandidateView * PFCol, std::vector<size_t> & indices, int n_pf_cands, int n_knn, bool debug = false);

  };
}  // namespace abcnet

#endif  // ABCNetMakeInputs_h
