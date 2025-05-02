//========================================================================================
// Class:      ABCNetProducer                     ---------------------------------------
//----------------------------------------------------------------------------------------
/**\class ABCNetProducer ABCNetProducer.cc PileupAlgos/plugins/ABCNetProducer.cc
------------------------------------------------------------------------------------------
 Description: This class produces ABCNet weights to be used to mitigate pileup  ---
 -----------------------------------------------------------------------------------------
 Implementation:                                                                       ---
     This EDProducer is meant to be used with CMSSW >= 10_6_26                         ---
*/
//========================================================================================
// Authors:  Fabio Iemmi (IHEP)                                      ---------------------
//         Created:  MON, 23 May 2022 10:00:28 GMT  ---------------------------------------
//========================================================================================

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h" // for edm::ParameterSet
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //for pat
#include "DataFormats/Candidate/interface/Candidate.h" // for reco
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"// for puppi candidates
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h" // to use TensorFlow
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"

#include "ABCNetProducer.h"
#include "ABCNetMakeInputs.h"

#include <limits>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <nlohmann/json.hpp>

using namespace abcnet;

//////////////////////////////////////////////////////////////////////////////////
// Defining methods --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////

// constructors
ABCNetProducer::ABCNetProducer(const edm::ParameterSet& iConfig, const ABCNetTFCache* cache):
  tokenPFCandidates_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candName"))),
  session_(nullptr),
  input_tensor_name_1(iConfig.getParameter<std::string>("input_tensor_name_1")),
  input_tensor_name_2(iConfig.getParameter<std::string>("input_tensor_name_2")),
  output_tensor_name_(iConfig.getParameter<std::string>("output_tensor_name")),
  n_pf_cands_(iConfig.getParameter<int>("n_pf_cands")),
  n_feats_(iConfig.getParameter<int>("n_feats")),
  n_knns_(iConfig.getParameter<int>("n_knns")),
  debug_(iConfig.getParameter<bool>("debug"))
{
  //parse data from preprocessing JSON file
  std::ifstream ifs(iConfig.getParameter<edm::FileInPath>("preprocess_json").fullPath());
  nlohmann::json js = nlohmann::json::parse(ifs);
  js.at("input_names").get_to(input_names_); //store names of input features from JSON inside input_names_
  const auto &var_info_pset = js.at("var_info"); //get the var_info block from JSON
  
  for (const auto & input_name : input_names_) {
    const auto &var_pset = var_info_pset.at(input_name);
    double upper_bound = var_pset.at("upper_bound");
    double lower_bound = var_pset.at("lower_bound");
    double norm_factor = var_pset.at("norm_factor");
    double replace_nan_value = var_pset.at("replace_nan_value");
    double pad_value = var_pset.at("pad_value");
    auto &prep_params = prep_info_map_[input_name];
    prep_params.var_info_map[input_name] = PreprocessParams::VarInfo(upper_bound,
								     lower_bound,
								     norm_factor,
								     replace_nan_value,
								     pad_value
								     );
  }
  
  //create the TF session using the meta graph from the cache
  session_ = tensorflow::createSession(cache->graphDef);

  // Produce a ValueMap of floats linking each PF candidate with its ABCNet weight
  ABCNetOut_ = produces<edm::ValueMap<float>> ();

};

// destructor
ABCNetProducer::~ABCNetProducer() {
};

void ABCNetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("candName", edm::InputTag("packedPFCandidates"));
  desc.add<edm::FileInPath>("graph_path", edm::FileInPath("CommonTools/PileupAlgos/data/AttentionBasedPileupRejectionModel_Run2_KDTree_dyn0pad_test.pb"));
  desc.add<edm::FileInPath>("preprocess_json", edm::FileInPath("CommonTools/PileupAlgos/data/preprocessing_info.json"));
  desc.add<std::string>("input_tensor_name_1", "input_1");
  desc.add<std::string>("input_tensor_name_2", "input_2");
  desc.add<std::string>("output_tensor_name", "Identity");
  desc.add<int>("n_pf_cands", 4000);
  desc.add<int>("n_feats", 19);
  desc.add<int>("n_knns", 20);
  desc.add<bool>("debug", false);

  descriptions.add("abc", desc);
};

std::unique_ptr<ABCNetTFCache> ABCNetProducer::initializeGlobalCache(const edm::ParameterSet& iConfig) {
  tensorflow::setLogging("3");
  // get the pb file
  std::string pbFile = iConfig.getParameter<edm::FileInPath>("graph_path").fullPath();
  // load the graph definition and save it in the cache
  ABCNetTFCache* cache = new ABCNetTFCache();
  cache->graphDef = tensorflow::loadGraphDef(pbFile);
  return std::unique_ptr<ABCNetTFCache>(cache);
};

void ABCNetProducer::globalEndJob(const ABCNetTFCache* cache) {
  if (cache->graphDef != nullptr) {
    delete cache->graphDef;
  }
};

std::vector<float> ABCNetProducer::minmax_scale(const std::vector<float> &input,
						float upper_bound,
						float lower_bound,
						float norm_factor,
						float pad_value,
						float replace_nan_value
						) {
  std::vector<float> out(input.size(), pad_value);
  for (unsigned i = 0; i < input.size(); ++i) {
    //first of all, check if input is inf/nan and replace it if needed
    auto val = std::isfinite(input[i]) ? input[i] : replace_nan_value;
    val = std::isnan(input[i]) ? replace_nan_value : input[i];
    //second of all, perform the actual min-max scaling
    val = (val>upper_bound) ? upper_bound : val;
    val = (val<lower_bound) ? lower_bound : val;
    out[i] = val/norm_factor; //need to check this value
  }
  return out;
};

void ABCNetProducer::preprocess(std::unordered_map<std::string, std::vector<float>> & featureMap, bool debug) {

  for (const auto & input_name : input_names_) {
    const auto & preprocessing_params = prep_info_map_[input_name].var_info_map[input_name];
    featureMap[input_name] = minmax_scale(featureMap[input_name],
					  preprocessing_params.upper_bound,
					  preprocessing_params.lower_bound,
					  preprocessing_params.norm_factor,
					  preprocessing_params.pad_value,
					  preprocessing_params.replace_nan_value
					  );
  }
  
  if(debug) {
    std::cout << "*** NOW CHECKING THE SANITY OF THE PREPROCESSING STEP ***" << std::endl;
    for(auto const & feat : featureMap) {
      bool issue = false;
      std::cout << feat.first << std::endl;
      if (feat.first == "PFCandEta" || feat.first == "PFCandPhi" || feat.first == "PFCandLogPt" || feat.first == "PFCandLogE" || feat.first == "PFCandCharge") {
	std::cout << "This feature is expected to take values outside [-1;1] interval. Skipping..." << std::endl;
	continue;
      } 
      for (auto const & element : feat.second) {
	if (element < -1 || element > 1) {
	  std::cout << "*** WARNING *** Found element with value outside [-1;1] interval" << std::endl;
	  issue = true;
	  break;
	}
      }
      if (!issue) std::cout << "No issues with input feature detected" << std::endl;
    }
  }
  
};

void ABCNetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Get PFCandidate Collection
  edm::Handle<reco::CandidateView> PFCandidates;
  iEvent.getByToken(tokenPFCandidates_, PFCandidates);
  const reco::CandidateView *pfCol = PFCandidates.product();
  
  //make feature map and KNNs, preprocess features
  std::vector<size_t> indices; 
  auto [features, KNNs] = ABCNetMakeInputs::makeFeatureMap(pfCol, indices, n_pf_cands_, n_knns_, debug_);
  ABCNetProducer::preprocess(features, debug_);
  
  //fill the input tensors
  tensorflow::Tensor inputs (tensorflow::DT_FLOAT, { 1, (int)features["PFCandEta"].size(), n_feats_ });
  inputs.flat<float>().setZero();
  for (int j = 0; j < n_feats_; j++) {
    for (int i = 0; i < (int)features["PFCandEta"].size(); i++) {
      inputs.tensor<float,3>()(0,i,j) = float(features[input_names_.at(j)].at(i)); //looks suboptimal; is there way of filling the tensor avoiding nested loops?
    }
  }
  tensorflow::Tensor knn_indices (tensorflow::DT_FLOAT, { 1, (int)features["PFCandEta"].size(), n_knns_ });
  knn_indices.flat<float>().setZero();
  for (int i = 0; i < (int)features["PFCandEta"].size(); i++) {
    for (int j = 0; j < n_knns_; j++){
      knn_indices.tensor<float,3>()(0,i,j) = float(KNNs.at(i*n_knns_+j));
    }
  }
  if (debug_) {
    std::cout << "PRINTING NETWORK INPUTS" << std::endl;
    for (int i = 0; i < 3; i++) {
      std::cout << "Features: \n";
    for (int j = 0; j < n_feats_; j++) {
      std::cout << inputs.tensor<float,3>()(0, i, j) << ", ";
    }
    std::cout << "\n";
    std::cout << "KNN indices: \n";
    for (int j = 0; j < n_knns_; j++) {
      std::cout << knn_indices.tensor<float,3>()(0, i, j) << ", ";
    }
    std::cout << "\n";
    }
  }
  std::vector<tensorflow::Tensor> outputs;
  tensorflow::run(session_, { { input_tensor_name_1, inputs }, { input_tensor_name_2, knn_indices } }, { output_tensor_name_ }, &outputs);
  if (debug_) {
    std::cout << "PRINTING NETWORK OUTPUTS" << std::endl;
    for (int i = 0; i < n_pf_cands_; i++) std::cout << outputs.at(0).tensor<float,3>()(0, i, n_feats_) << " ";
  
  }
  
  //initialize container for ABCNet weights and fill it
  std::vector<float> weights;
  
  int PFCounter = 0;
  for(auto const& aPF : *pfCol) {
    const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&aPF);
    float abcweight = -1.;
    if(lPack == nullptr) { 
      // throw error
      throw edm::Exception(edm::errors::LogicError,"ABCNetProducer: cannot get weights since inputs are not PackedCandidates");
    }
    else{
      if (indices.at(PFCounter) >= (unsigned)n_pf_cands_) { //if the particle wasn't considered in evaluation, assign 0 weight
	abcweight = 0.0;
      }
      else abcweight = outputs.at(0).tensor<float,3>()(0, indices.at(PFCounter), n_feats_);
      //abcweight = (abcweight > 0.01) ? abcweight : 0.0; //set a lower threshold to ABCNet weights? to be tested
    }
    weights.push_back(abcweight);
    PFCounter++;
  }

  edm::ValueMap<float> ABCNetOut;
  edm::ValueMap<float>::Filler ABCNetFiller(ABCNetOut);
  ABCNetFiller.insert(PFCandidates,weights.begin(),weights.end());
  ABCNetFiller.fill();
  iEvent.emplace(ABCNetOut_, ABCNetOut);

}
//define this as a plug-in
DEFINE_FWK_MODULE(ABCNetProducer);
