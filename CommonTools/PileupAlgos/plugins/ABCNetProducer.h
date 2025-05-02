#ifndef ABCNetProducer_h
#define ABCNetProducer_h 

// struct to hold preprocessing parameters
struct PreprocessParams {
  struct VarInfo {
    VarInfo() {}
    VarInfo(float upper_bound, float lower_bound, float norm_factor, float replace_nan_value, float pad_value)
      : upper_bound(upper_bound),
        lower_bound(lower_bound),
	norm_factor(norm_factor),
	replace_nan_value(replace_nan_value),
	pad_value(pad_value) {}
    float upper_bound = 100000;
    float lower_bound = -100000;
    float norm_factor = 1;
    float replace_nan_value = 0;
    float pad_value = 0;
  };

  std::vector<std::string> var_names;
  std::unordered_map<std::string, VarInfo> var_info_map;

  VarInfo info(const std::string &name) const { return var_info_map.at(name); }

};

struct ABCNetTFCache {
  ABCNetTFCache() : graphDef(nullptr) {}
  std::atomic<tensorflow::GraphDef*> graphDef;
};

///////////////////////////////////////////////////////////////////////////////////
// Declaring class ----------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////////

class ABCNetProducer : public edm::stream::EDProducer<edm::GlobalCache<ABCNetTFCache>> {

public:
  explicit ABCNetProducer(const edm::ParameterSet&, const ABCNetTFCache*);   
  ~ABCNetProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions &);

  static std::unique_ptr<ABCNetTFCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ABCNetTFCache*);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;
  std::vector<float> minmax_scale(const std::vector<float> &input,
				  float upper_bound = 100000,
				  float lower_bound = 100000,
				  float norm_factor = 1,
				  float pad_value = 0,
				  float replace_nan_value = 0
				  );
  void preprocess(std::unordered_map<std::string, std::vector<float>> &featureMap, bool debug);
  // tokens
  edm::EDGetTokenT<reco::CandidateView> tokenPFCandidates_;
  std::vector<std::string> input_names_; //names of the input features. Ordering matters!
  std::unordered_map<std::string, PreprocessParams> prep_info_map_;  //preprocessing info for each input feature
  //session for TF evaluation
  tensorflow::Session* session_;
  std::string input_tensor_name_1;
  std::string input_tensor_name_2;
  std::string output_tensor_name_;
  int n_pf_cands_;
  int n_feats_;
  int n_knns_;
  bool debug_;
  
  edm::EDPutTokenT<edm::ValueMap<float>> ABCNetOut_;
};

#endif //ABCNetProducer_h
