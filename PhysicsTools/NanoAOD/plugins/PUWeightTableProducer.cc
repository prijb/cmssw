#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/TriggerResults.h"

// Gen event info (pt hat)
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// Standard library includes
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <functional>

// Include the QCD weight calculator
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "PhysicsTools/NanoAOD/interface/QCDWeightCalc.h"
//#include "PhysicsTools/NanoAOD/interface/QCDReweightInfo.h"

// Make a table
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class PUWeightTableProducer : public edm::global::EDProducer<> {
public:
  explicit PUWeightTableProducer(const edm::ParameterSet& iConfig);
  ~PUWeightTableProducer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // The produce method: note the EDProducer for a global module now has an extra StreamID argument.
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

private:
  // Tokens to retrieve inputs from the event
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfoToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

  // Members for the QCD weight calculation
  float bxFreq = 30E6;
  edm::FileInPath qcdWeightFile;
  QCDWeightCalc qcdWeightCalc;

  bool muFilter, emFilter;
};

PUWeightTableProducer::PUWeightTableProducer(const edm::ParameterSet& iConfig)
: genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
  pileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupInfo"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  qcdWeightFile(iConfig.getParameter<edm::FileInPath>("qcdWeightFile")),
  qcdWeightCalc(iConfig.getParameter<edm::FileInPath>("qcdWeightFile"), bxFreq),
  muFilter(iConfig.getParameter<bool>("muFilter")),
  emFilter(iConfig.getParameter<bool>("emFilter"))
{
  // Declare that this producer will output a QCDReweightInfo product
  //produces<float>();
  produces<nanoaod::FlatTable>("puWeight");
}

void PUWeightTableProducer::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  // Create a new output product
  auto output = std::make_unique<float>(1.0);
  //output = std::make_unique<float>();

  // Retrieve the GenEventInfo product
  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(genEventInfoToken_, genEventInfo);
  
  // Retrieve the pileup information
  edm::Handle<std::vector<PileupSummaryInfo>> pileupInfo;
  iEvent.getByToken(pileupInfoToken_, pileupInfo);
  
  // Get the hard process scale (pT hat) from GenEventInfo
  float hardPtHat = genEventInfo->qScale();
  //std::cout << "Hard process pT hat: " << hardPtHat << std::endl;

  // Loop over pileup info and select only in-time (bunch crossing 0)
  std::vector<PileupSummaryInfo> pileupInfoIntime;
  for (const auto& pu : *pileupInfo) {
    if (pu.getBunchCrossing() == 0) {
      pileupInfoIntime.push_back(pu);
    }
  }
  
  // Extract the pileup pT hats (if available) and sort them in descending order
  std::vector<float> puPtHats;
  if (!pileupInfoIntime.empty()) {
    puPtHats = pileupInfoIntime[0].getPU_pT_hats();
  }
  std::sort(puPtHats.begin(), puPtHats.end(), std::greater<float>());
  
  // Here, trigger information is not used so we simply set these booleans to false
  bool pass_em = emFilter;
  bool pass_mu = muFilter;
  
  // Compute the new weight using your QCDWeightCalc
  //std::cout << "Calculating event weight..." << std::endl;
  float calcweight = qcdWeightCalc.weight(hardPtHat, puPtHats, pass_em, pass_mu);
  *output = calcweight;

  // DEBUG
  //std::cout << "Event weight produced with filter options (em: " << pass_em << ", mu: " << pass_mu << ") = " << *output << std::endl;

  // Place the product into the event so that downstream modules can use it
  //iEvent.put(std::move(output));

  auto table = std::make_unique<nanoaod::FlatTable>(1, "puWeight", true);
  table->addColumnValue<float>("weight", *output, "Event weight calculated from QCDWeightCalc", 14);
  iEvent.put(std::move(table), "puWeight");
}

void PUWeightTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  //desc.add<edm::InputTag>("genEventInfo", edm::InputTag("generator"));
  //desc.add<edm::InputTag>("PileupInfo", edm::InputTag("slimmedAddPileupInfo"));
  //desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  //desc.add<edm::FileInPath>("qcdWeightFile", edm::InputTag("PhysicsTools/NanoAOD/data/qcdWeights.json"));
  
  desc.setUnknown();
  descriptions.addDefault(desc);  
  //descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PUWeightTableProducer);
