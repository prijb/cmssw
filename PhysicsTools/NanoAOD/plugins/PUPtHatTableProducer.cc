// Adds the PUPtHats to the nanoAOD
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include <iostream>
#include <algorithm>

class PUPtHatTableProducer : public edm::global::EDProducer<> {
public:
    PUPtHatTableProducer(edm::ParameterSet const& params)
        : puTag_(consumes<std::vector<PileupSummaryInfo>>(params.getParameter<edm::InputTag>("pusrc"))),
          genTag_(consumes<GenEventInfoProduct>(params.getParameter<edm::InputTag>("gensrc"))){
    //produces<nanoaod::FlatTable>();
    produces<nanoaod::FlatTable>("pu");
    produces<nanoaod::FlatTable>("gen");
    }

    ~PUPtHatTableProducer() override {};

    void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override {
        edm::Handle<std::vector<PileupSummaryInfo>> puInfoHandle;
        iEvent.getByToken(puTag_, puInfoHandle);

        edm::Handle<GenEventInfoProduct> genInfoHandle;
        iEvent.getByToken(genTag_, genInfoHandle);

        // Get the pileup and gen pT hat
        float hardPtHat = genInfoHandle.isValid() ? genInfoHandle->qScale() : -1.0;

        std::vector<PileupSummaryInfo> puInfoIntime;
        std::vector<float> puPtHats;
        for (const auto& puInfo : *puInfoHandle){
            if (puInfo.getBunchCrossing() == 0){
                puInfoIntime.push_back(puInfo);
            }
        }
        if (!puInfoIntime.empty()){
            puPtHats = puInfoIntime[0].getPU_pT_hats();
        }
        std::sort(puPtHats.begin(), puPtHats.end(), std::greater<float>());
        
        // DEBUG
        //std::cout << "Number of PU pT hats: " << puPtHats.size() << ", PU pT hats: ";
        //for (const auto& ptHat : puPtHats){
        //    std::cout << ptHat << " ";
        //}

        // Fill the table
        auto puTab = std::make_unique<nanoaod::FlatTable>(puPtHats.size(), "PileupPtHat", false);
        puTab->addColumn<float>("puPtHats", puPtHats, "pT hats of the pileup interactions", 14);

        // Also add a generator pT hat table
        auto genTab = std::make_unique<nanoaod::FlatTable>(1, "GenPtHat", true);
        genTab->addColumnValue<float>("hardPtHat", hardPtHat, "pT hat of the hard interaction", 14);

        //iEvent.put(std::move(puTab));
        //iEvent.put(std::move(genTab));
        iEvent.put(std::move(puTab),  "pu");
        iEvent.put(std::move(genTab), "gen");
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
        edm::ParameterSetDescription desc;
        desc.add<edm::InputTag>("pusrc", edm::InputTag("slimmedAddPileupInfo"));
        desc.add<edm::InputTag>("gensrc", edm::InputTag("generator"));
        descriptions.addWithDefaultLabel(desc);
    }

    protected:
        const edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puTag_;
        const edm::EDGetTokenT<GenEventInfoProduct> genTag_;

};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PUPtHatTableProducer);