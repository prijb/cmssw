#include "PhysicsTools/NanoAOD/interface/QCDWeightCalc.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <algorithm>
#include <boost/property_tree/json_parser.hpp>
#include <stdexcept>
#include <iostream>

QCDWeightCalc::PtBinnedSample::PtBinnedSample(const boost::property_tree::ptree& sampleInfo)
{
    minPt = sampleInfo.get<float>("min_pt");
    maxPt = sampleInfo.get<float>("max_pt");
    xsec = sampleInfo.get<float>("xsec");
    nrIncl = sampleInfo.get<float>("nr_inclusive",0.);
    nrEm = sampleInfo.get<float>("nr_em",0.);
    emFiltEff = sampleInfo.get<float>("em_filt_eff",0.);
    emMuFiltEff = sampleInfo.get<float>("em_mu_filt_eff",0.);
    nrMu = sampleInfo.get<float>("nr_mu",0.);
    muFiltEff = sampleInfo.get<float>("mu_filt_eff",0.);
    muEmFiltEff = sampleInfo.get<float>("mu_em_filt_eff",0.);
    nrEmNoMuExpect = 0.;
    nrEmNoMuActual = 0.;
    nrMuNoEmExpect = 0.;
    nrMuNoEmActual = 0.;
    nrEmMuExpect = 0.;
    nrEmMuActual = 0.;
    
    std::cout <<"minPt "<<minPt<<" nr em "<<nrEm<<" nrIncl "<<nrIncl<<std::endl;
}

void QCDWeightCalc::PtBinnedSample::setEnrichedCounts(float nrMinBias,float minBiasXSec)
{
    auto mbInclCount = [this,nrMinBias,minBiasXSec](const float filtEff){
        float nrMB = nrMinBias*filtEff*this->xsec/minBiasXSec;
        float nrIncl = this->nrIncl*filtEff;
        return nrMB+nrIncl;
    };
    
    nrEmNoMuExpect = mbInclCount(emFiltEff*(1-emMuFiltEff));
    nrMuNoEmExpect = mbInclCount(muFiltEff*(1-muEmFiltEff));
    nrEmMuExpect = mbInclCount(emFiltEff*emMuFiltEff);
    
    nrEmNoMuActual = nrEmNoMuExpect + nrEm*(1-emMuFiltEff);
    nrMuNoEmActual = nrMuNoEmExpect + nrMu*(1-muEmFiltEff);
    nrEmMuActual = nrEmMuExpect + nrEm*emMuFiltEff + nrMu*muEmFiltEff;
}  

  

QCDWeightCalc::QCDWeightCalc(const edm::FileInPath& inputFile,float bxFreq):
    bxFreq_(bxFreq)
{
    std::string fullPath = inputFile.fullPath();
    boost::property_tree::ptree sampleData;
    boost::property_tree::read_json(fullPath,sampleData);
    const auto& qcdSamples = sampleData.get_child("v2").get_child("qcd");
    for(const auto& sample : qcdSamples){
        bins_.push_back(PtBinnedSample(sample.second));   
    }
    std::sort(bins_.begin(),bins_.end(),[](const auto& lhs,const auto& rhs){return lhs.minPt<rhs.minPt;});
    for(size_t binNr=1;binNr<bins_.size();binNr++){
        bins_[binNr].setEnrichedCounts(bins_[0].nrIncl,bins_[0].xsec);
    }

    //small validation check
    if(!bins_.empty()){
        if(bins_[0].minPt!=0){
        throw std::runtime_error("QCDWeightCalc: Error, no MinBias bin defined (as defined by minPt=0)");
        }
        auto dups = std::adjacent_find(bins_.begin(),bins_.end(),[](const auto& lhs,const auto& rhs){return lhs.minPt==rhs.minPt;});
        if(dups!=bins_.end()){
        throw std::runtime_error("QCDWeightCalc: Error, duplicated pt hat bins, first dupped minPt is "+std::to_string(dups->minPt));
        }
    }
        
}


size_t QCDWeightCalc::getBinNr(float ptHat)const
{
    //this function is written assuming
    //1) the number of bins is small (say <10)
    //2) that the vast majority of ptHats will be small and thus in bin 0 or 1
    //also note the minbias bin has an 0 - 9999 pt range we have to take into account it 
    //overlaps the other bins, more reason to treat it as special
    for(size_t binNr=1;binNr<bins_.size();binNr++){
        if(ptHat<bins_[binNr].minPt) return binNr-1;
        else if(ptHat<bins_[binNr].maxPt) return binNr;
    }  
    //not in a QCD binned sample, return 0 for minbias
    return 0;
}

float QCDWeightCalc::weight(float genPtHat,const std::vector<float>& puPtHats,bool passEm,bool passMu)const
{
    std::vector<float> binCounts(bins_.size()+1,0);
    for(const auto& ptHat : puPtHats){
        binCounts[getBinNr(ptHat)]++;
    }
    binCounts[getBinNr(genPtHat)]++;
    const float minBiasXsec = bins_[0].xsec;
    const float totCount = puPtHats.size()+1;//+1 for the genPtHat
    float expectEventsMC = 0;

    // DEBUG
    //float expectEventsMCunweighted = 0;
    //std::cout << "Number of bins (excluding overflow):  " << bins_.size() << std::endl;
    //std::cout << "Counts in each bin: "; 
    //for (size_t i = 0; i < binCounts.size(); i++) {
    //    std::cout << binCounts[i] << " ";
    //}
    //std::cout << std::endl;
    //std::cout << "Total count: " << totCount << std::endl;
    // END DEBUG

    for(size_t binNr=0;binNr<bins_.size();binNr++){
        
        float binFrac = binCounts[binNr]/totCount;
        float theoryFrac = bins_[binNr].xsec / minBiasXsec;
        //dont correct inclusively generated sample
        float probCorr = binNr!=0 ? binFrac / theoryFrac : 1.;
        expectEventsMC += bins_[binNr].nrIncl * probCorr;

        // DEBUG
        //expectEventsMCunweighted += bins_[binNr].nrIncl;
        //std::cout << "binNr: " << binNr << " binFrac: " << binFrac << " theoryFrac: " << theoryFrac << " probCorr: " << probCorr << " expectEventsMC from bin: " << (bins_[binNr].nrIncl * probCorr) << " expectEventsMC from bin (unweighted): " << (bins_[binNr].nrIncl)  << std::endl;
        // END DEBUG
    }

    // DEBUG
    //std::cout << "expectEventsMC: " << expectEventsMC << " expectEventsMC (unweighted): " << expectEventsMCunweighted << std::endl;
    // END DEBUG
    float weight = bxFreq_ / expectEventsMC;
    if(passEm || passMu){
        weight *= filtWeight(genPtHat,passEm,passMu);
    }
    return weight;  
}

float QCDWeightCalc::filtWeight(float genPtHat,bool passEm,bool passMu)const
{
    const auto& bin = bins_[getBinNr(genPtHat)];
    if(passEm && passMu){
        return bin.nrEmMuActual!=0 ? bin.nrEmMuExpect/bin.nrEmMuActual : 1.;
    }else if(passEm && !passMu){
        return bin.nrEmNoMuActual!=0 ? bin.nrEmNoMuExpect/bin.nrEmNoMuActual : 1.;
    }else if(!passEm && passMu){
        return bin.nrMuNoEmActual!=0 ? bin.nrMuNoEmExpect/bin.nrMuNoEmActual : 1.;
    }else{
        return 1.;
    }
}
