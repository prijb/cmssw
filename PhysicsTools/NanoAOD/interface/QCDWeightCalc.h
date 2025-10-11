#ifndef PHYSICSTOOLS_NANOAOD_QCDWEIGHTCALC_H
#define PHYSICSTOOLS_NANOAOD_QCDWEIGHTCALC_H

#include <string>
#include <vector>   
#include <boost/property_tree/ptree.hpp>
#include "FWCore/ParameterSet/interface/FileInPath.h"

class QCDWeightCalc {
public:
  struct PtBinnedSample {
    float minPt;
    float maxPt;
    float xsec;
    float nrIncl;
    float nrEm;
    float emFiltEff;
    float emMuFiltEff;
    float nrMu;
    float muFiltEff;
    float muEmFiltEff;
    float nrEmNoMuExpect;
    float nrEmNoMuActual;
    float nrMuNoEmExpect;
    float nrMuNoEmActual;
    float nrEmMuExpect;
    float nrEmMuActual;
  
    
    PtBinnedSample(const boost::property_tree::ptree& sampleInfo);
    void setEnrichedCounts(float nrMinBias,float minBiasXSec);
  };

  QCDWeightCalc(const edm::FileInPath& inputFile,float bxFreq=30E6);
  float weight(float genPtHat,const std::vector<float>& puPtHats,bool passEm,bool passMu)const;
  float filtWeight(float genPtHat,bool passEm,bool passMu)const;
  float operator()(float genPtHat,const std::vector<float>& puPtHats,bool passEm,bool passMu)const{
    return weight(genPtHat,puPtHats,passEm,passMu);
  }
private:
  size_t getBinNr(float ptHat)const;
  float bxFreq_;
  std::vector<PtBinnedSample> bins_;
  
};

#endif 

