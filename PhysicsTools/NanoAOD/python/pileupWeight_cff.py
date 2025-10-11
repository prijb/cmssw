import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.common_cff import *

pileupPtHatTable = cms.EDProducer("PUPtHatTableProducer",
    pusrc = cms.InputTag("slimmedAddPileupInfo"),
    gensrc = cms.InputTag("generator"),
)

pileupWeightTable = cms.EDProducer("PUWeightTableProducer",
    genEventInfo = cms.InputTag("generator"),
    PileupInfo = cms.InputTag("slimmedAddPileupInfo"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    qcdWeightFile = cms.FileInPath("PhysicsTools/NanoAOD/data/qcd_with_weights.json"),
    muFilter = cms.bool(False),
    emFilter = cms.bool(False),
) 

pileupWeightTask = cms.Task(pileupPtHatTable, pileupWeightTable)

def nanoPileupWeightCustomize(process):
    process.load("PhysicsTools.NanoAOD.pileupWeight_cff")
    process.nanoTableTaskCommon.add(process.pileupWeightTask)
    return process