import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("FlatCalib")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("HeavyIonsAnalysis.HiEvtPlaneCalib.evtplanecalibtree_cfi")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load("CondCore.CondDB.CondDB_cfi")
process.load("HeavyIonsAnalysis.Configuration.hfCoincFilter_cff")
process.load("HeavyIonsAnalysis.Configuration.analysisFilters_cff")
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v15', '')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
process.MessageLogger.cerr.FwkReport.reportEvery=1000

import FWCore.PythonUtilities.LumiList as LumiList
goodLumiSecs = LumiList.LumiList(filename = 'Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt').getCMSSWString().split(',')

#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring() 

#process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring())

process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/065400BB-5AB1-E611-A4F6-FA163E7F116F.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/0A4989E3-5AB1-E611-A792-02163E012A64.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/2649FC1A-53B1-E611-A563-FA163E8D77E4.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/2C089730-51B1-E611-9E7F-FA163ECE0FC7.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/30C43750-4EB1-E611-99EB-02163E01246E.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/3E424DD4-5AB1-E611-8A1F-02163E013451.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/40AC8EC7-5AB1-E611-AE5B-02163E0145D8.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/46A6A5C9-5AB1-E611-95E4-FA163ED7F38A.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/4A9B46EB-4EB1-E611-9727-02163E01475A.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/4C764BAE-4FB1-E611-ADEF-FA163E921868.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/521FDCBE-4FB1-E611-9199-02163E01410C.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/6272E4D6-4EB1-E611-981F-FA163E4CA5AC.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/62D5C1E6-4EB1-E611-8D3D-FA163E044471.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/64F6542E-51B1-E611-8E0C-FA163ED87989.root',
        'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/538/00000/686CE501-50B1-E611-A169-FA163EDC5F0C.root'
        ),
                             inputCommands=cms.untracked.vstring(
        'keep *',
        'drop *_hiEvtPlane_*_*'
        )
                             )



process.TFileService = cms.Service("TFileService",
    fileName = cms.string("calib.root")
)

from HeavyIonsAnalysis.Configuration.collisionEventSelection_cff import *


import HLTrigger.HLTfilters.hltHighLevel_cfi
process.minBias = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.minBias.HLTPaths = [
                "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_*_v*"
]

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.hiEvtPlane.loadDB = cms.bool(False)
process.p = cms.Path(process.collisionEventSelectionPA*process.minBias*process.hfCoincFilter3*process.hiEvtPlane* process.evtPlaneCalibTree )



                        

