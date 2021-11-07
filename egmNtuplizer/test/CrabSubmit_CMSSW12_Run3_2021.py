from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
from copy import deepcopy
import os
 
def submit(config):
    res = crabCommand('submit', config = config)
    #save crab config for the future
    with open(config.General.workArea + "/crab_" + config.General.requestName + "/crab_config.py", "w") as fi:
        fi.write(config.pythonise_())

samples = [
    
    ('/GJet_Pt-10to40_DoubleEMEnriched_TuneCP5_14TeV_Pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v1/MINIAODSIM', 'GJet_10to40'),
    
    ('/GJet_Pt-40toInf_DoubleEMEnriched_TuneCP5_14TeV_Pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v1/MINIAODSIM', 'GJet_40toInf'),
    
    ('/QCD_Pt-80to120_EMEnriched_TuneCP5_14TeV_pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v1/MINIAODSIM', 'QCD_80to120'),
    
    ('/QCD_Pt-120to170_EMEnriched_TuneCP5_14TeV_pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v2/MINIAODSIM', 'QCD_120to170'),
    
    ('/QCD_Pt-170to300_EMEnriched_TuneCP5_14TeV_pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v3/MINIAODSIM', 'QCD_170to300'),
    
    ('/QCD_Pt-300toInf_EMEnriched_TuneCP5_14TeV_pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v2/MINIAODSIM', 'QCD_300toInf'),
    
    ('/TauGun_Pt-15to500_14TeV-pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v2/MINIAODSIM', 'TauGun_15to500')

]

if __name__ == "__main__":
    for dataset, name in samples:
 
        conf = config()
        submitVersion = "ntuple_PFID_Run3Summer21_v1"
        mainOutputDir = '/store/group/phys_egamma/phazarik/%s' % submitVersion

        conf.General.workArea = 'crab_%s' % submitVersion
        conf.General.transferOutputs = True
        conf.JobType.pluginName  = 'Analysis'

        # Name of the CMSSW confuration file
        conf.JobType.psetName  = '/afs/cern.ch/user/p/phazarik/work/Egamma/CMSSW_12_0_0/src/EGMObjectDumper_Photons/egmNtuplizer/test/Run3_ConfFile_cfg.py'
        conf.JobType.allowUndistributedCMSSW = True
        conf.Data.allowNonValidInputDataset = True
        
        conf.Data.inputDBS = 'global'
        conf.Data.publication = False
        
        #conf.Data.publishDataName =
        conf.Site.storageSite = 'T2_CH_CERN'
        
        conf.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'mc')
        conf.Data.splitting     = 'FileBased'
        conf.Data.unitsPerJob   = 20
        conf.Data.allowNonValidInputDataset = True
        

        conf.General.requestName = name
        conf.Data.inputDataset = dataset
        
        submit(conf) 
