from CRABClient.UserUtilities import config
import sys
config = config()

submitVersion = "ntuples_Run3_ID_2023_v1"

mainOutputDir = '/store/group/phys_egamma/asroy/ID/Run3/egmNtuples/%s' % submitVersion

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = '/afs/cern.ch/work/a/asroy/public/EGammaWork/ID/EGMTree/CMSSW_11_1_0_pre8/src/EGMObjectDumper/egmNtuplizer/test/Run3_ConfFile_cfg.py'

config.JobType.sendExternalFolder = True
config.JobType.allowUndistributedCMSSW = True

config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 
config.Site.storageSite = 'T2_CH_CERN'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_%s' % submitVersion

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    ##### now submit DATA
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
    config.Data.splitting     = 'LumiBased'
    config.Data.totalUnits      = -1
#    config.Data.lumiMask      = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt' #UL2018  
    config.Data.unitsPerJob   = 40
#    config.Data.unitsPerJob   = 200


    ### For different DataSets
#    config.General.requestName  = 'EGamma_UL2018A'
#    config.Data.inputDataset    = '/EGamma/Run2018A-12Nov2019_UL2018-v2/AOD'
#    submit(config)


    config.General.requestName = 'job_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2023'
    config.Data.inputDataset   = '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/Run3Summer19MiniAOD-2023Scenario_106X_mcRun3_2023_realistic_v3-v1/MINIAODSIM'
    submit(config)    

    config.General.requestName = 'job_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2023'
    config.Data.inputDataset   = '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/Run3Summer19MiniAOD-2023Scenario_106X_mcRun3_2023_realistic_v3-v2/MINIAODSIM'
    submit(config)    

    config.General.requestName = 'job_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_Run3Summer19_2023'
    config.Data.inputDataset   = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_14TeV_Pythia8/Run3Summer19MiniAOD-2023Scenario_106X_mcRun3_2023_realistic_v3-v2/MINIAODSIM'
    submit(config)    

