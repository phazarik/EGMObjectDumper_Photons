import CRABClient
from CRABClient.UserUtilities import config #getUsernameFromSiteDB
import sys
config = config()

submitVersion = "ntuples_PFID_Run3Summer21_TauGun"
mainOutputDir = '/store/user/phazarik/%s' % submitVersion

# config.General.transferLogs = False

config.General.transferOutputs = True
config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = '/afs/cern.ch/user/p/phazarik/work/Egamma/CMSSW_11_1_0_pre8/src/EGMObjectDumper/egmNtuplizer/test/Run3_ConfFile_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 
config.Site.storageSite = 'T2_IN_TIFR'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from http.client import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_%s' % submitVersion

    def submit(config):
        try:
            crabdevCommand('submit', config = config)
        except HTTPException as hte:
            print("Failed submitting task: %s" % (hte.headers))
        except ClientException as cle:
            print("Failed submitting task: %s" % (cle))


    ##### submit MC
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'mc')
    config.Data.splitting     = 'FileBased'
    config.Data.unitsPerJob   = 20
    config.Data.allowNonValidInputDataset = True
    
    '''
    samples=[
        
        ('/TauGun_Pt-15to500_14TeV-pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v2/MINIAODSIM',
         'TauGun_Pt-15to500_14TeV-pythia8')
    ]
        
    for sample in samples:
        print(sample[0])
        config.Data.inputDataset=sample[0]
        config.General.requestName=sample[1]
        submit(config)

    '''
    config.Data.inputDataset = '/TauGun_Pt-15to500_14TeV-pythia8/Run3Summer21MiniAOD-FlatPU0to70_120X_mcRun3_2021_realistic_v5-v2/MINIAODSIM'
    config.General.requestName = 'TauGun_Pt-15to500_14TeV-pythia8'
    submit(config)
    
