#!/bin/bash                                                                                                              

InputFilePath=/eos/cms/store/group/phys_egamma/asroy/ID/Run3/egmNtuples/ntuples_Run3_ID_2021_v1/
InputFileName=egmTree_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_14TeV_Pythia8.root
OutputFile=output_PhotonTree_GJet_Pt-20toInf_MGG-40to80.root 

MaxEvent=1000
ReportEvery=100
EvtWeightXS=1.0
                                                  
echo "inside jobsubmit"
pwd
cd $CMSSW_BASE/src
pwd
source /cvmfs/cms.cern.ch/cmsset_default.sh
echo "sourced"
cmsenv
cd ${_CONDOR_SCRATCH_DIR}
pwd
hostname
make

echo "Running... ./photonTreeProducer.exe $InputFilePath$InputFileName $OutputFile $MaxEvent $ReportEvery $EvtWeightXS"
./photonTreeProducer.exe $InputFilePath$InputFileName $OutputFile $MaxEvent $ReportEvery $EvtWeightXS
