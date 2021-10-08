
InputFilePath="/eos/cms/store/group/phys_egamma/asroy/ID/Run3/egmNtuples/ntuples_Run3_ID_2021_v1/data/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/crab_job_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2021/201005_191500/0000/"

InputFileName=egmTree_50.root

OutputFile=output_PhotonTree.root

MaxEvent=1000
ReportEvery=100
EvtWeightXS=1.0

echo "Running... ./photonTreeProducer.exe $InputFilePath$InputFileName $OutputFile $MaxEvent $ReportEvery $EvtWeightXS"
./photonTreeProducer.exe $InputFilePath$InputFileName $OutputFile $MaxEvent $ReportEvery $EvtWeightXS
