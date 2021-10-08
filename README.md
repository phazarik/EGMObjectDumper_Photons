# EGMObjectDumper_Photons

This is designed to create flat trees from MINIAOD files. The output file contains useful variables of genparticles and photons. These are used to develop the PF-Photon ID using machine learning. The output flat trees has to be converted to a flattened csv file before feeding them to the neural network. This csv-maker can be found here - https://github.com/prachurjyacern/PFPhoton-ID-Trainer/blob/master/CSV_maker/root2csv.py 

Keep the 'EGMOBjectDumper_Photons' directory inside a CMSSW environment as : ```.../CMSSW_11_1_0_pre8/src/EGMObjectDumper_Photons```

The crab submission script and the Run3 Config file are located as : ```/EGMObjectDumper_Photons/egmNtuplizer/test/CrabSubmit_egmNtupleRun3_2021.py``` and  ```/EGMObjectDumper_Photons/egmNtuplizer/test/Run3_ConfFile_cfg.py```
