# Rivet
RIVET Analysis code to study single top event selection

## Installation

First, create a personal fork of the https://gitlab.cern.ch/cms-gen/Rivet repository
```bash
cmsrel CMSSW_11_1_0
cd CMSSW_11_1_0/src
cmsenv

git clone ssh://git@gitlab.cern.ch:7999/${USER}/Rivet.git
cd Rivet
git remote add cms-gen ssh://git@gitlab.cern.ch:7999/cms-gen/Rivet.git
git pull cms-gen master

#Clone the code from github
git clone --depth 1 https://github.com/WbWbX/Rivet.git TOP_NEW
cd TOP_NEW
git filter-branch --prune-empty --subdirectory-filter TOP HEAD
cd ../
mv TOP_NEW TOP

#setup and compile
source rivetSetup.sh
scram b -j
```

## running the code

Three options exist of running the code: 

1. Gridpack:
```bash
cmsRun ${CMSSW_BASE}/src/Rivet/TOP/test/runRivetAnalyzer_wbwbx_cfg.py maxEvents=50 seed=123451 gridpack=/eos/cms/store/cmst3/group/top/WbWb/gridpacks/ST_4f_w_lo_G1p0_slc7_amd64_gcc630_CMSSW_9_3_16_tarball.tar.xz
```

2. GEN file:
```bash
cmsRun ${CMSSW_BASE}/src/Rivet/TOP/test/runRivetAnalyzer_wbwbx_cfg.py yodafile=test.yoda inputFiles=PATH_TO_GEN isGEN=True
```

3. MINIAOD file:
```bash
cmsRun ${CMSSW_BASE}/src/Rivet/TOP/test/runRivetAnalyzer_wbwbx_cfg.py yodafile=test.yoda inputFiles=/store/mc/RunIISummer19UL17MiniAOD/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/230000/77CDCEDC-4D6A-614D-9D31-11A255B2012E.root isAOD=True
```