import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import glob
#original version from https://gitlab.cern.ch/psilva/Rivet/-/blob/master/TOP/test/runRivetAnalyzer_wbwb_cfg.py

options = VarParsing.VarParsing('python')
options.register('yodafile', 'test.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of yoda output file")
options.register('isAOD', False, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "set to True if running on AOD")
options.register('isGEN', False, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "set to True if running on GEN")
options.register('gridpack', 
                 '/eos/cms/store/cmst3/group/top/WbWb/gridpacks/ST_4f_w_lo_slc7_amd64_gcc630_CMSSW_9_3_16_tarball.tar.xz',
                 VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, 
                 "gridpack to use")
options.register('seed', 123456789, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "seed to use")

options.parseArguments()
print(options)
process = cms.Process("runRivetAnalysis")

# import of standard configurations
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(250)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

if options.isGEN or options.isAOD:
  process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(options.inputFiles))
else:
  process.source = cms.Source("EmptySource")
  process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(25)

  #configure the LHE producer (for gridpack run)
  process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    nEvents = cms.untracked.uint32(options.maxEvents),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh'),
    numberOfParameters = cms.uint32(1),
    args = cms.vstring(options.gridpack)
  )

  #hadronization without matching (ickkw=0 in datacards)
  from Configuration.Generator.Pythia8CommonSettings_cfi import *
  from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
  from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *
  process.generator = cms.EDFilter("Pythia8HadronizerFilter",
                                 maxEventsToPrint = cms.untracked.int32(1),
                                 pythiaPylistVerbosity = cms.untracked.int32(1),
                                 filterEfficiency = cms.untracked.double(1.0),
                                 pythiaHepMCVerbosity = cms.untracked.bool(False),
                                 comEnergy = cms.double(13000.),
                                 PythiaParameters = cms.PSet(
                                     pythia8CommonSettingsBlock,
                                     pythia8CP5SettingsBlock,
                                     pythia8PSweightsSettingsBlock,
                                     parameterSets = cms.vstring('pythia8CommonSettings',
                                                                 'pythia8CP5Settings',
                                                                 'pythia8PSweightsSettings',
                                     )
                                 )
  )
  
  
  from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
  randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
  process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=cms.untracked.int32(options.seed)
  randSvc.populate()
  print 'Seed is',process.RandomNumberGeneratorService.externalLHEProducer.initialSeed

if options.isAOD:
  process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
  process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("prunedGenParticles"),
    genEventInfo = cms.InputTag("generator"),
    signalParticlePdgIds = cms.vint32(),
  )
  process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
  
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_WbWb_internal')
process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.useGENweights = cms.bool(True)

if options.isGEN:
  process.p = cms.Path(process.rivetAnalyzer)
elif options.isAOD:
  process.p = cms.Path(process.generator*process.rivetAnalyzer)
else:
  process.p = cms.Path(process.externalLHEProducer*process.generator*process.rivetAnalyzer)


