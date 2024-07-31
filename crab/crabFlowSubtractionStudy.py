from WMCore.Configuration import Configuration
config = Configuration()

from datetime import date
today = date.today().strftime('%Y-%m-%d')

card = 'cardJetBackground.input'
jobTag='jetBackgroundAnalysis_defaultFlow_' + today
outputFile=jobTag+'.root'
fileLocation='2'  # Locations: 0 = Purdue, 1 = CERN, 2 = Vanderbilt, 3 = Search with xrootd

inputList='pythiaHydjet2018_miniAODforest.txt'

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card, 'output='+outputFile, 'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','jetBackgroundAnalysis.tar.gz',card]
config.JobType.outputFiles = [outputFile]
config.JobType.maxJobRuntimeMin = 120
config.JobType.maxMemoryMB = 800

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'jetBackgroundHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T2_US_Vanderbilt'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

