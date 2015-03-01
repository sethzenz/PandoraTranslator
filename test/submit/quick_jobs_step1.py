from math import ceil
from os import popen,F_OK,access

energy = 50
shortSampleName = "PhoGun"
stepName = "step1"
runN = 0
submit = True
randomSeedMult = 12345
randomSeedOffset = 42

fullName = "%s_%s%i" % (stepName,shortSampleName,energy)
basefile = "%s_cfi_GEN_SIM.py" % shortSampleName
cmsDir="/afs/cern.ch/work/s/sethzenz/lindsey3/CMSSW_6_2_0_SLHC23_patch2/src"
runDir="/afs/cern.ch/work/s/sethzenz/public/%s_run%i" % (fullName,runN)
#filelistfile = "inputlist_%s.txt" % shortSampleName
result = popen("mkdir -p %s" % runDir).read()

theVomsPath = popen("voms-proxy-info -path").read().strip()

if theVomsPath.count("/tmp") or len(theVomsPath) < 6:
  raise Exception,"voms-proxy-info -path is returning %s" % theVomsPath

eventsperjob = 100
queue = "1nh"

maxjobs = 100

#filelist = [ l.strip().strip(",").strip("'") for l in open(filelistfile,'r').read().split("\n") if len(l) > 20 ]
#print filelist
#print len(filelist)

njobs = maxjobs #len(filelist)

basetext = open(basefile,'r').read()

i = 0

#infilename = "file:/afs/cern.ch/work/s/sethzenz/lindsey3/CMSSW_6_2_0_SLHC23_patch2/src/test/Pho_step1_10k.root"
infilename = "unused"
while True:
  outfilename = "%s_%.3d.root" % (fullName,i)
  configname = "%s/config_%.3d.py" % (runDir,i)
  scriptname = "%s/run_%.3d.sh" % (runDir,i)
  outname = "%s/log_%.3d.stdout" % (runDir,i)
  errname = "%s/log_%.3d.stderr" % (runDir,i)
  if access(scriptname,F_OK):
    print "Skipping job %i because %s exists" % (i,scriptname)
    i += 1
    if i >= maxjobs:
      break
    continue
  newtext = basetext.replace("MAXEVENTS",str(eventsperjob)).replace("SKIPEVENTS",str(eventsperjob*i)).replace("OUTFILE",outfilename).replace("INPUTFILE",infilename).replace("YOURSEED",str(randomSeedOffset+randomSeedMult*i)).replace("ENERGY",str(energy))
  config = open(configname,'w')
  config.write(newtext)
  config.close()
  script = open(scriptname,'w')
  script.write("#!/bin/sh\n\n")
  script.write("export X509_USER_PROXY=%s\n\n" % theVomsPath)
  script.write("cd %s\n" % cmsDir)
  script.write("eval `scramv1 runtime -sh`\n\n")
  script.write("cd -\n")
  script.write("date\n")
  script.write("cmsRun %s\n" % configname)
  script.write("date\n")
  script.write("mv %s %s\n" % (outfilename,runDir))
  script.write("date\n")
  script.close()
  result = popen("chmod a+x %s"% scriptname).read()
  cmd = "bsub -q %s -o %s -e %s %s" % (queue,outname,errname,scriptname)
  print cmd
  if submit:
    result = popen(cmd).read()
    print result
  i += 1
  if i >= maxjobs:
    break
