from math import ceil
from os import popen,F_OK,access

energy = 20
shortSampleName = "KLGun"
stepName = "step2"
runN = 0
submit = True
inDir = "/afs/cern.ch/user/s/sethzenz/work/public/step1_%s%i_run0" % (shortSampleName,energy)
print inDir
fullName = "%s_%s%i" % (stepName,shortSampleName,energy)
basefile = "step2_DIGI_L1_DIGI2RAW.py"
cmsDir="/afs/cern.ch/work/s/sethzenz/lindsey3/CMSSW_6_2_0_SLHC23_patch2/src"
runDir="/afs/cern.ch/work/s/sethzenz/public/%s_run%i" % (fullName,runN)
result = popen("mkdir -p %s" % runDir).read()
theVomsPath = popen("voms-proxy-info -path").read().strip()

if theVomsPath.count("/tmp") or len(theVomsPath) < 6:
  raise Exception,"voms-proxy-info -path is returning %s" % theVomsPath

eventsperjob = -1
queue = "1nh"

maxjobs = 500

#filelist = [ l.strip().strip(",").strip("'") for l in open(filelistfile,'r').read().split("\n") if len(l) > 20 ]
#print filelist
#print len(filelist)

from os import listdir
filelist = sorted([ "file:%s/%s" % (inDir,l) for l in listdir(inDir) if l.count(".root") ]) # no gaurenteed mapping of numbers, sorry

njobs = len(filelist)

basetext = open(basefile,'r').read()

i = 0

for infilename in filelist:
  outfilename = "%s_%i.root" % (fullName,i)
  configname = "%s/config_%i.py" % (runDir,i)
  scriptname = "%s/run_%i.sh" % (runDir,i)
  outname = "%s/log_%i.stdout" % (runDir,i)
  errname = "%s/log_%i.stderr" % (runDir,i)
  if access(scriptname,F_OK):
    print "Skipping job %i because %s exists" % (i,scriptname)
    i += 1
    if i >= maxjobs:
      break
    continue
  newtext = basetext.replace("MAXEVENTS",str(eventsperjob)).replace("SKIPEVENTS","0").replace("OUTFILE",outfilename).replace("INPUTFILE",infilename)
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
