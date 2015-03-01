from ROOT import *

inFile = "hamster_Pho10.root"
outDir = "/afs/cern.ch/user/s/sethzenz/www/cms/quick_28Feb/Pho10"

from os import popen
cmd = "mkdir -p %s; cp ~/www/cms/JetPlots_Feb28/.htaccess %s" % (outDir,outDir)
result = popen(cmd).read()
print result

c1 = TCanvas()

tf = TFile(inFile)
tt = tf.Get("MyTree")
for b in tt.GetListOfBranches():
  name = b.GetName()
  tt.Draw(name,"abs(%s)<999999999"% name)
  c1.SaveAs("%s/%s.png" % (outDir,name))
  
