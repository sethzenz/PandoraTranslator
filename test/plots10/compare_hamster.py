# compare_hamster.py
#
# Dr. Seth Zenz AMInstP
#
# Written at Geneva this twenty eighth day of February in the Year of Our Lord Two thousand and fifteen 
#
# This is a python file for comparing hamsters.
#

inFilePho = "../Pho10/hamster.root"
inFileKL = "../KL10/hamster.root"
outDir = "/afs/cern.ch/user/s/sethzenz/www/cms/quick_28Feb/compare10"
energy = 10.

# If less than 7.5, need to look up more from
# https://github.com/lgray/LCContent/blob/master/src/LCPlugins/LCParticleIdPlugins.cc#L21-L44
if  energy < 7.5:
    raise Exception,"Missing info, see script comment"

if energy < 15:
    mipCut = 0.3
else:
    mipCut = 0.4

if energy < 40.:
    rmsCut = 40.
else:
    rmsCut = 90.

# Not plotted???
#    m_minCosAngle(0.3f),
#     m_highRadLengths(40.f),


cutLines = {"totalElectromagneticEnergy":[], # ???
            "mipFraction":[mipCut],
            "dCosR":[0.95],
            "clusterRms":[rmsCut],
            "innerLayerRadLengths":[10.],
            "energyAboveHighRadLengthsFrac":[0.04],
            "energyAboveHighRadLengths":[], # ???
            "radial90":[40.],
            "nRadiationLengths90":[4.0,30.0],
            "showerMaxRadLengths":[0.,25.],
            "showerProfileStart":[],
            "showerProfileChi2":[]
            }

hNBminMaxNBins = {#"totalElectromagneticEnergy":[],
            "mipFraction":[68,0.,0.68],
            "dCosR":[46,0.785,1.015],
            "clusterRms":[58,0.,58],
            "innerLayerRadLengths":[98,0.,98.],
            "energyAboveHighRadLengthsFrac":[110,0.,1.1],
#            "energyAboveHighRadLengths":[],
            "radial90":[137,0.,137.],
            "nRadiationLengths90":[81,0.,81.],
            "showerMaxRadLengths":[81,0.,81.],
            "showerProfileStart":[137,0.,137.],
            "showerProfileChi2":[112,0.,1.12]
            }


from ROOT import *
from os import popen
gStyle.SetOptStat(0)
cmd = "mkdir -p %s; cp ~/www/cms/JetPlots_Feb28/.htaccess %s" % (outDir,outDir)
result = popen(cmd).read()
print result

c1 = TCanvas("c1","c1",600,500)

tfp = TFile(inFilePho)
ttp = tfp.Get("MyTree")

tfk = TFile(inFileKL)
ttk = tfk.Get("MyTree")


for b in ttp.GetListOfBranches():
  name = b.GetName()
  print name
  if not hNBminMaxNBins.has_key(name):
      continue
  phoname = "pho%s" % name
  klname = "kL%s" % name
  hp = TH1F(phoname,name,hNBminMaxNBins[name][0],hNBminMaxNBins[name][1],hNBminMaxNBins[name][2])
  ttp.Project(phoname,name)
  hk = TH1F(klname,name,hNBminMaxNBins[name][0],hNBminMaxNBins[name][1],hNBminMaxNBins[name][2])
  ttk.Project(klname,name)
  hp.SetLineColor(2)
  hp.SetLineWidth(2)
  hk.SetLineWidth(2)
  hp.Scale(1./hp.Integral())
  hk.Scale(1./hk.Integral())
  hp.Draw()
  hk.Draw("same")
  print hp.Integral()
  print hk.Integral()
#  ymin = hp.GetYaxis().GetBinLowEdge(0)
#  ymax = hp.GetYaxis().GetBinLowEdge(hp.GetYaxis().GetNbins()+1)
#  ymin = hp.GetMinimum()
  ymin = 0
  ymax = 1.1*hp.GetMaximum()
  hp.SetMaximum(ymax)
  hp.Draw()
  hk.Draw("same")
  if cutLines.has_key(name):
    for num in cutLines[name]:
      print name,num  ,ymin,ymax
      l = TLine(num,ymin,num,ymax)
      l.SetLineWidth(2)
      l.SetLineColor(kBlack)
      l.SetLineStyle(2)
      l.Draw()
  else:
      print name,"No lines!"
  
  c1.SaveAs("%s/%s.png" % (outDir,name))
  
