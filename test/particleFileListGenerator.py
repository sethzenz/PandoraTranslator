# Depends on having eos mounted
# Fix the path accordingly
# I do this:
#
# eosmount -u ~/eos
# In home directory I have this softlink: store -> /afs/cern.ch/user/s/sethzenz/eos/cms/store

MIPs = [
"/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_%i.root" % i for i in range(1,101)
]

collnames = ["MIPs"]

patterns = {"Pho":"/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_%i_%s.root",
            "KL":"/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single130-FixE_CMSSW_6_2_0_SLHC23_patch2/Events_130_%i_%s.root"}
for ene in [5,10,20,40,50,75,100,125,175,250,400,500]:
  for key in patterns.keys():
    exec "%s%i = ['%s' %% i for i in range(1,101)]" % (key,ene,patterns[key]%(ene,"%i"))
    exec 'collnames += ["%s%i"]' % (key,ene)

#Pho5 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_5_%i.root" % i for i in range(1,101)]
#Pho10 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_10_%i.root" % i for i in range(1,101)]
#Pho20 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_20_%i.root" % i for i in range(1,101)]
#Pho40 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_40_%i.root" % i for i in range(1,101)]
#Pho50 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_50_%i.root" % i for i in range(1,101)]
#Pho75 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_75_%i.root" % i for i in range(1,101)]
#Pho100 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_100_%i.root" % i for i in range(1,101)]
#Pho250 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_250_%i.root" % i for i in range(1,101)]
#Pho500 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_500_%i.root" % i for i in range(1,101)]

from os import access,F_OK

outf = open("particleFileLists.py",'w')

for collname in collnames:
  exec "coll = %s" % collname
  print collname,len(coll),
  newlist = []
  outf.write('%s = [\n' % collname)
  n = 0
  for item in coll:
    if access(item,F_OK):
      outf.write('"%s",\n' % (item.replace("/afs/cern.ch/user/s/sethzenz/","/")))
      n += 1
  print n
  outf.write("]\n\n")

outf.close()    
