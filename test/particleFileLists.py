MIPs = [
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_1.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_10.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_2.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_3.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_4.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_5.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_6.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_7.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_8.root",
"/store/cmst3/group/hgcal/CMSSW/Single13_CMSSW_6_2_0_SLHC23_patch2/Events_13_10_9.root"
]

Pho5 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_5_%i.root" % i for i in range(1,101)]
Pho10 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_10_%i.root" % i for i in range(1,101)]
Pho20 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_20_%i.root" % i for i in range(1,101)]
Pho40 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_40_%i.root" % i for i in range(1,101)]
Pho50 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_50_%i.root" % i for i in range(1,101)]
Pho75 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_75_%i.root" % i for i in range(1,101)]
Pho100 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_100_%i.root" % i for i in range(1,101)]
Pho250 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_250_%i.root" % i for i in range(1,101)]
Pho500 = ["/afs/cern.ch/user/s/sethzenz/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC23_patch1/RECO-PU0/Events_22_500_%i.root" % i for i in range(1,101)]

from os import access,F_OK

for coll in [Pho5,Pho10,Pho20,Pho40,Pho50,Pho75,Pho100,Pho250,Pho500]:
  for item in coll:
    if access(item,F_OK):
      print item  
