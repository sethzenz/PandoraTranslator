from os import listdir

SZKL10dir = "/afs/cern.ch/work/s/sethzenz/public/step2_KLGun_run1"
SZPho10dir = "/afs/cern.ch/work/s/sethzenz/public/step2_PhoGun_run4"
SZKL20dir = "/afs/cern.ch/work/s/sethzenz/public/step2_KLGun20_run0"
SZPho20dir = "/afs/cern.ch/work/s/sethzenz/public/step2_PhoGun20_run0"
SZKL50dir = "/afs/cern.ch/work/s/sethzenz/public/step2_KLGun50_run0"
SZPho50dir = "/afs/cern.ch/work/s/sethzenz/public/step2_PhoGun50_run0"
SZKL100dir = "/afs/cern.ch/work/s/sethzenz/public/step2_KLGun100_run0"
SZPho100dir = "/afs/cern.ch/work/s/sethzenz/public/step2_PhoGun100_run0"


SZKL10 = ["file:%s/%s" % (SZKL10dir,x) for x in listdir(SZKL10dir) if x.count("root")]
SZPho10 = ["file:%s/%s" % (SZPho10dir,x) for x in listdir(SZPho10dir) if x.count("root")]
SZKL20 = ["file:%s/%s" % (SZKL20dir,x) for x in listdir(SZKL20dir) if x.count("root")]
SZPho20 = ["file:%s/%s" % (SZPho20dir,x) for x in listdir(SZPho20dir) if x.count("root")]
SZKL50 = ["file:%s/%s" % (SZKL50dir,x) for x in listdir(SZKL50dir) if x.count("root")]
SZPho50 = ["file:%s/%s" % (SZPho50dir,x) for x in listdir(SZPho50dir) if x.count("root")]
SZKL100 = ["file:%s/%s" % (SZKL100dir,x) for x in listdir(SZKL100dir) if x.count("root")]
SZPho100 = ["file:%s/%s" % (SZPho100dir,x) for x in listdir(SZPho100dir) if x.count("root")]


