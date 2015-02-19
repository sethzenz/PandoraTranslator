for collname in ["Pho10","Pho20","Pho40","Pho75","Pho100","Pho250","KL10","KL20","KL40","KL75","KL100","KL250"]:
  thetext = open("runPandora_cfg.py",'r').read().replace("Pho10",collname)
  thefile = open("runPandora_cfg_%s.py" % collname,'w')
  thefile.write(thetext)
  thefile.close()
  thecmd = "nohup cmsRun runPandora_cfg_%s.py >& output_%s.txt &" % (collname,collname)
  print thecmd
