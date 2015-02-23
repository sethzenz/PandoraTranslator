#!/bin/bash

if [ -z "$CMSSW_BASE" ]; then
    echo "you need to run \cmsenv\""
    return 1
fi  

if [ ! -f $CMSSW_BASE/src/.git/HEAD ];
then
  echo "CMSSW area appears not to be set up correctly. Have you run git cms-init before getting this repo?"
  echo "You might need to read the README carefully and start over!"
  return 1
fi

NFILES=`ls -1 ${CMSSW_BASE}/src | wc -l`
if [ ! ${NFILES} = "1" ]
then
  echo "CMSSW area appears to have extra files already. Start over and read README carefully."
  echo "You can remove this condition from the setup script if you wish, but proceed with caution!"
  return 1
fi

startDir=`pwd`
cd ${CMSSW_BASE}
git clone https://github.com/lgray/PandoraPFA.git ./PandoraPFA
git clone https://github.com/lgray/PandoraSDK.git ./PandoraPFA/PandoraSDK
git clone https://github.com/lgray/PandoraMonitoring.git ./PandoraPFA/PandoraMonitoring
git clone https://github.com/lgray/LCContent.git ./PandoraPFA/LCContent
export PANDORA_DIR=${CMSSW_BASE}/PandoraPFA
mkdir -p $PANDORA_DIR/lib
cd PandoraPFA/
MONITORING=1 make -j 9

cp ${CMSSW_BASE}/src/HGCal/PandoraTranslator/scripts/pandorapfanew_interface.xml ${CMSSW_BASE}/config/toolbox/${SCRAM_ARCH}/tools/selected
scram setup pandorapfanew_interface

TEST=`scram tool list | grep pandorapfa | wc -l`

if [ "$TEST" -eq "0" ]; then
    echo "pandora pfa was not successfullly installed :-("
    cd $startDir
    exit 1
fi

echo "you can now link to the pandora libraries through scram!"
cd $startDir