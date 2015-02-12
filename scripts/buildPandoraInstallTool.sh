#!/bin/bash

if [ -z "$CMSSW_BASE" ]; then
    echo "you need to run \cmsenv\""
    exit 1
fi  

startDir=`pwd`
cd ${CMSSW_BASE}/src
git clone https://github.com/lgray/PandoraPFA.git ./PandoraPFA
git clone https://github.com/lgray/PandoraSDK.git ./PandoraPFA/PandoraSDK
git clone https://github.com/lgray/PandoraMonitoring.git ./PandoraPFA/PandoraMonitoring
git clone https://github.com/lgray/LCContent.git ./PandoraPFA/LCContent
export PANDORA_DIR=${CMSSW_BASE}/src/PandoraPFA
mkdir -p $PANDORA_DIR/lib
cd PandoraPFA/
make -j 9

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