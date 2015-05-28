#!/bin/bash
#
echo SWITCH OFF DISPLAY
export DISPLAY=
echo DISPLAY = $DISPLAY

myDIR='/afs/cern.ch/user/m/mdonega/work/CMSSW_7_4_0/'
cd $myDIR
echo $PWD
echo cmsenv
eval `scramv1 runtime -sh`

eos='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'

# to export the grid autentification
export X509_USER_PROXY=~/x509up_u9880

# not to overwrite the same dirs/files 
JOBID=$1
outputroot="DoublePhotonGun_pT250-5000_${JOBID}_test.root"

# move to a tmp directory on the node
cd /tmp/mdonega/
echo $PWD
mkdir work_${JOBID}/
cd work_${JOBID}/
echo $PWD

echo  cmsRun 
cmsRun $myDIR/EGM-RunIIWinter15GS-00004_1_cfg.py output=$outputroot

echo copying output and log files
myEOS='/eos/cms/store/user/donega/tmp/'
$eos cp $outputroot $myEOS/

