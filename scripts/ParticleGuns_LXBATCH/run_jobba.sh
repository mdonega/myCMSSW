#!/bin/bash
#
queue='1nh'  #'8nm'
echo  submit jobs to queue $queue

for i in `seq 1 1`; 
do 
    echo Submitting job $i 
    bsub -q $queue -u mauro.donega@cern.ch -oo output_$i.txt jobba.sh $i;
done




