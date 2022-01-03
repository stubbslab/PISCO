#!/bin/bash

nexps=1

exp=$1
basefitsname=$2
objname=$3
comment=$4
#basefitsname='ThisIsBaseFitsName'

date=ut190412 

filenamesuffixes=('A') #the file name suffixes that will distinguish the file names 
raoffsets=() #the RA offsets 
decoffsets=() #The Decl offsets 
for (( c=0; c<=$nexps-1; c++ ))
do
    #command="source ~/.tcshrc; checkifreadoutisdone" 
    #while true; do
    #    readoutdone=$(ssh pisco@piscotest2 $command)
    #    if [[ $readoutdone -eq 1 ]]; then break; fi
    #    echo "waiting for previous exposure to readout"
    #    sleep 1
    #done
    echo "Starting exposure" $c

    fitsname="$basefitsname"_"${filenamesuffixes[c]}"
    echo $fitsname

    #command="source ~/.tcshrc; cd /data1/obs/pisco/'$date'/; exp_acq '$exp' fname='$fitsname' objname='$objname' comment='$comment' VERBOSE=0 kickonread=1"

    #ssh pisco@piscotest2 $command

    echo "done with exposure "$c 
    if (( $c < $nexps-1 ))
    then
        echo "offsetting telescope by " "${raoffsets[c]}" "${decoffsets[c]}"
    #	tcscoffset magellan2 "${raoffsets[c]}" "${decoffsets[c]}"
    fi 
    
done

echo "done with final exposure, reading out, can slew now."
