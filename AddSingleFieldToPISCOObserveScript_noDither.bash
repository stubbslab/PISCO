#!/bin/bash

repeat_dither=1 #The number of observations after which the dither will repeat

fileToWriteTo=$1
nexps=$2
exp=$3
basefitsname=$4
objname=$5
comment=$6
#basefitsname='ThisIsBaseFitsName'

date=ut191001

filenamesuffixes=("A") #the file name suffixes that will distinguish the file names
raoffsets=(0) #the RA offsets
decoffsets=(0) #The Decl offsets
total_ra_offset=0
total_dec_offset=0

echo "nexps is: " $nexps
echo "exp is: " $exp
echo "basefitsname is: " $basefitsname
echo "objname is: " $objname
echo "comment is: " $comment
for (( c=0; c<=$nexps-1; c++ ))
do
    echo 'command="source ~/.tcshrc; checkifreadoutisdone"' >> $fileToWriteTo
    echo 'while true; do' >> $fileToWriteTo
    echo '    readoutdone=$(ssh pisco@piscotest2 $command)'>> $fileToWriteTo
    echo '    if [[ $readoutdone -eq 1 ]]; then break; fi'>> $fileToWriteTo
    echo '    echo "waiting for previous exposure to readout"'>> $fileToWriteTo
#    echo '    break'>> $fileToWriteTo
    echo '    sleep 2'>> $fileToWriteTo
    echo 'done'>> $fileToWriteTo
    echo 'echo "Pausing 10s for escape"' >> $fileToWriteTo
    echo 'sleep 10' >> $fileToWriteTo
    echo 'echo Starting exposure ' $(( c + 1)) >> $fileToWriteTo

    fitsname=$basefitsname"_"${filenamesuffixes[$(($c % $repeat_dither))]}
    echo 'echo fitsname is '$fitsname >> $fileToWriteTo
    echo 'command="source ~/.tcshrc; cd /data1/obs/pisco/'$date'/; exp_acq '$exp' fname='$fitsname' objname='$objname' comment='$comment' VERBOSE=0 kickonread=1" ' >> $fileToWriteTo

    echo 'ssh pisco@piscotest2 $command' >> $fileToWriteTo

    echo 'echo done with exposure '$(( c + 1)) >> $fileToWriteTo
    if (( $c < $nexps-1 ))
    then
        total_ra_offset=$(($total_ra_offset + ${raoffsets[$(($c % $repeat_dither))]}))
        total_dec_offset=$(($total_dec_offset + ${decoffsets[$(($c % $repeat_dither))]}))
        echo 'echo offsetting telescope by ' "${raoffsets[$(($c % $repeat_dither))]}" "${decoffsets[$(($c % $repeat_dither))]}"  >> $fileToWriteTo
        echo 'tcscoffset magellan2 ' "${raoffsets[$(($c % $repeat_dither))]}" "${decoffsets[$(($c % $repeat_dither))]}" >> $fileToWriteTo
	    echo 'echo Pausing for 62 seconds to allow partial readout' >> $fileToWriteTo
        echo 'sleep 62' >> $fileToWriteTo #pause for a length to allow partial readout
    fi 
    

done
echo 'echo Recentering to undo dithering' >> $fileToWriteTo
echo 'echo tcscoffset magellan2 ' $(($total_ra_offset * -1)) ' ' $(($total_dec_offset * -1)) ' ' >> $fileToWriteTo
echo 'echo "done with final exposure, reading out, can slew now."' >> $fileToWriteTo
