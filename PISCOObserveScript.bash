#!/bin/bash
echo -e "\a"
echo -e "\a"
echo -e "\a"
read -p "Request that telescope operator move to object number  13 , por favor.  Press enter when ready." 
echo Starting observation of  CJ082823.94-060409.8
command="source ~/.tcshrc; checkifreadoutisdone"
while true; do
    readoutdone=$(ssh pisco@piscotest2 $command)
    if [[ $readoutdone -eq 1 ]]; then break; fi
    echo "waiting for previous exposure to readout"
    sleep 2
done
echo "Pausing 10s for escape"
sleep 10
echo Starting exposure  1
echo fitsname is CJ082823.94-060409.8_A
command="source ~/.tcshrc; cd /data1/obs/pisco/ut200229/; exp_acq 200 fname=CJ082823.94-060409.8_A objname=CJ082823.94-060409.8 comment=BleemField VERBOSE=0 kickonread=1" 
ssh pisco@piscotest2 $command
echo done with exposure 1
echo offsetting telescope by  10 7
tcscoffset magellan2  10 7
echo Pausing for 62 seconds to allow partial readout
sleep 62
command="source ~/.tcshrc; checkifreadoutisdone"
while true; do
    readoutdone=$(ssh pisco@piscotest2 $command)
    if [[ $readoutdone -eq 1 ]]; then break; fi
    echo "waiting for previous exposure to readout"
    sleep 2
done
echo "Pausing 10s for escape"
sleep 10
echo Starting exposure  2
echo fitsname is CJ082823.94-060409.8_B
command="source ~/.tcshrc; cd /data1/obs/pisco/ut200229/; exp_acq 200 fname=CJ082823.94-060409.8_B objname=CJ082823.94-060409.8 comment=BleemField VERBOSE=0 kickonread=1" 
ssh pisco@piscotest2 $command
echo done with exposure 2
echo Recentering to undo dithering
echo tcscoffset magellan2  -10   -7  
echo "done with final exposure, reading out, can slew now."
echo Done observing  CJ082823.94-060409.8
echo "-------------------------------------" 
echo -e "\a"
echo -e "\a"
echo -e "\a"
read -p "Request that telescope operator move to object number  14 , por favor.  Press enter when ready." 
echo Starting observation of  CJ092719.01-105357.1
command="source ~/.tcshrc; checkifreadoutisdone"
while true; do
    readoutdone=$(ssh pisco@piscotest2 $command)
    if [[ $readoutdone -eq 1 ]]; then break; fi
    echo "waiting for previous exposure to readout"
    sleep 2
done
echo "Pausing 10s for escape"
sleep 10
echo Starting exposure  1
echo fitsname is CJ092719.01-105357.1_A
command="source ~/.tcshrc; cd /data1/obs/pisco/ut200229/; exp_acq 200 fname=CJ092719.01-105357.1_A objname=CJ092719.01-105357.1 comment=BleemField VERBOSE=0 kickonread=1" 
ssh pisco@piscotest2 $command
echo done with exposure 1
echo offsetting telescope by  10 7
tcscoffset magellan2  10 7
echo Pausing for 62 seconds to allow partial readout
sleep 62
command="source ~/.tcshrc; checkifreadoutisdone"
while true; do
    readoutdone=$(ssh pisco@piscotest2 $command)
    if [[ $readoutdone -eq 1 ]]; then break; fi
    echo "waiting for previous exposure to readout"
    sleep 2
done
echo "Pausing 10s for escape"
sleep 10
echo Starting exposure  2
echo fitsname is CJ092719.01-105357.1_B
command="source ~/.tcshrc; cd /data1/obs/pisco/ut200229/; exp_acq 200 fname=CJ092719.01-105357.1_B objname=CJ092719.01-105357.1 comment=BleemField VERBOSE=0 kickonread=1" 
ssh pisco@piscotest2 $command
echo done with exposure 2
echo Recentering to undo dithering
echo tcscoffset magellan2  -10   -7  
echo "done with final exposure, reading out, can slew now."
echo Done observing  CJ092719.01-105357.1
echo "-------------------------------------" 
echo -e "\a" 
echo -e "\a" 
echo -e "\a" 
echo Done with all observations from this set. 
