#!/bin/bash
#Running this script will produce a scipt with a name specified by the [fileToWriteTo] variable.
# That script can be executed from the same directory and will walk PISCO through a predetermined
#   set of commands.  If something goes wrong, the user can just start that script again from
#   the point of failure by deleting everything above it.
fileToWriteTo="PISCOObserveScript.bash"
echo '#!/bin/bash' > $fileToWriteTo # This will overwrite the previous content of fileToWriteTo
chmod 744 PISCOObserveScript.bash

#Astronomers should typically change this stuff:
n_objs=2
nexps=(2 2) #Number of times each object is to be observed
master_nums=(13 14) #Catalogue number of object (to be given to telescope operator)
exp_times=(200 200) #Exposure times for each object
basefitsnames=("CJ082823.94-060409.8" "CJ092719.01-105357.1") #The base name of the .fits data file, along with a suffix and a sequence number to be added by the script
obj_names=("CJ082823.94-060409.8" "CJ092719.01-105357.1") #The object names as they will appear in the header
comments=(Gladders Gladders) #Comments as they will appear in the header

#End stuff astronomers should typically change

for ((i=0; i<=$n_objs-1; i++))
do
    echo 'echo -e "\a"' >> $fileToWriteTo
    echo 'echo -e "\a"' >> $fileToWriteTo
    echo 'echo -e "\a"' >> $fileToWriteTo
    echo 'read -p "Request that telescope operator move to object number ' ${master_nums[i]} ', por favor.  Press enter when ready." ' >> $fileToWriteTo
    echo 'echo Starting observation of ' ${obj_names[i]} >> $fileToWriteTo
    bash AddSingleFieldToPISCOObserveScript.bash "$fileToWriteTo" "${nexps[i]}" "${exp_times[i]}" "${basefitsnames[i]}" "${obj_names[i]}" "${comments[i]}"
    echo 'echo Done observing ' "${obj_names[i]}" >> $fileToWriteTo
    echo 'echo "-------------------------------------" '  >> $fileToWriteTo
done
echo 'echo -e "\a" '  >> $fileToWriteTo
echo 'echo -e "\a" '  >> $fileToWriteTo
echo 'echo -e "\a" '  >> $fileToWriteTo

echo 'echo Done with all observations from this set. '  >> $fileToWriteTo
