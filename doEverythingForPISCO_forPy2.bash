#!/bin/bash
#This is a comment
#user must still define lists of bias, flat[FILTER], science images
echo Beginning rough reduction of raw PISCO images in list $1
python PISCO_Reduce_V1.7_2x2_forPy2.py
echo Finished rough reduction of raw PISCO images in list $1
#echo Hello world1!
science_images_list=$(cat $1)
this_dir="$(pwd)/"
echo $this_dir
#echo Hello world2!
declare -a filters=('g' 'r' 'i' 'z')
for science_image in $science_images_list
do
    echo $science_image
    python updateFileWithImageNumberAndTime_py2.py $this_dir $science_image 
    for filter in "${filters[@]}"
    do
        bash doAnalysisOfSingleFilterPISCOImage_py2.bash $science_image $filter ${2:-1} &
    done
    sleep 10

    echo Waiting for analysis of g band image
    while [ ! -f $this_dir/g_done.txt ]
    do
	sleep 5
    done
    echo Done waiting for analysis of g band image! 
    echo Waiting for analysis of r band image 
    while [ ! -f $this_dir/r_done.txt ]
    do
	sleep 5
    done
    echo Done waiting for analysis of r band image! 
    echo Waiting for analysis of i band image 
    while [ ! -f $this_dir/i_done.txt ]
    do
	sleep 5
    done
    echo Done waiting for analysis of i band image! 
    echo Waiting for analysis of z band image 
    while [ ! -f $this_dir/z_done.txt ]
    do
	sleep 5
    done
    echo Done waiting for analysis of z band image! 
    
    #echo Hello world3!
    science_cat_root="crc_proc_${science_image//.fits}_"
    #echo $science_cat_root
    echo Beginning generation of region files and ellipticity plots for subimages of image $science_image
    python AnalyzePISCOData_forPy2.py $this_dir $science_cat_root
    echo Finished generation of region files and ellipticity plots for subimages of image $science_image
    echo You can now examine ellipticity and seeing profiles of this image 
done
python plotSeeingFunction.py
rm $this_dir/g_done.txt
rm $this_dir/r_done.txt
rm $this_dir/i_done.txt
rm $this_dir/z_done.txt
