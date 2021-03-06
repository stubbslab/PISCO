#!/bin/bash
#This is a comment
#user must still define lists of bias, flat[FILTER], science images
echo Beginning rough reduction of raw PISCO images in list $1 in directory $2
yes_or_no_on_PISCO_process=${3:-1}
yes_or_no_crc_correction=${4:-1}
rm_crc_images=${5:-1}
echo yes_or_no_on_PISCO_process $yes_or_no_on_PISCO_process
echo yes_or_no_crc_correction $yes_or_no_crc_correction
echo rm_crc_images $rm_crc_images
if [ $yes_or_no_on_PISCO_process -eq 1 ]
    then python3 PISCO_Reduce_V1.7_2x2_forPy3.py
    echo Finished rough reduction of raw PISCO images in list $1 in directory $2 
fi
#echo Hello world1!
target_dir=$2
science_images_list=$(cat $2/$1)
echo science_images_list is $science_images_list
echo Looking for images in directory $target_dir
#echo Hello world2!
declare -a filters=('g' 'r' 'i' 'z')
seeing_file=seeing.csv
seeing_image=seeing.png
for science_image in $science_images_list
do
    echo $science_image
    python3 updateFileWithImageNumberAndTime_py3.py $target_dir $science_image 
    for filter in "${filters[@]}"
    do
        bash doAnalysisOfSingleFilterPISCOImage_py3.bash $target_dir $science_image $filter $yes_or_no_crc_correction $rm_crc_images &
    done
    sleep 10
    echo Waiting for analysis of all bands... 
    wait
    echo All bands analyzed!  Continuing... 

    #echo Waiting for analysis of g band image
    #g_done_file=g_done.txt
    #while [ ! -f $target_dir$g_done_file ]
    #do
    #	sleep 5
    #done
    #echo Done waiting for analysis of g band image!
    #rm $target_dir$g_done_file
    #echo Waiting for analysis of r band image
    #r_done_file=r_done.txt
    #while [ ! -f $target_dir$r_done_file ]
    #do
    #	sleep 5
    #done
    #echo Done waiting for analysis of r band image!
    #rm $target_dir$r_done_file
    #echo Waiting for analysis of i band image
    #i_done_file=i_done.txt
    #while [ ! -f $target_dir$i_done_file ]
    #do
    #	sleep 5
    #done
    #echo Done waiting for analysis of i band image!
    #rm $target_dir$i_done_file
    #echo Waiting for analysis of z band image
    #z_done_file=z_done.txt
    #while [ ! -f $target_dir$z_done_file  ]
    #do
    #	sleep 5
    #done
    #echo Done waiting for analysis of z band image!
    #rm $target_dir$z_done_file 
    
    #echo Hello world3!
    if [ $yes_or_no_crc_correction -eq 1 ]
    then 
        science_cat_root="crc_proc_${science_image//.fits}_"
    else
	science_cat_root="proc_${science_image//.fits}_"
    fi
    #echo $science_cat_root
    echo Beginning generation of region files and ellipticity plots for subimages of image $science_image
    python3 AnalyzePISCOData_forPy3.py $target_dir $science_cat_root $seeing_file
    echo Finished generation of region files and ellipticity plots for subimages of image $science_image
    echo You can now examine ellipticity and seeing profiles of this image
    
done
python3 plotSeeingFunction.py $seeing_file $seeing_image $target_dir

#echo Hello world4!

