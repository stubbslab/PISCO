#!/bin/bash
#This is a comment
#user must still define lists of bias, flat[FILTER], science images
#python3 PISCO_Reduce_V1.7_2x2.py
echo Hello world1!
science_images_list=$(cat $1)
this_dir="$(pwd)/"
echo $this_dir
echo Hello world2!
declare -a filters=('g' 'r' 'i' 'z')
#declare -a filters=('g')
for science_image in $science_images_list
do
    echo $science_image
    for filter in "${filters[@]}"
    do 
        proc_science_image="proc_${science_image//.fits}_$filter.fits"
        echo $proc_science_image
        #python CosmicRayCleanPISCOData.py $this_dir $proc_science_image
        crc_proc_science_image="crc_$proc_science_image"
        echo $crc_proc_science_image
	echo "${crc_proc_science_image//.fits}.cat"
	science_catalog="${crc_proc_science_image//.fits}.cat"
        sex $crc_proc_science_image -CATALOG_NAME $science_catalog
	#echo "${science_catalog//$filter.cat}"
        #python AnalyzePISCOData.py $this_dir "${science_catalog//$filter.cat}" $filter
    done
    echo Hello world3!
    science_cat_root="crc_proc_${science_image//.fits}_"
    echo $science_cat_root
    python AnalyzePISCOData.py $this_dir $science_cat_root 
done
echo Hello world4!

