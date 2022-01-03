#!/bin/bash
science_image_root=$1
filter=$2
yes_or_no_crc=$3
echo yes_or_no_crc is $yes_or_no_crc
this_dir="$(pwd)/"
done_file_name=$filter"_done.txt"
touch $done_file_name
rm $done_file_name 
proc_science_image="proc_${science_image_root//.fits}_$filter.fits"
echo Beginning rough reduction of raw PISCO images in list $proc_science_image
crc_proc_science_image="crc_$proc_science_image"
if [ $yes_or_no_crc -eq 1 ]
then
    echo Beginning cosmic ray cleaning of reduced PISCO image $proc_science_image
    python2 CosmicRayCleanPISCOData_Py2.py $this_dir $proc_science_image
    echo Finished cosmic ray cleaning of reduced PISCO image $proc_science_image
    #echo "${crc_proc_science_image//.fits}.cat"
    science_catalog="${crc_proc_science_image//.fits}.cat"
    echo Beginning source extractor analysis of image $crc_proc_science_image
    sex $crc_proc_science_image -CATALOG_NAME $science_catalog
    echo Finished source extractor analysis of image $crc_proc_science_image
    rm $crc_proc_science_image
else
    cp $proc_science_image $crc_proc_science_image
    science_catalog="${proc_science_image//.fits}.cat"
    echo Beginning source extractor analysis of image $crc_proc_science_image
    sex $proc_science_image -CATALOG_NAME $science_catalog
    echo Finished source extractor analysis of image $proc_science_image
fi
touch $done_file_name

