
#!/bin/bash
target_dir=$1
echo target_dir is $target_dir
science_image_root=$2
filter=$3
yes_or_no_crc=$4
rm_crc_images=$5
echo yes_or_no_crc is $yes_or_no_crc
#this_dir="$(pwd)/"
proc_science_image="proc_${science_image_root//.fits}_$filter.fits"
echo Beginning rough reduction of raw PISCO images in list $proc_science_image
crc_proc_science_image="crc_$proc_science_image"
if [ $yes_or_no_crc -eq 1 ]
then
    echo Beginning cosmic ray cleaning of reduced PISCO image $proc_science_image
    python3 CosmicRayCleanPISCOData_py3.py $target_dir $proc_science_image
    echo Finished cosmic ray cleaning of reduced PISCO image $proc_science_image
    #echo "${crc_proc_science_image//.fits}.cat"
    science_catalog="${crc_proc_science_image//.fits}.cat"
    image_to_sextract=$crc_proc_science_image
else
    image_to_sextract=$crc_proc_science_image
    #cp $target_dir$proc_science_image $target_dir$crc_proc_science_image
    science_catalog="${proc_science_image//.fits}.cat"
fi

echo Beginning source extractor analysis of image $image_to_sextract
echo Looking for catalog $target_dir$science_catalog
sex $target_dir$image_to_sextract -CATALOG_NAME $target_dir$science_catalog
echo Finished source extractor analysis of image $target_dir$image_to_sextract 


if [ $yes_or_no_crc -eq 1 ]
then
    if [ $rm_crc_images -eq 1 ]
    then
       rm $target_dir$crc_proc_science_image
    fi
fi 

echo Generating object pixel coordinate lists from cataloge $science_catalog 
#python3 generateSortedStellarPositionsFromSexCat.py $science_catalog $target_dir 0
python3 generateSortedStellarPositionsFromSexCat.py $science_catalog $target_dir 1 
