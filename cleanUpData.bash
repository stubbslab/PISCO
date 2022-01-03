#!/bin/bash
#Assumes that data has been effectively cleaned up
echo Beginning process of cleaning extra files from analysis list $1 in directory $2
science_list_file=$1
target_dir=$2
master_catalog_name=$3

science_images_list=$(cat $target_dir$science_list_file)
#master_cat_file_name=$3
#crc_corrected=${4:-1}

declare -a filters=('g' 'r' 'i' 'z')
ref_filter=$filters
ref_image=$science_images_list
catalog_dir=catalog_data

echo Moving catalogue data to its own directry: $target_dir$catalog_dir
if [ ! -d "$target_dir$catalog_dir" ]; then
    mkdir $target_dir$catalog_dir
fi

mv $target_dir$master_catalog_name $target_dir$catalog_dir 

for science_image in $science_images_list
do
    
    echo Working on cleaning extra data catalogues of $science_image

    for filter in "${filters[@]}"
    do
	target_image="proc_${science_image//.fits}_$filter".fits
	echo Removing file $target_dir$target_image 
	rm $target_dir$target_image
	
	wcs_fits_file="crc_proc_${science_image//.fits}_$filter"_wcs.fits
	echo Removing file $target_dir$wcs_fits_file
	rm $target_dir$wcs_fits_file

	obj_comp_file_prefix="crc_proc_${science_image//.fits}_$filter"_wcs_VS_
	echo Moving comparison files beginning with $obj_comp_file_prefix to $target_dir$catalog_dir
	mv $target_dir$obj_comp_file_prefix*.txt $target_dir$catalog_dir

	wcs_text_file="crc_proc_${science_image//.fits}_$filter"_wcs.txt
	echo Removing file $target_dir$wcs_text_file
	rm $target_dir$wcs_text_file

	object_position_txt_file1="crc_proc_${science_image//.fits}_$filter"_positions.txt
	object_position_txt_file2="crc_proc_${science_image//.fits}_$filter"_positions_no_numbers.txt
	echo Removing file $target_dir$object_position_txt_file1
	rm $target_dir$object_position_txt_file1
	echo Removing file $target_dir$object_position_txt_file2
	rm $target_dir$object_position_txt_file2

	cat_file="crc_proc_${science_image//.fits}_$filter".cat
	echo Moving catalogue file $target_dir$cat_file to catalog directory $target_dir$catalog_dir/
	mv $target_dir$cat_file $target_dir$catalog_dir/

    done

done


	
