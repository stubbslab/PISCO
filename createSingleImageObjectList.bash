#!/bin/bash
#This is a comment
#Assumes that astrometry has been done 
echo Beginning process of creating single exposure catalogue 
science_image=$1
target_dir=$2
#science_images_list=$(cat $target_dir$science_list_file)
single_master_cat_file_name=$3
unified_cat_file_name=$4
echo unified_cat_file_names is $unified_cat_file_name
crc_corrected=${5:-1}

echo crc_corrected is $crc_corrected 

declare -a filters=('g' 'r' 'i' 'z')
ref_filter=$filters
ref_image=$science_images_list 
match_tol_in_deg=0.001 #tolerance, in pixels 
colmerge_dir=/Users/sasha/Documents/colmerge/colmerge2

echo Working on catalogues of $science_image

for filter in "${filters[@]}"
do
    if [ $crc_corrected -eq 1 ]
    then
	file_prefix=crc_proc_
    else
	file_prefix=proc_
    fi
    image_file_prefix=proc_
    #target_cat="$file_prefix${science_image//.fits}_$filter.cat"
    #target_positions_noNum="$file_prefix${science_image//.fits}_"$filter"_positions_no_numbers.txt"
    target_positions="$file_prefix${science_image//.fits}_"$filter"_positions.txt"
    ref_wcs_fits_prefix=wcs_proc_4C16.49_122_
    ref_wcs_fits=$ref_wcs_fits_prefix$filter.fits
    target_fake_wcs_cols="$file_prefix${science_image//.fits}_$filter"_wcs.txt
    target_image="$image_file_prefix${science_image//.fits}_$filter".fits
    ref_dir=/Users/sasha/Desktop/PISCORapidAnalysisToolkit/referenceWCSSolution/

    echo About to try creating columns for creation of single image object list 
    python3 createColumnsOfSuperImageCoordinates.py $target_dir $ref_dir $target_positions $ref_wcs_fits $target_fake_wcs_cols 
    #echo Removing file $target_dir$target_image  

    #if [ "$filter" == "$ref_filter" ]
    for ref_filter in "${filters[@]}"
    do
	if [ "$ref_filter" == "$filter" ]
	then
	    echo This reference filter $ref_filter is the current filter $filter so we break.
	    break
	else
	    echo Comparing objects in this filter $filter to those in reference filter $ref_filter
	    ref_wcs_cols="$file_prefix${science_image//.fits}_$ref_filter"_wcs.txt
            colmerge_file="$file_prefix${science_image//.fits}_"$ref_filter"_wcs_VS_"$file_prefix"${science_image//.fits}_"$filter".txt"
	    $colmerge_dir/colmerge 2,3 $target_dir$ref_wcs_cols 2,3 $target_dir$target_fake_wcs_cols -tol $match_tol_in_deg > $target_dir$colmerge_file
	fi
    done 
	    

    for prev_science_image in $science_images_list
    do
	if [ "$prev_science_image" == "$science_image" ]
	then
	    echo This is science_image $prev_science_image which is the one that I am already considering.  Break the loop.
	    break
	else
	    echo This is science_image $prev_science_image which is not the one that I am already considering.
	    echo Comparing $ref_wcs_cols to $target_wcs_cols 
	    ref_wcs_cols="$file_prefix${prev_science_image//.fits}_$filter"_wcs.txt
	    colmerge_file="$file_prefix${prev_science_image//.fits}_"$filter"_wcs_VS_"$file_prefix"${science_image//.fits}_"$filter".txt"
	    $colmerge_dir/colmerge 2,3 $target_dir$ref_wcs_cols 2,3 $target_dir$target_wcs_cols -tol $match_tol_in_deg > $target_dir$colmerge_file 
	fi
    echo Finished one filter! 
    done 	  		  
done
    echo Finished one filter!

    echo Creating unified  catalogue object...

temp_file_list="${science_image//.fits}.list"
echo $science_image > $target_dir$temp_file_list
prefix_cat_file_is_missing=""
python3 createMasterObjectList.py $target_dir $temp_file_list $file_prefix $single_master_cat_file_name $prefix_cat_file_is_missing
python3 createSingleObjectUnifiedCatalog.py $target_dir $target_dir $single_master_cat_file_name $unified_cat_file_name 
rm  $target_dir$temp_file_list  
#python3 createUnifiedCatalogueFile.py $target_dir $science_image_list $file_prefix $unified_cat_file_name 
#echo Cleaning up data...
#bash cleanUpData.bash $science_list_file $target_dir $master_cat_file_name

