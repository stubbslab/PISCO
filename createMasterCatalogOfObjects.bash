#!/bin/bash
#This is a comment
#user must still define lists of bias, flat[FILTER], science images
echo Beginning process of unifying catalogues in list $1 in directory $2
science_list_file=$1
target_dir=$2
science_images_list=$(cat $2$science_list_file)
crc_corrected=${3:-1}
do_astrometry=${4:-1}

echo crc_corrected is $crc_corrected 

declare -a filters=('g' 'r' 'i' 'z')
ref_filter=$filters
ref_image=$science_images_list 
match_tol_in_deg=0.001
colmerge_dir=/Users/sasha/Documents/colmerge/colmerge2
sashas_api_key=jcvzlxulnpbuksip
astrometry_dir=/Users/sasha/Documents/astrometry.net/net/client

for science_image in $science_images_list
do
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
	target_cat="$file_prefix${science_image//.fits}_$filter.cat"
	target_positions_noNum="$file_prefix${science_image//.fits}_"$filter"_positions_no_numbers.txt"
	target_positions="$file_prefix${science_image//.fits}_"$filter"_positions.txt"
	target_wcs_fits="$file_prefix${science_image//.fits}_$filter"_wcs.fits
	target_wcs_cols="$file_prefix${science_image//.fits}_$filter"_wcs.txt
	target_image="$image_file_prefix${science_image//.fits}_$filter".fits

	echo target_positions_noNum $target_positions_noNum
	echo target_positions $target_positions
	
	python3 generateSortedStellarPositionsFromSexCat.py $target_cat $target_dir 0
	python3 generateSortedStellarPositionsFromSexCat.py $target_cat $target_dir 1

	if [ $do_astrometry -eq 1 ]
	then
	    python3 $astrometry_dir/client.py --apikey=$sashas_api_key --upload=$target_dir/$target_positions_noNum --wcs=$target_dir$target_wcs_fits --scale-units='degwidth' --scale-lower=0.05 --scale-upper=0.15 --downsample=2
	    echo IIIIIIIIIIIIIIIIIIIIIII DID ASTROMETRY OF OBJECTS $target_positions_noNum IIIIIIIIIIIIIIIIIIIIIII
	fi 

	echo $target_dir/ $target_positions $target_wcs_fits $target_wcs_cols $target_image
	python3 createColumnsOfWCS.py $target_dir/ $target_positions $target_wcs_fits $target_wcs_cols $target_image

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
	        echo colmerge_file is $colmerge_file
	        $colmerge_dir/colmerge 2,3 $target_dir/$ref_wcs_cols 2,3 $target_dir/$target_wcs_cols -tol $match_tol_in_deg > $target_dir/$colmerge_file 
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
	        echo colmerge_file is $colmerge_file
	        $colmerge_dir/colmerge 2,3 $target_dir/$ref_wcs_cols 2,3 $target_dir/$target_wcs_cols -tol $match_tol_in_deg > $target_dir/$colmerge_file 
	    fi
	done 
	  		  
    done
    echo Finished one filter!
    
done
echo Finished with all! 

#python3 createMasterObjectList.py $target_dir $science_list_file $file_prefix 
