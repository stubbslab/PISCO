#!/bin/bash
#This is a comment
#user must still define lists of bias, flat[FILTER], science images
echo Beginning process of unifying catalogues in list $1 in directory $2
science_list_file=$1
target_dir=$2
science_images_list=$(cat $2/$science_list_file)
crc_corrected=${3:-1}
do_astrometry_local=${4:-0}
do_astrometry_API=${5:-1}

echo crc_corrected is $crc_corrected

declare -a filters=('g' 'r' 'i' 'z')
ref_filter=$filters
ref_image=$science_images_list
api_key=htmlavmvxjouwzzn # !!!!! THIS MUST BE CHANGED BY OTHER USERS !!!!!
astrometry_dir=/Users/sasha/Documents/astrometry.net/net/client ### THIS LINE SHOULD BE UPDATED TO CONFORM WITH THE USER'S ARCHITECTURE (most likely, leave the 'astrometry.net/net/client' part alone ###

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

	python generateSortedStellarPositionsFromSexCat.py $target_cat $target_dir 0
	python generateSortedStellarPositionsFromSexCat.py $target_cat $target_dir 1



	if [ $do_astrometry_API -eq 1 ]
	then
	    python $astrometry_dir/client.py --apikey=$api_key --upload=$target_dir$target_positions_noNum --wcs=$target_dir$target_wcs_fits --scale-units='degwidth' --scale-lower=0.05 --scale-upper=0.15 --downsample=2
	    echo IIIIIIIIIIIIIIIIIIIIIII DID ASTROMETRY OF OBJECTS $target_positions_noNum IIIIIIIIIIIIIIIIIIIIIII
  elif [ $do_astrometry_local -eq 1 ]
  then
      # python run local installation of astrometry.net api
      # something like: solve-field --scale-low 10 demo/apod4.xyls --width 719 --height 507 --overwrite
      python writeSortedStellarPositionsToFitsBinaryTable.py $target_cat $target_dir 0
	    echo IIIIIIIIIIIIIIIIIIIIIII Do astrometry for objects in $target_positions_noNum via local installtion  IIIIIIIIIIIIIIIIIIIIIII
	fi

    done
    echo Finished one filter!

done
echo Finished astrometry of all!

#python3 createMasterObjectList.py $target_dir $science_list_file $file_prefix
