PISCO Reduction Code - Taylor
Updated 2018/06/21

----------------------------------------
Prerequisites: 
Python 3
AstroPy (package)
NumPy (package)
----------------------------------------

----------------------------------------
To use:

1) Place script (PISCO_Reduce_V1.7_2x2.py) in a directory with relevant science, bias, and flatfield images. 
 
2) Generate simple text files containing one filename per line for:
	a) Bias
	b) Flatfield images for g-band
	c) Flatfield images for r-band
	d) Flatfield images for i-band
	e) Flatfield images for z-band
	f) Science images
The recommended method is to navigate to the directory in a terminal and call "ls bias*.fits > bias.list" or similar for each image type.  The '.list' file extension is unimportant, and any extension can be used. The names 'bias.list', 'flatg.list','flatr.list','flati.list','flatz.list','sci.list' are recommended but not required.  

3) Run the script by calling 'python PISCO_Reduce_V1.7_2x2.py' from the terminal in the same directory as the script and follow the prompts. 

4) At the prompt: 'Create combined bias? (y/n): ' input 'y' (without quotes) to produce a combined bias image.  The resulting FITS file will be named 'BIAS.fits' and will contain 4 images, one for each band. This image is necessary for the next few reduction steps, and a raw bias image from PISCO will not work.  If this step has already been run, entering 'n' (or literally anything but 'y') at the prompt will skip this step.
	a) If you answered 'y': At the prompt enter the previously generated bias list filenames with extension then press <return>. This will produce 'BIAS.list', median combining all images listed in bias.list.

5) Repeat this process for the flats, noting that 'BIAS.fits' must be present in the working directory (if you skipped the bias creation step.  If you didn't skip the bias creation step, then 'BIAS.fits' will be in the directory automatically).  This will produce 'FLAT.list', median combining all images listed in each list file for each band.

6) At the prompt: 'Enter science list filename: ' input 'sci.list' (without quotes) or whatever you chose for the science image list.  This will produce four FITS files (one for each filter) for each image in the science image list.  The output images will be named 'proc_[original filename]_[filter].fits'.  The header of each image will be the same as the original science image. 

7) The script will close automatically when completed.
----------------------------------------