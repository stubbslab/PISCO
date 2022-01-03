import sys 
import sextractorObject_py2 as sextractorObject
import numpy as np
import csv 

if __name__ == "__main__":
    cat_dir, cat_root = sys.argv[1:]
    #print 'args = ' + str(args) 
    #pisco_dir = '/Users/sasha/Documents/Harvard/physics/stubbs/PISCO/2018_06_20/'
    #cat_root = 'proc_SPT-CL0959-2801_35_'
    print 'Looking for catalogues of form ' + cat_dir + cat_root + '[filter].cat'

    filter_strs = ['g','r','i','z']
    #filter_strs = [image_filter]
    n_filters = len(filter_strs) 
    sex_cats = [cat_dir + cat_root + filter_str + '.cat' for filter_str in filter_strs]

    print 'Reading in sextractor catalogs...'
    PySex_objects_by_filter = [sextractorObject.PySex(cat) for cat in sex_cats]

    mag_cuts = [-np.inf, np.inf] #increase mag_cuts[0] to remove bright stars.  Decrease mag_cuts[1] to remove faint stars.  
    elong_thresh = 1.4
    star_gal_thresh = 0.7
    low_spread_funct_thresh = -0.01 
    n_brightest_for_fwhm = 15
    frac_stars_to_trim_for_fwhm = 0.2
    #0     The object has neighbors, bright and close enough to significantly bias the photometry, or bad pixels (more than 10% of the integrated area affected).
    #1     The object was originally blended with another one.
    #2     At least one pixel of the object is saturated (or very close to).
    #3     The object is truncates (to close to an image boundary).
    #4     Object's aperture data are incomplete or corrupted.
    #5     Object's isophotal data are incomplete or corrupted.
    #6     A memory overflow occurred during deblending.
    #7     A memory overflow occurred during extraction.
    flags_to_use = [2 ** flag_num for flag_num in [2, 3, 4, 5, 6, 7]]

    print 'Trimming object lists...'
    [PySex_objects.trimObjectList(cut_by_mag = 1, mag_cuts = mag_cuts, cut_by_flags = 1, flag_cuts = flags_to_use, cut_by_elong = 1, elong_thresh = elong_thresh, cut_by_star_gal = 1, star_gal_thresh = star_gal_thresh, cut_by_low_spread_funct = 0, low_spread_funct_thresh = low_spread_funct_thresh )
     for PySex_objects in PySex_objects_by_filter]
    print ('Computing FWHMs for best stars...')
    [PySex_objects.measureFWHM(n_brightest_for_fwhm, frac_to_trim = frac_stars_to_trim_for_fwhm) for PySex_objects in PySex_objects_by_filter]
    fwhms = [PySex_objects.fwhm_med for PySex_objects in PySex_objects_by_filter]
    stds = [PySex_objects.fwhm_std for PySex_objects in PySex_objects_by_filter]

    print 'Generating ellipse region files...'
    [PySex_objects_by_filter[i].generateRegionFile(file_name = cat_dir + cat_root + filter_strs[i] + '_ellipse.reg', color = 'green',  region_type = 'ellipse') for i in range(n_filters)]
    print 'Generating circle region files...'
    [PySex_objects_by_filter[i].generateRegionFile(file_name = cat_dir + cat_root + filter_strs[i] + '_circle.reg', color = 'red',  region_type = 'circle') for i in range(n_filters)]
    print ('Generating FWHM star region files...')
    [PySex_objects_by_filter[i].generateRegionFile(file_name = cat_dir + cat_root + filter_strs[i] + '_fwhm_sources.reg', color = 'blue',  region_type = 'fwhm') for i in range(n_filters)]
    print 'Generating ellipticity plots...'
    sextractorObject.makeGroupedEllipticityPlots(PySex_objects_by_filter, show = 0, save = 1, save_dir = cat_dir, file_name = cat_root + 'ellipticity_contoursAndArrows_plots.png', figsize = [1500.0 / 100.0, 2400.0 / 100.0], plot_contours = 1, plot_arrows = 1)

    pixel_to_arcsec = 0.2
    fwhms = [fwhm * 0.2 for fwhm in fwhms]
    stds = [std * 0.2 for std in stds]

    print 'fwhms = ' + str(fwhms)
    print ('stds = ' + str(stds))

    print ('fwhms = ' + str(fwhms))
    print ('stds = ' + str(stds))
    fwhms_and_stds = []
    for i in range(len(fwhms)):
        fwhms_and_stds = fwhms_and_stds + [fwhms[i], stds[i]]
                                                                                                  
     
    inputs = open("seeing.csv")
    all_lines = inputs.readlines()
    last_line = all_lines[-1]
    last_line = [float(elem) for elem in last_line.split(',')]
    all_lines.pop(len(all_lines)-1)  # removes last line
    inputs.close()
    print 'last_line = ' + str(last_line)
    new_last_line = last_line + fwhms_and_stds
    print 'new_last_line = ' + str(new_last_line)
    with open('seeing.csv', "w") as out:
        for line in all_lines:
            out.write(line.strip() + "\n")
    
    with open('seeing.csv', 'a') as fd:
        fd.write(','.join([str(elem) for elem in new_last_line]) + '\n' )
            
    
    #return fwhms 
    

    
