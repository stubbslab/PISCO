import sextractorObject_py3 as SextractorObject
import numpy as np
import math
import matplotlib.pyplot as plt
import csv

def loadAdvancedSextractorCat(cat_file, cat_direcotry = ''):
    return -1

def findBestDither(ref_xs, ref_ys, shift_xs, shift_ys, shift_fwhms, search_radius_scaling):
    shift_range = [[min(ref_xs) - max(shift_xs), max(ref_xs) - min(shift_xs)],
                   [min(ref_ys) - max(shift_ys), max(ref_ys) - min(shift_ys)]]
    for x_shift in range(shift_range[0][0], shift_range[0][1]):
        for y_shift in range(shift_range[1][0], shift_range[1][1]):
            shifted_xs = [x + x_shift for x in shift_xs]
            shifted_ys = [y + y_shift for y in shift_ys]
            for i in range(len(shifted_xs)):
                shifted_x = shifted_xs[i]
                shifted_y = shifted_ys[i]
                shifted_sqr_dists = (np.array(ref_xs) - shifted_xs) ** 2.0 + (np.array(ref_ys) - shifted_ys) ** 2.0
                possible_matches = [i for i in range(len(shifted_sqr_dists)) if shifted_sqr_dists[i] < (search_radius_scaling * shift_fwhms[i] / 2.0) ** 2.0]
    
def approximateShiftFunction(ref_cat, shift_cat, portion_to_compute_shift, x_key_str, y_key_str, fwhm_key_str):
    ref_xs = ref_cat.trimmed_dict[x_key_str]
    ref_ys = ref_cat.trimmed_dict[y_key_str]
    ref_fwhms = ref_cat.trimmed_dict[fwhm_key_str]
    all_shift_xs = shift_cat.trimmed_dict[x_key_str]
    all_shift_ys = shift_cat.trimmed_dict[y_key_str]
    all_shift_fwhms = shift_cat.trimmed_dict[fwhm_key_str]
    full_shift_x_range = [min(all_shift_xs), max(all_shift_xs)]
    full_shift_t_range = [min(all_shift_xs), max(all_shift_xs)]
    used_shift_x_range = [full_shift_x_range[0] + (full_shift_x_range[1] - full_shift_x_range[0]) * portion_to_compute_shift[0][0],
                          full_shift_x_range[0] + (full_shift_x_range[1] - full_shift_x_range[0]) * portion_to_compute_shift[0][1]]
    used_shift_y_range = [full_shift_y_range[0] + (full_shift_y_range[1] - full_shift_y_range[0]) * portion_to_compute_shift[1][0],
                          full_shift_y_range[0] + (full_shift_y_range[1] - full_shift_y_range[0]) * portion_to_compute_shift[1][1]]
    used_shift_xs = all_shift_xs[:]
    used_shift_ys = all_shift_ys[:]
    used_shift_fwhms = all_shift_fwhms[:] 
    for i in range(len(all_shift_xs)):
        x = used_shift_xs[i]
        y = used_shift_ys[i]
        if not(x < used_shift_x_range[1] and x > used_shift_x_range[0] and y < used_shift_y_range[1] and y > used_shift_y_range[0]):
            used_shift_xs[i] = -1
            used_shift_ys[i] = -1
            used_shift_fwhms[i] = -1
    used_shift_xs = [x for x in used_shift_xs if x > 0.0]
    used_shift_ys = [y for y in used_shift_ys if y > 0.0]
    used_shift_fwhms = [fw for fw in used_shift_fwhms if fw > 0.0]

    #Dither shift around and see where it fits best.
    best_shift = findBestDither()
    

def mergeAdvancedSextractorCats(cat_files, cat_directory = '', x_key_str = 'X_IMAGE', y_key_str = 'Y_IMAGE', fwhm_key_str = 'FWHM_IMAGE', 
                                shift_functs = None, portion_to_compute_shift = [[0.1, 0.9], [0.1, 0.9]], search_radius_scaling = 4.0):

    sexCatObjs = [SextractorObject.PySex(sextractor_file = cat_directory + cat_file)
                  for cat_file in cat_files]

    if shift_functs is None:
        shift_functs = [None for cat in cat_files[1:]]
        
    #if no shift function is given, the current default is to look for the best shift for a
    # subset of stars near the center of the frame.  Then we apply a small warping on top of that.
    for i in range(1, cat_files):
        shift_funct = shift_functs[i]
        if shift_funct is None:
            ref_cat_index = i-1
            shift_cat_index = i
            shift_funct = approximateShiftFunction(sexCatObjs[ref_cat_index], sexCatObjs[shift_cat_index], portion_to_compute_shift, x_key_str, y_key_str, fwhm_key_str, search_radius_scaling)
            
    if x_shift_functs is None:
        

    
    
    return 0 


