import MasterCatalogueObject as MCO 
import sys
import os.path 
import numpy as np

def measureColorGradient(filters, master_cat_files_for_fit, master_cat_directories_for_fit, color_gradient_file_name, stack_save_dir, stack_save_prefix, fit_degrees, extra_cat_prefix = 'crc_'):

    stackedObjectPositions = [0 for i in range(len(master_cat_files_for_fit))]
    #color_gradient_file_name = '/Users/sasha/Documents/Harvard/physics/stubbs/PISCO/stackedColorSlopes/ColorGradientFromStack1.npy'

    #fit_degrees = [1 for i in range(16)] + [2]

    pos_ref_filter = filters[0] 

    print ('filters to measure difference between = ' + str(filters))
    print ('master catalogue files from which we will measured 2d fit = ' + str(master_cat_files_for_fit))
    print ('directories where master catalogues will be found = ' + str(master_cat_directories_for_fit))
    print ('color gradient will be saved to ' + color_gradient_file_name)
    print ('stack_save_dir = ' + str(stack_save_dir))
    print ('stack_save_prefix = ' + str(stack_save_prefix)) 
    print ('fit_degress = ' + str(fit_degrees) ) 

    for i in range(len(master_cat_files_for_fit)):
        target_dir = master_cat_directories_for_fit[i] 
        master_cat_file = master_cat_files_for_fit[i]
        stack_file = stack_save_prefix + master_cat_file[0:-5] + '.npy'
        print ('Looking for file ' + stack_save_dir + stack_file + ' ... ') 
        if os.path.isfile(stack_save_dir + stack_file):
            print ('Found it.  Loading stack from file ' + stack_save_dir + stack_file)
            stackedObjectPositions[i] = np.load(stack_save_dir + stack_file)  
        else: 
            print ('Working on master catalogue ' + target_dir + master_cat_file)
            print ('Loading master catalogue...') 
            master_cat = MCO.MasterCatalogue(target_dir, master_cat_file, extra_cat_prefix = extra_cat_prefix)
            print ('Measuring color slopes...') 
            obj_position_functions = master_cat.FitColor(filt_diffs = filters, pos_ref_filter = pos_ref_filter, fit_degree = fit_degrees[i])
            print ('Stacking measured slopes (which we will save) ...') 
            stackedObjectPositions[i] = MCO.stackObjectPositions([obj_position_functions], filt_diffs = filters, pos_ref_filter = pos_ref_filter, save = 1, save_dir = stack_save_dir, save_file_name = stack_file)

    print ('Measuring color gradient from single-image stacks...') 
    colorGradient = MCO.measureColorGradientFromDirectionalDerivatives([stack[0] for stack in stackedObjectPositions], [stack[1] for stack in stackedObjectPositions], [stack[2] for stack in stackedObjectPositions])
    np.save(stack_save_dir + color_gradient_file_name, colorGradient)
    print ('Measuring 2d color gradient function ...') 
    best_fit_color_gradient = MCO.fitColorGradientTo2DPolynomial(*colorGradient)

    print ('best_fit_color_gradient = ' + str(best_fit_color_gradient))

    np.save(stack_save_dir + 'poly_fit_' + color_gradient_file_name, best_fit_color_gradient) 

    return best_fit_color_gradient 
    

#To invoke, type: python MeasurePISOCTwoDColorFit.py [filter1] [filter2] [MASTERCatalogue_a] [MASTERCatalogue_b] [MASTERCatalogue_c] ... [MASTERCatalogue_z] [ColorGradientFileName]
# Note that [MASTERCatalogue_a] etc. must include the FULL PATH to the catalogue (that is because the master catalogues are often stored in different directories
if __name__=="__main__":

    print ('sys.argv = ' + str(sys.argv))
    print ('len(sys.argv) = ' + str(len(sys.argv)) )
    filters = list(map(str, sys.argv[1].strip('[]').split(',')))
    filters = [filt.strip() for filt in filters]
    n_catalogues = len(sys.argv[3:-1]) // 2
    print ('n_catalogues = ' + str(n_catalogues))
    master_cat_files_for_fit = list(map(str, sys.argv[2].strip('[]').split(',')))
    master_cat_files_fot_fit = [mcat_file.strip() for mcat_file in master_cat_files_for_fit ]
    master_cat_directories_for_fit = list(map(str, sys.argv[3].strip('[]').split(',')))
    master_cat_directories_for_fit = [mcat_dir.strip() for mcat_dir in master_cat_directories_for_fit  ]
    color_gradient_file_name = sys.argv[4]
    stack_save_dir = sys.argv[5]
    stack_save_prefix = sys.argv[6] 
    if len (sys.argv) > 7: 
        fit_degrees = list(map(int, sys.argv[7].strip('[]').split(',')))
    else:
        fit_degrees = [1 for i in range(16)] + [2]
    #master_cat_files_for_fit = sys.argv[3:3+n_catalogues]
    #master_cat_directories_for_fit = sys.argv[3+n_catalogues:-1]
    #color_gradient_file_name = sys.argv[-1] 
    #cat_file_suffix = '.mcat'

    best_fit_color_gradient = measureColorGradient(filters, master_cat_files_for_fit, master_cat_directories_for_fit, color_gradient_file_name, stack_save_dir, stack_save_prefix, fit_degrees) 
    

    print ('best_fit_color_gradient = ' + str(best_fit_color_gradient)) 
        
