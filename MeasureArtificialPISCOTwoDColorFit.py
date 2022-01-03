from MeasurePISCOTwoDColorFit import measureColorGradient
import sys
import MasterCatalogueObject as MCO
import sextractorObject_py3 as SO
import numpy as np
from shutil import copyfile 

def generateArtificialCatalogues(master_cat_files_to_adjust, new_master_cat_files, master_cat_directories_for_fit, filters_to_correct, poly_vars_by_filt, artitifical_cat_file_prefix, x_pos_key_str = 'X_IMAGE', y_pos_key_str = 'Y_IMAGE', mag_key_str = 'MAG_ISO'):
    
    for i in range(len(master_cat_files_to_adjust)):
        master_cat_file = master_cat_files_to_adjust[i]
        new_master_cat_file = new_master_cat_files[i]
        target_dir = master_cat_directories_for_fit[i]
        master_cat = MCO.MasterCatalogue(target_dir, master_cat_file)
        [copyfile(target_dir + existing_cat_file + '.cat', target_dir + artificial_cat_file_prefix + existing_cat_file + '.cat') for existing_cat_file in master_cat.image_file_names[1:]]
        image_roots = np.unique([file_name[0:-2] for file_name in master_cat.image_file_names[1:] ] ).tolist()
        for j in range(len(filters_to_correct)):
            filter_to_correct = filters_to_correct[j] 
            cat_files_to_adjust = [root + '_' + filter_to_correct + '.cat' for root in image_roots]
            color_scaling, a0, ax, ay, axy, ax2, ay2 = poly_vars_by_filt[j] 
            #First, I need to create new (normal) catalogues for the desired filter with the adjusted magnitudes
            for cat_file in cat_files_to_adjust:
                cat_to_adjust = SO.PySex(target_dir + cat_file)
                print ('Changing catalogue ' + cat_file + ' and will save new catalogue as ' + artitifical_cat_file_prefix + cat_file)
                #print ('[color_scaling, a0, ax, ay, axy, ax2, ay2] = ' + str([color_scaling, a0, ax, ay, axy, ax2, ay2])) 
                lines_to_write = cat_to_adjust.opening_lines[:]
                old_mags = cat_to_adjust.full_dict[mag_key_str][:]
                xs = cat_to_adjust.full_dict[x_pos_key_str][:]
                ys = cat_to_adjust.full_dict[y_pos_key_str][:]
                new_mags = np.array(old_mags) * color_scaling + a0 + ax * np.array(xs) + ay * np.array(ys) + axy * np.array(xs) * np.array(ys) + ax2 * np.array(xs) ** 2.0 + ay2 * np.array(ys) ** 2.0
                lines_to_write = lines_to_write + [' '.join([str(cat_to_adjust.full_dict[key][i]) if not(key in [mag_key_str]) else str(new_mags[i]) for key in cat_to_adjust.full_dict.keys() ]) for i in range(len(cat_to_adjust.full_dict[mag_key_str])) ]
                with open(target_dir + artificial_cat_file_prefix + cat_file, 'w') as f:
                    for line in lines_to_write:
                        f.write(line + '\n')

        #Now I need to create new symbolic links to the other catalogues.  They are unchanged

        #Now I need to make new master catalogues with the names of the catalogues of the adjusting filter adjusted
        updated_cat_files = [file_name if not (file_name + '.cat' in cat_files_to_adjust) else artificial_cat_file_prefix + file_name for file_name in master_cat.image_file_names[1:]]

        with open(target_dir + new_master_cat_file, 'w') as f:
            f.write(master_cat.image_file_names[0] + ' ' + ' '.join([artificial_cat_file_prefix + file_name for file_name in master_cat.image_file_names[1:]]) + '\n')
            for master_obj_line in master_cat.master_objs:
                f.write(' '.join([str(obj_number) for obj_number in master_obj_line]) + '\n')


        print ('Finished creating new, artificial master catalogue: ' + new_master_cat_file)

    return 1         
    

if __name__=="__main__":
    print ('sys.argv = ' + str(sys.argv))
    print ('len(sys.argv) = ' + str(len(sys.argv)) )
    #It always adjusts the second filter of the two input filters 
    filters = list(map(str, sys.argv[1].strip('[]').split(',')))
    filters = [filt.strip() for filt in filters]
    n_catalogues = len(sys.argv[3:-1]) // 2
    print ('n_catalogues = ' + str(n_catalogues))
    master_cat_files_for_fit = list(map(str, sys.argv[2].strip('[]').split(',')))
    master_cat_files_fot_fit = [mcat_file.strip() for mcat_file in master_cat_files_for_fit ]
    artificial_master_cat_files = list(map(str, sys.argv[3].strip('[]').split(',')))
    artificial_master_cat_files = [mcat_file.strip() for mcat_file in artificial_master_cat_files ]
    master_cat_directories_for_fit = list(map(str, sys.argv[4].strip('[]').split(',')))
    master_cat_directories_for_fit = [mcat_dir.strip() for mcat_dir in master_cat_directories_for_fit  ]
    color_gradient_file_name = sys.argv[5]
    color_scaling_1, a0_1, ax_1, ay_1, axy_1, ax2_1, ay2_1 = list(map(float, sys.argv[6].strip('[]').split(',')))
    color_scaling_2, a0_2, ax_2, ay_2, axy_2, ax2_2, ay2_2 = list(map(float, sys.argv[7].strip('[]').split(',')))
    stack_save_dir = sys.argv[8]
    art_stack_save_prefix = sys.argv[9]
    artificial_cat_file_prefix = sys.argv[10] 
    if len (sys.argv) > 11: 
        fit_degrees = list(map(int, sys.argv[11].strip('[]').split(',')))
    else:
        fit_degrees = [1 for i in range(16)] + [2]

    update_master_cat_files_for_fit = generateArtificialCatalogues(master_cat_files_for_fit, artificial_master_cat_files, master_cat_directories_for_fit, filters, [[color_scaling_1, a0_1, ax_1, ay_1, axy_1, ax2_1, ay2_1], [color_scaling_2, a0_2, ax_2, ay_2, axy_2, ax2_2, ay2_2]], artificial_cat_file_prefix)

    
    artificial_color_gradient = measureColorGradient(filters, artificial_master_cat_files, master_cat_directories_for_fit, color_gradient_file_name, stack_save_dir, art_stack_save_prefix, fit_degrees, extra_cat_prefix = artificial_cat_file_prefix + 'crc_')

    print ('artificial_color_gradient = ' + str(artificial_color_gradient)) 
