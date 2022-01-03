import numpy as np
import MasterCatalogueObject as MCO
import sys 

if __name__=="__main__":
    target_dir, imgs_dir, single_master_cat_file_name, unified_cat_file_name  = sys.argv[1:]
    print ('[target_dir, imgs_dir, single_master_cat_file_name, unified_cat_file_name] = ' + str([target_dir, single_master_cat_file_name, unified_cat_file_name]))

    filters = ['g','r','i','z']
    master_cat = MCO.MasterCatalogue(target_dir, single_master_cat_file_name, imgs_dir = target_dir)

    unified_file = open(target_dir + unified_cat_file_name, 'w')
    unified_objects = sorted(list(master_cat.master_val_dict.keys())) 
    image_name = list(master_cat.master_val_dict[unified_objects[0]].keys())[0]
    #object_ordered_keys = ['g','r','i','z','z-i','r-i','g-r',']

    unified_file.write('Unified Object Number, g, g Err, r, r Err, i, i Err, z, zerr, i-z, i-z Errr, r-i r-i Err, g-r, g-r Err, g col, g row, r col, r row, i col, i row, z col, z row, g flag, r flag, i flag, z flag \n')
    for obj in unified_objects:
        obj_dict = master_cat.master_val_dict[obj][image_name]
        mags = []
        for filt in filters:
            mags = mags + obj_dict[filt] 
        colors = obj_dict['i-z'] + obj_dict['r-i'] + obj_dict['g-r'] 
        positions = obj_dict['xy_g'] + obj_dict['xy_r'] + obj_dict['xy_i'] + obj_dict['xy_z']
        flags = [obj_dict['flag'][filt] for filt in filters]
        vals_to_write = mags + colors + positions + flags
        unified_file.write(','.join([str(elem) for elem in vals_to_write]) + '\n')
 
    unified_file.close() 

    print ('Unified catalog created. ') 
