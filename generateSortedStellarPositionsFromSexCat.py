import numpy as np
import sextractorObject_py3 as SO
import sys
#import AstrometryClient 

if __name__ == "__main__":
    target_cat, target_dir, include_obj_nums = sys.argv[1:]
    include_obj_nums = int(include_obj_nums) 
    #print ('Here include_obj_nums = ' + str(include_obj_nums) )
    CatObject = SO.PySex(target_dir + target_cat)
    if include_obj_nums:
        print ('Here 1')
        CatObject.generateFileForAstrometry(include_obj_nums = include_obj_nums, file_suffix = '_positions')
    else:
        print ('Here 2' )
        CatObject.generateFileForAstrometry(include_obj_nums = include_obj_nums, file_suffix = '_positions_no_numbers')
    
