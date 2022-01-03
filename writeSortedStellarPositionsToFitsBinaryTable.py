import sextractorObject_py3 as SO
import sys

if __name__=="__main__":
    args = sys.argv[1:]
    print ('args = ' + str(args))
    target_cat, data_dir, verbose, binary_table_suffix  = sys.argv[1:]
    verbose = int(verbose)
    #print ('Here include_obj_nums = ' + str(include_obj_nums) )
    CatObject = SO.PySex(target_dir + target_cat, verbose = verbose)

    CatObject.generateFileForAstrometry(include_obj_nums = 0, include_fluxes = 1, file_suffix = binary_table_suffix, save_as_text = 0)
