import numpy as np
import math
import sys
import sextractorObject_py3 as SO 

class MasterCatalogue:

    def generateRegionFile(file_names = 'all', cat_file_addendum = '.cat',
                           obj_identifier_key_str = 'NUMBER', fwhm_key_str = 'FWHM_IMAGE',
                           x_pos_str = 'X_IMAGE', y_pos_str = 'Y_IMAGE',
                           reg_addendum = '_master_label.reg'):
        
        if file_names is 'all':
            file_names = self.image_file_names

        for file_name in file_names:
            reg_file_name = file_name + reg_addendum
            reg_file_lines = ['# Region file format: DS9 version 4.1' ,
                              'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1', 'image']
            file_master_obj_index = self.image_file_names.index(file_name)
            origPyObj = SO.PySex(self.target_dir + file_name + cat_file_addendum)

            local_obj_names = origPyObj.full_dict[obj_identifier_key_str]
            x_vals = origPyObj.full_dict[x_pos_str]
            y_vals = origPyObj.full_dict[y_pos_str]
            fwhm_vals = origPyObj.full_dict[fwhm_key_str]
            orig_obj_dict = {}
            for i in range(len(local_obj_names)):
                orig_obj_dict[local_obj_names[i]] = [x_vals[i], y_vals[i], fwhm_vals[i]]

            
                

            for master_obj in self.master_objs:
                local_obj_number = master_obj[file_master_obj_index]
                if local_obj_number >=1:
                    master_obj_number = master_obj[0]
                    orig_obj_vals = orig_obj_dict[local_obj_number]
                    text_addendum = '# text={M' + str(master_obj_number) + ', L' + str(local_obj_number) + ') } textangle=30'
                    reg_file_lines = reg_file_lines + ['circle(' + ','.join([str(val) for val in orig_obj_vals]) + ') ' + text_addendum]

            with open (target_dir + reg_file_name, 'w') as f:
                for line in reg_file_lines:
                    f.write(line + "\n")

            print ('Just wrote to file ' + target_dir + reg_file_name)

        

    def init(self, target_dir, master_cat_file_name):

        self.target_dir = target_dir
        with open(target_dir + master_cat_file_name, 'r') as f:
            lines = [ [entry for entry in line.rstrip('\n').split(' ')] for line in f]
            self.image_file_names = lines[0]
            self.master_objs = [[int(elem) for elem in master_obj] for master_obj in lines[1:]]

        
            

    
