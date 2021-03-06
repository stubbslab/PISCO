import sys
import numpy as np
import math
import sextractorObject_py3 as SO
import os.path
from cantrips import safeSortOneListByAnother
from cantrips import getElementsOfOneListNotInAnother
from cantrips import RemoveDuplicatesFromListWithOrder
import time
import MasterCatalogueObject as MCO 

def UpdateMasterCatWithUnmatchedObjs(master_cat_array, unmatched_new_objs, new_obj_index):
    #print ('new_obj_index = ' + str(new_obj_index)) 
    new_master_cat_index = master_cat_array[-1][0]
    n_cols = len(master_cat_array[0]) 
    appended_objs = [[new_master_cat_index + 1 + i]
                      + [-1 for elem in range(1, new_obj_index)]
                      + [unmatched_new_objs[i]]
                      + [-1 for elem in range(new_obj_index + 1, n_cols)] for i in range(len(unmatched_new_objs))]
    master_cat_array = master_cat_array + appended_objs

    return master_cat_array

def GetMasterIndexFromLocalIndex(image_name_index, local_index, master_cat_array, verbose = 0):
    if verbose: 
        print ('local_index = ' + str(local_index))
        print ('image_name_index = ' + str(image_name_index)) 
    for master_obj in master_cat_array[1:]:
        #if local_index == 192:
            #print ('[master_obj[0], image_name_index, master_obj[image_name_index]]  = ' + str([master_obj[0], image_name_index, master_obj[image_name_index]]))
        if master_obj[image_name_index] == local_index:
            if verbose: print ('master_obj where local index found = ' + str(master_obj) + ' Returning ' + str(master_obj[0])) 
            return master_obj[0]

    #print ('local_index ' + str(local_index) + ' for image index ' + str(image_name_index) + ', which corresponds to image ' + str(master_cat_array[0][image_name_index]) + ' does not seem to have been assigned a master object number. Returning placeholder value of -1. ')
    return -1 

def PickBestMasterObject(master_obj_matches, verbose = 0):
    if verbose: print ('master_obj_matches = ' + str(master_obj_matches)) 
    master_obj_matches.reverse()
    if verbose: print ('master_obj_matches = ' + str(master_obj_matches)) 
    master_obj_match_candidates = RemoveDuplicatesFromListWithOrder(master_obj_matches)
    if verbose: print ('master_obj_match_candidates = ' + str(master_obj_match_candidates)) 
    max_n_matches = -1
    best_master_obj = -1
    for master_obj in master_obj_match_candidates:
        n_occurences = len([elem for elem in master_obj_matches if elem == master_obj])
        if n_occurences > max_n_matches:
            max_n_matches = n_occurences
            best_master_obj = master_obj

    return best_master_obj

#priorities look like: [number of matches, index of closest match, master object number]
def GetMasterObjectPrioritiesFromNumberAndPosition(master_obj_matches, image_indeces_of_match):
    highest_priority = 0
    master_obj_priorities_dict = {} 
    for i in range(len(master_obj_matches)):
        master_obj = master_obj_matches[i] 
        if master_obj in master_obj_priorities_dict.keys():
            master_obj_priorities_dict[master_obj] = [master_obj_priorities_dict[master_obj][0] + 1, image_indeces_of_match[i], master_obj]
        else:
            master_obj_priorities_dict[master_obj] = [1, image_indeces_of_match[i], master_obj]
        if master_obj_priorities_dict[master_obj][0] > highest_priority:
            highest_priority = master_obj_priorities_dict[master_obj][0]

    master_obj_priorities = list(reversed(sorted([master_obj_priorities_dict[master_obj] for master_obj in master_obj_priorities_dict.keys()])))

    return master_obj_priorities

def UpdateMasterCatWithMatchedObjsOfSamePriority(objs_of_current_priority, included_new_objs, included_master_objs, master_cat_array, new_image_index):

    #print ('objs_of_current_priority = ' + str(objs_of_current_priority)) 
    for i in range(len(objs_of_current_priority[0])):
        local_obj = objs_of_current_priority[0][i]
        verbosity = (local_obj % 100 == 0 or local_obj == 192)
        verbosity = 0 
        master_obj = objs_of_current_priority[1][i]
        master_obj_index = master_obj
        if verbosity: print ('[local_obj, master_obj, master_obj_index, new_image_index] = ' + str([local_obj, master_obj, master_obj_index, new_image_index]))
        #verbosity = (local_obj == 192 and new_image_index == 3) or (master_obj == 154 and new_image_index == 3)
        #if verbosity: print ('[local_obj, master_obj, master_obj_index] = ' + str([local_obj, master_obj, master_obj_index]))
        if not(local_obj in included_new_objs) and not(master_obj in included_master_objs):
            if verbosity: print ('Adding local_obj ' + str(local_obj) + ' to master_cat ') 
            master_cat_array[master_obj_index][new_image_index] = local_obj
            included_new_objs = included_new_objs + [local_obj]
            included_master_objs = included_master_objs + [master_obj]
        else:
            if verbosity: print ('Not adding local_obj ' + str(local_obj) + ' to master_cat because either it was there or because the master object that it was going to match to was already included. ') 

    return [included_new_objs, included_master_objs, master_cat_array]
    
    

#def UpdateMasterCatWithMatchedObjs(obj_matches, master_cat_array, ref_image_index, new_image_index, new_obj_numbers):
def UpdateMasterCatWithMatchedObjs(all_new_obj_matches, master_cat_array, new_image_index):

    #print ('obj_matches = ' + str(obj_matches))
    #print ('[len(obj_matches[0], len(obj_matches[1]] = ' + str([len(obj_matches[0]), len(obj_matches[1])]) )
    #Sort and extend the list of matches to make the process of looking for existing entries faster
    obj_matches = [[], []]
    included_new_objs = []
    included_master_objs = [] 
    
    for key in all_new_obj_matches.keys():
        #print ('all_new_obj_matches[key] = ' + str(all_new_obj_matches[key]))
        verbosity = (key % 100 == 0 or key == 192)
        verbosity = 0
        master_obj_matches = [ GetMasterIndexFromLocalIndex(individual_obj_match[0], individual_obj_match[1], master_cat_array, verbose = verbosity) for individual_obj_match in all_new_obj_matches[key] ]
        master_obj_priorities = GetMasterObjectPrioritiesFromNumberAndPosition(master_obj_matches, [individual_obj_matches[0] for individual_obj_matches in all_new_obj_matches[key]])
            
        if -1 in master_obj_matches:
            print ('A -1 in master_obj_matches!!! That is not good. ') 
            print ('[key, all_new_obj_matches[key]] = ' + str([key, all_new_obj_matches[key]]))
            
            
        obj_matches[0] = obj_matches[0] + [key]
        #best_master_match = PickBestMasterObject(master_obj_matches, verbose = verbosity)
        obj_matches[1] = obj_matches[1] + [master_obj_priorities]

        if verbosity:
            print ('[new_image_index, key, all_new_obj_matches[key], master_obj_matches, master_obj_priorities] = ' + str([new_image_index, key, all_new_obj_matches[key], master_obj_matches, master_obj_priorities]))

    #if 192 in obj_matches[0]: print ('[obj_matches[0][obj_matches[0].index(192)], obj_matches[1][obj_matches[0].index(192)], all_new_obj_matches[192]] = ' + str([obj_matches[0][obj_matches[0].index(192)], obj_matches[1][obj_matches[0].index(192)], all_new_obj_matches[192]]))

    #print ('obj_matches = ' + str(obj_matches)) 
    objects_to_add = obj_matches
    ordered_priorities_for_all_objs = obj_matches[1]
    used_indeces = [-1 for priority in ordered_priorities_for_all_objs]

    SuperiorOrEqualPriority = lambda priority1, priority2: (priority1[0] > priority2[0] or (priority1[0]  ==  priority2[0] and priority1[1]  >=  priority2[1])) 
    prev_max_priority = [np.inf, np.inf]
    objs_still_to_match = 1 
    while objs_still_to_match:
        max_priority = [0, 0, 0]
        possible_next_priorities = [ordered_priorities_for_all_objs[i][used_indeces[i]+1] for i in range(len(used_indeces)) if used_indeces[i] + 1< len(ordered_priorities_for_all_objs[i])]
        if len(possible_next_priorities) == 0:
            #print ('Here 1')
            objs_still_to_match = 0
        else: 
            max_priority = max(possible_next_priorities)
            #print ('max_priority = ' + str(max_priority))
            objs_of_current_priority = [[], []]
            for obj_index in range(len(objects_to_add[0])):
                new_index = used_indeces[obj_index] + 1
                if new_index < len(ordered_priorities_for_all_objs[obj_index]):
                    ordered_priority = ordered_priorities_for_all_objs[obj_index][new_index]
                else:
                    ordered_priority = [0, 0, 0]
                if SuperiorOrEqualPriority(ordered_priority[0:2], max_priority[0:2]):
                    objs_of_current_priority[0] =  objs_of_current_priority[0] + [objects_to_add[0][obj_index]]
                    verbosity = (objs_of_current_priority[0][-1] % 100 == 0 or objs_of_current_priority[0][-1] == 192)
                    verbosity = 0
                    objs_of_current_priority[1] =  objs_of_current_priority[1] + [ordered_priority[2]]
                    if verbosity: print ('[objs_of_current_priority[0][-1], objs_of_current_priority[1][-1]] = ' + str([objs_of_current_priority[0][-1], objs_of_current_priority[1][-1]]))
                    used_indeces[obj_index] = new_index
        
            #print ('Adding objects of same priority to master catalogue.' )
            included_new_objs, included_master_objs, master_cat_array = UpdateMasterCatWithMatchedObjsOfSamePriority(objs_of_current_priority, included_new_objs, included_master_objs, master_cat_array, new_image_index)
            prev_max_priority = max_priority[:]                                                 
        

    
        #if local_obj == 192:
        #    print ('[local_obj, master_obj_index, new_image_index, master_cat_array[master_obj_index][new_image_index] ] = ' + str([local_obj, master_obj_index, new_image_index, master_cat_array[master_obj_index][new_image_index] ]))
        
    #sorted_obj_matches = safeSortOneListByAnother(obj_matches[0], [obj_matches[0], obj_matches[1]])
    #if len(sorted_obj_matches[0]) > 0: 
    #    max_matched_ref_obj = sorted_obj_matches[0][-1]
    #else:
    #    max_matched_ref_obj = 0 
    #extended_matched_ref_objs = list(range(1, max_matched_ref_obj + 1))
    #extended_matched_new_objs = [-1 for elem in extended_matched_ref_objs]
 
    #n_cols = len(master_cat_array[0])
    #for i in range(len(sorted_obj_matches[0])):
    #    ref_obj = sorted_obj_matches[0][i] 
    #    extended_matched_new_objs[ref_obj-1] = sorted_obj_matches[1][i]

    #for master_obj_index in range(1, len(master_cat_array)):
    #    ref_obj_identity = master_cat_array[master_obj_index][ref_image_index]
    #    #print ('ref_obj_identity = ' + str(ref_obj_identity)) 
    #    if ref_obj_identity >= 1 and ref_obj_identity <= max_matched_ref_obj:
    #        new_obj_identity = extended_matched_new_objs[ref_obj_identity-1]
    #        if new_obj_identity >=1: 
    #            #print ('Adding new_obj ' + str(new_obj_identity) + ' at master_obj_index ' + str(master_obj_index) + ' which has a ref_obj_identity of ' + str(ref_obj_identity) )
    #            master_cat_array[master_obj_index][new_image_index] = new_obj_identity

    return [included_new_objs, master_cat_array ]

def ReadObjectMatchesFromCompFile(target_dir, comp_file):
    first_match_marker = '|'
    duplicate_marker = 'D' 
    obj_number_index = 0
    x_pos_index = 1
    y_pos_index = 2
    split_index = 3
    first_obj_numbers = []
    second_obj_numbers = []
    first_obj_positions = [] 
    second_obj_positions = []
    pairs_with_duplicate_matches = []
    bad_objects = []

    with open(target_dir + comp_file, 'r') as f:
        all_lines = [ line.rstrip('\n') for line in f]
    for line in all_lines:
        
        line = [ elem for elem in line.split(' ') if not (elem in [' ', '']) ]
        first_object = line[0:split_index]
        splitter = line[split_index] 
        second_object = line[split_index+1:]
        first_obj_number = int(first_object[obj_number_index])
        first_obj_pos = [float(first_object[x_pos_index]), float(first_object[y_pos_index])] 
        second_obj_number = int(second_object[obj_number_index])
        second_obj_pos = [float(second_object[x_pos_index]), float(second_object[y_pos_index])]
        indeces_of_existing_pairs = []

        if splitter in [first_match_marker, duplicate_marker]:
            if first_obj_number in first_obj_numbers:
                
                indeces_of_existing_pairs = indeces_of_existing_pairs + [first_obj_numbers.index(first_obj_number)] 
            if second_obj_number in second_obj_numbers:
                
                indeces_of_existing_pairs = indeces_of_existing_pairs + [second_obj_numbers.index(second_obj_number)]
                
            if len(indeces_of_existing_pairs) == 0:
                first_obj_numbers = first_obj_numbers + [first_obj_number]
                first_obj_positions = first_obj_positions + [first_obj_pos]
                second_obj_numbers = second_obj_numbers + [second_obj_number]
                second_obj_positions = second_obj_positions + [second_obj_pos]
        #elif splitter is duplicate_marker:
        #    #This means that the first object in this line has already been paired.
        #    #So we add this to our list of duplicate matches, figure out which of the possible second objects
        #    # is closer to the first object, and move on with our lives.
        #    indeces_of_existing_pairs = indeces_of_existing_pairs + [first_obj_numbers.index(first_obj_number)] 
        else:
            print ('Objects separated by unknown separation indicator: ' + splitter + '  Do not know how to handle, so not including in set of matches. ') 
            bad_object = ' '.join(line)
            bad_objects = bad_objects + [bad_object]
        if len(indeces_of_existing_pairs) > 0:
            index_of_existing_pair = indeces_of_existing_pairs[0]
            new_pair = [first_obj_number, second_obj_number]
            old_pair = [first_obj_numbers[index_of_existing_pair], second_obj_numbers[index_of_existing_pair]]
            pairs_with_duplicate_matches = pairs_with_duplicate_matches + [old_pair, new_pair]
            old_sep = np.sqrt(np.sum((np.array(first_obj_positions[index_of_existing_pair]) -  np.array(second_obj_positions[index_of_existing_pair])) ** 2.0))
            new_sep = np.sqrt(np.sum((np.array(first_obj_pos) -  np.array(second_obj_pos)) ** 2.0))
            if new_sep < old_sep:
                #print ('first_obj_numbers[index_of_existing_pair] = ' + str(first_obj_numbers[index_of_existing_pair])) 
                first_obj_numbers[index_of_existing_pair] = first_obj_number
                #print ('first_obj_numbers[index_of_existing_pair] = ' + str(first_obj_numbers[index_of_existing_pair]))
                #print ('second_obj_numbers[index_of_existing_pair] = ' + str(second_obj_numbers[index_of_existing_pair])) 
                second_obj_numbers[index_of_existing_pair] = second_obj_number
                #print ('second_obj_numbers[index_of_existing_pair] = ' + str(second_obj_numbers[index_of_existing_pair])) 
                first_obj_positions[index_of_existing_pair] = first_obj_pos
                second_obj_positions[index_of_existing_pair] = second_obj_pos
                if len(indeces_of_existing_pairs) > 1:
                    index_of_second_overwritten_pair = indeces_of_existing_pairs[1]
                    first_obj_numbers = removeListElement(first_obj_numbers, index_of_second_overwritten_pair)
                    second_obj_numbers = removeListElement(second_obj_numbers, index_of_second_overwritten_pair)
                    first_obj_positions = removeListElement(first_obj_positions, index_of_second_overwritten_pair)
                    second_obj_positions = removeListElement(second_obj_positions, index_of_second_overwritten_pair)
            
            #print ('Object ' + str(bad_object) + ' of first catalogue in comp file ' + str(target_dir + comp_file) + ' is problematic (likely a duplicate match). ')

    #print ( 'We had the following pairs that contain at least one object that was matched more than once: ' + str(pairs_with_duplicate_matches) )
    #print ('We had the following objects with now known splitter key: ' + str(bad_objects) )
    
    return [first_obj_numbers, second_obj_numbers]

def addSextractorCatFilesToMasterCat(sex_cat_file_set, target_dir, master_cat_array, new_image_num, master_file_name):

    if os.path.isfile(master_file_name):
        prev_master_cat_entries = [ [entry for entry in line.rstrip('\n').split(' ')] for line in f]
        image_file_names = prev_master_cat_entries[0]
        prev_master_entries = [[int(elem) for elem in master_obj] for master_obj in prev_master_cat_entries[1:]]
    else:
        prev_master_entries = []
    #print ('master_cat_array = ' + str(master_cat_array))
    #print ('np.shape(master_cat_array) = ' + str(np.shape(master_cat_array) ))
    cat_file_part_to_cut = '.cat'
    comp_str = '_wcs_VS_'
    comp_suffix = '.txt' 
    SexPyObjs = [SO.PySex(target_dir + sex_cat_file) for sex_cat_file in sex_cat_file_set]
    n_filters = len(SexPyObjs)
    
    for filter_index in range(n_filters):
        SexPyObj = SexPyObjs[filter_index]
        sex_cat_file = sex_cat_file_set[filter_index] 
        new_object_numbers = SexPyObj.full_dict[obj_num_key_str]
        new_object_numbers = [int(num) for num in new_object_numbers]
        print ('Working on sex_cat_file = ' + str(sex_cat_file))
        if len(master_cat_array) == 1:
            master_cat_array = [master_cat_array[0]] + [[-1 for i in range(len(object_files) * 4 + 1)] for obj in new_object_numbers]
            
            for i in range(1,len(new_object_numbers)+1):
                master_cat_array[i][0] = new_object_numbers[i-1]
                master_cat_array[i][1] = new_object_numbers[i-1]
            print ('Created first entry in our master catalogue.' )
            #print ('master_cat_array = ' + str(master_cat_array))
        else:
            new_image_index = new_image_num * n_filters + filter_index + 1
            prev_images_in_filter_indeces = [nth_image * n_filters + filter_index + 1 for nth_image in range(new_image_num)]
            prev_filters_of_this_image_indeces = [new_image_num * n_filters + prev_filter_index + 1 for prev_filter_index in range(filter_index)]
            comp_image_indeces = prev_images_in_filter_indeces + prev_filters_of_this_image_indeces
            #print ('comp_image_indeces = ' + str(comp_image_indeces))
            all_new_obj_matches = {} 
            #unmatched_new_objs = new_object_numbers[:]
            for comp_image_index in comp_image_indeces:
                comp_image = master_cat_array[0][comp_image_index]
                new_image = sex_cat_file[0:-len(cat_file_part_to_cut)]
                standard_comp_file = comp_image + comp_str + new_image + comp_suffix
                inverse_comp_file = new_image + comp_str + comp_image + comp_suffix
                if os.path.isfile(target_dir + standard_comp_file):
                    print ('Found standard comp file ' + target_dir + standard_comp_file)
                    ref_image_obj_matches, new_image_obj_matches = ReadObjectMatchesFromCompFile(target_dir, standard_comp_file)
                elif os.path.isfile(target_dir + inverse_comp_file):
                    print ('Found inverse comp file ' + target_dir + inverse_comp_file)
                    new_image_obj_matches, ref_image_obj_matches = ReadObjectMatchesFromCompFile(target_dir, inverse_comp_file)
                else:
                    print ('Could not find comp file ' + target_dir + standard_comp_file + ' nor comp_file ' + target_dir + inverse_comp_file)
                    print ('Will not update master object array. ')
                    new_image_obj_matches = []
                    ref_image_obj_matches = []
                for matched_obj_index in range(len(new_image_obj_matches)):
                    new_matched_obj = new_image_obj_matches[matched_obj_index ]
                    ref_matched_obj = ref_image_obj_matches[matched_obj_index]
                    all_new_obj_matches[new_matched_obj] = (all_new_obj_matches[new_matched_obj] + [[comp_image_index, ref_matched_obj]] if new_matched_obj in all_new_obj_matches.keys() else [[comp_image_index, ref_matched_obj]])
                    
                #obj_matches = [ref_image_obj_matches, new_image_obj_matches]
                #master_cat_array = UpdateMasterCatWithMatchedObjs(obj_matches, master_cat_array, comp_image_index, new_image_index, new_object_numbers)
                #update master list of matches with these matches
                
            #unmatched_new_objs = getElementsOfOneListNotInAnother(new_object_numbers, list(all_new_obj_matches.keys()))
            #print ('all_new_obj_matches = ' + str(all_new_obj_matches)) 
            matched_new_objs, master_cat_array = UpdateMasterCatWithMatchedObjs(all_new_obj_matches, master_cat_array, new_image_index)
            unmatched_new_objs = getElementsOfOneListNotInAnother(new_object_numbers, matched_new_objs)
            #print ('The following objects went unmatched: ' + str(unmatched_new_objs))
            #print ('unmatched_new_objs = ' + str(unmatched_new_objs)) 
            master_cat_array = UpdateMasterCatWithUnmatchedObjs(master_cat_array, unmatched_new_objs, new_image_index)
            #print ('master_cat_array = ' + str(master_cat_array)) 
                
                    
                #obj_matches = ReadObjectMatchesFromCompFile(target_dir, comp_file)

        #print ('master_cat_array = ' + str(master_cat_array))

    return master_cat_array 
        

if __name__=="__main__":
    target_dir, object_list_file, file_prefix, master_cat_name, missing_cat_prefix = sys.argv[1:]
    print ('[target_dir, object_list_file, file_prefix, master_cat_name] =' + str([target_dir, object_list_file, file_prefix, master_cat_name])) 
    #master_file_name = 'ObjectCrossReferenes.mcat'

    filters = ['g','r','i','z']

    with open(target_dir + object_list_file, 'r') as f:
        object_files = [ line.rstrip('\n') for line in f]

    print ('object_files = ' + str(object_files)) 
    sextractor_master_cat = target_dir + master_cat_name

    file_name_part_to_cut = '.fits'
    master_cat_array = [[file_prefix + object_file[0:-len(file_name_part_to_cut)] + '_' + str(filt) for filt in filters] for object_file in object_files]
    #print ('master_cat_array = ' + str(master_cat_array))
    master_cat_array =  [['master_object_index'] + [j for i in master_cat_array for j in i]]
    #print ('master_cat_array = ' + str(master_cat_array))
    #print ('object_files = ' + str(object_files))

    
    sex_cat_files = [[file_prefix + object_file[0:-len(file_name_part_to_cut)] + '_' + str(filt) + '.cat'
                      for filt in filters]
                     for object_file in object_files] 

    obj_num_key_str = 'NUMBER'
    for image_num in range(len(sex_cat_files)):
        start = time.time()
        #print ('master_cat_array = ' + str(master_cat_array))
        sex_cat_file_set = sex_cat_files[image_num]                                                                                   
        master_cat_array = addSextractorCatFilesToMasterCat(sex_cat_file_set, target_dir, master_cat_array, image_num, master_cat_name)
        end = time.time()
        #SexPyObjs = [SO.PyObj(target_dir + sex_cat_file) for sex_cat_file in sex_cat_file_set]
        #[addSextractorCatObjToMasterCat(SexPyObj, sex_cat_files, target_dir, master_cat, image_num) for SexPyObj in SexPyObjs]
        print ('Took = ' + str(end - start) + 's for image_num = ' + str(image_num)) 

    with open(sextractor_master_cat, 'w') as f:
        for line in master_cat_array:
            f.write(' '.join([str(elem) for elem in line]) + "\n")
    #print ('master_cat_array = ' + str(master_cat_array))

    #I have not yet moved the objects around, so I need to specify that the directory with the images is the same as the directory with the master catalogue
    # (normally, the images sit one directory higher than the directory with the catalogues) 
    master_cat = MCO.MasterCatalogue(target_dir, master_cat_name, imgs_dir = target_dir, missing_cat_prefix, missing_cat_prefix = missing_cat_prefix )
    master_cat.generateRegionFiles() 
        
                
        
        

    

    
