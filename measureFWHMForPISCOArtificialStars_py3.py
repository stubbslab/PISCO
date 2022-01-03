import numpy as np
import csv
import sys
import sextractorObject_py3 as SO
import matplotlib.pyplot as plt 

if __name__ == "__main__":
    #cat_dir, images_list = sys.argv[1:]
    arguments = sys.argv[1:]
    cat_dir = arguments[0]
    images_list = arguments[1:]
    print ('arguments = ' + str(arguments))
    print ('images_list = ' + str(images_list)) 
    print (type(images_list) ) 
    focus_vals = [89, 90, 91, 92, 93, 94, 95]
    pix_scale = 0.1 
    filters = ['g','r','i','z']
    colors = ['green','red','black','blue']
    filter_names = ['SDSS-g','SDSS-r','SDSS-i','SDSS-z']
    plot_object_locations = [[0,0], [0,1], [1,0], [1,1]]
    plot_object_text = ["Upper left 'star'", "Lower left 'star'", "Middle 'star'", "Lower right 'star'"]
    fwhm_key_str = 'FWHM_IMAGE'
    elongation_key_str = 'ELONGATION'
    x_key_str = 'X_IMAGE'
    y_key_str = 'Y_IMAGE'

    figsize = [10, 8] 
    f, axarr = plt.subplots(2,2, figsize = figsize)
    
    count = 0
    for image in images_list:
        print ('image = ' + image)
        image_root = image[0:-len('.fits')]
        print ('image_root = ' + str(image_root)) 
        cats = ['proc_' + image_root + '_' + filter_str + '.cat' for filter_str in filters]
        print ('Creating sextractor objects to measure closest object to star positions...')
        sextractor_objects = [SO.PySex(sextractor_file = cat) for cat in cats]
        [sextractor_objects[i].generateRegionFile(file_name = cats[i][0:-len('.fits')] + filters[i] + '.reg') for i in range(len(cats))]
        target_pixels_dict = {'g':[[268.0, 3710.0], [234.0, 562.0], [1559.0, 2125.0], [2855.0, 538.0]],
                              'r':[[148.0, 3745.0], [161.0, 585.0], [1468.0, 2171.0], [2790.0, 600.0]],
                              'i':[[147.0, 3705.0], [119.0, 545.0], [1446.0, 2125.0], [2760.0, 530.0]],
                              'z':[[208.0, 3758.0], [193.0, 604.0], [1511.0, 2189.0], [2826.0, 595.0]] }

        closest_objects_dict = {}
        better_fwhms_dict = {} 

        for i in range(len(filters)):
            filter_str = filters[i]
            target_pixels = target_pixels_dict[filter_str]
            sextractor_object = sextractor_objects[i]
            closest_objects = [sextractor_object.findClosestObjectToPoint(pixel)
                               for pixel in target_pixels]
            closest_objects_dict[filter_str] = closest_objects
            better_fwhms_dict[filter_str]= [sextractor_object.computeRobustFWHMOfObject(sex_object['NUMBER'], 'proc_' + image_root + '_' + filter_str + '.fits', cat_dir) for sex_object in closest_objects]
            
        #print ('closest_objects_dict = ' + str(closest_objects_dict))

        n_objects = len(closest_objects_dict[filters[0]])
        for i in range(n_objects):
            #filter_str = filters[i]
            plot_object_location = plot_object_locations[i]
            #closest_objects = closest_objects_dict[filter_str]
            focus_val_to_use = focus_vals[count]
            scats = [axarr[plot_object_location[0], plot_object_location[1]].scatter([focus_val_to_use], pix_scale * np.array([better_fwhms_dict[filters[j]][i]]),
                                                                             marker = 'o', c = colors[j])
             for j in range(len(filters))]
            axarr[plot_object_location[0], plot_object_location[1]].text(89.0, 2.00, plot_object_text[i])
        #axarr[plot_filter_location[0], plot_filter_location[1]].scatter(focus_vals_to_use, [py_obj[elongation_key_str] for py_obj in closest_objects], marker = '^')
        count = count + 1
        #axarr[plot_object_location[0], plot_object_location[1]].legend(scats, colors)

    f.legend(scats, filter_names, ncol = 4) 
    f.subplots_adjust(hspace = 0.0, wspace = 0.0)
    y_lims = [0.0,0.9]
    axarr[0,0].set_xticks([])
    axarr[0,0].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    axarr[0,0].set_ylabel('FWHM (arcsec)')
    axarr[0,0].set_ylim(y_lims)
    axarr[0,1].set_xticks([])
    axarr[0,1].set_yticks([])
    axarr[0,1].set_ylim(y_lims)
    #axarr[1,0].set_xticks([0,400,800,1200])
    axarr[1,0].set_xlabel('focus value (mm)')
    #axarr[1,0].set_yticks([0,5,10,15,20])
    axarr[1,0].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    axarr[1,0].set_ylabel('FWHM (arcsec)')
    axarr[1,0].set_ylim(y_lims)
    #axarr[1,1].set_xticks([0,400,800,1200])
    axarr[1,1].set_xlabel('focus value (mm)')
    axarr[1,1].set_yticks([])
    axarr[1,1].set_ylim(y_lims)
    plt.show() 
