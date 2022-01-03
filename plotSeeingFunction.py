import matplotlib.pyplot as plt
import csv
import numpy as np
import sys 

if __name__ == "__main__":
    seeing_file_name, seeing_image_name, target_dir = sys.argv[1:] 
    with open(target_dir + seeing_file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        jds = []
        file_numbers = []
        fwhms = []
        stds = []
        
        for row in csv_reader:
            #print ('row = ' + str(row) )
            jds = jds + [float(row[0])]
            file_numbers = file_numbers + [int(float(row[1]))]
            fwhm_stds = row[2:]
            n_filters = int(len(fwhm_stds)/2) 
            fwhms = fwhms + [[float(fwhm_stds[2*i]) for i in range(n_filters)]]
            stds = stds + [[float(fwhm_stds[2*i+1]) for i in range(n_filters)]]

    figsize = (20, 10)
    f, axarr = plt.subplots(2,2, figsize = figsize)
    colors = ['g','r','black','b']
    filter_names = ['g','r','i','z']
    plt.subplots_adjust(wspace = 0.0, hspace = 0.0)
    #print ('fwhms = ' + str(fwhms))
    #print ('stds = ' + str(stds) )
    #print (np.array(stds))
    #print (np.array(fwhms))
    #print (np.shape(np.array(stds)))
    #print (np.shape(np.array(fwhms))) 
    #ax2 = axarr[0,0].twiny()
    #ax2.set_xticks([fwhm[0] for fwhm in fwhms])
    #ax2.set_xticklabels([fwhm[0] for fwhm in fwhms])
    #scats = [plt.scatter((jds), [fwhm[i] for fwhm in fwhms], c = colors[i]) for i in range(n_filters)]
    #[plt.errorbar((jds), [fwhm[i] for fwhm in fwhms], yerr = [std[i] for std in stds], ecolor = colors[i], fmt = 'none') for i in range(n_filters)]
    #plt.legend(scats, filter_names)
    print ('[fwhms, stds] = ' + str([fwhms, stds]))
    y_lims = [np.min(np.array(fwhms) - np.array(stds)) * 0.9, np.max(np.array(fwhms) + np.array(stds)) * 1.1]
    #plt.ylim(y_lims)
    y_extent = y_lims[1] - y_lims[0] 
    axarr[0,0].scatter(jds, [fwhm[0] for fwhm in fwhms], c = colors[0])
    axarr[0,0].errorbar(jds, [fwhm[0] for fwhm in fwhms], yerr = [std[0] for std in stds], fmt = 'none', ecolor = colors[0])
    axarr[0,0].set_xticks([])
    axarr[0,0].set_ylabel('Median seeing FWHM (arcsec)')
    [axarr[0,0].annotate(str(int(file_numbers[i])) + '_' + filter_names[0], xy = (jds[i], fwhms[i][0] - y_extent * 0.1)) for i in range(len(jds))]
    axarr[0,0].set_ylim(y_lims)
    axarr[0,1].scatter(jds, [fwhm[1] for fwhm in fwhms], c = colors[1])
    axarr[0,1].errorbar(jds, [fwhm[1] for fwhm in fwhms], yerr = [std[1] for std in stds], fmt = 'none', ecolor = colors[1])
    axarr[0,1].set_xticks([])
    axarr[0,1].set_yticks([])
    [axarr[0,1].annotate(str(int(file_numbers[i])) + '_' + filter_names[1], xy = (jds[i], fwhms[i][1] - y_extent * 0.1)) for i in range(len(jds))]
    axarr[0,1].set_ylim(y_lims)
    axarr[1,0].scatter(jds, [fwhm[2] for fwhm in fwhms], c = colors[2])
    axarr[1,0].errorbar(jds, [fwhm[2] for fwhm in fwhms], yerr = [std[2] for std in stds], fmt = 'none', ecolor = colors[2])
    axarr[1,0].set_xlabel('Julian time of observation (hrs)')
    axarr[1,0].set_ylabel('Median seeing FWHM (arcsec)')
    [axarr[1,0].annotate(str(int(file_numbers[i])) + '_' + filter_names[2], xy = (jds[i], fwhms[i][2] - y_extent * 0.1)) for i in range(len(jds))]
    axarr[1,0].set_ylim(y_lims)
    axarr[1,1].scatter(jds, [fwhm[3] for fwhm in fwhms], c = colors[3])
    axarr[1,1].errorbar(jds, [fwhm[3] for fwhm in fwhms], yerr = [std[3] for std in stds], fmt = 'none', ecolor = colors[3])
    axarr[1,1].set_xlabel('Julian time of observation (hrs)')
    axarr[1,1].set_yticks([])
    [axarr[1,1].annotate(str(int(file_numbers[i])) + '_' + filter_names[3], xy = (jds[i], fwhms[i][3] - y_extent * 0.1)) for i in range(len(jds))]
    axarr[1,1].set_ylim(y_lims)

    plt.savefig(target_dir + seeing_image_name)
    #plt.show() 
