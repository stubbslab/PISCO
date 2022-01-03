import cantrips as c
import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    filts = [r'$g$',r'$r$',r'$i$',r'$z$']
    n_filts = len(filts)
    exposure_keyword = 'EXPTIME'
    volt_keyword = 'COMMENT'
    adu_thresholds = [[1000, 50000], [1000, 50000], [1000, 20000], [1000, 12000]] #Should cut off either at saturation or smearing threshold
    adu_thresholds = [[-np.inf, np.inf], [-np.inf, np.inf], [-np.inf, np.inf], [-np.inf, np.inf]]
    measurement_sections = [[[[600, 1200], [600, 3800]], [[1600, 2400], [600, 3800]]],
                            [[[600, 1200], [600, 3800]], [[1600, 2400], [600, 3800]]],
                            [[[600, 1200], [600, 3800]], [[1600, 2400], [600, 3800]]],
                            [[[600, 1200], [600, 3800]], [[1600, 2400], [600, 3800]]] ]
    n_amps = len(measurement_sections[0][0])
    shared_prefix = 'proc_domeflat_'
    data_dir = '/Volumes/PISCO Passport/PISCOData_Raw/ut200301/'
    files = [[shared_prefix + str(i) + '_g.fits' for i in range(2, 64)], #+ [shared_prefix + '_iz_' + str(i) + '_g.fits' for i in range(537, 594)],
             [shared_prefix + str(i) + '_r.fits' for i in range(2, 64)], # range(537, 594) [shared_prefix + '_g_' + str(i) + '_r.fits' for i in range(188, 270)], #+ [shared_prefix + '_iz_' + str(i) + '_r.fits' for i in range(537, 594)],
             [shared_prefix + str(i) + '_i.fits' for i in range(2, 64)],# range(537, 594) [shared_prefix + '_g_' + str(i) + '_i.fits' for i in range(188, 270)],
             #[shared_prefix + '_iz_' + str(i) + '_z.fits' for i in range(537, 594)] # + [shared_prefix + '_g_' + str(i) + '_z.fits' for i in range(188, 270)]
             [shared_prefix + str(i) + '_z.fits' for i in list(range(2, 64)) + list(range(2, 64))] #The 3s images are weird for the z band, so we cut them out
             ]
    #print ([ [ int(np.transpose(c.readInDataFromFitsFile(file, data_dir)[0])[100:101,200:201])
    #             for file in files[i] ]  for i in range(len(filts)) ] )
    means = [{}  for i in range(n_filts) ]
    stds = [{}  for i in range(n_filts) ]

    filt_strs = ['g','r','i','z']
    img_pairs_by_exptime_and_voltage = [[[shared_prefix + str(i*2) + '_' + filt_str + '.fits', shared_prefix + str(i*2+1) + '_' + filt_str + '.fits'] for i in range(1, 32)] for filt_str in filt_strs]
    all_filt_stds = [ [[] for img_pair_num in range(len(img_pairs))] for img_pairs in img_pairs_by_exptime_and_voltage ]
    all_filt_means = [ [[] for img_pair_num in  range(len(img_pairs))] for img_pairs in img_pairs_by_exptime_and_voltage ]
    print ('[np.shape(all_filt_means), np.shape(all_filt_stds)] = ' + str([np.shape(all_filt_means), np.shape(all_filt_stds)] ) )
    print ('[len(all_filt_means), len(all_filt_means[0])] = ' + str([len(all_filt_means), len(all_filt_means[0])] ))

    for filt_index in range(len(filt_strs)):
        single_band_img_pairs = img_pairs_by_exptime_and_voltage[filt_index]
        print ('single_band_img_pairs = ' + str(single_band_img_pairs))
        for img_pair_num in range(len(single_band_img_pairs)):
            img_pair = single_band_img_pairs[img_pair_num]
            img_data_pairs = [c.readInDataFromFitsFile(img_file, data_dir) for img_file in img_pair]
            data_arrays = [np.transpose(img_data[0]) for img_data in img_data_pairs]
            headers = [img_data[1] for img_data in img_data_pairs]
            diff_array = data_arrays[1] - data_arrays[0]
            all_filt_stds[filt_index][img_pair_num] = [0.0 for amp_num in range(n_amps)]
            all_filt_means[filt_index][img_pair_num] = [0.0 for amp_num in range(n_amps)]
            for amp_num in range(n_amps):
                meas_section = measurement_sections[filt_index][amp_num]
                mean_left = np.mean(data_arrays[0][meas_section[0][0]:meas_section[0][1],meas_section[1][0]:meas_section[1][1]])
                mean_right = np.mean(data_arrays[1][meas_section[0][0]:meas_section[0][1],meas_section[1][0]:meas_section[1][1]])
                mean = np.mean([mean_left, mean_right])
                std = np.std(diff_array[meas_section[0][0]:meas_section[0][1],meas_section[1][0]:meas_section[1][1]])
                print ('[img_pair, filt_index, amp_num, meas_section, mean_left, mean_right, mean, std] = ' + str([img_pair, filt_index, amp_num, meas_section, mean_left, mean_right, mean, std]))
                all_filt_means[filt_index][img_pair_num][amp_num] = mean
                all_filt_stds[filt_index][img_pair_num][amp_num]  = std
    print ('all_filt_means = ' + str(all_filt_means))
    print ('all_filt_stds = ' + str(all_filt_stds))

    #means = [ [ np.mean(np.transpose(c.readInDataFromFitsFile(file, data_dir)[0])[measurement_sections[i][0][0]:measurement_sections[i][0][1],measurement_sections[i][1][0]:measurement_sections[i][1][1]])
    #             for file in files[i] ]  for i in range(len(filts)) ]
    #stds = [ [ np.std(np.transpose(c.readInDataFromFitsFile(file, data_dir)[0])[measurement_sections[i][0][0]:measurement_sections[i][0][1],measurement_sections[i][1][0]:measurement_sections[i][1][1]])
    #         for file in files[i] ]  for i in range(len(filts)) ]

    figsize = [8,8]
    f, axarr = plt.subplots(n_filts, n_amps, figsize = figsize)
    for filt_num in range(n_filts):
        for amp_num in range(n_amps):
            single_filt_means = all_filt_means[filt_num]
            single_filt_stds = all_filt_stds[filt_num]
            xs = [single_pair_means[amp_num] for single_pair_means in single_filt_means if (single_pair_means[amp_num] > adu_thresholds[filt_num][0] and single_pair_means[amp_num] < adu_thresholds[filt_num][1])]
            ys = [single_filt_stds[img_pair_num][amp_num] ** 2.0 / 2.0 for img_pair_num in range(len(single_filt_stds)) if (single_filt_means[img_pair_num][amp_num] > adu_thresholds[filt_num][0]  and single_filt_means[img_pair_num][amp_num] < adu_thresholds[filt_num][1])]
            poly_fit = np.polyfit(xs, ys, 1)
            funct = np.poly1d(poly_fit)
            axarr[filt_num, amp_num].scatter(xs, ys )
            x_lims = [np.min(xs) * 0.9, np.max(xs) * 1.1]
            axarr[filt_num, amp_num].plot([x_lims[0]] + sorted(xs) + [x_lims[1]] , funct([x_lims[0]] + sorted(xs) + [x_lims[1]]), c = 'r')
            axarr[filt_num, amp_num].set_xlim(x_lims)
            axarr[filt_num, amp_num].set_ylabel(r'$\sigma_{ADU}^2$')
            axarr[filt_num, amp_num].set_xlabel('Mean ADU')
            axarr[filt_num, amp_num].set_title(r'$\sigma^2$ vs mean for filt ' + filts[filt_num] + (' right' if amp_num else ' left') + ' amp')
            axarr[filt_num, amp_num].text((x_lims[1] - x_lims[0])/2, 8.0 * (np.max(ys) - np.min(ys)) / 9.0, r'$G=$slope$=$' + str(c.round_to_n(poly_fit[0], 3)), horizontalalignment = 'center')
        #axarr[i%2, i//2].scatter([means[i][j] for j in range(len(means[i])) if means[i][j] < adu_threshold] , [means[i][j] / (stds[i][j] ** 1.0 ) for j in range(len(means[i])) if means[i][j] < adu_threshold] )
    plt.tight_layout()
    plt.show()
