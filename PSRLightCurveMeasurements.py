import cantrips as c
import SashasAstronomyTools as atools
import photutils
import numpy as np
import datetime
import matplotlib.pyplot as plt
import MasterCatalogueObject as mco

#Images of PSR J1311-3430 taken on 2020/03/02

def getImageFiles(wcs = 0, proc = 0):
    img_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/PISCO/Feb_2020/ut200301/'
    image_roots = ['PSRJ1311-3430_' + suffix for suffix in ['A_169', 'A_178', 'A_185', 'A_191', 'A_200', 'A_203', 'A_212', 'A_216',
                                                            'B_170', 'B_179', 'B_186', 'B_192', 'B_201', 'B_204', 'B_217',
                                                            'C_171', 'C_180', 'C_187', 'C_193', 'C_205', 'C_218', #C_202 is bad; rotated
                                                            'D_172', 'D_181', 'D_188', 'D_194', 'D_206', 'D_219',
                                                            'E_173', 'E_182', 'E_189', 'E_195', 'E_207', 'E_220',
                                                            'F_174', 'F_183', 'F_190', 'F_208', 'F_221',
                                                            'G_175', 'G_184', 'G_209', 'G_222',
                                                            'H_176', 'H_210', 'H_223', 'I_177', 'I_211'] ]
    filt_strs = ['g','r','i','z']
    if wcs:
        prefix = 'wcs_proc_'
        suffixes = ['_' + filt_str + '.fits' for filt_str in filt_strs]
    elif proc:
        prefix = 'proc_'
        suffixes = ['_' + filt_str + '.fits' for filt_str in filt_strs]
    else:
        prefix = ''
        suffixes = ['.fits']
    processed_files = c.flattenListOfLists([[prefix + img + suffix for suffix in suffixes] for img in image_roots])

    return [img_dir, [proc_file for proc_file in processed_files ]]

def measureSeeingInPSRImages(pixel_thresholds_by_band = [35000, 35000, 15000, 12000]):
    guess_seeing_in_arcsec = 1.0
    psr_dir, psr_files = getImageFiles(wcs = 0, proc = 1)

    psr_fwhms = [-1.0 for psr_file in psr_files]


    #for i in range(len(PSR_files)):
    for i in range(len(psr_files)):
        psr_file = psr_files[i]
        good_pix_thresh = pixel_thresholds_by_band[i % len(pixel_thresholds_by_band)]
        print ('good_pix_thresh = ' + str(good_pix_thresh))
        star_fwhms = atools.getSeeingInImage(psr_file, '', guess_seeing_in_arcsec,
                                             desired_n_stars = 20, min_n_iters_before_completion = 3, show_fit = 0, good_pix_thresh = pixel_thresholds_by_band[i % len(pixel_thresholds_by_band)])
        #star_fwhms = [star_fwhm for star_fwhm in star_fwhms if not(star_fwhm is None)]
        psr_fwhms[i] = star_fwhms
        print ('star_fwhms = ' + str(star_fwhms))

    return psr_fwhms

#/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/PISCO/Feb_2020/ut200301/PSRJ1311-3430.mcat
def measureMedianADUOfBackgroundStars(mcat_file, mcat_dir = '', meas_stars = []):
    psr_mcat = mco.MasterCatalogue(mcat_dir, mcat_file)


    return  []

def measureADUOfPSR(star_aperture_in_fwhm = 1.5, bg_annulus_in_fwhm = [2.5, 3.0], #centroid_finding_in_fwhm = 5.0, centroid_bg_in_fwhm  = [6.0, 7.0],
                    guess_seeing_in_arcsec = 1.0, n_stars_for_seeing = 20,
                    centroiding_star_fwhm = 4.5, bg_centroiding_fwhm = [5.0,5.5],
                    target_coord = ["13:11:45.724 hours", "-34:30:30.34 degrees"], pixel_thresholds_by_band = [35000, 35000, 15000, 12000],
                    pixel_scale = 0.1, gains = {'g':1.88, 'r':1.88, 'i':1.84, 'z':1.82}, plot_light_curves = 0,
                    show_fig = 0, save_fig = 0,
                    save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/PISCO/plots/', fig_name = 'J1311-3430_uncal_LC_ut200308_PISCO.pdf'):
    psr_dir, psr_files = getImageFiles(wcs = 1)

    target_pixel_coords = [[0.0, 0.0] for psr_file in psr_files]
    filt_strs = ['g','r','i','z']
    light_curves = {filt_str:[] for filt_str in filt_strs}
    for i in range(len(psr_files)):
        psr_file = psr_files[i]
        print ('Measuring seeing of psr_file ' + psr_dir + psr_file)
        band = psr_file[-6]
        gain = gains[band]

        star_fwhms = atools.getSeeingInImage(psr_file, psr_dir, guess_seeing_in_arcsec,
                                             desired_n_stars =  n_stars_for_seeing, min_n_iters_before_completion = 3, show_fit = 0,
                                             good_pix_thresh = pixel_thresholds_by_band[i % len(pixel_thresholds_by_band)], verbose = 0)
        seeing = np.median(star_fwhms) * pixel_scale
        print ('Seeing (in arcseconds) is ' + str(seeing) )
        data, header = c.readInDataFromFitsFile(psr_file, '')
        exp_time = float(header['EXPTIME'])
        exposure_date = header['DATEOBS']
        #shutter_open_str = float(header['UT-OPEN'])
        shutter_close_str = header['UT-CLOSE']
        #start_date_time = datetime.datetime.strptime(exposure_date + ' ' + shutter_open_str, '%Y-%m-%d %H:%M:%S.%f')
        end_time = datetime.datetime.strptime(exposure_date + ' ' + shutter_close_str, '%Y-%m-%d %H:%M:%S.%f')


        x, y = atools.getPixelFromWCSHeader(header, target_coord, ra_dec_in_deg = 0)
        centroid_fit = atools.measure1dPSFOfStar(x, y, psr_file, '', seeing,
                                                 fit_radius_scaling = centroiding_star_fwhm, bg_radii_scalings = bg_centroiding_fwhm ,
                                                 show_fit = 1)
        if centroid_fit is None:
            adj_x, adj_y = [x,y]
        else:
            adj_x, adj_y = [centroid_fit[0], centroid_fit[1]]
        source_aperture = photutils.CircularAperture([adj_x, adj_y], seeing / pixel_scale * star_aperture_in_fwhm)
        bg_annulus = photutils.CircularAnnulus([adj_x, adj_y], seeing / pixel_scale * bg_annulus_in_fwhm[0], seeing * bg_annulus_in_fwhm[1] / pixel_scale * star_aperture_in_fwhm)

        phot_table = photutils.aperture_photometry(data, [source_aperture, bg_annulus])

        bg_rate = phot_table['aperture_sum_1'][0] / bg_annulus.area
        source_adu = phot_table['aperture_sum_0'][0] - bg_rate * source_aperture.area

        source_adu_rate = source_adu / exp_time
        source_e_rate = source_adu_rate / gain
        source_and_bg_and_bias_uncertainty = np.sqrt(phot_table['aperture_sum_0'][0] + header['ENOISE'])
        source_e_rate_uncertainty = source_and_bg_and_bias_uncertainty / exp_time

        print ("[adj_x, adj_y, seeing / pixel_scale, source_e_rate] = "
                + str([adj_x, adj_y, seeing / pixel_scale, source_e_rate] ))

        light_curves[band] = light_curves[band] + [[end_time, source_e_rate, source_e_rate_uncertainty, exp_time, seeing]]

    if show_fig or save_fig:
        colors = ['g','r','magenta','cyan']
        plt.rc('font', family='serif')
        plt.rc('text', usetex=True)
        f, axarr = plt.subplots(2,4, sharex = True, figsize = (20, 6))
        plt.subplots_adjust(hspace = 0.0)
        for ind in range(len(filt_strs)):
            filt = filt_strs[ind]
            light_curve = light_curves[filt]
            time_floats = [elem[0].hour * 3600 + elem[0].minute * 60 + elem[0].second for elem in light_curve]
            flux = [elem[1] for elem in light_curve]
            errs = [elem[2] for elem in light_curve]
            exp_times = [int(elem[3]) for elem in light_curve]
            seeing  = [float(elem[4]) for elem in light_curve]
            scat_300 = axarr[1][ind].scatter([time_floats[obs_num] for obs_num in range(len(time_floats)) if exp_times[obs_num] == 300], [flux[obs_num] for obs_num in range(len(flux)) if exp_times[obs_num] == 300], c = colors[ind], marker = 'o')
            scat_180 = axarr[1][ind].scatter([time_floats[obs_num] for obs_num in range(len(time_floats)) if exp_times[obs_num] == 180], [flux[obs_num] for obs_num in range(len(flux)) if exp_times[obs_num] == 180], c = colors[ind], marker = 'x')
            axarr[1][ind].errorbar(time_floats, flux, yerr = errs, c = colors[ind], fmt = 'none')
            axarr[1][ind].legend([scat_300, scat_180], ['300s exp', '180s exp'])
            axarr[0][ind].scatter([time_floats[obs_num] for obs_num in range(len(time_floats)) if exp_times[obs_num] == 300], [seeing[obs_num] for obs_num in range(len(seeing)) if exp_times[obs_num] == 300], c = colors[ind], marker = 'o')
            axarr[0][ind].scatter([time_floats[obs_num] for obs_num in range(len(time_floats)) if exp_times[obs_num] == 180], [seeing[obs_num] for obs_num in range(len(seeing)) if exp_times[obs_num] == 180], c = colors[ind], marker = 'x')
            axarr[0][ind].set_ylabel('Seeing FWHM (arcsec)')
            axarr[1][ind].set_ylabel('J1311-3430 e-/s')
            axarr[1][ind].set_xlabel('UT Time (s) on date 2020/03/02')
            axarr[0][ind].set_title('PISCO ' + filt_strs[ind] )
            axarr[1][ind].set_yscale('log')
        plt.tight_layout()
        if save_fig:
            plt.savefig(save_dir + fig_name)
        if show_fig:
            plt.show()

    return light_curves

#I want to use my developed infrastructure that cross indexes objects using SExtractror.
# I measure the relative number of ADU of the PSR to those objects.
# I should be able to select only objects shared between all targets
def measureRelativeMagOfPSR():
    return 1


class PSRLightCurve
