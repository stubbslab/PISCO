#Defines a light curve by measuring photometry in a series of images_in_master_dict
# Ideally,
import SharedPISCOCommands as spc
import numpy as np
import matplotlib.pyplot as plt
import cantrips as c
import SashasAstronomyTools as atools
import datetime
import photutils
import MasterCatalogueObject as mco


class PISCOLightCurve:

    def measureLightCurve(self):
        return 1

    def getImageSeeing(self, centroiding_fwhm = 4.5, bg_centroiding_fwhm = [5.0, 5.5], #centroid_finding_in_fwhm = 5.0, centroid_bg_in_fwhm  = [6.0, 7.0],
                       guess_seeing_in_arcsec = 1.0, n_stars_for_seeing = 20, min_iters_for_measuring_fwhm = 3,
                       show_seeing_fits = 0):
        img_fwhms = [0.0 for obj_file in self.obj_files]
        for file_num in range(len(self.obj_files)):
            obj_file = self.obj_files[file_num]
            suffix = self.pisco_commands.getSuffix()
            img_band = self.img_bands[file_num]

            star_fwhms = atools.getSeeingInImage(obj_file, self.target_dir, guess_seeing_in_arcsec,
                                                 desired_n_stars =  n_stars_for_seeing, min_n_iters_before_completion = min_iters_for_measuring_fwhm,
                                                 show_fit = show_seeing_fits, fit_radius_scaling = centroiding_fwhm, bg_radii_scalings = bg_centroiding_fwhm,
                                                 good_pix_thresh = self.pixel_thresholds_by_band[img_band], verbose = 0)
            seeing = np.median(star_fwhms) * self.pixel_scale
            print ('seeing = ' + str(seeing))
            img_fwhms[file_num] = seeing
        self.seeings = img_fwhms
        return 1

    def measureADUOfTarget(self, star_aperture_in_fwhm = 1.5, bg_annulus_in_fwhm = [2.5, 3.0], centroiding_fwhm_scaling = 4.5, bg_centroiding_scaling = [5.0, 5.5], show_source_fit = 1):
        target_pixel_coords = [[0.0, 0.0] for obj_file in self.obj_files]
        #calculates reference time with respect to first image in file\

        for i in range(len(self.obj_files)):
            obj_file = self.obj_files[i]
            suffix = self.pisco_commands.getSuffix()
            img_band = self.img_bands[i]
            gain = self.gains_by_band[img_band]
            seeing = self.seeings[i]
            target_coord = self.target_coords[i]
            print ('Seeing (in arcseconds) is ' + str(seeing) )
            data, header = c.readInDataFromFitsFile(obj_file, self.target_dir)
            exp_time = float(header[self.pisco_commands.getExpTimeKeyword()])
            exposure_date_str =  header[self.pisco_commands.getDateObsKeyword()]#header['DATEOBS']
            exposure_start_str =  header[self.pisco_commands.getExpStartKeyword()] #shutter_open_str = float(header['UT-OPEN'])
            exposure_end_str =  header[self.pisco_commands.getExpEndKeyword()] # shutter_close_str = header['UT-CLOSE']
            #start_date_time = datetime.datetime.strptime(exposure_date + ' ' + shutter_open_str, '%Y-%m-%d %H:%M:%S.%f')
            start_time = datetime.datetime.strptime(exposure_date_str + ' ' + exposure_start_str, '%Y-%m-%d %H:%M:%S.%f')
            print ('start_time = ' + str(start_time))
            end_time = datetime.datetime.strptime(exposure_date_str + ' ' + exposure_end_str, '%Y-%m-%d %H:%M:%S.%f')
            print ('end_time = ' + str(end_time))
            if i == 0:
                delta_t = 0.0
                ref_time = start_time
            else:
                delta_t = start_time - ref_time

            print ('delta_t = '  + str(delta_t))
            if not(type(delta_t)) == float:
                delta_t = delta_t.total_seconds()

            print ('self.delta_ts = ' + str(self.delta_ts))


            x, y = atools.getPixelFromWCSHeader(header, target_coord, ra_dec_in_deg = 0)
            centroid_fit = atools.measure1dPSFOfStar(x, y, obj_file, self.target_dir, seeing,
                                                     fit_radius_scaling = centroiding_fwhm_scaling, bg_radii_scalings = bg_centroiding_scaling,
                                                     show_fit = show_source_fit)
            if centroid_fit is None:
                adj_x, adj_y = [x,y]
            else:
                adj_x, adj_y = [centroid_fit[0], centroid_fit[1]]
            source_aperture = photutils.CircularAperture([adj_x, adj_y], seeing / self.pixel_scale * star_aperture_in_fwhm)
            bg_annulus = photutils.CircularAnnulus([adj_x, adj_y], seeing / self.pixel_scale * bg_annulus_in_fwhm[0], seeing * bg_annulus_in_fwhm[1] / self.pixel_scale * star_aperture_in_fwhm)

            phot_table = photutils.aperture_photometry(data, [source_aperture, bg_annulus])

            bg_rate = phot_table['aperture_sum_1'][0] / bg_annulus.area
            source_adu = phot_table['aperture_sum_0'][0] - bg_rate * source_aperture.area

            source_adu_rate = source_adu / exp_time
            source_e_rate = source_adu_rate / gain
            source_and_bg_and_bias_uncertainty = np.sqrt(phot_table['aperture_sum_0'][0] + header['ENOISE'])
            source_e_rate_uncertainty = source_and_bg_and_bias_uncertainty / exp_time

            print ("[adj_x, adj_y, seeing / self.pixel_scale, delta_t, source_e_rate] = " + str([adj_x, adj_y, seeing / self.pixel_scale, delta_t, source_e_rate] ))

            self.source_xs[i] = adj_x
            self.source_ys[i] = adj_y
            self.source_fwhms[i] = seeing / self.pixel_scale
            self.delta_ts[i] = delta_t
            self.source_e_rates[i] = source_e_rate
            self.source_e_errs[i] = np.sqrt(source_e_rate)

            #self.target_xs, self.target_ys, self.target_fwhms, self.target_e_rates, self.target_e_errs = [[-1 for target_file in self.target_files] for i in range(5)
        #light_curves[band] = light_curves[band] + [[end_time, source_e_rate, source_e_rate_uncertainty, exp_time, seeing]]

        return 1

    def determineImageBands(self):
        img_bands = ['' for i in range(len(self.obj_files))]
        for i in range(len(self.obj_files)):
            obj_file = self.obj_files[i]
            suffix = self.pisco_commands.getSuffix()
            img_band = obj_file[-(len(suffix) + 1)]
            img_bands[i] = img_band
        return img_bands

    def readInMasterCatalogue(self, new_mcat_file = None, specific_imgs_dir = None):
        if not(new_mcat_file == None):
            self.mcat_file = new_mcat_file
        if not(specific_imgs_dir == None):
            imgs_dir = specific_imgs_dir
        else:
            imgs_dir = self.target_dir
        self.target_mcat = mco.MasterCatalogue(self.target_dir, self.mcat_file, imgs_dir = imgs_dir)
        return 1

    #We want to use the master catalogue object to identify the best background objects for reference.
    # "Good" means the most star like objects that are also: (a) measured in every image, (b) no warning flags.
    # The function will take either the self.n_stars_for_relative_photometry or the number of good objects, whichever is lower.

    def identifyReferenceObjects(self, missing_prefix = 'wcs_', missing_suffix = '.fits', filt_suffix_funct = lambda filt_str: '_' + filt_str):
        good_object_keyword_prefix = self.pisco_commands.getGoodObjectKeywordPrefix()
        star_gal_keyword_prefix = self.pisco_commands.getStarVsGalKeywordPrefix()
        object_peak_flux_keyword = self.pisco_commands.getObjectMaxFluxKeywordPrefix()
        mcat_keyword_band_suffix = self.pisco_commands.getKeywordBandSuffix()

        mcat_img_files = self.target_mcat.image_file_names[1:]
        img_files_to_use = c.intersection(self.obj_files, [missing_prefix + mcat_img_file + missing_suffix for mcat_img_file in mcat_img_files])
        print ('img_files_to_use are ' + str(img_files_to_use))
        for obj_file in self.obj_files:
            if not(obj_file in img_files_to_use):
                print ('Warning: file ' + obj_file + ' not found in master catalogue object. ')
        img_file_keys = list(set([file[len(missing_prefix):-(len(filt_suffix_funct(self.filt_strs[0])) + len(missing_suffix)) ] for file in img_files_to_use ]))
        mcat_obj_ids = list(self.target_mcat.master_val_dict.keys())
        well_detected_objects = mcat_obj_ids[:]
        for filt_str in self.filt_strs:
            well_detected_objects_in_filt = [obj_id for obj_id in mcat_obj_ids if np.all([self.target_mcat.master_val_dict[obj_id][file_key][good_object_keyword_prefix + mcat_keyword_band_suffix + filt_str] for file_key in img_file_keys ]) ]
            peak_threshold = self.pixel_thresholds_by_band[filt_str]
            print ('peak_threshold = ' + str(peak_threshold))
            objs_below_threshold_in_filt = [obj_id for obj_id in mcat_obj_ids if np.all([self.target_mcat.master_val_dict[obj_id][file_key][object_peak_flux_keyword + mcat_keyword_band_suffix + filt_str] < peak_threshold for file_key in img_file_keys ]) ]
            well_detected_objects = c.intersection(c.intersection(well_detected_objects_in_filt, objs_below_threshold_in_filt), well_detected_objects)
        average_star_classes = [np.median(c.flattenListOfLists([[self.target_mcat.master_val_dict[obj_id][file_key][star_gal_keyword_prefix + mcat_keyword_band_suffix + filt_str] for file_key in img_file_keys] for filt_str in self.filt_strs])) for obj_id in well_detected_objects]
        sorted_star_classes, sorted_well_detected_objects = c.safeSortOneListByAnother(average_star_classes, [average_star_classes, well_detected_objects])
        self.sorted_star_classes = c.niceReverse(sorted_star_classes)
        self.sorted_well_detected_objects = c.niceReverse(sorted_well_detected_objects)
        print ('self.sorted_star_classes = ' + str(self.sorted_star_classes))
        print ('self.sorted_well_detected_objects = ' + str(self.sorted_well_detected_objects))

        return 1


    def measureReferenceADUs(self, n_target_ref_objs = 10, star_aperture_in_fwhm = 1.5, bg_annulus_in_fwhm = [2.5, 3.0], centroiding_fwhm_scaling = 4.5, bg_centroiding_scaling = [5.0, 5.5], show_source_fit = 0,
                             missing_prefix = 'wcs_', missing_suffix = '.fits', filt_suffix_funct = lambda filt_str: '_' + filt_str, verbose = 0):
            obj_ADUs = [[] for obj_file in self.obj_files]
            star_position_keyword_prefix = self.pisco_commands.getStarPositionKeywordPrefix()
            mcat_keyword_band_suffix = self.pisco_commands.getKeywordBandSuffix()
            n_ref_objs = min(n_target_ref_objs, len(self.sorted_well_detected_objects))
            ref_objs = self.sorted_well_detected_objects[0:n_ref_objs]
            print ('ref_objs are ' + str(ref_objs))
            weighted_mean_e_rates = [0.0 for obj_file in self.obj_files]
            #calculates reference time with respect to first image in file\

            for i in range(len(self.obj_files)):
                obj_file = self.obj_files[i]
                suffix = self.pisco_commands.getSuffix()
                img_band = self.img_bands[i]
                gain = self.gains_by_band[img_band]
                seeing = self.seeings[i]
                if verbose: print ('In measureReferenceADUs, [obj_file, suffix, img_band, gain, seeing] = ' + str([obj_file, suffix, img_band, gain, seeing]) )
                data, header = c.readInDataFromFitsFile(obj_file, self.target_dir)
                exp_time = float(header[self.pisco_commands.getExpTimeKeyword()])
                file_root = obj_file[len(missing_prefix):-(len(filt_suffix_funct(img_band)) + len(missing_suffix)) ]
                ref_obj_e_rates = [-1.0 for ref_obj in ref_objs]
                ref_obj_e_rate_errs = [-1.0 for ref_obj in ref_objs]
                for ref_obj_num in range(n_ref_objs):
                    obj_id = ref_objs[ref_obj_num]

                    x, y = self.target_mcat.master_val_dict[obj_id][file_root][star_position_keyword_prefix  + mcat_keyword_band_suffix + img_band]
                    centroid_fit = atools.measure1dPSFOfStar(x, y, obj_file, self.target_dir, seeing,
                                                             fit_radius_scaling = centroiding_fwhm_scaling, bg_radii_scalings = bg_centroiding_scaling,
                                                             show_fit = show_source_fit, verbose = verbose)
                    if centroid_fit is None:
                        adj_x, adj_y = [x,y]
                    else:
                        adj_x, adj_y = [centroid_fit[0], centroid_fit[1]]

                    print ('I adjusted centroid of master object ' + str(obj_id) + ' in file ' + str(obj_file) + ' from ' + str([x,y]) + ' to ' + str([adj_x, adj_y]))
                    source_aperture = photutils.CircularAperture([adj_x, adj_y], seeing / self.pixel_scale * star_aperture_in_fwhm)
                    bg_annulus = photutils.CircularAnnulus([adj_x, adj_y], seeing / self.pixel_scale * bg_annulus_in_fwhm[0], seeing * bg_annulus_in_fwhm[1] / self.pixel_scale * star_aperture_in_fwhm)
                    phot_table = photutils.aperture_photometry(data, [source_aperture, bg_annulus])
                    bg_rate = phot_table['aperture_sum_1'][0] / bg_annulus.area
                    source_adu = phot_table['aperture_sum_0'][0] - bg_rate * source_aperture.area

                    source_adu_rate = source_adu / exp_time
                    source_e_rate = source_adu_rate / gain
                    source_and_bg_and_bias_uncertainty = np.sqrt(phot_table['aperture_sum_0'][0] + header['ENOISE'])
                    source_e_rate_uncertainty = source_and_bg_and_bias_uncertainty / exp_time
                    ref_obj_e_rates[ref_obj_num] = source_e_rate
                    ref_obj_e_rate_errs[ref_obj_num] = source_e_rate_uncertainty

                weighted_mean_e_rates[i] = c.weighted_mean(ref_obj_e_rates, ref_obj_e_rate_errs)
            self.ref_e_rates = weighted_mean_e_rates

            return 1


    def __init__(self, obj_files, target_coords, mcat_file,
                 target_dir = '', coords_type = 'wcs',
                 star_aperture_in_fwhm = 1.5, bg_annulus_in_fwhm = [2.5, 3.0], #centroid_finding_in_fwhm = 5.0, centroid_bg_in_fwhm  = [6.0, 7.0],
                 guess_seeing_in_arcsec = 1.0, n_stars_for_seeing = 20, min_iters_for_measuring_fwhm = 3,
                 centroiding_fwhm = 4.5, bg_centroiding_fwhm = [5.0,5.5],
                 n_stars_for_relative_photometry = 5,
                 show_seeing_fits = 0, show_source_fit = 1, show_fig = 0, save_fig = 0,
                 save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/PISCO/plots/', fig_name = 'J1311-3430_uncal_LC_ut200308_PISCO.pdf'):

        self.target_dir = target_dir
        self.mcat_file = mcat_file
        self.obj_files = obj_files
        self.pisco_commands = spc.CommandHolder()
        pixel_thresholds_by_band = self.pisco_commands.getPixelThresholds()
        self.filt_strs = self.pisco_commands.getFilters()
        n_filts = len(self.filt_strs)
        self.pixel_thresholds_by_band = {self.filt_strs[i]:pixel_thresholds_by_band[i] for i in range(n_filts) }
        gains_by_amp = self.pisco_commands.getAmplifierGains()
        self.gains_by_band = {self.filt_strs[i]:gains_by_amp[i*2] for i in range(n_filts) }
        self.pixel_scale = self.pisco_commands.getPixelScale()
        self.n_stars_for_relative_photometry = n_stars_for_relative_photometry

        self.img_bands = self.determineImageBands()

        #If less coordinates are given than the number of images, assume that
        # we are looking at an object at the same location.
        if type(target_coords[0]) in [str, float, int]:
            self.target_coords = [target_coords for obj_file in obj_files]
        else:
            self.target_coords= target_pixel_coords

        self.seeings = [0.0 for obj_file in self.obj_files]
        self.getImageSeeing(guess_seeing_in_arcsec = guess_seeing_in_arcsec, n_stars_for_seeing = n_stars_for_seeing, min_iters_for_measuring_fwhm = min_iters_for_measuring_fwhm,
                                           centroiding_fwhm = centroiding_fwhm, bg_centroiding_fwhm = bg_centroiding_fwhm, show_seeing_fits = show_seeing_fits)

        #print('self.seeings = ' + str(self.seeings))

        self.source_xs, self.source_ys, self.source_fwhms, self.delta_ts, self.source_e_rates, self.source_e_errs = [[-1 for target_file in self.obj_files] for i in range(6)]
        #measured_targets = self.measureADUOfTarget(star_aperture_in_fwhm = star_aperture_in_fwhm, bg_annulus_in_fwhm = bg_annulus_in_fwhm, show_source_fit = 1)#centroid_finding_in_fwhm = 5.0, centroid_bg_in_fwhm  = [6.0, 7.0],):

        #print ('[self.img_bands, self.source_xs, self.source_ys, self.source_fwhms, self.source_e_rates, self.source_e_errs ] = ' + str([self.img_bands, self.source_xs, self.source_ys, self.source_fwhms, self.source_e_rates, self.source_e_errs ]))
