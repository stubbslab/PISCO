import fitsAnalysisScripts as fas
import cantrips as c
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def measureFWHMFromCatFile(file_name, n_cat_keys = 19, fwhm_index = 16):
    columns = c.readInColumnsToList(file_name, n_ignore = n_cat_keys, delimiter = ' ', verbose = 0)
    fwhms = [float(elem) for elem in columns[fwhm_index-1] ]
    return fwhms

if __name__=="__main__":
    obj_center = SkyCoord("12:06:50.41 +02:01:42.46", unit=(u.hourangle, u.deg))
    data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/PISCO/Apr_2019/SDSS-J120650/'
    shared_prefix = 'proc_SDSS-J12065_ut1904130_'
    file_prefixes = [shared_prefix + str(img_num) + '_' for img_num in range(72, 77)]
    g_files = ['wcs_' + prefix + 'g' + '.fits' for prefix in file_prefixes]
    r_files = ['wcs_' + prefix + 'r' + '.fits' for prefix in file_prefixes]
    i_files = ['wcs_' + prefix + 'i' + '.fits' for prefix in file_prefixes]
    z_files = ['wcs_' + prefix + 'z' + '.fits' for prefix in file_prefixes]
    g_cat_files = [prefix + 'g' + '.cat' for prefix in file_prefixes]
    r_cat_files = [prefix + 'r' + '.cat' for prefix in file_prefixes]
    i_cat_files = [prefix + 'i' + '.cat' for prefix in file_prefixes]
    z_cat_files = [prefix + 'z' + '.cat' for prefix in file_prefixes]

    g_fwhms = [np.median(measureFWHMFromCatFile(data_dir + g_cat_file, n_cat_keys =  19, fwhm_index = 16)) for g_cat_file in g_cat_files]
    r_fwhms = [np.median(measureFWHMFromCatFile(data_dir + r_cat_file, n_cat_keys =  19, fwhm_index = 16)) for r_cat_file in r_cat_files]
    i_fwhms = [np.median(measureFWHMFromCatFile(data_dir + i_cat_file, n_cat_keys =  19, fwhm_index = 16)) for i_cat_file in i_cat_files]
    z_fwhms = [np.median(measureFWHMFromCatFile(data_dir + z_cat_file, n_cat_keys =  19, fwhm_index = 16)) for z_cat_file in z_cat_files]

    g_center = [1663, 1540]
    r_center = [1664, 1568]
    i_center = [1607, 1529]
    z_center = [1651, 1565]
    g_rstar_nfwhm = 3.0
    g_rinner_nfwhm = g_rstar_nfwhm + 2.0
    g_router_nfwhm = g_rinner_nfwhm + 1.0
    r_rstar_nfwhm = g_rstar_nfwhm
    r_rinner_nfwhm = g_rinner_nfwhm
    r_router_nfwhm = g_router_nfwhm
    i_rstar_nfwhm = g_rstar_nfwhm
    i_rinner_nfwhm = g_rinner_nfwhm
    i_router_nfwhm = g_router_nfwhm
    z_rstar_nfwhm = g_rstar_nfwhm
    z_rinner_nfwhm = g_rinner_nfwhm
    z_router_nfwhm = g_router_nfwhm

    g_data = [c.readInDataFromFitsFile(g_file, data_dir) for g_file in g_files]
    g_data_arrays = [g_datum[0] for g_datum in g_data]
    g_headers = [g_datum[1] for g_datum in g_data]
    r_data = [c.readInDataFromFitsFile(r_file, data_dir) for r_file in r_files]
    r_data_arrays = [r_datum[0] for r_datum in r_data]
    r_headers = [r_datum[1] for r_datum in r_data]
    i_data = [c.readInDataFromFitsFile(i_file, data_dir) for i_file in i_files]
    i_data_arrays = [i_datum[0] for i_datum in i_data]
    i_headers = [i_datum[1] for i_datum in i_data]
    z_data = [c.readInDataFromFitsFile(z_file, data_dir) for z_file in z_files]
    z_data_arrays = [z_datum[0] for z_datum in z_data]
    z_headers = [z_datum[1] for z_datum in z_data]

    g_wcs = [fas.loadWCSFromHeader(header) for header in g_headers]
    r_wcs = [fas.loadWCSFromHeader(header) for header in r_headers]
    i_wcs = [fas.loadWCSFromHeader(header) for header in i_headers]
    z_wcs = [fas.loadWCSFromHeader(header) for header in z_headers]

    g_stars = [fas.measureAnnulusAperture(g_data_arrays[i], g_wcs[i].all_world2pix([[obj_center.ra.deg, obj_center.dec.deg]], 0)[0], g_rstar_nfwhm * g_fwhms[i], g_rinner_nfwhm * g_fwhms[i], g_router_nfwhm * g_fwhms[i],
                 method='subpixel', subpixels = 5) for i in range(len(g_data_arrays))]
    r_stars = [fas.measureAnnulusAperture(r_data_arrays[i], r_wcs[i].all_world2pix([[obj_center.ra.deg, obj_center.dec.deg]], 0)[0], r_rstar_nfwhm * r_fwhms[i], r_rinner_nfwhm * r_fwhms[i], r_router_nfwhm * r_fwhms[i],
                 method='subpixel', subpixels = 5) for i in range(len(r_data_arrays))]
    i_stars = [fas.measureAnnulusAperture(i_data_arrays[i], i_wcs[i].all_world2pix([[obj_center.ra.deg, obj_center.dec.deg]], 0)[0], i_rstar_nfwhm * i_fwhms[i], i_rinner_nfwhm * i_fwhms[i], i_router_nfwhm * i_fwhms[i],
                 method='subpixel', subpixels = 5) for i in range(len(i_data_arrays))]
    z_stars = [fas.measureAnnulusAperture(z_data_arrays[i], z_wcs[i].all_world2pix([[obj_center.ra.deg, obj_center.dec.deg]], 0)[0], z_rstar_nfwhm * z_fwhms[i], z_rinner_nfwhm * z_fwhms[i], z_router_nfwhm * z_fwhms[i],
                 method='subpixel', subpixels = 5) for i in range(len(z_data_arrays))]

    print ('[g_star[0] - g_star[1] / g_star[3] * g_star[2] for g_star in g_stars] = ' + str([g_star[0] - g_star[1] / g_star[3] * g_star[2] for g_star in g_stars]))
    print ('[r_star[0] - r_star[1] / r_star[3] * r_star[2] for r_star in r_stars] = ' + str([r_star[0] - r_star[1] / r_star[3] * r_star[2] for r_star in r_stars]))
    print ('[i_star[0] - i_star[1] / i_star[3] * i_star[2] for i_star in i_stars] = ' + str([i_star[0] - i_star[1] / i_star[3] * i_star[2] for i_star in i_stars]))
    print ('[z_star[0] - z_star[1] / z_star[3] * z_star[2] for z_star in z_stars] = ' + str([z_star[0] - z_star[1] / z_star[3] * z_star[2] for z_star in z_stars]))
