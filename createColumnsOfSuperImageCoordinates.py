import numpy as np
import sys
from cantrips import readInDataFromFitsFile
from cantrips import saveDataToFitsFile
from astropy.io import fits
from astropy import wcs
from AstronomicalParameterArchive import AstronomicalParameterArchive 

if __name__ == "__main__":
    args = sys.argv[1:]
    print ('args = ' + str(args))
    target_dir, ref_dir, target_positions_file, ref_wcs_fits_file, target_fake_wcs_txt_file = sys.argv[1:]
    #print ('target_positions_file = ' + str(target_positions_file)) 

    #astro_arch = AstronomicalParameterArchive ()
    #deg_to_rad = astro_arch.getDegToRad()
    #super_image_dim = [1590.0003, 2185.498]
    #super_image_center = [(dim + 2.0) / 2.0 - 1.0 for dim in super_image_dim]
    #single_filter_dim = [1499, 2139]
    #super_image_pix_coords = {'g':[776.67 - 1.0, 1062.04 - 1.0], 'r':[729.75 - 1.0, 1084.75 - 1.0], 'i':[719.52 - 1.0, 1057.657 - 1.0], 'z':[750.46 - 1.0, 1086.93 - 1.0]}
    #super_to_single_angle = {'g':359.132 * deg_to_rad, 'r':0.0 * deg_to_rad, 'i':359.393 * deg_to_rad, 'z':359.61 * deg_to_rad}
    #filters = super_to_single_angle.keys() 
    #single_to_super_angle = { filt:2.0 * np.pi - super_to_single_angle[filt] for filt in filters }
    #single_to_super_shift = { filt: [super_image_center[0] - (np.cos(single_to_super_angle[filt]) * super_image_pix_coords[filt][0] + np.sin(single_to_super_angle[filt]) * super_image_pix_coords[filt][1]),
    #                                      super_image_center[1] - (-np.sin(single_to_super_angle[filt]) * super_image_pix_coords[filt][0] + np.cos(single_to_super_angle[filt]) * super_image_pix_coords[filt][1])]
    #                                for filt in filters }
    #super_image_coord_trans = lambda x, y, center, angle: [np.cos(angle) * x + np.sin(angle) * y + center[0], -np.sin(angle) * x + np.cos(angle) * y + center[1]]
    #hard_coded_coordinate_transforms = { filt: lambda xs, ys, filt = filt: super_image_coord_trans(np.array(xs), np.array(ys), single_to_super_shift[filt], single_to_super_angle[filt]) for filt in filters }

    wcs_prefix = 'wcs_'
    
    print ('Trying to open file ' + str(target_dir + target_positions_file))
    with open(target_dir + target_positions_file, 'r') as f:
        lines = [line for line in f]
        lines = [ [float(pos) for pos in line.rstrip('\n').split(' ')] for line in lines]
        #print ('lines[0] = ' + str(lines[0])) 
        obj_nums = [line[0] for line in lines]
        x_y_positions = np.array([[line[1], line[2]] for line in lines], np.float_) 

    keywords_to_copy = ['WCSAXES', 'CTYPE1','CTYPE2','EQUINOX','LONPOLE','LATPOLE','CRVAL1',
                        'CRVAL2','CRPIX1','CRPIX2','CUNIT1','CUNIT2','CD1_1','CD1_2','CD2_1','CD2_2',
                        'IMAGEW','IMAGEH',
                        'A_ORDER','A_0_0','A_0_1','A_0_2','A_1_0','A_1_1','A_2_0',
                        'B_ORDER','B_0_0','B_0_1','B_0_2','B_1_0','B_1_1','B_2_0',
                        'AP_ORDER','AP_0_0','AP_0_1','AP_0_2','AP_1_0','AP_1_1','AP_2_0',
                        'BP_ORDER','BP_0_0','BP_0_1','BP_0_2','BP_1_0','BP_1_1','BP_2_0']

    ref_image_header = ref_dir + ref_wcs_fits_file 
    ref_hdulist = fits.open(ref_image_header)
    wcs_object = wcs.WCS(ref_hdulist[0].header)
    fake_RA_Dec_positions = wcs_object.wcs_pix2world(x_y_positions, 1)
    
    with open(target_dir + target_fake_wcs_txt_file, 'w') as f:
        for i in range(len(fake_RA_Dec_positions)):
            ra_dec = fake_RA_Dec_positions[i]
            obj_num = obj_nums[i] 
            f.write(' '.join([str(int(obj_num))] + [str(elem) for elem in ra_dec]) + '\n')
        
        

    #Do the analysis
