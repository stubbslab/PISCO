import numpy as np
import sys
import sextractorObject_py3 as SO
from cantrips import safeSortOneListByAnother
from cantrips import readInDataFromFitsFile
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sextractorObject_py3 import checkSexFlag
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib
from cantrips import flattenListOfLists
from AstronomicalParameterArchive import AstronomicalParameterArchive
import scipy.optimize as optimize
from randomSortData import randomSortListWithNans
import random

def function_to_minimize(params, flat_x_mesh, flat_y_mesh, flat_x_grad, flat_x_grad_err, flat_y_grad, flat_y_grad_err):
    normalized_diff_between_flat_arrays_with_nans = lambda expected, measured, errs: np.nansum(((expected - measured) / errs) ** 2.0)
    ax, ay, axy, ax2, ay2 = params
    grad_x_sum_of_sqrs = normalized_diff_between_flat_arrays_with_nans (ax + axy * np.array(flat_y_mesh) + 2.0 * ax2 * np.array(flat_x_mesh), flat_x_grad, flat_x_grad_err)
    grad_y_sum_of_sqrs = normalized_diff_between_flat_arrays_with_nans (ay + axy * np.array(flat_x_mesh) + 2.0 * ay2 * np.array(flat_y_mesh), flat_y_grad, flat_y_grad_err)
    sum_of_sqrs = grad_x_sum_of_sqrs + grad_y_sum_of_sqrs
    print ('For [ax, ay, axy, ax2, ay2] = ' + str(params), ' np.log(sum_of_sqrs) = ' + str(np.log(sum_of_sqrs)) )
    return np.log(sum_of_sqrs)

#Currently, we only fit the gradient to a 2nd order, color polynomaial:
# (f(x,y) = a0 + ax * x + ay * y + axy xy + ax2 x**2 + ay2 y ** 2 => grad f (x, y) = (ax + axy * y + 2 * ax2 * x, ay + axy * x + 2 * ay2 * y)
def fitColorGradientTo2DPolynomial(x_grad_img, y_grad_img, x_grad_errs_img, y_grad_errs_img,
                                   param_scalings = [10.0 ** -0.0, 10.0 ** -0.0, 10.0 ** -0.0, 10.0 ** -0.0, 10.0 ** -0.0], lower_lim_frac = 0.05, upper_lim_frac = 0.95, random_arrange_data = 0, randomize_around_zero = 0, randomize_bootstrap = 0 ):
    #randomization allows you to determine the statistcal range of the polynomial parameters
    if random_arrange_data:
        flat_x_grad = x_grad_img.flatten()
        flat_x_grad_errs = x_grad_errs_img.flatten()
        flat_y_grad = y_grad_img.flatten()
        flat_y_grad_errs = y_grad_errs_img.flatten()

        print ('randomizing x grads...')
        new_flat_x_grad, new_flat_x_grad_errs = randomSortListWithNans(flat_x_grad, flat_x_grad_errs)
        print ('randomizing y grads...')
        new_flat_y_grad, new_flat_y_grad_errs = randomSortListWithNans(flat_y_grad, flat_y_grad_errs)

        print ('Generating new (randomized) x_grad img...')
        x_grad_img = np.reshape(new_flat_x_grad, np.shape(x_grad_img))
        x_grad_errs_img = np.reshape( new_flat_x_grad_errs, np.shape(x_grad_errs_img) )
        print ('Generating new (randomized) y_grad img...')
        y_grad_img = np.reshape(new_flat_y_grad, np.shape(y_grad_img))
        y_grad_errs_img = np.reshape( new_flat_y_grad_errs, np.shape(y_grad_errs_img) )
    elif randomize_around_zero:
        x_ratios_img = x_grad_img / x_grad_errs_img
        x_ratios_centered_img = x_ratios_img - np.mediannan(x_ratios_img)
        y_ratios_img = y_grad_imt / y_grad_errs_img
        y_ratios_centered_img = y_ratios_img - np.mediannan(y_ratios_img)

        rand_x_imgs = np.random.normal(0.0, x_grad_errs_img)
        rand_y_imgs = np.random.normal(0.0, y_grad_errs_img)
        x_grad_img = rand_x_imgs
        y_grad_img = rand_y_imgs

    n_ys, n_xs = np.shape(x_grad_img)
    x_mesh, y_mesh = np.meshgrid(range(n_xs), range(n_ys))

    print ('Trimming partial x derivatives...')
    flat_x_grad_img = x_grad_img.flatten()
    flat_x_grad_errs_img = x_grad_errs_img.flatten()
    full_sorted_flat_x_grad = sorted([ flat_x_grad_img[i] / flat_x_grad_errs_img[i] for i in range(len(x_grad_img.flatten())) if not(np.isnan( flat_x_grad_img[i] )) ])
    lower_x_grad, upper_x_grad = [ full_sorted_flat_x_grad[int(lower_lim_frac * len(full_sorted_flat_x_grad))], full_sorted_flat_x_grad[int(upper_lim_frac * len(full_sorted_flat_x_grad))] ]
    scaled_x_grad_img = (x_grad_img / x_grad_errs_img )[:]
    trimmed_x_grad_img = x_grad_img[:]
    trimmed_x_grad_img[scaled_x_grad_img < lower_x_grad] = np.nan
    trimmed_x_grad_img[scaled_x_grad_img > upper_x_grad] = np.nan
    trimmed_x_grad_errs_img = x_grad_errs_img[:]
    trimmed_x_grad_errs_img[np.isnan(trimmed_x_grad_img)] = np.nan

    print ('Trimming partial y derivatives...')
    flat_y_grad_img = y_grad_img.flatten()
    flat_y_grad_errs_img = y_grad_errs_img.flatten()
    full_sorted_flat_y_grad = sorted([ flat_y_grad_img[i] / flat_y_grad_errs_img[i] for i in range(len(y_grad_img.flatten())) if not(np.isnan( flat_y_grad_img[i] )) ])
    lower_y_grad, upper_y_grad = [ full_sorted_flat_y_grad[int(lower_lim_frac * len(full_sorted_flat_y_grad))], full_sorted_flat_y_grad[int(upper_lim_frac * len(full_sorted_flat_y_grad))] ]
    scaled_y_grad_img = (y_grad_img / y_grad_errs_img )[:]
    trimmed_y_grad_img = y_grad_img[:]
    trimmed_y_grad_img[scaled_y_grad_img < lower_x_grad] = np.nan
    trimmed_y_grad_img[scaled_y_grad_img > upper_x_grad] = np.nan
    trimmed_y_grad_errs_img = y_grad_errs_img[:]
    trimmed_y_grad_errs_img[np.isnan(trimmed_y_grad_img)] = np.nan

    flat_x_mesh = x_mesh.flatten()
    flat_y_mesh = y_mesh.flatten()
    flat_x_grad = trimmed_x_grad_img.flatten()
    flat_x_grad_err = trimmed_x_grad_errs_img.flatten()
    flat_y_grad = trimmed_y_grad_img.flatten()
    flat_y_grad_err = trimmed_y_grad_errs_img.flatten()

    if randomize_bootstrap:
        print ('Performing bootstrap...')
        trimmed_flat_x_grad = [ elem for elem in flat_x_grad if not(np.isnan(elem)) ]
        trimmed_flat_x_grad_err = [ elem for elem in flat_x_grad_err if not(np.isnan(elem)) ]
        trimmed_flat_y_grad = [ elem for elem in flat_y_grad if not(np.isnan(elem)) ]
        trimmed_flat_y_grad_err = [ elem for elem in flat_y_grad_err if not(np.isnan(elem)) ]
        trimmed_flat_x_mesh = [ flat_x_mesh[i] for i in range(len(flat_x_mesh)) if not(np.isnan(flat_x_grad[i])) ]
        trimmed_flat_y_mesh = [ flat_y_mesh[i] for i in range(len(flat_y_mesh)) if not(np.isnan(flat_x_grad[i])) ]
        for j in range(len(flat_x_mesh)):
            rand_index = random.randint(0, len(trimmed_flat_x_grad) - 1)
            if j % 100000 == 0: print ( 'For index j = ' + str(j) + ', random index is ' + str(rand_index) )
            if not(np.isnan(flat_x_grad[j])): flat_x_mesh[j] = trimmed_flat_x_mesh[rand_index]
            if not(np.isnan(flat_x_grad[j])): flat_y_mesh[j] = trimmed_flat_y_mesh[rand_index]
            if not(np.isnan(flat_x_grad[j])): flat_x_grad[j] = trimmed_flat_x_grad[rand_index]
            if not(np.isnan(flat_x_grad_err[j])): flat_x_grad_err[j] = trimmed_flat_x_grad_err[rand_index]
            if not(np.isnan(flat_y_grad[j])): flat_y_grad[j] = trimmed_flat_y_grad[rand_index]
            if not(np.isnan(flat_y_grad_err[j])): flat_y_grad_err[j] = trimmed_flat_y_grad_err[rand_index]

    #plt.hist([elem for elem in flat_x_grad / flat_x_grad_err if not(np.isnan(elem)) ], bins = 1001)
    #plt.title('Hist of trimmed x gradients to be fit')
    #plt.show()
    #plt.hist([elem for elem in flat_y_grad / flat_y_grad_err if not(np.isnan(elem)) ], bins = 1001)
    #plt.title('Hist of trimmed y gradients to be fit')
    #plt.show()


    #function_to_minimize = lambda ax, ay, axy, ax2, ay2: ( normalized_diff_between_flat_arrays_with_nans (ax + axy * np.array(flat_y_mesh) + 2.0 * ax2 * np.array(flat_x_mesh), flat_x_grad, flat_x_grad_err) +
    #                                                                        normalized_diff_between_flat_arrays_with_nans (ay + axy * np.array(flat_x_mesh) + 2.0 * ay2 * np.array(flat_y_mesh), flat_y_grad, flat_y_grad_err) )
    local_function_to_minimize = lambda params: function_to_minimize([params[i] * param_scalings[i] for i in range(len(params))], flat_x_mesh, flat_y_mesh, flat_x_grad, flat_x_grad_err, flat_y_grad, flat_y_grad_err)

    print ('Working on minimization... ')
    min_results = optimize.minimize(local_function_to_minimize, [0.0, 0.0, 0.0, 0.0, 0.0])

    return min_results

def computePartialDerivFromDirectionalDeriv(directD1, directD2, vec1, vec2, direction_index):
    if direction_index == 0: perp_index = 1
    else: perp_index = 0
    mask_section = (1 - np.isnan(directD1)) * (1 - np.isnan(directD2))
    mask_section = 1 - mask_section
    numerator = np.nanprod([directD1, vec2[perp_index]], axis = 0) - np.nanprod([directD2, vec1[perp_index]], axis = 0)
    numerator[mask_section == 1] = np.nan

    denominator = np.nanprod([vec1[direction_index], vec2[perp_index]], axis = 0) - np.nanprod([vec2[direction_index], vec1[perp_index]], axis = 0)
    denominator[mask_section == 1] = np.nan

    partial = np.nanprod([numerator, 1.0 / denominator], axis = 0)
    partial[mask_section == 1] = np.nan

    return partial

def computePartialDerivErrFromDirectionalDerivErr(directD1_err, directD2_err, vec1, vec2, direction_index):
    if direction_index == 0: perp_index = 1
    else: perp_index = 0
    mask_section = (1 - np.isnan(directD1_err)) * (1 - np.isnan(directD2_err))
    mask_section = 1 - mask_section
    numerator = (np.nanprod([directD1_err, vec2[perp_index]], axis = 0)) ** 2.0 + (np.nanprod([directD2_err, vec1[perp_index]], axis = 0)) ** 2.0
    numerator[mask_section == 1] = np.nan

    denominator = (np.nanprod([vec1[direction_index], vec2[perp_index]], axis = 0) - np.nanprod([vec2[direction_index], vec1[perp_index]], axis = 0)) ** 2.0
    denominator[mask_section == 1] = np.nan

    partial_err = np.sqrt(np.nanprod([numerator, 1.0 / denominator], axis = 0))
    partial_err[mask_section == 1] = np.nan

    return partial_err

def UpdateDerivs(current_partial, current_partial_errs, new_partial, new_partial_errs, additional_masks = [1, 1]):


    stack_mask_components = [np.isnan(new_partial), np.isnan(current_partial)]
    stack_mask = 1
    print ('Updating mask... ' )
    for i in range(len(stack_mask_components)):
        stack_mask_component = stack_mask_components[i]
        additional_mask = additional_masks[i]
        stack_mask = stack_mask * stack_mask_component * additional_mask

    print ('Updating numerator...' )
    weighted_current_partial = np.nanprod([current_partial, 1.0 / (current_partial_errs ** 2.0) ], axis = 0)
    weighted_new_partial = np.nanprod([new_partial, 1.0 / (new_partial_errs ** 2.0) ], axis = 0)
    weighted_current_partial[np.isnan(current_partial)] = np.nan
    weighted_new_partial[np.isnan(new_partial)] = np.nan
    partial_stack_numerator = np.nansum([weighted_current_partial, weighted_new_partial], axis = 0)
    print ('Updating denominator...' )
    partial_stack_denominator = np.nansum([1.0 / (current_partial_errs ** 2.0), 1.0 / (new_partial_errs ** 2.0)], axis = 0)
    partial_stack_numerator[stack_mask == 1] = np.nan
    partial_stack_denominator[stack_mask == 1] = np.nan
    print ('Generating dx...')
    current_partial = np.nanprod([partial_stack_numerator, 1.0 / partial_stack_denominator], axis = 0)
    current_partial[stack_mask == 1] = np.nan
    current_partial_errs = np.sqrt(1.0 / partial_stack_denominator)

    return [current_partial, current_partial_errs]

def measureColorGradientFromDirectionalDerivatives(color_slopes, color_stds, color_directions):
    if not(len(color_slopes) == len(color_stds)) or not(len(color_stds) == len(color_directions)):
        print ('Number of images of colors slopes, color stds, and color directions must be the same.  They are not.  Retuning 1. ')
        return 1
    else:
        n_stacks = len(color_slopes)
        n_paired_stacks = int( (n_stacks -1) ** 2.0 / 2 + (n_stacks -1) / 2 )
        img_size = np.shape(color_slopes[0])
        partial_xs = [0 for i in range(n_stacks)]
        partial_x_errs = [0 for i in range(n_stacks)]
        partial_ys = [0 for i in range(n_stacks)]
        partial_y_errs = [0 for i in range(n_stacks)]
        n_stacks_in = 0
        current_dx = None
        current_dy = None
        current_dx_err = None
        current_dy_err = None
        min_angle_sep = 0.0 * np.pi / 4.0
        for i in range(n_stacks):
            color_slope1 = color_slopes[i]
            color_std1 = color_stds[i]
            color_directions1 = color_directions[i]
            vec1 = [1.0 / np.sqrt(1.0 + color_directions1 ** 2.0) * (np.zeros(np.shape(color_directions1)) + 1.0) ,
                    1.0 / np.sqrt(1.0 + color_directions1 ** 2.0) * (color_directions1) ]
            for j in range(i+1, n_stacks):
                n_stacks_in = n_stacks_in + 1
                print ('Working on stack combination: ' + str([i,j]) + ', which is ' + str(n_stacks_in) + ' of ' + str(n_paired_stacks) )
                color_slope2 = color_slopes[j]
                color_std2 = color_stds[j]
                color_directions2 = color_directions[j]
                vec2 = [1.0 / np.sqrt(1.0 + color_directions2 ** 2.0) * (np.zeros(np.shape(color_directions2)) + 1.0) ,
                        1.0 / np.sqrt(1.0 + color_directions2 ** 2.0) * (color_directions2) ]
                sep_angles = np.arccos(abs(vec1[0] * vec2[0] + vec1[1] * vec2[1]))
                sep_angle_mask = (sep_angles < min_angle_sep) * np.isnan(sep_angles)


                partial_x = computePartialDerivFromDirectionalDeriv(color_slope1, color_slope2, vec1, vec2, 0)
                partial_x_err = computePartialDerivErrFromDirectionalDerivErr(color_std1, color_std2, vec1, vec2, 0)
                #partial_xs = partial_xs + [partial_x]
                #partial_x_errs = partial_x_errs + [partial_x_err]
                partial_y = computePartialDerivFromDirectionalDeriv(color_slope1, color_slope2, vec1, vec2, 1)
                partial_y_err = computePartialDerivErrFromDirectionalDerivErr(color_std1, color_std2, vec1, vec2, 1)
                #partial_ys = partial_ys + [partial_y]
                #partial_y_errs = partial_y_errs + [partial_y_err]

                if current_dx is None:
                    current_dx = partial_x
                    current_dx_err = partial_x_err
                else:
                    current_dx, current_dx_err = UpdateDerivs(current_dx, current_dx_err, partial_x, partial_x_err, additional_masks = [1, sep_angle_mask])
                if current_dy is None:
                    current_dy = partial_y
                    current_dy_err = partial_y_err
                else:
                    current_dy, current_dy_err = UpdateDerivs(current_dy, current_dy_err, partial_y, partial_y_err, additional_masks = [1, sep_angle_mask])


        print ('n_paired_stacks = ' + str(n_paired_stacks))
        stack_mask_components = [np.isnan(partial_xs[i]) for i in range(n_paired_stacks)]
        stack_mask = 1
        print ('Generating mask for pixels where less than 2 observations were made... ' )
        for stack_mask_component in stack_mask_components:
            stack_mask = stack_mask * stack_mask_component

        return [current_dx, current_dy, current_dx_err, current_dy_err]


def fitColorsTo1DFunction(master_cat, theoretical_fit_funct,
                          targ_objs = 'all', filt_diffs = ['g','r'], zero_to_first = 1, overall_names = 'all', zenith_key_str = 'ZD',
                          min_n_images_for_use = 10, err_estimation = 'std', obj_position_key_str_prefix = 'xy_', x_divide = 763, divide_buffer = 20 ):

    astro_arch = AstronomicalParameterArchive()
    deg_to_rad = astro_arch.getDegToRad()

    if targ_objs is 'all':
        targ_objs = list(master_cat.master_val_dict.keys())

    if overall_names is 'all':
        overall_names = list(master_cat.master_val_dict[targ_objs[0]].keys())

    overall_names = sorted(overall_names)

    master_cat.measureColor(overall_names = overall_names, filters = filt_diffs)
    for considered_filter in filt_diffs:
        master_cat.measurePositionInRefFilter(ref_filter = considered_filter, overall_names = overall_names)

    color_key_str = filt_diffs[0] + '-' + filt_diffs[1]

    good_objs = []
    for targ_obj in targ_objs:
        obs_quality_lists = [[master_cat.master_val_dict[targ_obj][name]['good_' + filt] for name in overall_names] for filt in filt_diffs]
        obs_quality_list = list(np.array(obs_quality_lists[0]) * np.array(obs_quality_lists[1]))
        n_good_objs = len([1 for obs_quality in obs_quality_list if obs_quality])
        if n_good_objs >= min_n_images_for_use:
            good_objs = good_objs + [targ_obj]

    print ('Have a total of ' + str(len(good_objs)) + ' good objects. ')

    colors_by_obj = []
    stds_by_obj = []
    pos_by_obj = []


    all_colors = [[] for obj in good_objs]
    all_stds = [-1 for obj in good_objs]
    all_positions = [[[-1, -1], [-1, -1]] for obj in good_objs]
    all_zeniths = [[] for obj in good_objs]
    all_n_atmospheres = [[] for obj in good_objs]
    print ('Pulling positions, zeniths, and colors and measuring stds and medians for each good object...' )
    for i in range(len(good_objs)):
        good_obj = good_objs[i]
        obs_quality_lists = [[master_cat.master_val_dict[good_obj][name]['good_' + filt] for name in overall_names] for filt in filt_diffs]
        obs_quality_list = list(np.array(obs_quality_lists[0]) * np.array(obs_quality_lists[1]))
        obj_positions = [[master_cat.master_val_dict[good_obj][overall_names[i]][obj_position_key_str_prefix + filt] for filt in filt_diffs] for i in range(len(overall_names))
                         if obs_quality_list[i] ]
        obj_colors = [master_cat.master_val_dict[good_obj][overall_names[i]][color_key_str][0] for i in range(len(overall_names))
                         if obs_quality_list[i] ]
        obj_zeniths = [float(master_cat.object_ind_quantities[overall_names[i]][zenith_key_str]) for i in range(len(overall_names))
                         if obs_quality_list[i] ]
        obj_n_atmospheres = [1.0 / np.cos(zenith * deg_to_rad ) for zenith in obj_zeniths]
        obj_fully_left = [filt_pair_position[0][0] <= x_divide - divide_buffer and filt_pair_position[1][0] <= x_divide - divide_buffer for filt_pair_position in obj_positions]
        obj_fully_right = [filt_pair_position[0][0] > x_divide + divide_buffer and filt_pair_position[1][0] > x_divide + divide_buffer for filt_pair_position in obj_positions]
        left_obj_colors = [obj_colors[i] for i in range(len(obj_colors)) if obj_fully_left[i] and not np.isnan(obj_colors[i])]
        right_obj_colors = [obj_colors[i] for i in range(len(obj_colors)) if obj_fully_right[i] and not np.isnan(obj_colors[i])]
        left_obj_mean = np.mean(left_obj_colors)
        right_obj_mean = np.mean(right_obj_colors)
        left_right_std = np.sqrt((np.sum([(color - left_obj_mean) ** 2.0 for color in left_obj_colors]) + np.sum([(color - right_obj_mean) ** 2.0 for color in right_obj_colors]))
                                / len(left_obj_colors + right_obj_colors))
        raw_std = np.std([color for color in obj_colors if not np.isnan(color) ])
        all_colors[i] = obj_colors
        all_stds[i] = [left_right_std, raw_std]
        all_positions[i] = obj_positions
        all_zeniths[i] = obj_zeniths
        all_n_atmospheres[i] = obj_n_atmospheres

    all_medians = [ np.median(obj_colors) for obj_colors in all_colors ]
    all_zerod_colors = [np.array(all_colors[i]) - all_medians[i] for i in range(len(all_colors)) ]
    print ('Done.')

    medianCorrection = lambda *args: flattenListOfLists([ (np.array(all_zerod_colors[i]) + np.median(theoretical_fit_funct(np.array(all_n_atmospheres[i]), *args))).tolist() for i in range(len(all_colors)) ])


    function_to_minimize = lambda xs, a, b: xs * a + b + medianCorrection(xs, a, b)
    function_to_minimize = lambda xs, *args: theoretical_fit_funct(xs, *args) + medianCorrection(*args)

    return all_n_atmospheres, all_positions, all_zerod_colors, all_stds, function_to_minimize

#def getInterpolated

def getInterpolatedColorSlopeAlongLine(line_xs, line_y, mean_x, mean_y, left_end, right_end, obj_lin_fit, color_fit, color_fit_err, line_inference_dist = 1.0, uncertainty_type = 'poly_covar' ):
    #Makes sure that the positions in question are formatted as 1-d numpy array
    if (type(line_xs) in [list]):
        line_xs = np.array(line_xs)
    elif not (type(line_xs) in [np.ndarray]):
        line_xs = np.array([line_xs])
    n_cols = len(line_xs)
    centered_xs = line_xs - mean_x
    centered_y = line_y - mean_y
    slope, intercept = obj_lin_fit
    closest_xs = ( (line_y + 1.0 / slope  * line_xs - intercept) / (slope + 1.0 / slope) )
    closest_ys = slope * closest_xs + intercept
    closest_dists_sqr = (((closest_xs - mean_x) ** 2.0 + (closest_ys - mean_y) ** 2.0) * (closest_xs - mean_x) / abs(closest_xs - mean_x))
    left_dist_sqr = ((left_end[0] - mean_x) ** 2.0 + (left_end[1] - mean_y) ** 2.0) * (left_end[0] - mean_x) / abs(left_end[0] - mean_x)
    right_dist_sqr = ((right_end[0] - mean_x) ** 2.0 + (right_end[1] - mean_y) ** 2.0) * (right_end[0] - mean_x) / abs(right_end[0] - mean_x)
    interp_xs = list(range(len(line_xs)))
    interp_ys = list(range(len(line_xs)))
    for i in range(n_cols):
        closest_dist_sqr = closest_dists_sqr[i]
        closest_x = closest_xs[i]
        closest_y = closest_ys[i]
        interp_xs[i] =  (left_end[0] if closest_dist_sqr < left_dist_sqr else right_end[0] if closest_dist_sqr > right_dist_sqr else closest_x )
        interp_ys[i] =  (left_end[1] if closest_dist_sqr < left_dist_sqr else right_end[1] if closest_dist_sqr > right_dist_sqr else closest_y )

    point_to_line_dists = [ (interp_xs[i] - line_xs[i] ) ** 2.0 + (interp_ys[i] - line_y ) ** 2.0 for i in range(n_cols) ]
    close_enough = [ point_to_line_dist <= line_inference_dist ** 2.0 for point_to_line_dist in point_to_line_dists ]
    dists_along_line = [ np.sqrt((interp_xs[i] - mean_x) ** 2.0 + (interp_ys[i] - mean_y) ** 2.0) * (interp_xs[i] - mean_x) / abs(interp_xs[i] - mean_x) if close_enough[i] else np.nan for i in range(len(close_enough)) ]

    fit_degree = len(color_fit) - 1

    interp_color_slopes = [ sum([(fit_degree - i) * dists_along_line[j] ** (fit_degree - i - 1) * color_fit[i] for i in range(fit_degree + 1)]) if close_enough[j] else np.nan for j in range(len(dists_along_line)) ]
    if uncertainty_type in ['poly_covar']:
        interp_color_slope_errs = [ sum([(fit_degree - i) * dists_along_line[j] ** (fit_degree - i - 1) * color_fit_err[i] for i in range(fit_degree + 1)]) if close_enough[j] else np.nan for j in range(len(dists_along_line)) ]
    else:
        interp_color_slope_errs = [ color_fit_err for dist in dists ]
    #interp_color_slope_errs = [ -1 if close_enough[j] else np.nan for j in range(len(dists_along_line)) ]
    interp_colors = [ sum([ dists_along_line[j] ** (fit_degree - i) * color_fit[i] for i in range(fit_degree+1)]) if close_enough[j] else np.nan for j in range(len(dists_along_line)) ]


    #return interp_color_slopes, interp_color_slope_errs, slope
    return interp_color_slopes, interp_color_slope_errs, slope


def MeasureFitSlopeOverObject(obj_position_function, size = [1500, 2140], line_inference_dist = 1.0, uncertainty_type = 'poly_covar' ):
    left_end_index = 2
    right_end_index = 3
    x_lims = [int(obj_position_function[2][0]), int(obj_position_function[3][0] + 1)]
    x_lims = sorted(x_lims)
    x_lims = [int(max(1.0, x_lims[0] - line_inference_dist)), int(min(size[0], x_lims[1] + line_inference_dist + 1))]
    y_lims = [int(obj_position_function[2][1]), int(obj_position_function[3][1] + 1)]
    y_lims = sorted(y_lims)
    y_lims = [int(max(1.0, y_lims[0] - line_inference_dist)), int(min(size[1], y_lims[1] + line_inference_dist + 1))]
    color_slope_frame = np.empty(size)
    color_slope_frame[:] = np.nan
    color_slope_err_frame = np.empty(size)
    color_slope_err_frame[:] = np.nan

    for y in range(y_lims[0], y_lims[1] + 1):
        new_line_slopes, new_line_slope_errs, direction = getInterpolatedColorSlopeAlongLine(list(range(x_lims[0], x_lims[1] + 1)), y, *obj_position_function, line_inference_dist = line_inference_dist, uncertainty_type = uncertainty_type)
        color_slope_frame[x_lims[0]-1:x_lims[1], y-1] = new_line_slopes
        color_slope_err_frame[x_lims[0]-1:x_lims[1], y-1] = new_line_slope_errs

    return [color_slope_frame, color_slope_err_frame, direction]

def stackObjectPositions(obj_position_functions_set, filt_diffs = ['g','r'], pos_ref_filter = 'g', save = 0, save_dir = '/Users/sasha/Documents/Harvard/physics/stubbs/PISCO/stackedColorSlopes/', save_file_name = 'default.npy'):

    color_key_str = filt_diffs[0] + '-' + filt_diffs[1]

    first_obj_color_slope_data = MeasureFitSlopeOverObject(obj_position_functions_set[0][list(obj_position_functions_set[0].keys())[0]], line_inference_dist = 5.0 )
    current_color_slope_combine = first_obj_color_slope_data[0]
    if type(first_obj_color_slope_data[1]) in [list, np.ndarray]:
        current_weights_img = (first_obj_color_slope_data[1]) ** (-2.0)
    else:
        current_weights_img = ((np.zeros(np.shape(current_color_slope_combine)) + first_obj_color_slope_data[1]) ** (-2.0)) * (1 - np.isnan(current_color_slope_combine))
        current_weights_img[current_weights_img == 0.0] = np.nan
    if type (first_obj_color_slope_data[2]) in [list, np.ndarray]:
        current_direction_img = first_obj_color_slope_data[2]
    else:
        current_direction_img = np.zeros(np.shape(current_color_slope_combine)) + first_obj_color_slope_data[2]
        current_direction_img[np.isnan(current_color_slope_combine)] = np.nan

    for obj_position_functions in obj_position_functions_set:
        good_objs = list(obj_position_functions.keys())
        for i in range(1, len(good_objs)):
            good_obj = good_objs[i]
            print ('Adding object ' + str(good_obj) + ' (' + str(i) + ' of ' + str(len(good_objs)) + ') to the stack. ')
            obj_color_slope_data = MeasureFitSlopeOverObject(obj_position_functions[good_obj], line_inference_dist = 5.0 )
            new_color_slope = obj_color_slope_data[0]
            #interpolated_color_slope[i] = obj_color_slope_data[0]
            if type(obj_color_slope_data[1]) in [list, np.ndarray]:
                new_weights_img = (obj_color_slope_data[1]) ** (-2.0)
            else:
                new_weights_img = ((np.zeros(np.shape(new_color_slope)) + obj_color_slope_data[1]) ** (-2.0) )  * (1 - np.isnan(new_color_slope))
                new_weights_img[new_weights_img == 0.0] = np.nan

            if type(obj_color_slope_data[2]) in [list, np.ndarray]:
                new_direction_img = obj_color_slope_data[2]
            else:
                new_direction_img = np.zeros(np.shape(new_color_slope)) + obj_color_slope_data[2]
                new_direction_img[np.isnan(new_color_slope)] = np.nan

            new_combined_weights = np.nansum([current_weights_img, new_weights_img], axis = 0)
            new_combined_weights[new_combined_weights == 0.0] = np.nan

            scaled_current_color_slope_img = np.nanprod([ current_color_slope_combine, current_weights_img, 1.0 / new_combined_weights ], axis = 0)
            scaled_current_color_slope_img[np.isnan(current_color_slope_combine)] = np.nan #Any spot where the original had no values should still have no values
            scaled_new_color_slope_img = np.nanprod([ new_color_slope, new_weights_img, 1.0 / new_combined_weights ], axis = 0)
            scaled_new_color_slope_img[np.isnan(new_color_slope)] = np.nan #Any spot where the original had no values should still have no values
            current_color_slope_combine = np.nansum([scaled_current_color_slope_img, scaled_new_color_slope_img], axis = 0)
            current_color_slope_combine[ np.isnan(scaled_current_color_slope_img) * np.isnan(scaled_new_color_slope_img) ] = np.nan

            scaled_current_direction_img = np.nanprod([ current_direction_img, current_weights_img, 1.0 / new_combined_weights ], axis = 0)
            scaled_current_direction_img[np.isnan(current_direction_img)] = np.nan #Any spot where the original had no values should still have no values
            scaled_new_direction_img = np.nanprod([ new_direction_img, new_weights_img, 1.0 / new_combined_weights ], axis = 0)
            scaled_new_direction_img[np.isnan(new_direction_img)] = np.nan #Any spot where the original had no values should still have no values
            current_direction_img = np.nansum([scaled_new_direction_img, scaled_current_direction_img], axis = 0)
            current_direction_img[ np.isnan(scaled_current_direction_img) * np.isnan(scaled_new_direction_img) ] = np.nan

            current_weights_img = new_combined_weights[:]

    if save:
        np.save(save_dir + save_file_name, np.array([current_color_slope_combine, 1.0 / np.sqrt(current_weights_img), current_direction_img]))

    return [current_color_slope_combine, 1.0 / np.sqrt(current_weights_img), current_direction_img]

class MasterCatalogue:

    def MeasureFitSlopeOverObjects(self, line_inference_dist = 1.0 ):
        for obj in self.master_objects:
            single_obj = MeasureFitSlopeOverObject(obj_position_functions[16], line_inference_dist = line_inference_dist)
        groupObject = combine_median_of_slopes(single_objs)
        return

    def FitColor(self, targ_objs = 'all', overall_names = 'all', filt_diffs = ['g','r'], pos_ref_filter = 'g',
                 obj_position_key_str_prefix = 'xy_', obj_fit_prefix = 'obj_line_movement_', fit_degree = 2,
                 x_divide = 763, divide_buffer = 20):
        min_n_images_for_use = fit_degree + 3 + 1
        if targ_objs is 'all':
            targ_objs = list(self.master_val_dict.keys())

        if overall_names is 'all':
            overall_names = list(self.master_val_dict[targ_objs[0]].keys())

        overall_names = sorted(overall_names)

        self.measureColor(overall_names = overall_names, filters = filt_diffs)
        for considered_filter in filt_diffs:
            self.measurePositionInRefFilter(ref_filter = considered_filter, overall_names = overall_names)

        color_key_str = filt_diffs[0] + '-' + filt_diffs[1]

        good_objs = []
        good_objs_usable_images_list = []

        for targ_obj in targ_objs:
            obs_quality_lists = [[self.master_val_dict[targ_obj][name]['good_' + filt] for name in overall_names] for filt in filt_diffs]
            obs_quality_list = list(np.array(obs_quality_lists[0]) * np.array(obs_quality_lists[1]))
            obj_positions = [[self.master_val_dict[targ_obj][name][obj_position_key_str_prefix + filt] for filt in filt_diffs]  for name in overall_names]
            obj_fully_left = [filt_pair_position[0][0] <= x_divide - divide_buffer and filt_pair_position[1][0] <= x_divide - divide_buffer for filt_pair_position in obj_positions]
            obj_fully_right = [filt_pair_position[0][0] > x_divide + divide_buffer and filt_pair_position[1][0] > x_divide + divide_buffer for filt_pair_position in obj_positions]
            obj_in_good_position_list = [obj_fully_left[i] or obj_fully_right[i] for i in range(len(overall_names)) ]
            obs_usable_list = list(np.array(obs_quality_list) * np.array(obj_in_good_position_list))
            n_good_objs = len([1 for usable in obs_usable_list if usable])
            if n_good_objs >= min_n_images_for_use:
                #print ('Adding obj ' + str(targ_obj) + ' to list of good_objs.')
                good_objs = good_objs + [targ_obj]
                good_objs_usable_images_list = good_objs_usable_images_list + [ obs_usable_list ]

        print ('good_objs = ' + str(good_objs))

        obj_positions_key_str = obj_position_key_str_prefix + pos_ref_filter
        obj_position_functions = {}
        for n_obj in range(len(good_objs)):
            good_obj = good_objs[n_obj]
            obj_usable_images_list = good_objs_usable_images_list[n_obj]
            n_frames = len(overall_names)
            obj_positions = [ self.master_val_dict[good_obj][overall_names[i]][obj_positions_key_str] for i in range(n_frames) if obj_usable_images_list[i] ]
            obj_xs, obj_ys = [ [pos[0] for pos in obj_positions if pos], [pos[1] for pos in obj_positions] ]
            if len(obj_xs) > fit_degree + 3:
                min_x, max_x = [min(obj_xs), max(obj_xs)]
                left_end = [min_x, obj_ys[np.argmin(obj_xs)]]
                right_end = [max_x, obj_ys[np.argmax(obj_xs)]]
                mean_x, mean_y = [np.mean(obj_xs), np.mean(obj_ys) ]

                centered_xs, centered_ys = [[x - mean_x for x in obj_xs], [y - mean_y for y in obj_ys]]
                obj_lin_fit = list(np.polyfit(obj_xs, obj_ys, 1))
                obj_fit_xs = obj_xs
                obj_fit_ys = [obj_lin_fit[0] * x + obj_lin_fit[1] for x in obj_fit_xs]
                dists_along_line = [ np.sqrt((obj_fit_xs[i] - mean_x) ** 2.0 + (obj_fit_ys[i] - (obj_lin_fit[0] * mean_x + obj_lin_fit[1])) ** 2.0) * (obj_fit_xs[i] - mean_x) / abs(obj_fit_xs[i] - mean_x) for i in range(len(centered_xs)) ]
                colors_along_line = [ self.master_val_dict[good_obj][overall_names[i]][color_key_str][0] for i in range(n_frames) if obj_usable_images_list[i] ]

                color_fit, color_cov = np.polyfit(dists_along_line, colors_along_line, fit_degree, cov = True)

                color_fit_errs = [color_cov[i,i] for i in range(fit_degree + 1)]
                direct_color_std = np.std(colors_along_line)
                interp_color_funct_params = [mean_x, mean_y, left_end, right_end, obj_lin_fit, color_fit]

                interpolated_left_end = [ left_end[0], obj_lin_fit[0] * left_end[0] + obj_lin_fit[1] ]
                interpolated_right_end = [ right_end[0], obj_lin_fit[0] * right_end[0] + obj_lin_fit[1] ]
                obj_position_functions[good_obj] = [mean_x, mean_y, interpolated_left_end, interpolated_right_end, obj_lin_fit, color_fit, color_fit_errs]
                #obj_position_functions[good_obj] = [mean_x, mean_y, interpolated_left_end, interpolated_right_end, obj_lin_fit, color_fit, direct_color_std]
            else:
                obj_position_functions[good_obj] = np.nan

        self.obj_position_functions_of_color[color_key_str] = obj_position_functions

        return obj_position_functions

    def MeasureDitherDirection(self, targ_objs = 'all', overall_names = 'all', ref_filter = 'g', obj_position_key_str_prefix = 'xy_', obj_fit_prefix = 'obj_line_movement_'):
        if targ_objs is 'all':
            targ_objs = list(self.master_val_dict.keys())

        if overall_names is 'all':
            overall_names = list(self.master_val_dict[targ_objs[0]].keys())

        overall_names = sorted(overall_names)

        self.measurePositionInRefFilter(ref_filter = ref_filter, overall_names = overall_names)
        obj_positions_key_str = obj_position_key_str_prefix + ref_filter
        obj_position_functions = {}
        for targ_obj in targ_objs:
            obj_positions = [self.master_val_dict[targ_obj][overall_name][obj_positions_key_str] for overall_name in overall_names]
            obj_positions = [pos for pos in obj_positions if not(np.isnan(pos[0])) and not(np.isnan(pos[1])) ]
            obj_xs, obj_ys = [ [pos[0] for pos in obj_positions if pos], [pos[1] for pos in obj_positions] ]
            print ('targ_obj = ' + str(targ_obj))
            print ('obj_xs = ' + str(obj_xs))
            print ('obj_ys = ' + str(obj_ys))
            if len(obj_xs) > 1 and len(obj_ys) > 1:
                obj_lin_fit = list(np.polyfit(obj_xs, obj_ys, 1))
                obj_position_functions[targ_obj] = obj_lin_fit + [min(obj_xs), max(obj_xs)]
            else:
                obj_position_functions[targ_obj] = np.nan
        self.addFrameIndependentQuantities(master_objs = 'all', quantity_key_strs = [obj_fit_prefix + ref_filter], quantity_values_by_obj = [obj_position_functions])

        return 1


    def PlotPositionsAcrossImage(self, targ_objs = 'all', filt_diffs = ['g','r'], zero_to_first = 1, color_lims = [-0.5, 0.5], x_lims = [0, 1500], y_lims = [0, 2140],
                                 overall_names = 'all', plot_type = 'plane', max_objs_to_plot = 'all', color_scheme = 'scaled', zenith_key_str = 'ZD',
                                 min_n_images_for_use = 10, n_best_objs = 20, how_to_determine_best = 'std', obj_position_key_str_prefix = 'xy_', x_divide = 763, divide_buffer = 20,
                                 single_obj_colors = ['blue','red','orange','green','yellow','cyan','purple','grey','black','pink']):

        if targ_objs is 'all':
            targ_objs = list(self.master_val_dict.keys())

        if overall_names is 'all':
            overall_names = list(self.master_val_dict[targ_objs[0]].keys())

        overall_names = sorted(overall_names)

        self.measureColor(overall_names = overall_names, filters = filt_diffs)
        for considered_filter in filt_diffs:
            self.measurePositionInRefFilter(ref_filter = considered_filter, overall_names = overall_names)

        color_key_str = filt_diffs[0] + '-' + filt_diffs[1]

        good_objs = []

        for targ_obj in targ_objs:
            obs_quality_lists = [[self.master_val_dict[targ_obj][name]['good_' + filt] for name in overall_names] for filt in filt_diffs]
            obs_quality_list = list(np.array(obs_quality_lists[0]) * np.array(obs_quality_lists[1]))
            n_good_objs = len([1 for obs_quality in obs_quality_list if obs_quality])
            if n_good_objs >= min_n_images_for_use:
                good_objs = good_objs + [targ_obj]

        stds = [-1.0 for good_obj in good_objs]
        median_color_errs = [-1.0 for good_obj in good_objs]
        for i in range(len(good_objs)):
            good_obj = good_objs[i]
            obj_colors = [self.master_val_dict[good_obj][image_name][color_key_str][0] for image_name in overall_names]
            #obj_zeniths = [self.object_ind_quantities[image_name][zenith_key_str] for image_name in overall_names]
            obj_color_errs = [self.master_val_dict[good_obj][image_name][color_key_str][1] for image_name in overall_names]
            #Because of the lack of the amplifier dividing the image in two, we can only safely measured stds for those observations where
            # an object is on the same side of that divide in both filters
            obj_positions = [[self.master_val_dict[good_obj][image_name][obj_position_key_str_prefix + filt] for filt in filt_diffs]  for image_name in overall_names]
            obj_fully_left = [filt_pair_position[0][0] <= x_divide - divide_buffer and filt_pair_position[1][0] <= x_divide - divide_buffer for filt_pair_position in obj_positions]
            obj_fully_right = [filt_pair_position[0][0] > x_divide + divide_buffer and filt_pair_position[1][0] > x_divide + divide_buffer for filt_pair_position in obj_positions]
            left_obj_colors = [obj_colors[i] for i in range(len(obj_colors)) if obj_fully_left[i] and not np.isnan(obj_colors[i])]
            right_obj_colors = [obj_colors[i] for i in range(len(obj_colors)) if obj_fully_right[i] and not np.isnan(obj_colors[i])]
            left_obj_mean = np.mean(left_obj_colors)
            right_obj_mean = np.mean(right_obj_colors)
            left_right_std = np.sqrt((np.sum([(color - left_obj_mean) ** 2.0 for color in left_obj_colors]) + np.sum([(color - right_obj_mean) ** 2.0 for color in right_obj_colors]))
                                / len(left_obj_colors + right_obj_colors))
            raw_std = np.std([color for color in obj_colors if not np.isnan(color) ])
            stds[i] = [left_right_std, raw_std]
            median_color_err = np.median(obj_color_errs)
            median_color_errs[i] = median_color_err

        if how_to_determine_best.lower() in ['std', 'stds']:
            print ('Determining best images by inherrent standard deviation. ')
            sorted_stds, sorted_color_errs, sorted_good_objs = safeSortOneListByAnother([std_pair[0] for std_pair in stds], [[std_pair[0] for std_pair in stds], median_color_errs, good_objs])
        else:
            print ('Determining best images by reported magnitude errors. ')
            sorted_stds, sorted_color_errs, sorted_good_objs = safeSortOneListByAnother(median_color_errs, [[std_pair[0] for std_pair in stds], median_color_errs, good_objs])

        best_stds = [0 for i in range(n_best_objs)]
        best_errs = [0 for i in range(n_best_objs)]
        best_objs = [0 for i in range(n_best_objs)]
        for i in range(n_best_objs):
            best_stds[i] = sorted_stds[i]
            best_objs[i] = sorted_good_objs[i]
            best_errs[i] = sorted_color_errs[i]

        print ('good_objs = ' + str(good_objs))
        print ('best_objs = ' + str(best_objs))
        if max_objs_to_plot is 'all': max_objs_to_plot = len(best_objs)

        if plot_type is 'plot3D':
            fig = plt.figure()
            #norm = colors.Normalize(vmin=color_lims[0], vmax=color_lims[1])
            ax = fig.add_subplot(111, projection = '3d')
            ax.set_zlim(color_lims)
        else:
            fig, ax = plt.subplots()

        scats = []
        for i in range(len(best_objs[0:max_objs_to_plot])):
            best_obj = best_objs[0:max_objs_to_plot][i]
            obs_quality_lists = [[self.master_val_dict[best_obj][name]['good_' + filt] for name in overall_names] for filt in filt_diffs]
            obs_quality_list = list(np.array(obs_quality_lists[0]) * np.array(obs_quality_lists[1]))

            xs = [self.master_val_dict[best_obj][overall_names[i]]['xy_' + self.pos_filter][0] for i in range(len(overall_names))
                  if obs_quality_list[i]]
            ys = [self.master_val_dict[best_obj][overall_names[i]]['xy_' + self.pos_filter][1] for i in range(len(overall_names))
                  if obs_quality_list[i]]
            colors = [self.master_val_dict[best_obj][overall_names[i]][color_key_str][0] for i in range(len(overall_names))
                  if obs_quality_list[i]]
            zeniths = [self.object_ind_quantities[overall_names[i]][zenith_key_str] for i in range(len(overall_names))
                  if obs_quality_list[i]]

            if zero_to_first:
                color_zero_point = np.median(colors)
                colors = [color - color_zero_point for color in colors]

            if plot_type is 'plot3D':
                ax.scatter(xs, ys, colors)
            elif plot_type is 'plotVsX':
                cm = plt.cm.get_cmap('RdYlBu')
                if color_scheme in ['scaled']:
                    sc = plt.scatter(xs, colors, c=colors, vmin=color_lims[0], vmax=color_lims[1], cmap=cm)
                else:
                    sc = plt.scatter(xs, colors, c= single_obj_colors[i%len(single_obj_colors)])
                scats = scats + [sc]
                if i == len(best_objs[0:max_objs_to_plot]) - 1:
                    if color_scheme in ['scaled']:
                        cbar = plt.colorbar(sc)
                        cbar.set_label('g-r, \n shifted so median = 0')
                    else:
                        plt.legend(scats, best_objs )
                    plt.xlabel('CCD row (pix) for ' + self.pos_filter + '-filter')
                    plt.ylabel('g-r, \n shifted so median = 0')
                    plt.xlim(x_lims)
                    plt.ylim(color_lims)

            elif plot_type is 'plotVsY':
                cm = plt.cm.get_cmap('RdYlBu')
                if color_scheme in ['scaled']:
                    sc = plt.scatter(ys, colors, c=colors, vmin=color_lims[0], vmax=color_lims[1], cmap=cm)
                else:
                    sc = plt.scatter(ys, colors, c= single_obj_colors[i%len(single_obj_colors)] )
                scats = scats + [sc]
                if i == len(best_objs[0:max_objs_to_plot]) - 1:
                    if color_scheme in ['scaled']:
                        cbar = plt.colorbar(sc)
                        cbar.set_label('g-r, \n shifted so median = 0')
                    else:
                        print ('len(scats) = ' + str(len(scats)))
                        print ('good_objs = ' + str(best_objs) )
                        plt.legend(scats, best_objs )
                    plt.xlabel('CCD col (pix) for ' + self.pos_filter + '-filter')
                    plt.ylabel('g-r, \n shifted so median = 0')
                    plt.xlim(y_lims)
                    plt.ylim(color_lims)

            elif plot_type is 'plotVsZenith':
                cm = plt.cm.get_cmap('RdYlBu')
                if color_scheme in ['scaled']:
                    sc = plt.scatter(zeniths, colors, c = colors, vmin=color_lims[0], vmax=color_lims[1], cmap=cm)
                else:
                    sc = plt.scatter(zeniths, colors, c= single_obj_colors[i%len(single_obj_colors)] )
                scats = scats + [sc]
                if i == len(best_objs[0:max_objs_to_plot]) - 1:
                    if color_scheme in ['scaled']:
                        cbar = plt.colorbar(sc)
                        cbar.set_label('g-r, \n shifted so median = 0')
                    else:
                        print ('len(scats) = ' + str(len(scats)))
                        print ('good_objs = ' + str(best_objs) )
                        plt.legend(scats, best_objs )
                    plt.xlabel('Zenith angle')
                    plt.ylabel('g-r, \n shifted so median = 0')
                    #plt.xlim(y_lims)
                    plt.ylim(color_lims)

            else:
                #cmap = matplotlib.cm.get_cmap('viridis')
                #normalize = matplotlib.colors.Normalize(color_lims[0], color_lims[1])
                #colors = [cmap(normalize(value)) for value in colors]
                #ax[0,0].scatter(xs, ys, c = colors)

                #sc = plt.scatter(xs, ys, c = colors, cmap = 'seismic')
                #if good_obj == good_objs[-1]: plt.colorbar(sc)

                cm = plt.cm.get_cmap('RdYlBu')
                xy = range(20)
                z = xy
                sc = plt.scatter(xs, ys, c=colors, vmin=color_lims[0], vmax=color_lims[1], cmap=cm)
                if best_obj == best_objs[-1]:
                    cbar = plt.colorbar(sc)
                    cbar.set_label('g-r, \n shifted so median = 0')
                    plt.xlabel('CCD col (pix) for ' + self.pos_filter + '-filter')
                    plt.ylabel('CCD row (pix) for ' + self.pos_filter + '-filter')
                    plt.xlim(x_lims)
                    plt.ylim(y_lims)

                #im = ax.scatter(xs, ys, c=colors, cmap=plt.cm.jet)
                #plt.colorbar(scat)
                #cax, _ = matplotlib.colorbar.make_axes(ax)
                #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
                #colors = [color_lims[0] if color < color_lims[0] else color_lims[1] if color > color_lims[1] else color for color in colors ]

                #scat = plt.scatter(xs, ys, c = colors, vmin = cmap=plt.cm.coolwarm)
                #plt.colorbar(scat)

        #if not(plot3D):
        #    fig.colorbar(im, ax=ax)

        #plt.xlabel('CCD col (pix) for ' + self.pos_filter + '-filter')
        #plt.ylabel('CCD row (pix) for ' + self.pos_filter + '-filter')
        #plt.xlim(x_lims)
        #plt.ylim(y_lims)
        plt.show()



    def DetermineGoodSources(self, filters = ['g','r','i','z'], flags_to_cut_on = [1, 2, 4, 8, 16, 32, 64, 128], overall_names = 'all', mag_mcat_key_str = 'MAG_ISO'):
        master_objs = list(self.master_val_dict.keys())
        objs_usable = [1 for obj in master_objs]

        master_dict_good_obj_key_str = 'good_'

        if overall_names is 'all':
            overall_names = self.image_file_names[1:]
            suffix_length = 2
            overall_names = list(set([name[0:-suffix_length] for name in overall_names]))

        for filt in filters:
            self.measureMags(filt = filt, overall_names = overall_names)

        for filt in filters:
            self.PrepMasterValDict(overall_names, [master_dict_good_obj_key_str + filt ])
        self.measurePositionInRefFilter(ref_filter = self.pos_filter, overall_names = overall_names)
        self.measureFlags(overall_names = overall_names)

        for i in range(len(master_objs)):
            master_obj = master_objs[i]
            images = list(self.master_val_dict[master_obj].keys())
            flags = [[self.master_val_dict[master_obj][image]['flag'][filt] for filt in filters] for image in images]
            flagged_bad_by_filter = [[any(checkSexFlag(flags_to_cut_on, flag)) for flag in flag_filter_set] for flag_filter_set in flags]
            flagged_bad_images = [any(flagged_bad_set) for flagged_bad_set in flagged_bad_by_filter]
            mags = [[self.master_val_dict[master_obj][image][mag_mcat_key_str + '_' + filt] for filt in filters] for image in images]
            bad_detections_by_filter = [[np.isnan(mag) for mag in mag_filter_set] for mag_filter_set in mags]
            bad_detections = [any(bad_mag_set) for bad_mag_set in bad_detections_by_filter]
            overall_bad_by_filter = [ [ bad_detections_by_filter[j][k] or flagged_bad_by_filter[j][k] for k in range(len(filters)) ] for j in range(len(bad_detections_by_filter)) ]
            overall_bad_object = [ bad_detections[j] or flagged_bad_images[j] for j in range(len(bad_detections)) ]
            for j in range(len(images)):
                for k in range(len(filters)):
                    filt = filters[k]
                    self.master_val_dict[master_obj][images[j]][master_dict_good_obj_key_str + filt] = not(overall_bad_by_filter[j][k])


    def Plot2DColor(self, pos_var = 'x', filters = ['g','r'], targ_objs =  [155, 181, 173, 511, 438, 4281, 616],
                          color_cycle = ['blue','red','green','yellow','orange','cyan','black','purple','grey'],
                          zero_start = 1, x_lims = [0, 1500], color_lims = [-1.5, 1.5]):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        color_lims = color_lims
        for i in range(len(targ_objs)):
            targ_obj = targ_objs[i]
            color = color_cycle[i%len(color_cycle)]
            ordered_images = self.getObservationOrder()
            gmrs = []
            gmr_errs = []
            xs = []
            ys = []
            res_in_order =  [self.master_val_dict[targ_obj][image] for image in ordered_images]
            #print (self.master_val_dict[targ_obj])
            for j in range(len(res_in_order )):
                color_str = str(filters[0]) + '-' + str(filters[1])
                if not(np.isnan(res_in_order[j]['g-r'][0])) and not(np.isnan(res_in_order[j]['xy_g'][0])) and not(np.isnan(res_in_order[j]['xy_g'][1])):
                    gmrs = gmrs + [ res_in_order[j]['g-r'][0] ]
                    gmr_errs = gmr_errs + [res_in_order[j]['g-r'][1]]
                    xs = xs + [res_in_order[j]['xy_g'][0]]
                    ys = ys + [res_in_order[j]['xy_g'][1]]
                else:
                    print ('For [targ_obj, j] = ' + str([targ_obj, j]))
                    print ("[(res_in_order[j]['g-r'][0]), (res_in_order[j]['xy_g'][0]), (res_in_order[j]['xy_g'][1])] = " + str([(res_in_order[j]['g-r'][0]), (res_in_order[j]['xy_g'][0]), (res_in_order[j]['xy_g'][1])]))
            if pos_var in ['x','xs','X','XS','Xs']:
                positions = xs
                pos_lims = x_lims
            else:
                positions = ys
                pos_lims = y_lims

            if zero_start:
                gmr_zero = gmrs[0]
                gmrs = [gmr - gmr_zero for gmr in gmrs]

            print ('gmrs[0] = ' + str(gmrs[0]))

            ax.scatter(positions, gmrs, color = color  )
            ax.set_xlim(pos_lims)
            ax.set_ylim(color_lims)

        plt.show()

    def Plot3DColor(self, filters = ['g','r'], targ_objs =  [155, 181, 173, 511, 722, 438, 4281, 616],
                          color_cycle = ['blue','red','green','yellow','orange','cyan','black','purple','grey'],
                          zero_start = 1, x_lims = [0, 1500], y_lims = [0, 2140], color_lims = [-0.5, 0.5]):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        color_lims = color_lims
        for i in range(len(targ_objs)):
            targ_obj = targ_objs[i]
            color = color_cycle[i%len(color_cycle)]
            ordered_images = self.getObservationOrder()
            xs = []
            ys = []
            gmrs = []
            gmr_errs = []
            res_in_order =  [self.master_val_dict[targ_obj][image] for image in ordered_images]
            for j in range(len(res_in_order )):
                color_str = str(filters[0]) + '-' + str(filters[1])
                gmrs = gmrs + [res_in_order[j]['g-r'][0]]
                gmr_errs = gmr_errs + [res_in_order[j]['g-r'][1]]
                xs = xs + [res_in_order[j]['xy_g'][0]]
                ys = ys + [res_in_order[j]['xy_g'][1]]

            if zero_start:
                gmr_zero = gmrs[0]
                gmrs = [gmr - gmr_zero for gmr in gmrs]

            print ('gmrs[0] = ' + str(gmrs[0]))

            ax.scatter(xs, ys, gmrs, color = color  )
            ax.set_xlim(x_lims)
            ax.set_ylim(y_lims)
            ax.set_zlim(color_lims)

        plt.show()


    def getObservationOrder(self):
        str_root = 'crc_proc_4C16.49_'
        all_images = list(self.master_val_dict[1].keys())
        order_terms = [int(image[len(str_root):]) for image in all_images]
        #print ('all_images = ' + str(all_images))
        #print ('order_terms = ' + str(order_terms))
        ordered_terms, ordered_images = safeSortOneListByAnother(order_terms, [order_terms, all_images])

        return ordered_images


    def measurePositionInRefFilter(self, overall_names = 'all', ref_filter = 'g',
                                   obj_identifier_key_str = 'NUMBER',
                                    x_pos_key_str = 'X_IMAGE', y_pos_key_str = 'Y_IMAGE'):
        self.pos_ref_filter = ref_filter
        master_dict_pos_key_str = 'xy' + '_' + ref_filter

        filter_suffix = '_' + ref_filter
        if overall_names is 'all':
            overall_names = self.image_file_names[1:]
            suffix_length = len(filter_suffix)
            overall_names = list(set([name[0:-suffix_length] for name in overall_names]))

        if np.all([np.all([master_dict_pos_key_str in self.master_val_dict[obj][name].keys() for name in overall_names]) for obj in self.master_val_dict.keys()]):
            print ('keyword ' + master_dict_pos_key_str + ' already in master dictionary, so nothing to update here. ')
            return 1

        self.PrepMasterValDict(overall_names, [master_dict_pos_key_str])

        full_names = [overall_name + filter_suffix for overall_name in overall_names]
        self.UpdateObjDictsIfNeeded([x_pos_key_str, y_pos_key_str], full_names)

        file_indeces = [self.image_file_names.index(file_name) for file_name in full_names]

        for master_obj in self.master_objs:
            for i in range(len(full_names)):
                ref_file_index = file_indeces[i]
                overall_name = overall_names[i]
                full_name = full_names[i]
                local_obj_val = master_obj[ref_file_index]
                if local_obj_val >= 1:
                    x_val = self.image_dicts[self.image_file_names[ref_file_index]][master_obj[ref_file_index]][x_pos_key_str]
                    y_val = self.image_dicts[self.image_file_names[ref_file_index]][master_obj[ref_file_index]][y_pos_key_str]
                else:
                    x_val = np.nan
                    y_val = np.nan

                self.master_val_dict[master_obj[0]][overall_name][master_dict_pos_key_str] = [x_val, y_val]

        return 1

    def measureFlags(self, overall_names = 'all', filters = ['g','r','i','z'],
                     obj_identifier_key_str = 'NUMBER', flag_key_str = 'FLAGS'):

        master_dict_flag_key_str = 'flag'

        filter_suffixes = ['_' + filt for filt in filters ]
        if overall_names is 'all':
            overall_names = self.image_file_names[1:]
            suffix_length = len(filter_suffixes[0])
            overall_names = list(set([name[0:-suffix_length] for name in overall_names]))
        self.PrepMasterValDict(overall_names, [master_dict_flag_key_str])

        all_file_names_by_filter = [[overall_name + filter_suffix for filter_suffix in filter_suffixes] for overall_name in overall_names]
        all_file_names = [name for names_by_filter in all_file_names_by_filter for name in names_by_filter]

        self.UpdateObjDictsIfNeeded([flag_key_str], all_file_names)

        file_indeces = [[self.image_file_names.index(file_name) for file_name in name_in_filters] for name_in_filters in all_file_names_by_filter]

        for master_obj in self.master_objs:
            #print ('Determining flags for master obj number ' + str(master_obj[0]) )
            for i in range(len(overall_names)):
                color_set_indeces = file_indeces[i]
                #print ('color_set_indeces = ' + str(color_set_indeces))
                overall_name = overall_names[i]
                local_obj_vals = [master_obj[index] for index in color_set_indeces]
                flags_by_filter_dict = {}
                for filt_index in range(len(local_obj_vals)):
                    filt = filters[filt_index]
                    local_obj_val = local_obj_vals[filt_index]
                    if local_obj_val >= 1:
                        index = color_set_indeces[filt_index]
                        #print ('flag_key_str = ' + str(flag_key_str))
                        #print ('master_obj[index] = ' + str(master_obj[index]))
                        #print ('self.image_file_names[index]  = ' + str(self.image_file_names[index] ))
                        #print ('self.image_dicts[self.image_file_names[index]][master_obj[index]][flag_key_str] = ' + str(self.image_dicts[self.image_file_names[index]][master_obj[index]][flag_key_str]))
                        flag_val = self.image_dicts[self.image_file_names[index]][master_obj[index]][flag_key_str]
                    else:
                        flag_val = np.nan
                    flags_by_filter_dict[filt] = flag_val

                self.master_val_dict[master_obj[0]][overall_name][master_dict_flag_key_str] = flags_by_filter_dict

        return 0

    #def measureCatalogueParameter(self, param_key_str,
    #                              overall_names = 'all', filt = 'g',
    #                              obj_identifier_key_str = 'NUMBER', fwhm_key_str = 'FWHM_IMAGE',
    #                              ref_filter_index = 0):
    #
    #    filter_suffix = '_' + filt
    #    master_dict_param_key_str = param_key_str + filter_suffix
    #    if overall_names is 'all':
    #        overall_names = self.image_file_names[1:]
    #        suffix_length = len(filter_suffix)
    #        overall_names = list(set([name[0:-suffix_length] for name in overall_names]))
    #
    #    valid_key_str = self.PrepMasterValDict(overall_names, [master_dict_param_key_str])
    #    if not(valid_key_str):
    #        print ('Key ' + obj_identifier_key_str + ' not found in catalogue file. Not adding to master catalogue. ')
    #        return 1
    #
    #    all_file_names = [overall_name + filter_suffix for overall_name in overall_names]
    #
    #    self.UpdateObjDictsIfNeeded([fwhm_key_str, param_key_str], all_file_names)
    #
    #    file_indeces = [self.image_file_names.index(file_name) for file_name in all_file_names]
    #
    #    for master_obj in self.master_objs:
    #        #print ('Computing ' + filt + '-magnitude for master obj number ' + str(master_obj[0]) )
    #        for i in range(len(all_file_names)):
    #            color_index = file_indeces[i]
    #            overall_name = overall_names[i]
    #            local_obj_val = master_obj[color_index]
    #            if obj is valid:
    #                param_val = self.image_dicts[self.image_file_names[color_index]][master_obj[color_index]][obj_key_str]
    #                param_val = np.nan
    #
    #            self.master_val_dict[master_obj[0]][overall_name][master_dict_param_key_str] = param_val
    #    print ('master value dictionary is updated. ')
    #
    #    return 0

    def measureMags(self, overall_names = 'all', filt = 'g',
                             obj_identifier_key_str = 'NUMBER', fwhm_key_str = 'FWHM_IMAGE',
                             mag_key_str = 'MAG_ISO', #, mag_err_key_str = 'MAGERR_ISO',
                             ref_filter_index = 0):

        return self.measureCatalogueParameter(mag_key_str,
                                              overall_names = overall_names, filt = filt,
                                              obj_single_cat_key_str = obj_identifier_key_str, fwhm_key_str = fwhm_key_str,
                                              ref_filter_index = ref_filter_index)

    def measureEllipticities(self, overall_names = 'all', filt = 'g',
                             obj_identifier_key_str = 'NUMBER', fwhm_key_str = 'FWHM_IMAGE',
                             ell_key_str = 'ELLIPTICITY',
                             ref_filter_index = 0):

        return self.measureCatalogueParameter(ell_key_str,
                                              overall_names = overall_names, filt = filt,
                                              obj_identifier_key_str = obj_identifier_key_str, fwhm_key_str = fwhm_key_str,
                                              ref_filter_index = ref_filter_index)


    def measureCatalogueParameter(self, param_key_str,
                                  overall_names = 'all', filt = 'g',
                                  obj_single_cat_key_str = 'NUMBER', fwhm_key_str = 'FWHM_IMAGE',
                                  ref_filter_index = 0):

        #master_dict_mag_key_str = filt

        filter_suffix = '_' + filt
        master_dict_param_key_str = param_key_str + filter_suffix
        if overall_names is 'all':
            overall_names = self.image_file_names[1:]
            suffix_length = len(filter_suffix)
            overall_names = list(set([name[0:-suffix_length] for name in overall_names]))
        all_file_names = [overall_name + filter_suffix for overall_name in overall_names]
        dict_updated_successfully = self.UpdateObjDictsIfNeeded([fwhm_key_str, param_key_str], all_file_names)

        if not(dict_updated_successfully):
            print ('Key ' + param_key_str + ' not found in catalogue file. Not adding to master catalogue. ')
            return 1

        self.PrepMasterValDict(overall_names, [master_dict_param_key_str])

        file_indeces = [self.image_file_names.index(file_name) for file_name in all_file_names]

        print ('Updating field ' + str(param_key_str))
        for master_obj in self.master_objs:
            #print ('Computing ' + filt + '-magnitude for master obj number ' + str(master_obj[0]) )
            for i in range(len(all_file_names)):
                color_index = file_indeces[i]
                overall_name = overall_names[i]
                local_obj_val = master_obj[color_index] #The catalogue name of the index in the sextractor catalogue for this object.  This will be -1 if the object is not detected in that image
                if local_obj_val >= 1:
                    param_val = self.image_dicts[self.image_file_names[color_index]][master_obj[color_index]][param_key_str]
                else:
                    param_val = np.nan

                self.master_val_dict[master_obj[0]][overall_name][master_dict_param_key_str] = param_val
        print ('master value dictionary is updated. ')

        return 0

    def measureColor(self, overall_names = 'all', filters = ['g','r'],
                     obj_identifier_key_str = 'NUMBER', fwhm_key_str = 'FWHM_IMAGE',
                     mag_key_str = 'MAG_ISO', mag_err_key_str = 'MAGERR_ISO',
                     ref_filter_index = 0):

        master_dict_color_key_str = filters[0] + '-' + filters[1]

        filter_suffixes = ['_' + filt for filt in filters ]
        if overall_names is 'all':
            overall_names = self.image_file_names[1:]
            suffix_length = len(filter_suffixes[0])
            overall_names = list(set([name[0:-suffix_length] for name in overall_names]))
        self.PrepMasterValDict(overall_names, [master_dict_color_key_str])

        all_file_names_by_filter = [[overall_name + filter_suffix for filter_suffix in filter_suffixes] for overall_name in overall_names]
        all_file_names = [name for names_by_filter in all_file_names_by_filter for name in names_by_filter]

        self.UpdateObjDictsIfNeeded([fwhm_key_str, mag_key_str, mag_err_key_str], all_file_names)

        file_indeces = [[self.image_file_names.index(file_name) for file_name in name_in_filters] for name_in_filters in all_file_names_by_filter]

        for master_obj in self.master_objs:
            #print ('Computing color for master obj number ' + str(master_obj[0]) )
            for i in range(len(overall_names)):
                color_pair_indeces = file_indeces[i]
                overall_name = overall_names[i]
                local_obj_vals = [master_obj[index] for index in color_pair_indeces]
                if all([local_obj_val >=1 for local_obj_val in local_obj_vals]):
                    mag_vals = [self.image_dicts[self.image_file_names[index]][master_obj[index]][mag_key_str] for index in color_pair_indeces]
                    mag_errs = [self.image_dicts[self.image_file_names[index]][master_obj[index]][mag_err_key_str] for index in color_pair_indeces]
                    color_val = mag_vals[0] - mag_vals[1]
                    color_err = np.sqrt(sum([err ** 2.0 for err in mag_errs]))
                else:
                    color_val = np.nan
                    color_err = np.nan

                self.master_val_dict[master_obj[0]][overall_name][master_dict_color_key_str] = [color_val, color_err]

    def PrepMasterValDict(self, overall_names, used_key_strs):
        obj_nums = list(self.master_val_dict.keys())
        existing_names = self.images_in_master_dict
        existing_key_strs = self.keywords_in_master_dict

        for key_str in used_key_strs:
            if not(key_str in existing_key_strs):
                for obj_num in obj_nums:
                    for name in existing_names:
                        self.master_val_dict[obj_num][name][key_str] = 'UNASSIGNED'

        self.keywords_in_master_dict = list(set(self.keywords_in_master_dict).union(set(used_key_strs)))
        for overall_name in overall_names:
            if not(overall_name in existing_names):
                for obj_num in obj_nums:
                   self.master_val_dict[obj_num][overall_name] = {key_str:'UNASSIGNED' for key_str in self.keywords_in_master_dict}

        self.images_in_master_dict = list(set(self.images_in_master_dict).union(set(overall_names)))


    def generateRegionFiles(self, file_names = 'all', cat_file_addendum = '.cat',
                           obj_identifier_key_str = 'NUMBER', fwhm_key_str = 'FWHM_IMAGE',
                           x_pos_key_str = 'X_IMAGE', y_pos_key_str = 'Y_IMAGE',
                           reg_addendum = '_mcat_objs.reg'):
        if file_names is 'all':
            file_names = self.image_file_names[1:]

        self.UpdateObjDictsIfNeeded([x_pos_key_str, y_pos_key_str, fwhm_key_str], file_names)

        for file_name in file_names:
            reg_file_name = file_name + '_' + self.master_cat_file_name[0:-len('.mcat')] + reg_addendum
            reg_file_lines = ['# Region file format: DS9 version 4.1' ,
                              'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1', 'image']
            file_master_obj_index = self.image_file_names.index(file_name)

            orig_obj_dict = self.image_dicts[file_name]

            for master_obj in self.master_objs:
                local_obj_number = master_obj[file_master_obj_index]
                if local_obj_number >=1:
                    master_obj_number = master_obj[0]
                    orig_obj_vals_dict = orig_obj_dict[local_obj_number]
                    x = orig_obj_vals_dict[x_pos_key_str]
                    y = orig_obj_vals_dict[y_pos_key_str]
                    fwhm = orig_obj_vals_dict[fwhm_key_str]
                    orig_obj_vals = [x,y,fwhm]
                    text_addendum = '# text(' + str(x) + ',' + str(y + fwhm) + ') text={M' + str(master_obj_number) + ', L' + str(local_obj_number) + '}'
                    reg_file_lines = reg_file_lines + ['circle(' + ','.join([str(val) for val in orig_obj_vals]) + ') ', text_addendum]

            with open (self.target_dir + reg_file_name, 'w') as f:
                for line in reg_file_lines:
                    f.write(line + "\n")

        print ('Just wrote region files for master cat file ' + self.master_cat_file_name)

    def UpdateObjDictsIfNeeded(self, key_strs, file_names):
        all_key_strs_included = all([key in self.used_key_strs for key in key_strs])
        all_files_included = all([file_name in self.loaded_files for file_name in file_names])
        successfully_updated = 1
        if not (all_key_strs_included and all_files_included):
            updated_key_strs = list(set(self.used_key_strs).union(set(key_strs)))
            updated_file_names = list(set(self.loaded_files).union(set(file_names)))
            successfully_updated = self.CreateObjDicts(updated_key_strs, file_names = file_names, new_keys_to_add = not(all_key_strs_included))
        return successfully_updated


    def CreateObjDicts(self, key_strs, file_names = 'all', obj_identifier_key_str = 'NUMBER', cat_file_addendum = '.cat', new_keys_to_add = 1):
        if file_names is 'all':
            file_names = self.image_file_names[1:]

        #print ('file_names = ' + str(file_names))
        for file_name in file_names:
            if not(file_name in self.loaded_files) or new_keys_to_add:
                origPyObj = SO.PySex(self.target_dir + file_name + cat_file_addendum)

                local_obj_names = origPyObj.full_dict[obj_identifier_key_str]

                orig_obj_dict = {}
                for i in range(len(local_obj_names)):
                    if np.all(np.array([key_str in origPyObj.full_dict.keys() for key_str in key_strs])):
                        orig_obj_dict[int(local_obj_names[i])] = {key_str:origPyObj.full_dict[key_str][i] for key_str in key_strs}
                    else:
                        print ('Parameter not found in catalogue.')
                        return 0
                self.image_dicts[file_name] = orig_obj_dict

        self.used_key_strs = key_strs
        self.loaded_files = file_names

        return 1

    def addObjectIndependentQuantities(self, imgs_dir = None, overall_names = 'all', header_key_strs = [], prefix_to_remove = 'crc_', prefix_to_add = 'wcs_', file_suffix = '_g.fits'):
        if imgs_dir is None:
            target_dir_components = self.target_dir.split('/')[1:-1]
            imgs_dir = '/' + '/'.join(target_dir_components[0:-1]) + '/'

        if overall_names is 'all':
            overall_names = list(self.master_val_dict[self.master_objs[0][0]].keys())

        print ('overall_names = ' + str(overall_names ))
        for overall_name in overall_names:
            if not(overall_name) in self.object_ind_quantities.keys(): self.object_ind_quantities[overall_name] = {}
            if len(header_key_strs) > 0:
                if overall_name[0:len(prefix_to_remove)] in [prefix_to_remove]:
                    data, header = readInDataFromFitsFile(prefix_to_add + overall_name[len(prefix_to_remove):] + file_suffix, imgs_dir)
                else:
                    print ('overall_name = ' + str(overall_name))
                    data, header = readInDataFromFitsFile(prefix_to_add + overall_name + file_suffix, imgs_dir)
                for header_key_str in header_key_strs:
                    self.object_ind_quantities[overall_name][header_key_str] = header[header_key_str]
        return 1

    def addFrameIndependentQuantities(self, master_objs = 'all', quantity_key_strs = [], quantity_values_by_obj = []):

        if master_objs is 'all':
            master_objs = [master_obj[0] for master_obj in self.master_objs]

        #print ('master_objs = ' + str(master_objs))
        for master_obj in master_objs:
            if not(master_obj) in self.frame_ind_quantities.keys(): self.frame_ind_quantities[master_obj] = {}
            for i in range(len(quantity_key_strs)):
                quantity_key_str = quantity_key_strs[i]
                quantity_value_by_obj = quantity_values_by_obj[i][master_obj]
                self.frame_ind_quantities[master_obj][quantity_key_str] = quantity_value_by_obj

        return 1

    def __init__(self, target_dir, master_cat_file_name,
                 start_key_strs = ['NUMBER', 'FWHM_IMAGE', 'X_IMAGE', 'Y_IMAGE', 'MAG_ISO', 'FLAGS', 'MAGERR_ISO'],
                 files_to_load = 'all', pos_filter = 'r', filters_to_consider = ['g','r','i','z'], overall_names ='all',
                 flags_to_cut_on = [1, 2, 4, 8, 16, 32, 64, 128], imgs_dir = None,
                 extra_cat_prefix = 'crc_', missing_cat_prefix = 'wcs_', params_to_add_key_str = ['MAG_ISO', 'MAGERR_ISO', 'ELLIPTICITY', 'CLASS_STAR', 'FLUX_MAX'], verbose = 0):

        #print ('[extra_cat_prefix, missing_cat_prefix] = ' + str([extra_cat_prefix, missing_cat_prefix]))
        self.target_dir = target_dir
        self.master_cat_file_name = master_cat_file_name
        with open(target_dir + master_cat_file_name, 'r') as f:
            lines = [ [entry for entry in line.rstrip('\n').split(' ')] for line in f]
            self.image_file_names = lines[0]
            self.master_objs = [[int(elem) for elem in master_obj] for master_obj in lines[1:]]

        self.pos_filter = pos_filter

        self.start_key_strs = start_key_strs
        self.loaded_files = []
        self.image_dicts = {}
        self.CreateObjDicts(self.start_key_strs, file_names = files_to_load)

        self.keywords_in_master_dict = []
        self.images_in_master_dict = []
        self.master_val_dict = {master_obj[0]:{} for master_obj in self.master_objs}

        print ('Measuring magnitudes of master objects ... ')
        for param_key_str in params_to_add_key_str:
            for filt in filters_to_consider:
                self.measureCatalogueParameter(param_key_str,
                                               overall_names = overall_names, filt = filt)
        print ('Measuring colors of master objects... ')
        self.measureColor(filters = ['i','z'], overall_names = overall_names)
        self.measureColor(filters = ['r','i'], overall_names = overall_names)
        self.measureColor(filters = ['g','r'], overall_names = overall_names)
        print ('Measuring positions of master objects in filter: ' + pos_filter)
        for filt in filters_to_consider:
            self.measurePositionInRefFilter(ref_filter = filt, overall_names = overall_names)
        print ('Loading flags of all objects: ' + pos_filter)
        self.measureFlags(overall_names = overall_names)
        print ('Picking which objects are usable, based on flags and existence in each image...')
        self.DetermineGoodSources(overall_names = overall_names, filters = filters_to_consider, flags_to_cut_on = flags_to_cut_on)

        self.object_ind_quantities = {}
        self.addObjectIndependentQuantities(imgs_dir = imgs_dir, header_key_strs = ['ZD'], prefix_to_remove = extra_cat_prefix, prefix_to_add = missing_cat_prefix)

        self.frame_ind_quantities = {}
        self.obj_position_functions_of_color = {}
        self.addFrameIndependentQuantities(master_objs = 'all', quantity_key_strs = [], quantity_values_by_obj = [])








        print ('Done')
