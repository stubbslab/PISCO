from astropy.io import fits

def readInFitsHeaderAndArray(load_dir, file_name):
    hdul = fits.open(load_dir + file_name)
    header = hdul[0].header
    data_array = hdul[0].data
    return [data_array, header] 
