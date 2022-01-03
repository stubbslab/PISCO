import sys
import readInFitsFile_py2 as readInFitsFile
from astropy.io import fits
import datetime 

def convertToHours(ut_time):
    time_array = [float(elem) for elem in ut_time.split(':')]
    time = time_array[0] + time_array[1] / 60.0 + time_array[2] / 3600.0
    return time
        

def readInFitsHeaderAndArray(load_dir, file_name):
    hdul = fits.open(load_dir + file_name)
    header = hdul[0].header
    data_array = hdul[0].data
    return [data_array, header] 

if __name__ == "__main__":
    image_dir, image_name = sys.argv[1:]
    image, header = readInFitsHeaderAndArray(image_dir, image_name)
    #print 'header = '
    #print header 
    shutter_open_keyword = 'UT-OPEN'
    shutter_close_keyword = 'UT-CLOSE'
    image_acq_time_keyword = 'TELUT'
    ut_date_obs_keyword = 'DATEOBS'
    ut_date = header[ut_date_obs_keyword]
    shutter_open_time = header[shutter_open_keyword]
    shutter_close_time = header[shutter_close_keyword]
    image_acq_ut_time = header[image_acq_time_keyword]
    print convertToHours(image_acq_ut_time) 

    #print 'ut_date = ' + str(ut_date)
    #print 'shutter_open_time = ' + str(shutter_open_time)
    #print 'shutter_close_time = ' + str(shutter_close_time) 
    
