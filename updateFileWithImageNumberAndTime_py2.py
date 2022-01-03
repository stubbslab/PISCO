import sys
from astropy.io import fits
import astropy.time
import dateutil.parser
import datetime
import os.path
import csv 

def getJulianTime(ut_time, ut_date):
    dt = dateutil.parser.parse(ut_date + ' ' + ut_time)
    time = astropy.time.Time(dt)
    julian_time = time.jd
    return julian_time

def readInFitsHeaderAndArray(load_dir, file_name):
    hdul = fits.open(load_dir + file_name)
    header = hdul[0].header
    data_array = hdul[0].data
    return [data_array, header]

def writeToCSVFile(time, image_number, file_dir, file_name):
    full_file_name = file_dir + file_name
    new_row = [str(time), str(image_number) + '\n']
    if os.path.isfile(full_file_name):
        existing_rows = []
        with open(full_file_name) as fin:
            existing_rows = [row.split(',') for row in fin]
        #print 'existing_rows = ' + str(existing_rows)
        fin.close()
        existing_image_numbers = [int(float(row[1])) for row in existing_rows]
        #print 'existing_image_numbers = ' + str(existing_image_numbers)
        #print 'image_number = ' + str(image_number)
        if int(image_number) in existing_image_numbers:
            print 'Old image number.  Replacing old entry and readding to end of file.'
            replace_row_index = existing_image_numbers.index(int(image_number))
            del existing_rows[replace_row_index]
            existing_rows = existing_rows + [new_row] 
        else:
            print 'New image number.  Adding to end of file. '
            existing_rows = existing_rows + [new_row]
        with open(full_file_name, 'w') as fout:
            for row in existing_rows: 
                fout.write(','.join(row))
        fout.close() 
        #print 'Wrote to file ' + str(full_file_name)
    else:
        with open(full_file_name, 'w') as fd:
            fd.write(','.join([str(time), str(image_number)]) + '\n' )
        print 'Created and wrote to file ' + str(full_file_name)
    return 1 
    

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
    time = getJulianTime(image_acq_ut_time, ut_date)
    print 'time = ' + str(time)

    broken_name = image_name.split('_')
    
    image_number = broken_name[-1].split('.')[0]
    print 'image_number = ' + str(image_number)

    writeToCSVFile(time, image_number, image_dir, 'seeing.csv')

    #print 'ut_date = ' + str(ut_date)
    #print 'shutter_open_time = ' + str(shutter_open_time)
    #print 'shutter_close_time = ' + str(shutter_close_time) 
    
