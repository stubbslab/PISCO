import cosmics_py3 as cosmics
import sys 

def CleanCosmics(image_dir, image_names, readnoise = 5.0, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0, maxiter = 2, new_image_prefix = 'crc_'):
    
    for image_name in image_names:
        #print ('Beginning cosmic ray cleaning for image ' + image_dir + image_name)
        image_array, image_header = cosmics.fromfits(image_dir + image_name)
        c = cosmics.cosmicsimage(image_array, readnoise = readnoise, sigclip = sigclip, sigfrac = sigfrac, objlim = objlim, verbose = False)
        c.run(maxiter = maxiter)
        image_header['CRCLEANED'] = 'Cosmic rays removed by cosmics.py'
        cosmics.tofits(image_dir + new_image_prefix + image_name, c.cleanarray, image_header)

    #print 'Done cleaning cosmic rays. '

    return 1

if __name__ == "__main__":
    image_dir, image_name = sys.argv[1:]
    CleanCosmics(image_dir, [image_name]) 
        
