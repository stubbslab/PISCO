import subprocess as sp
wcslist=open('wcs.list').read().splitlines()
for image in wcslist:
    sp.call(['solve-field','--downsample', '4','--no-plots', '--scale-units', 'arcsecperpix', '--scale-low', '0.2', '--scale-high', '0.24',image])
    imagename = image.replace(".fits", "")
    sp.call(['mv', imagename+'.new',imagename+'_WCS.fits'])
    sp.call(['rm', imagename+'-indx.xyls', imagename+'.axy', imagename+'.corr', imagename+'.match', imagename+'.rdls', imagename+'.solved', imagename+'.wcs'])
