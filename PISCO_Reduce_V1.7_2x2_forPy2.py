#Import necesary packages

from astropy.io import fits
import numpy as np

#Load Functions

def loadimage(name): #loads a PISCO Image into a 8-image 3D numpy array
    imagelist = []
    hdulist_temp = fits.open(name)
    for i in list(range(1,9)):
        image_data=hdulist_temp[i].data
        imagelist.append(image_data.astype(float))
    hdulist_temp.close() 
    imagearray = np.array(imagelist)
    return imagearray

def loadimageFB(name): #loads a Flatfield or Bias image created by PISCOBias(), PISCOFlatSeparate(), or PISCOFlatCombined()
    imagelist = []
    hdulist_temp = fits.open(name)
    for i in [0,1,2,3]:
        image_data=hdulist_temp[i].data
        imagelist.append(image_data.astype(float))
    hdulist_temp.close() 
    imagearray = np.array(imagelist)
    return imagearray

def PISCOStitch(rawimages): #loads the result of loadimage(), crops and combined the 8 raw PISCO images into 4 combined images returned as a 4-image 3D numpy array
    images=[]
    for i in list(range(0,8)):
        images.append(np.delete(rawimages[i],np.s_[773:],1))
    images=np.asarray(images)
    images[1]=np.fliplr(images[1])
    images[2]=np.fliplr(images[2])
    images[4]=np.flipud(images[4])
    images[5]=np.fliplr(np.flipud(images[5]))
    images[6]=np.fliplr(np.flipud(images[6]))
    images[7]=np.flipud(images[7])
    comimages=[]
    comimages.append(np.hstack((images[0],images[1])))
    comimages.append(np.hstack((images[3],images[2])))
    comimages.append(np.hstack((images[4],images[5])))
    comimages.append(np.hstack((images[7],images[6])))
    return np.asarray(comimages)

def PISCOImport(name): #combines the functions of loadimage() and PISCOStitch()
    temp=loadimage(name)
    return PISCOStitch(temp)

def PISCOBias(biaslistfile,outputbiasfilename): #loads a group of bias images named in a text file with one image name per line.  Bias files are processed, median combined, and written out as a 4 extension FITS file.
    biasimages=np.loadtxt(biaslistfile, dtype='str')
    allbias=[]
    for bias in biasimages:
        allbias.append(loadimage(bias))
    allbias=np.asarray(allbias)
    biasimages=PISCOStitch(np.median(allbias,axis=0))
    WriteFITSFB(biasimages,outputbiasfilename)

def PISCOFlatSeparate(flatlistfileG,flatlistfileR,flatlistfileI,flatlistfileZ,biasfile,outputflatfilename): #loads groups of flatfield images named in separate text files with one image name per line for each filter.  Also loads a biasfile created by PISCOBias() for bias subtraction of the flatfield images.  A combined flatfield image is output as a 4 extension FITS file.
    bandflats=[]
    for band in [flatlistfileG,flatlistfileR,flatlistfileI,flatlistfileZ]:
        flatimages=np.loadtxt(band, dtype='str')
        allflat=[]
        for flat in flatimages:
            allflat.append(PISCOStitch(loadimage(flat)))
        allflat=np.asarray(allflat)
        allflat=allflat-loadimageFB(biasfile)
        bandflats.append(np.median(allflat,axis=0))
    bandflats=np.asarray(bandflats)
    selectflats=np.asarray([bandflats[0,0,:,:],bandflats[1,1,:,:],bandflats[2,2,:,:],bandflats[3,3,:,:]])
    normarray=[]
    normarray.append(np.median(selectflats[0,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(selectflats[1,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(selectflats[2,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(selectflats[3,550:2500,20:1500],axis=(0,1)))
    normarray=np.asarray(normarray)
    WriteFITSFB(selectflats/normarray[:,np.newaxis,np.newaxis],outputflatfilename)

def PISCOFlatCombined(flatlistfile,biasfile,outputflatfilename): #Operates the same as PISCOFlatSeparate() but without the different lists for each filter.
    flatimages=np.loadtxt(flatlistfile, dtype='str')
    allflat=[]
    for flat in flatimages:
        allflat.append(PISCOStitch(loadimage(flat)))
    allflat=np.asarray(allflat)
    allflat=allflat-loadimageFB(biasfile)
    comboflat=np.median(allflat,axis=0)
    allflat=[]
    normarray=[]
    normarray.append(np.median(comboflat[0,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(comboflat[1,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(comboflat[2,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(comboflat[3,550:2500,20:1500],axis=(0,1)))
    normarray=np.asarray(normarray)
    WriteFITSFB(comboflat/normarray[:,np.newaxis,np.newaxis],outputflatfilename)
    
def WriteFITSFB(images,outfile): #Outputs a 4 extension FITS file suitable for processed Bias or Flatfield files as part of PISCOBias(), PISCOFlatSeparate(), or PISCOFlatCombined()
    new_hdul = fits.HDUList()
    one=fits.ImageHDU(images[0])
    two=fits.ImageHDU(images[1])
    three=fits.ImageHDU(images[2])
    four=fits.ImageHDU(images[3])
    new_hdul.append(one)
    new_hdul.append(two)
    new_hdul.append(three)
    new_hdul.append(four)
    new_hdul.writeto(outfile,clobber=True)
    
def PISCOProcBatch(scilist,biasfile,flatfile): #Processes scientific images loaded from a text file with one image name per line.  Also loads bias and flatfield files for science image processing.  Outputs 1-extension FITS files for each of the 4 color bands.  Images are saved in a 32-bit integer format for campatibility with Astrometrica.
    bias=loadimageFB(biasfile)
    flat=loadimageFB(flatfile)
    scilistimages=np.loadtxt(scilist, dtype='str', ndmin = 1)
    for sci in scilistimages:
        data=loadimage(sci)
        images=PISCOStitch(data)
        #If the images are not all the same size (most likely because of binning problems),
        #    then just assume no 0 bias and uniform flat.
        all_images_same_size = (np.all([(np.shape(image) == np.shape(bias)) and (np.shape(image)== np.shape(flat))
                                 for image in images])
                                and np.shape(bias) == np.shape(flat)) 
        if not(all_images_same_size):
            images = (images - 0.0) / 1.0
        else:
            images=(images-bias)/flat
        images[np.isnan(images)] = 0
        count=0
        lowy=[365,429,300,375]
        highy=[2505,2579,2440,2515]
        for band in ['g','r','i','z']:
            cropped=images[count,lowy[count]:highy[count],10:1510]
            file=fits.PrimaryHDU()
            file.header=fits.getheader(sci)
            file.data=cropped
            #file.header.rename_keyword('DATEOBS','DATE-OBS')
            #file.header.rename_keyword('TELUT','UT-START')
            #file.header.rename_keyword('UT-OPEN','UT-START')
            #file.header.rename_keyword('UT-CLOSE','UT-END')
            file.scale('int32')
            file.writeto('proc_'+sci[:-5]+'_'+band+'.fits',clobber=True)
            count=count+1
    return 'Done'

#Script starts here

choice1=input('Create combined bias? (y/n): ')
if choice1=='y':
    biaslist=input('Enter bias list filename: ')
    PISCOBias(biaslist,'BIAS.fits')
choice2=input('Create combined flat? (y/n): ')
if choice2=='y':
    flatglist=input('Enter g-band flat list filename: ')
    flatrlist=input('Enter r-band flat list filename: ')
    flatilist=input('Enter i-band flat list filename: ')
    flatzlist=input('Enter z-band flat list filename: ')
    PISCOFlatSeparate(flatglist,flatrlist,flatilist,flatzlist,'BIAS.fits','FLAT.fits')
scilist=input('Enter science list filename: ')
PISCOProcBatch(scilist,'BIAS.fits','FLAT.fits')

