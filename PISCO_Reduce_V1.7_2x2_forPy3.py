#Import necesary packages

from astropy.io import fits
import numpy as np

#Load Functions

def loadimage(name): #loads a PISCO Image into a 8-image 3D numpy array
    imagelist = []
    hdulist_temp = fits.open(name)
    for i in list(range(1,9)):
        image_data=hdulist_temp[i].data
        #print ('np.shape(image_data) = ' + str(np.shape(image_data)) ) 
        imagelist.append(image_data.astype(float))
    hdulist_temp.close() 
    imagearray = np.array(imagelist)
    #print ('np.shape(imagearray) = ' + str(np.shape(imagearray))) 
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
    bin_2x2_cut = 773
    #print (np.shape(rawimages[0]) )
    #print ('np.shape(rawimages[0])[0] = ' + str(np.shape(rawimages[0])[0])) 
    if np.shape(rawimages[0])[0] < 4000:
        bin_cut = bin_2x2_cut
    else:
        bin_cut = bin_2x2_cut * 2 
    for i in list(range(0,8)):
        images.append(np.delete(rawimages[i],np.s_[bin_cut:],1))
    images=np.asarray(images)
    #print ('(np.shape(images) ) = ' + str(np.shape(images) ))
    images[1]=np.flip(images[1],1)
    images[2]=np.flip(images[2],1)
    images[4]=np.flip(images[4],0)
    images[5]=np.flip(np.flip(images[5],0),1)
    images[6]=np.flip(np.flip(images[6],0),1)
    images[7]=np.flip(images[7],0)
    comimages=[]
    comimages.append(np.hstack((images[0],images[1])))
    comimages.append(np.hstack((images[3],images[2])))
    comimages.append(np.hstack((images[4],images[5])))
    comimages.append(np.hstack((images[7],images[6])))
    return np.asarray(comimages)

def PISCOImport(name): #combines the functions of loadimage() and PISCOStitch()
    temp=loadimage(name)
    return PISCOStitch(temp)

def PISCOBias(biaslistfile,outputbiasfilename, target_dir): #loads a group of bias images named in a text file with one image name per line.  Bias files are processed, median combined, and written out as a 4 extension FITS file.
    biasimages=np.loadtxt(target_dir + biaslistfile, dtype='str')
    allbias=[]
    for bias in biasimages:
        allbias.append(loadimage(target_dir + bias))
    allbias=np.asarray(allbias)
    biasimages=PISCOStitch(np.median(allbias,axis=0))
    WriteFITSFB(biasimages, target_dir + outputbiasfilename)

def PISCOFlatSeparate(flatlistfileG,flatlistfileR,flatlistfileI,flatlistfileZ,biasfile,outputflatfilename, target_dir): #loads groups of flatfield images named in separate text files with one image name per line for each filter.  Also loads a biasfile created by PISCOBias() for bias subtraction of the flatfield images.  A combined flatfield image is output as a 4 extension FITS file.
    bandflats=[]
    for band in [flatlistfileG,flatlistfileR,flatlistfileI,flatlistfileZ]:
        flatimages=np.loadtxt(target_dir + band, dtype='str')
        allflat=[]
        for flat in flatimages:
            allflat.append(PISCOStitch(loadimage(target_dir + flat)))
        allflat=np.asarray(allflat)
        allflat=allflat-loadimageFB(target_dir + biasfile)
        bandflats.append(np.median(allflat,axis=0))
    bandflats=np.asarray(bandflats)
    selectflats=np.asarray([bandflats[0,0,:,:],bandflats[1,1,:,:],bandflats[2,2,:,:],bandflats[3,3,:,:]])
    normarray=[]
    normarray.append(np.median(selectflats[0,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(selectflats[1,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(selectflats[2,550:2500,20:1500],axis=(0,1)))
    normarray.append(np.median(selectflats[3,550:2500,20:1500],axis=(0,1)))
    normarray=np.asarray(normarray)
    WriteFITSFB(selectflats/normarray[:,np.newaxis,np.newaxis],target_dir + outputflatfilename)

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
    new_hdul.writeto(outfile,overwrite=True)
    
def PISCOProcBatch(scilist,biasfile,flatfile, target_dir): #Processes scientific images loaded from a text file with one image name per line.  Also loads bias and flatfield files for science image processing.  Outputs 1-extension FITS files for each of the 4 color bands.  Images are saved in a 32-bit integer format for campatibility with Astrometrica.
    bias=loadimageFB(target_dir + biasfile)
    flat=loadimageFB(target_dir + flatfile)
    scilistimages=np.loadtxt(target_dir + scilist, dtype='str', ndmin = 1)
    for sci in scilistimages:
        data=loadimage(target_dir + sci)
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
        #print ('np.shape(images[0])  = ' + str(np.shape(images[0])))
        lowx_bin_2x2 = 10
        highx_bin_2x2 = 1510
        lowy_bin_2x2 = [365,429,300,375]
        highy_bin_2x2 = [2505,2569,2440,2515]
        if np.shape(images[0])[0] < 4000:
            lowy=[elem for elem in lowy_bin_2x2]
            highy=[elem for elem in highy_bin_2x2]
            lowx = lowx_bin_2x2
            highx = highx_bin_2x2
        else:
            lowy=[elem * 2 for elem in lowy_bin_2x2]
            highy=[elem * 2 for elem in highy_bin_2x2]
            lowx = lowx_bin_2x2 * 2
            highx = highx_bin_2x2 * 2
        #print ('lowx = ' + str(lowx))
        #print ('highx = ' + str(highx))
        for band in ['g','r','i','z']:
            #print ('np.shape(images)  = ' + str(np.shape(images) ))
            cropped=images[count,lowy[count]:highy[count],lowx:highx]
            #print ('np.shape(cropped)  = ' + str(np.shape(cropped) ))
            file=fits.PrimaryHDU()
            file.header=fits.getheader(target_dir + sci)
            file.data=cropped
            #file.header.rename_keyword('DATEOBS','DATE-OBS')
            #file.header.rename_keyword('TELUT','UT-START')
            #file.header.rename_keyword('UT-OPEN','UT-START')
            #file.header.rename_keyword('UT-CLOSE','UT-END')
            file.scale('int32')
            file.writeto(target_dir + 'proc_'+sci[:-5]+'_'+band+'.fits',overwrite=True)
            count=count+1
    return 'Done'

#Script starts here

target_dir=input("Enter directory, including trailing / , where images are stored (just hit return for the current directory): ")
print ('target_dir = ' + str(target_dir)) 

choice1=input('Create combined bias? (y/n): ')
if choice1=='y':
    biaslist=input('Enter bias list filename: ')
    PISCOBias(biaslist,'BIAS.fits', target_dir)
choice2=input('Create combined flat? (y/n): ')
if choice2=='y':
    flatglist=input('Enter g-band flat list filename: ')
    flatrlist=input('Enter r-band flat list filename: ')
    flatilist=input('Enter i-band flat list filename: ')
    flatzlist=input('Enter z-band flat list filename: ')
    PISCOFlatSeparate(flatglist,flatrlist,flatilist,flatzlist,'BIAS.fits','FLAT.fits', target_dir)
scilist=input('Enter science list filename: ')
PISCOProcBatch(scilist,'BIAS.fits','FLAT.fits', target_dir)

