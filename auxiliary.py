#!/usr/bin/python
# -*- coding: utf-8 -*-

# MAKE SEGMENTATION IMAGE 

smoothtype='gauss'

# from auxiliary import region,create_circular_mask,fit_gauss,replace_secondary_sources_with_bkg,calc_bg,rmsdiff,cartesian_product,makeGaussian,petrosian_radius,binary_detection_pawlik,background_pawlik,RMAX_calc,smoothed_mask_skysub,asymmetry_function,asymmetry_function2,asymmetry_simple,Cutout2D_mine,minimize_mu,asymmetry_verysimple_4bkg,M20_simple,asymmetry_function2_shape,asymmetry_function_bkg,bkg_search,asymmetry_bkg_simple
import os,sys,time
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy,math,glob,random

from gauss_fit_routines import gaussian,moments,fitgaussian

from skimage import data
from skimage.morphology import disk
from skimage.filters import rank
from skimage.morphology import erosion, dilation, opening, closing, white_tophat
from skimage.morphology import black_tophat, skeletonize, convex_hull_image
from astLib import astCalc,astCoords,astImages,astPlots,astStats,astWCS

import numpy as np
from numpy.linalg import norm
from astropy.nddata.utils import Cutout2D
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy.utils.data import download_file
from astropy.io.fits import getdata
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve

from photutils.aperture import CircularAperture,RectangularAperture
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from photutils.aperture import CircularAperture, CircularAnnulus,aperture_photometry
from photutils.aperture import ApertureStats

import skimage

#import fits
from scipy import signal,ndimage
#import matplotlib.animation as animation
import photutils
print('Photutils version =',photutils.__version__)
from photutils import detect_threshold, detect_sources, deblend_sources 
# from photutils import source_properties, properties_table  DEPRECATED
from photutils import EllipticalAperture
from photutils.morphology import data_properties
from photutils.datasets import make_4gaussians_image
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize 
from photutils.utils import make_random_cmap
from photutils.segmentation import SourceCatalog
from astropy.convolution import Box2DKernel,Gaussian2DKernel,convolve

from smoothing import smooth
#from segm_utils import segmentation
import scipy.ndimage as ndi
from astropy.visualization import LogStretch
from astropy.modeling import models
import photutils,statmorph,skimage
from rotate import rotate
from astropy.cosmology import WMAP9 as cosmo
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,FixedLocator
from scipy import signal,ndimage,stats



def import_table(homefolder,filename):
    import numpy as np
    table=np.genfromtxt(homefolder+filename,names=True,dtype=None)
    return table

def region(ra,dec,s,color,font,t2):
    reg= 'fk5;circle('+str(ra)+','+str(dec)+','+str(s)+'") #color='+color+' font='+str(font)+' text={'+str(t2)+'}'
    r = pyregion.parse(reg)

def create_circular_mask(h, w, center=None, radius=None):
    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y- center[1])**2)
    mask = dist_from_center <= radius
    return mask

def fit_gauss(lnspc,mean,std,amplitude) :
    #m, s = stats.norm.fit(ser) # get mean and standard deviation  
    pdf_g = amplitude*stats.norm.pdf(lnspc, mean, std)
    return pdf_g

def replace_secondary_sources_with_bkg(imagein33,maskgoodsource,segmentation_all33) :
  from photutils.datasets import make_noise_image
  import numpy as np
  from scipy.optimize import curve_fit
  check_fit_noise=False 
  show_initial_images=False

  # Create mask of secondary sources
  inverse_maskgoodsource=np.absolute(maskgoodsource-1)
  segmentation_all33[segmentation_all33>0.5]=1
  mask_secondary_sources=segmentation_all33*inverse_maskgoodsource
  inverse_segmentationall=np.absolute(segmentation_all33-1)

  flux_bkg_1Darray=imagein33[segmentation_all33==0]
  mean_guess3 = np.mean(flux_bkg_1Darray)
  std_guess3 = np.std(flux_bkg_1Darray)
  area_guess3 = np.sum(flux_bkg_1Darray)
  entries3, bin_edges3 = np.histogram(flux_bkg_1Darray, range=[min(flux_bkg_1Darray), max(flux_bkg_1Darray)], bins=100,density=True)  # range=[limitsbkg[0], limitsbkg[1]]
  # calculate binmiddles
  bin_middles3 = 0.5*(bin_edges3[1:] + bin_edges3[:-1])
  width3 = 0.8*(bin_edges3[1]-bin_edges3[0])
  # FIT NORMAL :
  norm_opt33, _ = curve_fit(fit_gauss, bin_middles3, entries3, p0=[mean_guess3,std_guess3,area_guess3])
  mean_opt33=norm_opt33[0] ; std_opt33=norm_opt33[1]
  print('Fitted mean and std (gaussian noise bkg) =',mean_opt33,std_opt33)
  synthetic_noise=make_noise_image(imagein33.shape, distribution='gaussian',mean=mean_opt33, stddev=std_opt33)


  good_source=imagein33*maskgoodsource
  secondary_sources_replaced=mask_secondary_sources*synthetic_noise
  background33=imagein33*inverse_segmentationall

  final_image33=good_source+secondary_sources_replaced+background33  

  if show_initial_images==True :
    f, ((ax1, ax2, ax3), (ax4, ax5,ax6)) = plt.subplots(2,3, sharey=True, sharex=True)
    imagestat1=calc_bg(imagein33)
    back1=imagestat1[0]   #'background'
    sigmaback1=imagestat1[1]
    norm1=matplotlib.colors.Normalize(vmin=back1-2*sigmaback1, vmax=back1+7*sigmaback1, clip=False)
    imagestat3=calc_bg(imagein33)
    #imagestat3=calc_bg(np.ravel(image_4_clumpiness))
    back3=imagestat3[0]   #'background'
    sigmaback3=imagestat3[1]
    norm3=matplotlib.colors.Normalize(vmin=back3-2*sigmaback3, vmax=back3+7*sigmaback3, clip=False)
    ax1.imshow(imagein33, origin='lower', norm=norm1, cmap='hot')
    ax1.set_title('Image original')
    ax2.imshow(maskgoodsource, origin='lower', cmap='hot', norm=norm1)
    ax2.set_title('Mask good source')
    # Immagine originale con maschera
    #ax3.imshow(maskgood*image_4_clumpiness, origin='lower', norm=norm1, cmap='hot')
    #ax3.set_title(str(IDgal))
    #plt.title(str(IDgal))
    #plt.colorbar()
    ax3.imshow(mask_secondary_sources, origin='lower', norm=norm3, cmap='hot')
    ax3.set_title('Mask secondary sources')
    ax4.imshow(synthetic_noise, origin='lower', norm=norm3, cmap='hot')
    ax4.set_title('Synthetic noise')
    ax5.imshow(imagein33*inverse_segmentationall, origin='lower', norm=norm3, cmap='hot')
    ax5.set_title('Only background\nfrom original image')
    ax6.imshow(final_image33, origin='lower', norm=norm3, cmap='hot')
    ax6.set_title('Final image after replacement')
    # Segmentation map (Deblended) on clumps
    plt.subplots_adjust(wspace=0.2, hspace=0.01)
    #plt.tight_layout()
    #plt.show()

  if check_fit_noise==True :
    plt.bar(bin_middles, entries, align='center', width=width, label = 'Normalised data Gauss', alpha=0.5,color='salmon')
    plt.plot(bin_middles, fit_gauss(bin_middles, *norm_opt33), color='red',ls='dashed', label='Normal fit')
    
    plt.axvline(norm_opt33[0],lw=2,color='k',ls='dashed')
    #plt.plot(bin_middles, fit_poisson(bin_middles, *pois_opt), 'r--', label='Poisson fit')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()
  return final_image33


def calc_bg(a):
    """ Estimate the background level and rms. """
    import numpy as np
    good = ~np.isnan(a)
    assert good.sum(), 'no good pixels!'
    # poor man's source detection...
    vmax = np.percentile(a[good], 70)
    c0 = a[good] < vmax
    temp = a[good][c0]
    bg = np.median(temp)
    # now find the rms in the background
    belowbg = temp[temp < bg]
    #print(len(belowbg))
    # remove lowest 2% to get rid of any outliers
    flo = np.percentile(belowbg, 5)
    belowbg = belowbg[belowbg > flo]
    rms = np.concatenate([belowbg, 2*bg - belowbg]).std()
    return bg, rms    # background and rms

# Altro modo per calcolare rms di un'immagine (non ho provato)


def rmsdiff(im1, im2):
    "Calculate the root-mean-square difference between two images"
    import ImageChops
    import math, operator
    h = ImageChops.difference(im1, im2).histogram()

    # calculate rms
    return math.sqrt(reduce(operator.add,map(lambda h, i: h*(i**2), h, range(256))) / (float(im1.size[0]) * im1.size[1]))

def cartesian_product(tup1, tup2):
  """Returns a tuple that is the Cartesian product of tup_1 and tup_2
  >>> X = (1, 2)
  >>> Y = (4, 5)
  >>> cartesian_product(X, Y)
  ((1, 4), (1, 5), (2, 4), (2, 5))
  """
  return tuple((t1, t2) for t1 in tup1 for t2 in tup2)

def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    if center is None:  x0 = y0 = size / 2
    else: x0 = center[0] ; y0 = center[1]
    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


# def create_circular_mask(h, w, center=None, radius=None) :
#     if center is None: # use the middle of the image
#         center = [int(w/2), int(h/2)]
#     if radius is None: # use the smallest distance between the center and image walls
#         radius = min(center[0], center[1], w-center[0], h-center[1])
#     Y, X = np.ogrid[:h, :w]
#     dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
#     mask = dist_from_center <= radius
#     return mask


def petrosian_radius(image_newX,centroid_singleX,radii_test) :
  from astropy.stats import sigma_clip

  #print('shape imagein =',image_newX.shape)
  mu_array_1=[]
  mu_array_2=[]
  for ik in np.arange(len(radii_test)) : 

    # METODO 1 : 
    aperture_DD = CircularAperture(centroid_singleX, r=radii_test[ik])
    annulus_aperture_DD = CircularAnnulus(centroid_singleX, r_in=radii_test[ik], r_out=radii_test[ik]+1)
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    aper_stats = ApertureStats(image_newX, aperture_DD, sigma_clip=None)
    #print('aper stats =',aper_stats.median)
    bkg_stats = ApertureStats(image_newX, annulus_aperture_DD, sigma_clip=sigclip)
    #print('bkg stats =',bkg_stats.median)
    #print('Mu = ',bkg_stats.median/aper_stats.median)
    mu_array_1.append(bkg_stats.median/aper_stats.median)
    # #phot_table = aperture_photometry(image_newX, aperture_DD)
    # aperstats1 = ApertureStats(image_newX, aperture_DD)
    # mean_internal = aperstats1.mean
    # print('mean inside =',mean_internal)
    

    # METODO 2 :

    circ_mask1= create_circular_mask(image_newX.shape[0], image_newX.shape[1], center=centroid_singleX, radius=radii_test[ik])
    circ_mask2= create_circular_mask(image_newX.shape[0], image_newX.shape[1], center=centroid_singleX, radius=radii_test[ik]+1.)
    circ_mask1=circ_mask1*1
    circ_mask2=circ_mask2*1

    #print('circ_mask1 =',circ_mask1*1)
    centro=image_newX*circ_mask1
    centro_1d=centro[circ_mask1==1]
    centro_avg=np.median(centro_1d)
    #print('centro-avg =',centro_avg)
    #plt.imshow(centro,origin='lower')
    #plt.show()
    #quit()
    
    annulus_mask=np.zeros(centro.shape)
    annulus_mask[ (circ_mask1==0) & (circ_mask2==1)]=1
    imageLL=annulus_mask*image_newX
    imageLLb=imageLL[annulus_mask==1]
    annulus_image=np.ravel(imageLLb)
    #annulus_image=np.ravel(annulus_image)
    #filtered_data = sigma_clip(annulus_image, sigma=3, maxiters=10)
    #print(annulus_image)
    annulus_avg=np.median(annulus_image)
    #print('annulus-avg =',annulus_avg)
    #print('Ratio =',annulus_avg/centro_avg)
    #plt.imshow(imageLL,origin='lower')
    #plt.show()

    mu_array_2.append(annulus_avg/centro_avg)
    # masks_DD = annulus_aperture_DD.to_mask(method='center')  # methods=exact or center
    # masks_DD.to_image(image_newX.shape)
    # plt.imshow(masks_DD)
    # plt.show()
    # quit()
    # data_weighted_DD = masks_DD.multiply(image_newX)
    # avg_DD=np.average(np.ravel(data_weighted_DD))
    # 
    # masks_JJ = aperture_DD.to_mask(method='center')  # methods=exact or center
    # masks_JJ.to_image(image_newX.shape)
    # plt.imshow(masks_JJ)
    # plt.show()
    # data_weighted_JJ = masks_JJ.multiply(image_newX)
    # avg_JJ=np.average(np.ravel(data_weighted_JJ))
  return mu_array_1,mu_array_2


def binary_detection_pawlik(imagein_smoothed_skysub,std_pawlik):
  def GetKey11(item):
    return item[1]
  # Binary detection mask creation :
  binary_pawlik_single=np.zeros(imagein_smoothed_skysub.shape)
  binary_pawlik_single[imagein_smoothed_skysub>std_pawlik]=1
  #binary_pawlik=np.zeros(imagein_smoothed_skysub.shape)
  imagein_smoothed_skysub_seg = photutils.detect_sources(imagein_smoothed_skysub, np.ones(imagein_smoothed.shape)*std_pawlik, 10,connectivity=8)
  # Take closest to the center
  #photutils.detect_sources(imagein, threshold1, npixelsV)
  propsMASK = source_properties(imagein_smoothed_skysub, imagein_smoothed_skysub_seg)
  distM=np.empty(len(imagein_smoothed_skysub_seg.labels))
  centr0=np.empty(len(imagein_smoothed_skysub_seg.labels))
  centr1=np.empty(len(imagein_smoothed_skysub_seg.labels))
  for l2 in imagein_smoothed_skysub_seg.labels :
      propMASK=propsMASK[l2-1]
      positionM = (propMASK.xcentroid.value, propMASK.ycentroid.value)
      positionP = (propMASK.maxval_xpos.value, propMASK.maxval_ypos.value)
      centr0[l2-1]= positionM[0]
      centr1[l2-1]= positionM[1]
      distM[l2-1] = math.hypot(positionM[0] - imagein_smoothed_skysub.shape[0]/2, positionM[1] - imagein_smoothed_skysub.shape[1]/2)
  mlabelsM= np.column_stack((imagein_smoothed_skysub_seg.labels,distM,centr0,centr1))
  #print 'len(imagein_smoothed_skysub_seg.labels)=',len(imagein_smoothed_skysub_seg.labels)
  
  # Qua sto prendendo la regione piu' vicina al centro dell'immagine !!!
  mgood2=sorted(mlabelsM,key=GetKey11)
  labelsingle=int(mgood2[0][0])
  centroid=[mgood2[0][2],mgood2[0][3]]
  #print('centroid =',centroid)
  #Rmax=mgood2[-1][1]
  #print('label and distance (Rmax of Pawlik+16) of chosen label =',labelsingle) #,Rmax
  imagein_smoothed_skysub_seg.keep_labels(labelsingle)
  binary_pawlik=imagein_smoothed_skysub_seg.data
  #binary_pawlik[binary_pawlik>0]=1
  binary_pawlik=1*binary_pawlik 
  return binary_pawlik,centroid


def background_pawlik(flux_histogram_bkg) :
  merit=[0,0]
  for i in np.arange(1000):
    #print len(flux_histogram_bkg)
    flux_histogram_bkg,low,upp=scipy.stats.sigmaclip(flux_histogram_bkg, low=3, high=3)
    bkg=2.5*np.median(flux_histogram_bkg)-1.5*np.average(flux_histogram_bkg)
    #print 'bkg Pawlik =',bkg
    #print 'std(bkg)',np.std(flux_histogram_bkg)
    #print 'difference =',np.abs(bkg-merit[-1])
    if ( np.abs(bkg-merit[-1]) < 1e-8) :
      #print 'Difference =',bkg-merit[i]
      #print 'Number of iterations =',i
      return bkg,np.std(flux_histogram_bkg)
    else :
      merit=np.append(merit,bkg)
    if (i==999):
      print('Reached maximum number of iterations. Returning last value of bkg !!')
      return bkg,np.std(flux_histogram_bkg)

def RMAX_calc(binary_pawlik,centroid):
  # Calcolo Rmax:
  distanze=[]
  whereC=np.where(binary_pawlik>0)
  numpix_pawlik=len(whereC[0])
  #print 'How many pixels has the pawlik galaxy =', numpix_pawlik
  pixel_search_pawlik=[  (whereC[0][ig],whereC[1][ig]) for ig in np.arange(len(whereC[0])) ]
  for i in np.arange(len(pixel_search_pawlik)):
    distX=math.hypot(pixel_search_pawlik[i][0] - centroid[0], pixel_search_pawlik[i][1] - centroid[1])
    distanze=np.append(distanze,distX)
  distanze=sorted(distanze,reverse=True)
  Rmax=distanze[0]
  # print('Rmax =',Rmax)
  radius=int(Rmax)
  return radius,Rmax

def smoothed_mask_skysub(imagein8,segmap8) :
  imagein_smoothed8= ndi.uniform_filter(imagein8, size=3)
  segmap_smoothed8= ndi.uniform_filter(segmap8, size=3)
  flux_histogram_bkg8=np.asarray(imagein_smoothed8[(segmap_smoothed8==0)])
  
  bkg_pawlik8,std_pawlik8=background_pawlik(flux_histogram_bkg8)
  # print('background for Pawlik binary mask (and std) 8 =',bkg_pawlik8,std_pawlik8)
  
  #imagein_skysub=imagein-bkg_pawlik
  imagein_smoothed_skysub8=imagein_smoothed8-np.ones(imagein_smoothed8.shape)*bkg_pawlik8
  return imagein_smoothed_skysub8,bkg_pawlik8,std_pawlik8
  

# def asymmetry_function
def asymmetry_function(xcen,ycen,image_smaller,maskgood_smaller,sizeK,maskgood,imagein_smoothed_skysub) :
  center_search_x=np.arange(xcen-1.2,xcen+1.2,0.1)
  center_search_y=np.arange(ycen-1.2,ycen+1.2,0.1)
  
  # Just to show an example with rotation around the center
  #f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, sharex=True,figsize=(12,6))
  #ax1.imshow(image_smaller,origin='lower', cmap='hot')

  #rotated_image_smaller=rotate(image_smaller, 180.0, xcen, ycen)
  #ax2.imshow(rotated_image_smaller,origin='lower', cmap='hot')
  
  #absVV=np.absolute(image_smaller)
  #sumabsVV=np.sum(absVV)
  #diffVV=np.absolute(image_smaller-rotated_image_smaller)
  #ax3.imshow(diffVV/sumabsVV,origin='lower', cmap='hot')
  #plt.show()

  # print('\n\nEntered asymmetry function !!')
  asymmetry_all_temp=[] ; asymmetry_all_temp2=[] ; asymmetry_all_temp3=[]
  #for (xc,yc) in good_pixel : # cartesian_product(whereC[0],whereC[1]) :   # cartesian_product(dimx,dimy) :
  for xc in center_search_x :
    for yc in center_search_y :
      # center = [int(xc),int(yc)]
      center = [xc,yc]
      # center = [100,90]
      
      # xc=80.3 ; yc=78.8
      image_smaller_abs=np.absolute(image_smaller)
      rotated_image_smaller=rotate(image_smaller, 180.0, xc, yc)
      differenceASY=np.absolute(image_smaller-rotated_image_smaller)
      maskJ= (image_smaller!=0) & (rotated_image_smaller!=0)
      maskJ=1*maskJ
      
      #print('Other center')
      # Background estimation:
      imagein_smoothed_skysub_reduced_bkg=bkg_search(sizeK,maskgood,imagein_smoothed_skysub)
      asymmetry0b=asymmetry_function_bkg(sizeK,sizeK,imagein_smoothed_skysub_reduced_bkg,maskgoodBKG_smaller)
      # maskgood_smaller serve solo per il primo modo di calcolare l'asimmetria, ma in realta' non serve per il metodo standard
      #print 'asymmetry background (in the two ways):',asymmetry0b,asymmetry_good_b
      asymmetry0=(np.sum(differenceASY*maskJ))/(2*np.sum(image_smaller_abs*maskJ))
      asymmetry1=(np.sum(differenceASY)-asymmetry0b)/(2*np.sum(image_smaller_abs))
      #print 'numeratori 1, 2 e 3 =',np.sum(differenceASY*maskJ),np.sum(differenceASY),asymmetry0b
      
      rotated_maskgood_smaller=rotate(maskgood_smaller, 180.0, xc, yc)
      maskfinal=rotated_maskgood_smaller+maskgood_smaller
      maskfinal=maskfinal-1
      maskfinal[maskfinal<0]=0
      image_smaller_abs_good=image_smaller_abs*maskfinal
      denominatore_comune=np.sum(image_smaller_abs_good)
      differenceASY_good=differenceASY*maskfinal
      asymmetry_good=(np.sum(differenceASY_good))/(2*denominatore_comune)
      #rotated_maskgoodBKG_smaller=rotate(maskgoodBKG_smaller, 180.0, xc, yc)
      #maskfinalBKG=rotated_maskgoodBKG_smaller+maskgoodBKG_smaller
      #maskfinalBKG=maskfinalBKG-1
      #maskfinalBKG[maskfinalBKG<1]=0
      #image_smaller_abs_bkg=image_smaller_abs*maskfinalBKG
      #differenceASY_bkg=differenceASY*maskfinalBKG
      #asymmetry_bkg=(np.sum(differenceASY_bkg))/(2*np.sum(image_smaller_abs))
      
      #print 'asymmetry good =',asymmetry_good
      ##print 'differenceASY =',np.sum(differenceASY)
      ##print 'differenceASY bkg =',np.sum(differenceASY_bkg)
      #print 'asymmetry0 =',asymmetry0
      ##print 'asymmetry bkg =',asymmetry_bkg
      ##print 'asymmetry final =',asymmetry0-asymmetry_bkg
      #f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, sharex=True,figsize=(12,6))
      ##limit1=int(len(imagein)/2-3.2*morphsingle.rhalf_ellip)
      ##limit2=int(len(imagein)/2+3.2*morphsingle.rhalf_ellip)
      #ax1.imshow(image_smaller, origin='lower', cmap='hot') 
      #plt.show()
      #ax2.imshow(rotated_image_smaller, origin='lower', cmap='hot') 
      #plt.show()
      #ax3.imshow(differenceASY, origin='lower', cmap='hot') 
      #plt.show()
      if ( (asymmetry0>0.) & (asymmetry1>0.) & (asymmetry_good>0.)) :
        asymmetry_all_temp=np.append(asymmetry_all_temp,asymmetry_good)
        asymmetry_all_temp2=np.append(asymmetry_all_temp2,asymmetry0)
        asymmetry_all_temp3=np.append(asymmetry_all_temp3,asymmetry1)
      ###### --------------- ######   
  AS3=sorted(asymmetry_all_temp)
  AS4=sorted(asymmetry_all_temp2)
  AS5=sorted(asymmetry_all_temp3)
  return AS3[0],AS4[0],AS5[0]



# def asymmetry_function2
def asymmetry_function2(xcen,ycen,image_smaller,maskgood_smaller,show_image) :
  center_search_x=np.arange(xcen-6,xcen+6,0.3)
  center_search_y=np.arange(ycen-6,ycen+6,0.3)
  
  print('\n\nEntered asymmetry_function2')
  #print 'Show plot first time !!'
  # Just to show an example with rotation around the center
  
  if (show_image==True):
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, sharex=True,figsize=(12,6))
    ax1.imshow(image_smaller,origin='lower', cmap='hot')
    rotated_image_smaller=rotate(image_smaller, 180.0, xcen, ycen)
    
    ax2.imshow(maskgood_smaller,origin='lower', cmap='hot')
    
    absVV=np.absolute(image_smaller)
    sumabsVV=np.sum(absVV)
    diffVV=np.absolute(image_smaller-rotated_image_smaller)
    ax3.imshow(diffVV/sumabsVV,origin='lower', cmap='hot')
    plt.show()
  

  def GetKey00(item):
    return item[0] # Distance from center image
  asymmetry_all_temp=[] ; asymmetry_all_temp2=[] ; asymmetry_all_temp3=[]
  area_all_temp=[] ; xc_temp=[] ; yc_temp=[] ; denominatore_temp=[]
  #for (xc,yc) in good_pixel : # cartesian_product(whereC[0],whereC[1]) :   # cartesian_product(dimx,dimy) :
  #image_smaller_abs=np.absolute(image_smaller)
  for xc in center_search_x :
    for yc in center_search_y :
      #center = [int(xc),int(yc)]
      center = [xc,yc]
      # center = [100,90]
      # print(cent)
      
      # xc=80.3 ; yc=78.8
      rotated_image_smaller=rotate(image_smaller, 180.0, xc, yc)
      differenceASY=np.absolute(image_smaller-rotated_image_smaller)

      rotated_maskgood_smaller=rotate(maskgood_smaller, 180.0, xc, yc)
      maskJ=(maskgood_smaller!=0) & (rotated_maskgood_smaller!=0)
      # Puoi scegliere se mettere o no questa maschera
      #maskJ= (image_smaller!=0) & (rotated_image_smaller!=0)
      maskJ=1*maskJ
      #plt.imshow(maskJ)
      #plt.show()

      image_smaller_abs=np.absolute(image_smaller)*maskJ

      #print('Other center')
      # Background estimation:
      #imagein_smoothed_skysub_reduced_bkg=bkg_search(sizeK,maskgood,imagein_smoothed_skysub)
      #asymmetry0b=asymmetry_function_bkg(sizeK,sizeK,imagein_smoothed_skysub_reduced_bkg,maskgoodBKG_smaller)
      ## maskgood_smaller serve solo per il primo modo di calcolare l'asimmetria, ma in realta' non serve per il metodo standard
      ##print 'asymmetry background (in the two ways):',asymmetry0b,asymmetry_good_b
      #asymmetry0=(np.sum(differenceASY*maskJ))/(2*np.sum(image_smaller_abs*maskJ))
      
      denominatore=(2.*np.sum(image_smaller_abs))
      #print 'denominatore =',denominatore
      asymmetry1=np.sum(differenceASY*maskJ)/denominatore # -asymmetry0b) #/(2*np.sum(image_smaller_abs))
      #print 'numeratori 1, 2 e 3 =',np.sum(differenceASY*maskJ),np.sum(differenceASY),asymmetry0b
      area=np.sum(maskJ)
      
      #rotated_maskgood_smaller=rotate(maskgood_smaller, 180.0, xc, yc)
      #maskfinal=rotated_maskgood_smaller+maskgood_smaller
      #maskfinal=maskfinal-1
      #maskfinal[maskfinal<0]=0
      #image_smaller_abs_good=image_smaller_abs*maskfinal
      #denominatore_comune=np.sum(image_smaller_abs_good)
      #differenceASY_good=differenceASY*maskfinal
      #asymmetry_good=(np.sum(differenceASY_good))/(2*denominatore_comune)
      ##rotated_maskgoodBKG_smaller=rotate(maskgoodBKG_smaller, 180.0, xc, yc)
      ##maskfinalBKG=rotated_maskgoodBKG_smaller+maskgoodBKG_smaller
      ##maskfinalBKG=maskfinalBKG-1
      ##maskfinalBKG[maskfinalBKG<1]=0
      ##image_smaller_abs_bkg=image_smaller_abs*maskfinalBKG
      ##differenceASY_bkg=differenceASY*maskfinalBKG
      ##asymmetry_bkg=(np.sum(differenceASY_bkg))/(2*np.sum(image_smaller_abs))
      
      ##print('asymmetry good =',asymmetry_good)
      ###print('differenceASY =',np.sum(differenceASY))
      ###print('differenceASY bkg =',np.sum(differenceASY_bkg))
      ##print('asymmetry0 =',asymmetry0)
      ###print('asymmetry bkg =',asymmetry_bkg)
      ###print('asymmetry final =',asymmetry0-asymmetry_bkg)
      ##f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, sharex=True,figsize=(12,6))
      ###limit1=int(len(imagein)/2-3.2*morphsingle.rhalf_ellip)
      ###limit2=int(len(imagein)/2+3.2*morphsingle.rhalf_ellip)
      ##ax1.imshow(image_smaller,origin='lower',cmap='hot') 
      ##plt.show()
      ##ax2.imshow(rotated_image_smaller,origin='lower',cmap='hot') 
      ##plt.show()
      ##ax3.imshow(differenceASY,origin='lower',cmap='hot') 
      ##plt.show()
      #if ( (asymmetry0>0.) & (asymmetry1>0.) & (asymmetry_good>0.)) :
      #  asymmetry_all_temp=np.append(asymmetry_all_temp,asymmetry_good)
      #  asymmetry_all_temp2=np.append(asymmetry_all_temp2,asymmetry0)
      #  asymmetry_all_temp3=np.append(asymmetry_all_temp3,asymmetry1)
      #print asymmetry1
      if (asymmetry1>=0.) :
        asymmetry_all_temp=np.append(asymmetry_all_temp,asymmetry1)
        area_all_temp=np.append(area_all_temp,area)
        xc_temp=np.append(xc_temp,xc)
        yc_temp=np.append(yc_temp,yc)
        denominatore_temp=np.append(denominatore_temp,denominatore)
        total00=np.column_stack((asymmetry_all_temp,area_all_temp,xc_temp,yc_temp,denominatore_temp))
      ###### ---------------
  
  print('Qua puo dare errore. Controlla binary_pawlik in caso. L errore e la')
  mgood00=sorted(total00,key=GetKey00)
  
  asym_metry=mgood00[0][0]
  area_best=mgood00[0][1]
  xc_best=mgood00[0][2]
  yc_best=mgood00[0][3]
  denominatore_best=mgood00[0][4]
  #labelnucleus2X=int(mgood00[0][0]) 
  #AS3=sorted(asymmetry_all_temp)
  #AS4=sorted(asymmetry_all_temp2)
  #AS5=sorted(asymmetry_all_temp3)
  
  if (show_image==True):
    rotatedBEST=rotate(image_smaller, 180.0, xc_best, yc_best)
    differenceBEST=np.absolute(image_smaller-rotatedBEST)
  
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, sharex=True,figsize=(12,6))
    ###limit1=int(len(imagein)/2-3.2*morphsingle.rhalf_ellip)
    ###limit2=int(len(imagein)/2+3.2*morphsingle.rhalf_ellip)
    ax1.imshow(image_smaller, origin='lower', cmap='hot') 
    ##plt.show()
    ax2.imshow(differenceBEST, origin='lower', cmap='hot') 
    ##plt.show()
    ax3.imshow(differenceBEST/denominatore_best, origin='lower', cmap='hot') 
    plt.show()  
  return asym_metry,area_best,xc_best,yc_best,denominatore_best #,AS4[0],AS5[0]

#f=23.8 
#print(round(f))
#print(int(np.round(f,0)))
#quit()

def Cutout2D_mine(image_smaller55,xy_center55,shape55) :
    x_center55=round(xy_center55[0])
    y_center55=round(xy_center55[1])
    dimlato=round(shape55/2.)
    tolerance=3
    dimlato_ok_x=dimlato-tolerance
    dimlato_ok_y=dimlato-tolerance
    new_image55=image_smaller55[y_center55-dimlato_ok_y:y_center55+dimlato_ok_y,x_center55-dimlato_ok_x:x_center55+dimlato_ok_x]
    return new_image55


def asymmetry_simple(image_smaller,maskgood_smaller,xc,yc) :
  # Simple function which does not subtract the background !!!
  
  # print(xc,yc)
  # print(image_smaller.shape[0]-1)
  # plt.imshow(image_smaller)
  # plt.show()

  
  cutout2 = Cutout2D(image_smaller, (xc, yc), image_smaller.shape[0])
  cutout2d=cutout2.data
  rotated_image_smaller=np.rot90(cutout2d, 2)
  differenceASY=np.absolute(cutout2d-rotated_image_smaller)
  image_smaller_abs=np.absolute(cutout2d)
  denominatorexx=2*np.sum(image_smaller_abs)
  asymmetry1xx=np.sum(differenceASY)/denominatorexx
  #cutout2d=Cutout2D_mine(image_smaller,(xc, yc), image_smaller.shape[0])
  # cutout3 = Cutout2D(maskgood_smaller, (xc, yc), maskgood_smaller.shape[0]-2)
  

  #rotated_image_smaller=rotate(image_smaller, 180.0, xc, yc)
  #print(cutout2d.shape)
  #print(type(cutout2d))
  #sys.exit()
  # rotated_image_smaller=np.rot90(cutout2d, 2)
  '''
  rotated_image_smaller=skimage.transform.rotate(image_smaller, 180.0, center=(yc, xc))
  differenceASY=np.absolute(image_smaller-rotated_image_smaller)
  image_smaller_abs=np.absolute(image_smaller)
  denominatorexx=2*np.sum(image_smaller_abs)
  asymmetry1xx=np.sum(differenceASY)/denominatorexx
  '''

  # rotated_maskgood_smaller=np.rot90(cutout3, 2)
  # maskJ=(cutout3d!=0) & (rotated_maskgood_smaller!=0)
  # # Puoi scegliere se mettere o no questa maschera
  # #maskJ= (image_smaller!=0) & (rotated_image_smaller!=0)
  # maskJ=1*maskJ
  # image_smaller_abs=np.absolute(cutout2d)*maskJ  
  #print(rotated_image_smaller)
  
  #  sys.exit()
  #print('Other center')
  # Background estimation:
  #imagein_smoothed_skysub_reduced_bkg=bkg_search(sizeK,maskgood,imagein_smoothed_skysub)
  #asymmetry0b=asymmetry_function_bkg(sizeK,sizeK,imagein_smoothed_skysub_reduced_bkg,maskgoodBKG_smaller)
  ## maskgood_smaller serve solo per il primo modo di calcolare l'asimmetria, ma in realta' non serve per il metodo standard
  ##print 'asymmetry background (in the two ways):',asymmetry0b,asymmetry_good_b
  #asymmetry0=(np.sum(differenceASY*maskJ))/(2*np.sum(image_smaller_abs*maskJ))
  
  #denominatore=(2.*np.sum(image_smaller_abs))
  #print(image_smaller_abs)
  #print('denominatore =',denominatore)
  #sys.exit()
  # asymmetry1=np.sum(differenceASY*maskJ)/denominatore
  #print('asymmetry1 =',asymmetry1xx)

  if plot_asymmetry==True :
    
    #print('xc and yc =',xc,yc)
    #plt.imshow(image_smaller,origin='lower')
    ##  #plt.imshow(rotated_maskgood_smaller)
    ##  #plt.imshow(image_smaller_abs)
    #plt.title('image_smaller',fontsize=24)
    #plt.show()
    
    plt.imshow(cutout2d,origin='lower')
    ##  #plt.imshow(rotated_maskgood_smaller)
    ##  #plt.imshow(image_smaller_abs)
    plt.title('cutout2d',fontsize=24)
    plt.show()
  
    plt.imshow(differenceASY,origin='lower')
    #  #plt.imshow(rotated_maskgood_smaller)
    #  #plt.imshow(image_smaller_abs)
    plt.title('differnce for ASYMM',fontsize=24)
    plt.show()
  return asymmetry1xx,denominatorexx


def minimize_mu(image_enter,xcenterA,ycenterA) :
  # print('Image-enter deve essere quadrata !!!!!!')
  mu_all=[] ; x_all=[]  ; y_all=[]
  for rt in xcenterA :
      for rw in ycenterA :
          flussoT=[]  ; coordinateT=[]
          for pt in np.arange(image_enter.shape[0]) : 
              for pw in np.arange(image_enter.shape[1]) :
                if (image_enter[pt,pw]>0) : 
                  flussoT.append(image_enter[pt,pw])
                  quantity00=(pt-rt)**2+(pw-rw)**2
                  coordinateT.append(quantity00)
          flussoT=np.array(flussoT) ; coordinateT=np.array(coordinateT)
          mu_X=sum(flussoT*coordinateT)
          mu_all.append(mu_X)
          x_all.append(rt)
          y_all.append(rw)
  
  x_all=np.array(x_all) ; y_all=np.array(y_all)
  mu_final_value=min(mu_all)
  where_minMU=np.where(mu_all<mu_final_value+0.00001)[0]
  
  x_center_mu=x_all[where_minMU]
  y_center_mu=y_all[where_minMU]
  return mu_final_value,x_center_mu,y_center_mu


def M20_simple(image_enterB,x_center_mm,y_center_mm) :
  flussoTb=[]  ; coordinateTb=[]
  for ptb in np.arange(image_enterB.shape[0]) : 
      for pwb in np.arange(image_enterB.shape[1]) :
        if (image_enterB[ptb,pwb]>0) : 
          flussoTb.append(image_enterB[ptb,pwb])
          quantity00b=(ptb-x_center_mm)**2+(pwb-y_center_mm)**2
          coordinateTb.append(quantity00b)
  
  flussoTb=np.array(flussoTb)
  flusso_parziale=flussoTb/sum(flussoTb)
  data = {'flux': list(flussoTb), 'partial': list(flusso_parziale), 'coord': coordinateTb}
  df = pd.DataFrame(data)
  df_sorted= df.sort_values('flux',ascending=False)

  F1=np.array(df_sorted['flux'])
  F2=np.array(df_sorted['partial'])
  F3=np.array(df_sorted['coord'])
  momenti=F1*F3
  # print('devono essere in ordine decrescente :')
  # print('brightest fluxes =',F1)
  
  parzialiX=[]
  for ppp in np.arange(len(F2)) :
      parzialiX=sum(F2[0:ppp])
      if parzialiX>=0.2 : 
          index_ok=ppp*1
          break
  mu_calculated=sum(momenti[0:index_ok])
  mu_calculated=mu_calculated[0]
  return mu_calculated


def asymmetry_verysimple_4bkg(image_smaller33,maskgood_smaller33,xc33,yc33) :
  # Simple function which does not subtract the background !!!
  
  
  rotated_image_smaller33b=np.rot90(image_smaller33, 2)
  differenceASY33b=np.absolute(image_smaller33-rotated_image_smaller33b)
  image_smaller_abs33b=np.absolute(image_smaller33)
  #denominatore=(2.*np.sum(image_smaller_abs))
  denominatore33b=2*np.sum(image_smaller_abs33b)   # sum(np.ravel(image_smaller_abs33))
  #print('denominatore =',denominatore33)
  # asymmetry1=np.sum(differenceASY*maskJ)/denominatore
  #print('sum difference =',sum(differenceASY33))
  #asymmetry133=sum(np.ravel(differenceASY33))/denominatore33
  asymmetry133b=np.sum(differenceASY33b)/denominatore33b
  
  '''
  rotated_image_smaller33b=skimage.transform.rotate(image_smaller33, 180.0, center=(yc33, xc33))
  differenceASY33b=np.absolute(image_smaller33-rotated_image_smaller33b)
  image_smaller_abs33b=np.absolute(image_smaller33)
  denominatore33b=2*np.sum(image_smaller_abs33b)   # sum(np.ravel(image_smaller_abs33))
  asymmetry133b=np.sum(differenceASY33b)/denominatore33b
  '''
  #print('asymmetry133 =',asymmetry133)
  # if asymmetry133<0 : print('negative asymm') ; sys.exit()
  return asymmetry133b,denominatore33b



def asymmetry_function2_shape(xcen,ycen,image_smaller,maskgood_smaller,show_image) :
  center_search_x=np.arange(xcen-4,xcen+4,0.2)
  center_search_y=np.arange(ycen-4,ycen+4,0.2)
  
  # print('\n\nEntered asymmetry_function2 for shape asymmetry')

  def GetKey00(item):
    return item[0] # Distance from center image
  asymmetry_all_temp=[] ; asymmetry_all_temp2=[] ; asymmetry_all_temp3=[]
  area_all_temp=[] ; xc_temp=[] ; yc_temp=[] ; denominatore_temp=[]
  #for (xc,yc) in good_pixel : # cartesian_product(whereC[0],whereC[1]) :   # cartesian_product(dimx,dimy) :
  image_smaller_abs=np.absolute(image_smaller)
  for xc in center_search_x :
    for yc in center_search_y :
      #center = [int(xc),int(yc)]
      center = [xc,yc]
      #xc=80.3 ; yc=78.8
      rotated_image_smaller=rotate(image_smaller, 180.0, xc, yc)
      #rotated_image_smaller[rotated_image_smaller>0]=1
      differenceASY=np.absolute(image_smaller-rotated_image_smaller)

      # Non divido per il denominatore cosi' magari e' un po' piu' veloce. Il denominatore tanto e' sempre uguale e lo uso per dividere fuori dal loop e addirittura fuori dalla funzione !!!
      #denominatore=(2.*np.sum(image_smaller_abs))
      #print 'denominatore sum =',np.sum(denominatore)
      #print 'numeratore sum =',np.sum(differenceASY)
      asymmetry1=np.sum(differenceASY) #/denominatore # -asymmetry0b) #/(2*np.sum(image_smaller_abs))
      #print 'numeratori 1, 2 e 3 =',np.sum(differenceASY*maskJ),np.sum(differenceASY),asymmetry0b
      area=np.sum(differenceASY)
      
      if (asymmetry1>=0.) :
        asymmetry_all_temp=np.append(asymmetry_all_temp,asymmetry1)
        area_all_temp=np.append(area_all_temp,area)
        xc_temp=np.append(xc_temp,xc)
        yc_temp=np.append(yc_temp,yc)
        #denominatore_temp=np.append(denominatore_temp,denominatore)
        total00=np.column_stack((asymmetry_all_temp,area_all_temp,xc_temp,yc_temp))
      ###### ---------------
  
  # print('Qua puo dare errore. Controlla binary_pawlik in caso. L errore e la')
  mgood00=sorted(total00,key=GetKey00)
  
  asym_metry=mgood00[0][0]
  area_best=mgood00[0][1]
  xc_best=mgood00[0][2]
  yc_best=mgood00[0][3]
  #denominatore_best=mgood00[0][4]
  #labelnucleus2X=int(mgood00[0][0]) 
  #AS3=sorted(asymmetry_all_temp)
  #AS4=sorted(asymmetry_all_temp2)
  #AS5=sorted(asymmetry_all_temp3)
  
  if (show_image_shape_asymmetry==True):
    rotatedBEST=rotate(image_smaller, 180.0, xc_best, yc_best)
    #rotatedBEST[rotatedBEST>0]=1
    differenceBEST=np.absolute(image_smaller-rotatedBEST)
  
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, sharex=True,figsize=(12,6))
    ###limit1=int(len(imagein)/2-3.2*morphsingle.rhalf_ellip)
    ###limit2=int(len(imagein)/2+3.2*morphsingle.rhalf_ellip)
    ax1.imshow(image_smaller, origin='lower', cmap='hot') 
    ##plt.show()
    ax2.imshow(rotatedBEST, origin='lower', cmap='hot') 
    ##plt.show()
    ax3.imshow(differenceBEST, origin='lower', cmap='hot') 
    plt.show()  

  return asym_metry,area_best,xc_best,yc_best #,AS4[0],AS5[0]



def asymmetry_function_bkg(xcen,ycen,image_smaller,maskgood_bkg) :
  center_search_x=np.arange(xcen-2,xcen+2,0.2)
  center_search_y=np.arange(ycen-2,ycen+2,0.2)
  # print('\n\nEntered asymmetry_function_bkg')
  
  asymmetry_all_temp=[] ; asymmetry_all_temp2=[]
  #for (xc,yc) in good_pixel : # cartesian_product(whereC[0],whereC[1]) :   # cartesian_product(dimx,dimy) :
  image_smaller_abs=np.absolute(image_smaller)
  for xc in center_search_x :
    for yc in center_search_y :
      #center = [int(xc),int(yc)]
      center = [xc,yc]
      #center = [100,90]
      #print cent
      
      rotated_image_smaller=rotate(image_smaller, 180.0, xc, yc)
      #rotated_maskgood_bkg=rotate(maskgood_bkg, 180.0, xc, yc)
      differenceASY=np.absolute(image_smaller-rotated_image_smaller)
      # Puoi scegliere se mettere o no questa maschera
      #maskJ= (rotated_maskgood_bkg==0) & (maskgood_bkg==0)        #(image_smaller!=0) & (rotated_image_smaller!=0)
      #maskJ=1*maskJ

      asymmetry1=np.sum(differenceASY)  #/(2*np.sum(image_smaller_abs)) #*maskJ) #/np.sum(maskJ)  # cosi' e' normalizzato al numero di pixels
      
      # maskJ= (image_smaller!=0) & (rotated_image_smaller!=0)
      # maskJ=1*maskJ
      # asymmetry0=(np.sum(differenceASY*maskJ)) # /(2*np.sum(image_smaller_abs*maskJ))
      # 
      # rotated_maskgood_smaller=rotate(maskgood_smaller, 180.0, xc, yc)
      # maskfinal=rotated_maskgood_smaller+maskgood_smaller
      # maskfinal=maskfinal-1
      # maskfinal[maskfinal<0]=0
      # image_smaller_abs_good=image_smaller_abs*maskfinal
      # denominatore_comune=np.sum(image_smaller_abs_good)
      # differenceASY_good=differenceASY*maskfinal
      # asymmetry_good=(np.sum(differenceASY_good))  #/(2*denominatore_comune)
      # if ( (asymmetry0>0.) & (asymmetry_good>0.)) :
      #   asymmetry_all_temp=np.append(asymmetry_all_temp,asymmetry_good)
      #   asymmetry_all_temp2=np.append(asymmetry_all_temp2,asymmetry0)
      
      asymmetry_all_temp=np.append(asymmetry_all_temp,asymmetry1)
      ###### ---------------
  AS3=sorted(asymmetry_all_temp)
  #AS4=sorted(asymmetry_all_temp2)
  # print('numerator for asymmetry background =',AS3[0])
  return AS3[0] #,AS4[0]


def bkg_search(sizeK,singlegood,maskgood,segm,imagein_smoothed_skysub,show_result):
  threshold=1e4
  # Quadro e' uguale a 1 dove posso cercare il background
  quadro=np.zeros(maskgood.shape) ; quadroL=len(quadro)/2 ; quadroQ=len(quadro)/2-(sizeK+3)
  quadro[quadroL-quadroQ:quadroL+quadroQ,quadroL-quadroQ:quadroL+quadroQ]=1
  
  radius_bkg=[sizeK,sizeK/2,sizeK/4]
  whereB=np.where( (maskgood==0) & (quadro==1)) # & (imagein_masked!=0) )
  #print 'How many pixels has the background =', len(whereB[0])
  #pixel_search=cartesian_product(whereB[0],whereB[1])
  pixel_searchB=[  (whereB[0][ig],whereB[1][ig]) for ig in np.arange(len(whereB[0])) ]
  #print 'pixel_search =',pixel_search
  #imagein_smoothed_skysub2=imagein_smoothed_skysub*1
  #imagein_smoothed_skysub2[maskgood==1]=threshold # cioe 1e6
  #pixel_fluxesB=np.array([ imagein_smoothed_skysub2[i] for i in pixel_searchB])
  #imagein_reduced_smoothed_skysub_bkg=(1-background_reduced)*imagein_reduced_smoothed_skysub
  #imagein_smoothed_skysub_bkg=(1-maskgood)*imagein_smoothed_skysub
  #find=False
  #plt.imshow(binary_pawlik, origin='lower', cmap='hot')
  #plt.show()
  ##plt.imshow(imagein_smoothed_skysub_bkg, origin='lower', cmap='hot')
  ##plt.show()
  find=False
  wp=0
  for radiusX in radius_bkg :
    for (yc,xc) in pixel_searchB :
      #hB, wB = maskgood.shape[:2]
      centerB = [int(xc),int(yc)]
      ##mask = create_circular_mask(hB, wB, center=centerB,radius=radiusX)
      #mask_quadrata= np.zeros(maskgood.shape)
      #mask_quadrata[centerB[1]-radiusX:centerB[1]+radiusX,centerB[0]-radiusX:centerB[0]+radiusX]=1
    
      # masked_bkg_all = imagein_smoothed_skysub2.copy()
      # masked_bkg_all[~mask] = 0 
      # masked_bkg=masked_bkg_all[centerB[1]-radiusX:centerB[1]+radiusX,centerB[0]-radiusX:centerB[0]+radiusX]
      
      segm_reduced=segm[centerB[1]-radiusX:centerB[1]+radiusX,centerB[0]-radiusX:centerB[0]+radiusX]
      #masked_bkg_all=imagein_smoothed_skysub2*mask_quadrata
      #masked_bkg=imagein_smoothed_skysub2[centerB[1]-radiusX:centerB[1]+radiusX,centerB[0]-radiusX:centerB[0]+radiusX]
      #maskgood_bkg=1-maskgood[centerB[1]-radiusX:centerB[1]+radiusX,centerB[0]-radiusX:centerB[0]+radiusX]   # Create inverted ma of maskgood (invert 0 and 1)
      #plt.imshow(masked_bkg, origin='lower', cmap='hot')
      #plt.show()
      #quit()
      if ( np.sum(segm_reduced)==0 ):
        find=True
        break
    if find :
      # print('Trovato skybox')
      imageinHH=imagein_smoothed_skysub*1
      maskHH=np.zeros(imagein_smoothed_skysub.shape)
      maskHH[centerB[1]-radiusX:centerB[1]+radiusX,centerB[0]-radiusX:centerB[0]+radiusX]=1

      image_out=imagein_smoothed_skysub[centerB[1]-radius_bkg[wp]:centerB[1]+radius_bkg[wp],centerB[0]-radius_bkg[wp]:centerB[0]+radius_bkg[wp]]
      mask_out=maskgood[centerB[1]-radius_bkg[wp]:centerB[1]+radius_bkg[wp],centerB[0]-radius_bkg[wp]:centerB[0]+radius_bkg[wp]]
      if (show_result==True):
       f, (ax1,ax2) = plt.subplots(1, 2,sharey=True,sharex=True)
       ax1.imshow(imageinHH, origin='lower', cmap='hot')
       ax2.imshow(imageinHH*maskHH, origin='lower', cmap='hot')
       #ax1.set_title('Original image') #plt.title(str(IDgal)) #if (show_all==True) :
       plt.show()
      break
    else :
      wp+=1
      print('Not found any good skybox')
      #print "Potresti, quando calcoli l'asimmetria del background, mettere una maschera alla fine a dove ci sono degli oggetti o parti di oggetti"
      #print " Cercare uno skybox piu' piccolo non credo vada bene, perche' dopo non e' consistente piu' con la size usata per l'asimmetria della galassia"
  return wp,mask_out,image_out



def asymmetry_bkg_simple(sizeK7_max,maskgood77,image77,Rmax77,which_calc,smoothradiusX745):
  # sizeK should be integer !!!
  # Usa uno skybox definito dai !! sizeK e' fissato nella simple version !!

  # THis function is very important, because if it does not find a good background, it tries with a smaller background box !!!
  
  # print('asymmetry-bkg-simple') 
  # print('shape input image =',image77.shape)
  # print('shape maskgood input =',maskgood77.shape)
  asymm_array23=[] ; smoothness_background745=[]
  # print('sizeK7_max =',sizeK7_max)
  assert sizeK7_max>0
  for sizeK7 in [round(sizeK7_max),round(sizeK7_max/2),round(sizeK7_max/3),round(sizeK7_max/4)] :

    centrix=np.arange(round(sizeK7/2+1),round(image77.shape[0]-sizeK7/2-1),1)
    centriy=np.arange(round(sizeK7/2+1),round(image77.shape[1]-sizeK7/2-1),1)
  
    #smoothradiusX745=Rmax77/4.
    factor_smooth745=1

    asymm_array23=[] ; smoothness_background745=[]
    for xcc in centrix :
      for ycc in centriy :
        cutout23 = Cutout2D(image77, (xcc, ycc), sizeK7)
        cutout23d=cutout23.data
        mask23 = Cutout2D(maskgood77, (xcc, ycc), sizeK7)
        mask23d=mask23.data
        
        #cutout23d = Cutout2D_mine(image77, (ycc, xcc), sizeK7)
        #mask23d = Cutout2D_mine(maskgood77, (ycc, xcc), sizeK7)
  
        if (mask23d!=[]) & (sum(np.ravel(mask23d))==0) :
          #print('Sum zero')
          if which_calc=='asymmetry' :
            asymm23,denom23=asymmetry_verysimple_4bkg(cutout23d,mask23d,xcc,ycc) # Ricorda che mask23d non viene usata !!!
            asymm_array23.append(asymm23)
            smoothness_background745.append(0)
  
          elif which_calc=='smoothness' :
            # Calculate smoothness background
            imagesmooth745=smooth(cutout23d,smoothtype,int(smoothradiusX745), factor_smooth745)  # Il 15 alla fine che vuol dire ??
            # Calculate the residual image
            residuals_orig745=cutout23d-imagesmooth745  # image_4_clumpiness e' immagine background subtracted
            residuals745=residuals_orig745*1
            residuals745[residuals745 < 0] = 0
            #plt.imshow(imagesmooth745)
            #plt.imshow(residuals745)
            #plt.show()
            
            padx=1
            #residuals745=np.absolute(np.ravel(residuals745[padx:-padx,padx:-padx]))
            #denomin_745=np.absolute(np.ravel(cutout23d[padx:-padx,padx:-padx]))
            residuals745=np.absolute(np.ravel(residuals745))
            denomin_745=np.absolute(np.ravel(cutout23d))
  
            _smoothness_background745=sum(residuals745)/sum(denomin_745)
            #print('\n\n Somme')
            #print(sum(residuals745))
            #print(sum(denomin_745))
            smoothness_background745.append(_smoothness_background745)
            asymm_array23.append(0)
  
    if (asymm_array23==[]) & (which_calc=='asymmetry') :
      print('Continue loop to smaller bkg region') 
    elif (asymm_array23!=[]) & (which_calc=='asymmetry') :
      asymm_array23=np.array(asymm_array23)
      #asymm_array23=asymm_array23[asymm_array23>0]
      median_asymme_backg=np.median(asymm_array23)
      median_smoothness_background745=0
      #if (which_calc=='asymmetry') : 
      # print('median asymmetry bkg =',median_asymme_backg)
      break

    
    if (smoothness_background745==[]) & (which_calc=='smoothness') :
      print('continue with smaller bkg size')
      #median_smoothness_background745=0
    elif (smoothness_background745!=[]) & (which_calc=='smoothness') : 
      #print(smoothness_background745)
      #print(min(smoothness_background745),max(smoothness_background745))
      #time.sleep(0.2)
      smoothness_background745=np.array(smoothness_background745)
      #print('\nsmoothness background =')
      #print(smoothness_background745)
      #quit()
      median_asymme_backg=0
      median_smoothness_background745=np.median(smoothness_background745)
      #if (which_calc=='smoothness') : 
      print('smoothness bkg =',median_smoothness_background745)
      break
  
  if smoothness_background745==[] : # per asimmetria non dovrebbero esserci problemi
    median_smoothness_background745=0

  if asymm_array23==[] : # per asimmetria non dovrebbero esserci problemi
    median_asymme_backg=0
  # print('which_calc =',which_calc)
  # print('median_smoothness_background745 =',median_smoothness_background745)
  print('best size background =',sizeK7)
  #quit()

  return median_asymme_backg,median_smoothness_background745





# -------------------------------------------------------------------------------



# MAKE SEGMENTATION IMAGE 

def make_segmentation(segmap_detection,segmentation_type,imagein,IDpixscaleFF,input_images_folder,IDgal,banda_obj,galaxy_index,test_segmap_index) : 
  save_segmap=True
  deblend_segmap=False

  # maskgood definition
  if segmentation_type=='petrosian_based' : # Dovrebbe essere Rmax based magari, non Petrosian based
    print('\nPetrosian radius')
    centroid_image=(imagein.shape[0]/2.,imagein.shape[1]/2.)
    print('centroid single =',centroid_image)
    radii_test_L=np.arange(1,10,0.2)
    mu_array_11,mu_array_22=petrosian_radius(imagein,centroid_image,radii_test_L)
    #print('Mu =',mu_array_11)
    #print('Mu 22 =',mu_array_22)
    mu_array_11=np.array(mu_array_11) ; mu_array_22=np.array(mu_array_22)

    #print(mu_array_11)
    print('-----------\n')
    #print(mu_array_22)
    try :
      petr_index=np.where(mu_array_22<=0.2)[0]
      petr_index_ok=petr_index[0]
      petrosian_radius=radii_test_L[petr_index_ok]
    except :
      petrosian_radius=20
    print('Petrosian radius =',petrosian_radius)
    print("Controlla che sia piu' basso di quello che ha calcolato statmorph !!!!")
    time.sleep(sleep_time)
    #petrosian_radius=petrosian_radius(imagein,centroid_single)

    # STEP 1) Create the segmentation map :
    sigma_conv=petrosian_radius/5.
    #kernel_1 = Gaussian2DKernel(x_stddev=FWHM_F090W/2.355/pixel_scale_SW,y_stddev=FWHM_F090W/2.355/pixel_scale_SW)
    #image_conv_1 = scipy_convolve(image1a, kernel_1, mode='same', method='direct')
    imagein_conv = scipy.ndimage.gaussian_filter(imagein, sigma_conv, mode='constant')

    # STEP 1) Flux at petrosian radius :

    #circ_mask1_original= create_circular_mask(imagein.shape[0], imagein.shape[1], center=centroid_image, radius=petrosian_radius)
    circ_mask1a= create_circular_mask(imagein.shape[0], imagein.shape[1], center=centroid_image, radius=petrosian_radius-0.5)
    circ_mask2a= create_circular_mask(imagein.shape[0], imagein.shape[1], center=centroid_image, radius=petrosian_radius+0.5)
    circ_mask1a=circ_mask1a*1
    circ_mask2a=circ_mask2a*1
    annulus_mask_b=np.zeros(imagein.shape)
    annulus_mask_b[ (circ_mask1a==0) & (circ_mask2a==1)]=1
    imageLLX=annulus_mask_b*imagein_conv
    imageLLbX=imageLLX[annulus_mask_b==1]
    annulus_image_b=np.ravel(imageLLbX)
    avg_Gini=np.median(annulus_image_b)
    #aperture_Gini = CircularAnnulus(centroid_single, r_in=petrosian_radius-0.5, r_out=petrosian_radius+0.5)
    #masks_Gini = aperture_Gini.to_mask(method='center')  # methods=exact or center
    #masks_Gini.to_image(imagein.shape)
    #data_weighted_Gini = masks_Gini.multiply(imagein)
    #data_weighted_Gini_1d = np.ravel(data_weighted_Gini)
    #data_weighted_Gini_1d = masks_Gini.get_values(imagein)
    #avg_Gini=np.average(data_weighted_Gini_1d)
    print('avg Gini =',avg_Gini)

    # STEP LESS THAN 10 SIGMA THAN NEIGHBOURING PIXELS : 
    mask_galaxy_neigh=np.zeros(imagein.shape)
    for i in np.arange(1,mask_galaxy_neigh.shape[0]-1,1) :
      for j in np.arange(1,mask_galaxy_neigh.shape[1]-1,1) :
        neigh_pixels=np.array([imagein_conv[i-1,j-1],imagein_conv[i,j-1],imagein_conv[i+1,j-1],imagein_conv[i+1,j],imagein_conv[i+1,j+1],imagein_conv[i,j+1],imagein_conv[i-1,j+1],imagein_conv[i-1,j] ])
        std_neigh=np.std(neigh_pixels) ; avg_neigh=np.median(neigh_pixels)
        if (imagein_conv[i,j]<10*std_neigh+avg_neigh) & (imagein_conv[i,j]>avg_neigh-10*std_neigh) :
          mask_galaxy_neigh[i,j]=1

    inner_parts_mask=create_circular_mask(imagein.shape[0], imagein.shape[1], center=centroid_image, radius=petrosian_radius)
    mask_galaxy=np.zeros(imagein.shape)
    mask_galaxy[ (imagein_conv >= avg_Gini) ]=1 # & (inner_parts_mask == 1)        
    final_mask_galaxy_segm=mask_galaxy_neigh*mask_galaxy*inner_parts_mask
    imagein_masked=final_mask_galaxy_segm*imagein
    #plt.imshow(mask_galaxy_segm*imagein,origin='lower')
    maskgood=final_mask_galaxy_segm*1

  # ----------------------------------------------------------------------------------

  elif segmentation_type=='square' :
    which_mask='mask_square'
    mask_source_radius=0.5 # arcseconds

    # Create a rectangular mask :
    # if IDdist_clumps[iyy]==1.5 : 
    #   sizebox=1.5*conv_kpc_to_arcsec/IDpixscale[iyy]
    # elif IDdist_clumps[iyy]==3 :
    #   sizebox=3*conv_kpc_to_arcsec/IDpixscale[iyy]
    # else : print('Errore stop') ; sys.exit()

    if IDdist_clumps[iyy]==size1 : 
      sizebox=size_square_box_small/IDpixscaleFF
    elif IDdist_clumps[iyy]==size2 :
      sizebox=size_square_box_large/IDpixscaleFF
    else : print('Errore stop') ; sys.exit()
    mask_square = np.zeros(imagein.shape)
    mask_square[int(positions[0]-sizebox):int(positions[0]+sizebox)+1,int(positions[1]-sizebox):int(positions[1]+sizebox)+1] = 1

    #mask_source_radius=IDdist_clumps[iyy]  # arcsec
    #sigma_clip = SigmaClip(sigma=3.)
    #apertureV = CircularAperture(positions, r=aperture_arcsec/pixel_scale_TT)
    #print('Pixel scale ID =',IDpixscale[iyy])
    mask_source = create_circular_mask(imagein.shape[0], imagein.shape[1], center=positions, radius=mask_source_radius/IDpixscaleFF)
    #print('Sum mask_source = ',np.sum(mask_source))

    # maskgood definition
    if which_mask=='mask_square' :
      maskgood=mask_square*1
    elif which_mask=='mask_source' :
      maskgood=mask_source*1      

  # ----------------------------------------------------------------------------------

  elif segmentation_type=='photutils' :
    print('Photutils based segmentation')
    firstsegmap_snr=2 ; firstsegmap_npixels=10  

    try : 
      threshold1 = photutils.detect_threshold(imagein, nsigma=firstsegmap_snr)
      #npixelsV = firstsegmap_npixels*1  # minimum number of connected pixels
      segm1 = photutils.detect_sources(imagein, threshold1, firstsegmap_npixels)
      # Although statmorph is designed to process all the sources labeled by the segmentation map, in this example we only focus on the main (largest) source found in the image.
      #areas=segm.areas ; #print segm.areas

      # Keep only the largest segment
      #label = np.argmax(segm.areas) + 1  # Prende la regione con area piu' grande
      #segmap = segm.data == label

      # CALCULATE STATISTICS ON BACKGROUND ALL (mi serve per le simulazioni)
      background_all1 = segm1.data < 0.5
      #print('Created segmap, now convert bool to int')
      background_all1=1*background_all1
      image_bkg1=imagein*background_all1
      background_only=image_bkg1[background_all1>0]
      imagestatS=calc_bg(background_only)
      backS=imagestatS[0]   #'background'
      sigmabackS=imagestatS[1]

      # Riapplica, stavolta con un background migliorato
      threshold2 = backS + (2.0 * sigmabackS)
      segmV = photutils.detect_sources(imagein, threshold1, firstsegmap_npixels)

      # norm=matplotlib.colors.Normalize(vmin=backS-2*sigmabackS, vmax=backS+20*sigmabackS, clip=False)
      # f, (ax1,ax2) = plt.subplots(1, 2,sharey=True,sharex=True)
      # ax1.imshow(imagein, origin='lower', norm=norm, cmap='hot')
      # ax2.imshow(segmV,origin='lower')
      # ax1.set_title('Original image')
      # #plt.colorbar()
      # #if (show_all==True) :
      # plt.show()
      # SOURCE
      sourceallV = segmV.data > 0.5
      #print('Created segmap, now convert bool to int')
      maskgood=1*sourceallV

    except :
      maskgood=np.zeros(imagein.shape)
    # Qui però il problema sta nell'isolare solo l'oggetto centrale !!!!

  # ----------------------------------------------------------------------------------

  elif segmentation_type=='Pawlik' :
    firstsegmap_snr=2 ; firstsegmap_npixels=10
    #try :
    imagein_smoothed= ndi.uniform_filter(imagein, size=3)
    threshold1 = photutils.detect_threshold(imagein_smoothed, nsigma=firstsegmap_snr)
    segm1 = photutils.detect_sources(imagein_smoothed, threshold1, firstsegmap_npixels)
    #plt.imshow(segm1.data)
    #plt.show()
    if np.sum(segm1.data)==0 :
      maskgood=np.zeros(imagein.shape)
      binary_pawlik_4=maskgood*1
      binary_pawlik_4_complete=maskgood*1
    else :
      # CALCULATE STATISTICS ON BACKGROUND ALL (mi serve per le simulazioni)
      background_all1 = segm1.data < 0.5
      #print('Created segmap, now convert bool to int')
      background_all1=1*background_all1
      image_bkg1=imagein_smoothed*background_all1
      background_only=image_bkg1[background_all1>0]
      imagestatS=calc_bg(background_only)
      backS=imagestatS[0]   #'background'
      sigmabackS=imagestatS[1]
      print('back and sigma =',backS,sigmabackS)
      # Riapplica, stavolta con un background migliorato
      threshold2 = backS + (2.0 * sigmabackS)
      segmV = photutils.detect_sources(imagein_smoothed, threshold1, firstsegmap_npixels)
      imagein_smoothed_skysub=imagein_smoothed-np.ones(imagein_smoothed.shape)*backS
      imagein_smoothed_skysub_debl=photutils.deblend_sources(imagein_smoothed, segmV, firstsegmap_npixels, kernel=None, labels=None, nlevels=32, contrast=0.001, mode='exponential', connectivity=8, relabel=True)

      # Save this segmentation image :

      if deblend_segmap==True :
        segmap_filename=input_images_folder+'ID_'+str(IDgal)+'_'+banda_obj+'_segm.fits'
        if save_segmap==True :
          astImages.saveFITS(segmap_filename,segmV)
        print('Saved segmentation map in '+segmap_filename)
        segmap_final=segmV.data
      else :
        segmap_filename=input_images_folder+'ID_'+str(IDgal)+'_'+banda_obj+'_segm.fits'
        if save_segmap==True :
          astImages.saveFITS(segmap_filename,imagein_smoothed_skysub_debl)
        print('Saved segmentation map in '+segmap_filename)
        segmap_final=imagein_smoothed_skysub_debl.data
      
      segmap_final_ones=segmap_final*1
      maskgood=segmap_final*1 
      segmap_final_ones[segmap_final_ones>0]=1
      binary_pawlik_4_complete=segmap_final_ones*1

      if segmap_detection=='automatic' :
      	# Look for the segmentation region closest to the center :
        
        center_x=round(maskgood.shape[0]/2)
        running_x=[center_x,center_x-1,center_x+1,center_x-2,center_x+2]
        center_y=round(maskgood.shape[1]/2)
        running_y=[center_y,center_y-1,center_y+1,center_y-2,center_y+2]
        
        #while True:
        indice_centrale=0
        for tj in running_x :
            for th in running_y :
              indice_centrale=maskgood[tj,th]
              if indice_centrale>0 :
                print('indice centrale =',indice_centrale)
                break

        if indice_centrale>0 :
          maskgood[maskgood!=indice_centrale]=0
          maskgood[maskgood==indice_centrale]=1
        else :
          maskgood=maskgood*0
        binary_pawlik_4=maskgood*1

      else :
        maskgood[maskgood!=galaxy_index]=0
        maskgood[maskgood==galaxy_index]=1
        binary_pawlik_4=maskgood*1

  else :
    print('Problem with segmentation - Exit')
    sys.exit()



  if segmentation_type!='Pawlik' :
    try :
      which_image_segmap_SNpixel='smoothed'  # normal or smoothed
      # Create Pawlik segmentation map :
      imagein_smoothed_skysub_4,bkg_pawlik_4,std_pawlik_4=smoothed_mask_skysub(imagein,maskgood)

      if which_image_segmap_SNpixel=='smoothed' :
        min_pixel_source=40   
        min_pixel_source_complete=10
        image_for_segmap=imagein_smoothed_skysub_4*1
        std_4_segmap=std_pawlik_4*1
        bkg_4_segmap=bkg_pawlik_4*1
        #print('std 4 SNpixel =',std_4_segmap)
        Nsigma=2
        threshold4 = np.ones(image_for_segmap.shape)*std_4_segmap*Nsigma

  

      elif which_image_segmap_SNpixel=='normal' :

        min_pixel_source=10
        image_for_segmap=imagein*1
        bkg_1d=imagein[maskgood==0]
        bkg_4_segmap,std_4_segmap=background_pawlik(bkg_1d)
        #print('std 4 SNpixel =',bkg_4_segmap,std_4_segmap)
        NSIGMA=2
        threshold4 = np.ones(image_for_segmap.shape)*(bkg_4_segmap+NSIGMA*std_4_segmap)

      # CHECK HERE :
      # See definition SN-pixel (see Lotz et al. 2004, formula 5 !!!) : 
      # Prima devi calcolare una buona segmentation map per la galassia, non un quadrato

      # SEGMENTATION MAP ONLY SOURCE : 

      imagein_smoothed_skysub_seg_4 = photutils.detect_sources(image_for_segmap, threshold4 , npixels=min_pixel_source,connectivity=4)
      #plt.imshow(imagein_smoothed_skysub_seg_4.data)
      #plt.show()
      imagein_smoothed_skysub_seg_4_debl=photutils.deblend_sources(image_for_segmap, imagein_smoothed_skysub_seg_4_complete, min_pixel_source_complete, kernel=None, labels=None, nlevels=32, contrast=0.001, mode='exponential', connectivity=8, relabel=True)
      #plt.imshow(imagein_smoothed_skysub_seg_4_complete_debl,origin='lower')
      #plt.show()

      # SEGMENTATION MAP COMPLETE :
      imagein_smoothed_skysub_seg_4_complete = photutils.detect_sources(image_for_segmap, threshold4 , npixels=min_pixel_source_complete,connectivity=4)

      try :
        binary_pawlik_4=imagein_smoothed_skysub_seg_4.data
        binary_pawlik_4=1*binary_pawlik_4
        binary_pawlik_4_complete=imagein_smoothed_skysub_seg_4_complete.data
        binary_pawlik_4_complete=1*binary_pawlik_4_complete
      except :
        binary_pawlik_4=np.zeros(image_for_segmap.shape)
        binary_pawlik_4_complete=np.zeros(image_for_segmap.shape)
        print('No source found binary_pawlik !!!')

    except :
      binary_pawlik_4=np.zeros(image_for_segmap.shape)
      binary_pawlik_4_complete=np.zeros(image_for_segmap.shape)

  return maskgood,binary_pawlik_4,binary_pawlik_4_complete,segmap_final

















