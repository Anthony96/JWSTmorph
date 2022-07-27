#!/usr/bin/python
# -*- coding: utf-8 -*-

###################################################################
__author__ = "Antonello Calabro"
__affiliation__= "INAF - Osservatorio Astronomico di Roma"
#__copyright__ = "Copyright 2007, The Cogent Project"
#__credits__ = []
__license__ = "GPL"
#__version__ = "1.0.1"
__maintainer__ = "Antonello Calabro"
__email__ = "antonello.calabro@inaf.it"
__status__ = "Production"
###################################################################

which_subset='TEST_jwst'    # it's just a name that you give to distinguish among different subsets

# INPUT CATALOG, INPUT AND OUTPUT FOLDERS
# START PROGRAM AND LOOP ON IMAGES  

# To convert to 'from_image_list' you need essentially to update two main parts 

calculate_parameters=True # Default is True
import os,sys,time,glob
import numpy as np
sleep_time=0.01
cwd=os.getcwd()+'/'

# __*****************************************************************__
#            PARAMETERS TO CHANGE BY THE USER                        #
# __*****************************************************************__

ID_list='catalog'      # 'catalog' or 'automatic'
segmap_detection_type='automatic'  # automatic or indices (in the second case you have to write all the indices of the segmentation region corresponding to the galaxy itself)
test_segmap_index=False  # Put True if you have to select the segmentation region of the galaxies
use_statmorph=True
make_images=False
show_figures_clumps=False
show_images_masks_smoothness=False
remove_nucleus=True # You can also put a file nucleus with all zeros and it does not remove anything
# N.B. If you are not interested in the clumpiness parameter, you can skip this step, just put a nuclei file with all zeros.

# SET THE PIXEL SCALE OF YOUR IMAGES 
pixel_scale_SW=0.031 # arcsec/pixel   # Short wavelength channel (from F090W to F200W)
pixel_scale_LW=0.031 # arcsec/pixel   # Long wavelength channel (from F277W to F444W)

# HERE DEFINE HOW MANY BANDS YOU HAVE FOR EACH OBJECT (these can be modified, but be sure that these keywords appear in the filenames of the images)
bande_all=['f090w','f115w','f150w','f200w','f277w','f356w','f444w']
IDpixscale=[pixel_scale_SW,pixel_scale_SW,pixel_scale_SW,pixel_scale_SW,pixel_scale_LW,pixel_scale_LW,pixel_scale_LW]
IDfwhm_banda=[0.034,0.040,0.050,0.066,0.091,0.115,0.145] # FWHM resolution in arcsec of your images ; For JWST they are taken from # Taken from https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-predicted-performance/nircam-point-spread-functions

segmentation_type='Pawlik'  #  'square', 'photutils', 'Pawlik'      # 'petrosian_based' or 'photutils' 
skybox_arcsec=0.5 ; # maximum skybox region for determining the asymmetry and smoothness of the background
size_psf=5 # This is only used by statmorph

# ---------------------------------------------------------------------------------

# CLUMPINESS PARAMETERS (you can skip if you don't care about clumpiness calculation)
clumps_method = 0   # 0=standard (Conselice 2003/4) - DO NOT CHANGE !!!
#npixels = 5  # minimum number of connected pixels
#npixels_original = 5
threshold_clumpiness=3 ; deblend_clumpmap=1 # threshold should be at least 3 - 5 is the best and suggested value 
# deblend_clumpmap - put True if you want to deblend the detected clumpy structures, which might be initially connected each other

# ---------------------------------------------------------------------------------

# __*******************************************************************************__
#          OTHER PARAMETERS NOT TO CHANGE (many might be useless)                 #
# __*******************************************************************************__
dil_times=0 ; cutstamp=False 
MASK=False ; semialtezza=100 ; semibase = 100 ; MASK2=True 
smoothtype='gauss' ; deblend=True ; pix_consecutive=3  
plots=False ; save_int=True ; plots1=False ; plots2=False ; saveoutput=True ; silent=True 
boxsize=8 ; numberclumps=False 
segmentation=True ; sigmasegm= 1.2 ; plotsNUC = False
# smoothtype available : gauss, median, average, boxcar, percentile, bilateral, normal, convolve
#show_all=True ;
# -----------------------------------------------------------------------------------------









# __*******************************************************************************__
#          # INPUT CATALOG, INPUT AND OUTPUT FOLDERS                              #
# __*******************************************************************************__

# Put all the images, catalog with ID and ascii table with the nuclei and the size of the central mask for the object
input_folder=cwd+"data/"
input_folder_catalog=cwd+"catalog/"
output_folder=cwd+'results/'
output_folder_segmap=cwd+'segmaps/'
filenameoutput=output_folder+"results_subset_"+which_subset+".txt"

  # IMPORTING APERTURE MASKS
  table_apertures=np.genfromtxt(input_folder_catalog+file_apertures,names=True,dtype=None)
  size_square_source_all=table_apertures['size']
  ID_size_all=table_apertures['ID']


if ID_list=='catalog' :
  file_catalog=input_folder_catalog+"subset1_catalog.txt" # there should be at least a column with ID, redshift, size
  file_apertures='subset1_boxsize.txt'
  table_catalog=np.genfromtxt(file_catalog,names=True,dtype=None)
  IDobj_all=table_catalog['ID'] # ; ra_all=table_catalog['RA'] ; dec_all=table_catalog['DEC']
  redshift_all=table_catalog['redshift']
  size_square_source_all=table_catalog['size']
  # RA and DEC not needed, cutouts should be already centered.
  # ID in the input catalog should be integers
  IDlista=IDobj_all*1
  if remove_nucleus==True :
    file_catalog_N=input_folder_catalog+"file_nuclei.txt" #
    table_catalog_N=np.genfromtxt(file_catalog_N,names=True,dtype=None)
    IDxx=table_catalog_N['ID']
    N1x=table_catalog_N['N1x']
    N1y=table_catalog_N['N1y']
    N2x=table_catalog_N['N2x']
    N2y=table_catalog_N['N2y']
else :
  all_images_x=glob.glob(input_folder+'*w.fits')
  all_images=sorted(all_images_x)
  IDlista=np.arange(len(all_images))
  for tt in all_images : print(tt)

# OPEN GALAXY INDICES SEGMENTATION MAPS : 
# IN_catalog=input_folder+"galaxy_segmap_indices.txt"
# indices_catalog=np.genfromtxt(IN_catalog,names=True,dtype=None)


# __*******************************************************************************__
#          # MORPHOLOGICAL PROPERTIES YOU WANT TO COMPUTE                             #
# __*******************************************************************************__

calculate_area_SNpixel=True 
replace_secondary_sources=True
calculate_Rmax=True
calculate_smoothness=True
calculate_asymmetry_center=True # contains also asymmetry calculation
calculate_shape_asymmetry_easier=True
calculate_concentration=True
calculate_gini=True
calculate_m20=True
calculate_clumpiness=True
calculate_shape_asymmetry_old=False
calculate_asymmetry_oldstyle=False 
calculate_m20_oldstyle=False
show_figures=False
show_figures_Gini=False
show_image_shape_asymmetry=False
plot_asymmetry=False

# ---------------- ---------------- ---------------- -------------- -------------

# --------------------------------------------------------------------------------------------
# Import necessary modules

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy,math,glob,random
import pandas as pd

import skimage
from skimage import data
from skimage.morphology import disk
from skimage.filters import rank
from skimage.morphology import erosion, dilation, opening, closing, white_tophat
from skimage.morphology import black_tophat, skeletonize, convex_hull_image
from astLib import astCalc,astCoords,astImages,astPlots,astStats,astWCS

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

from scipy import signal,ndimage
#import matplotlib.animation as animation
import photutils
#print('Photutils version =',photutils.__version__)
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
#redshift=6
#_conv_kpc_to_arcsec=cosmo.arcsec_per_kpc_proper(redshift)
#_conv_arcsec_to_kpc=1./_conv_kpc_to_arcsec
#conv_kpc_to_arcsec=_conv_kpc_to_arcsec.value
#conv_arcsec_to_kpc=_conv_arcsec_to_kpc.value
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,FixedLocator
from scipy import signal,ndimage,stats

from auxiliary_easy import import_table,region,create_circular_mask,fit_gauss,replace_secondary_sources_with_bkg,calc_bg,rmsdiff,cartesian_product,makeGaussian,petrosian_radius,binary_detection_pawlik,background_pawlik,RMAX_calc,smoothed_mask_skysub,asymmetry_function,asymmetry_function2,asymmetry_simple,Cutout2D_mine,minimize_mu,asymmetry_verysimple_4bkg,M20_simple,asymmetry_function2_shape,asymmetry_function_bkg,bkg_search,asymmetry_bkg_simple
from auxiliary_easy import make_segmentation











# __*******************************************************************************__
#                  # START PROGRAM AND LOOP ON IMAGES                             #
# __*******************************************************************************__

ID_all=[] ; Rmax_all=[] ; clumpiness_all=[] ; clumpiness2_all=[] ; clumpiness3_all=[] ; clumpiness4_all=[] ; clumpiness5_all=[] ; clumpiness6_all=[]
number_clumps=[] ; nc_nodebl=[] ; areagood=[] ; distancegood=[] ; photometry=[] ; peakgood=[]
sn_per_pixel=[] ; sn_per_pixel2=[] ; sn_per_pixel3=[] ; area=[] ; distlotz_all=[]
GiniM20_1=[] ; asymmetry_all1=[] ; asymmetry_all2=[] ; concentration_all=[] ; smoothness_all=[]
merger_all=[] ; gini_all=[] ; m20_all=[]
ellipticity_all=[] ; elongation_all=[] ; rpetro_all=[] ; rhalf_all=[]
shape_asymmetry_all=[] ; photoband_all=[] ; final_table=[]
sn_per_pixel_SM=[] ; asymmetry_SM=[] ; concentration_SM=[] ; smoothness_SM=[] ; gini_SM=[] ; m20_SM=[]
ellipticity_SM=[] ; elongation_SM=[] ; rpetro_SM=[] ; rhalf_SM=[] ; flag_SM=[]
centerax_all=[]
centeray_all=[]

if calculate_parameters==True :
  counter989=0
  
  # START LOOP ON IMAGES 

  #for iyy in lista_immagini_ID :  # THIS IS FOR FIGURES OF 6 in a row !!!!!
  for iyy in np.arange(len(IDlista)) :  # So I do them all
  #for iyy in [2,3] :
    
    IDgal=IDlista[iyy]
    print('\n\n\nID-image =',IDgal)
    print('processing ID-image '+str(iyy)+' ... ID = '+str(IDgal))
    time.sleep(sleep_time)
    count_image=0

    redshift_source=redshift_all[ikk]
    _conv_kpc_to_arcsec=cosmo.arcsec_per_kpc_proper(redshift) # cosmology defined above
    conv_kpc_to_arcsec=_conv_kpc_to_arcsec.value
    
    dd=np.where(IDxx==IDgal)[0]
    #print(N1x)
    N1_x=N1x[dd[0]]
    N1_y=N1y[dd[0]]
    N2_x=N2x[dd[0]]
    N2_y=N2y[dd[0]] 
    print('nuclei =',N1_x,N1_y,N2_x,N2_y)

    # aperture mask
    size_square_source=size_square_source_all[iyy]
    print('\nsize square source =',size_square_source)


    # INITIALIZE FIGURE :::
    if ( (test_segmap_index == False) & (make_images==True)) :
      fig66, ((axs)) = plt.subplots(1, 6, figsize=(18,3)) # , sharex=True) sharey=True,

    for jkk in np.arange(len(bande_all)) :
    #for jkk in [1,2] :
      
      # This is needed only if importing segmentation index
      # IDH=np.where(indices_catalog['ID']==IDgal)[0]
      # galaxy_index_7=indices_catalog['index'][IDH] # The segmentation region of the galaxy
      # galaxy_index=galaxy_index_7[jkk]
      # print('galaxy index =',galaxy_index)
      #quit()
      galaxy_index=1

      counter989+=1
      banda_obj=bande_all[jkk]
      banda_obj_index=jkk*1
      IDpixscale_banda=IDpixscale[jkk]

      print('\n\nLooking at band '+banda_obj)
      print('Open original image, segmentation map, and mask for the same galaxy')

      smoothradiusX22=round(1*conv_kpc_to_arcsec/IDpixscale_banda)  # per smoothness and clumpiness
      smoothradius_4segmap=round(1*conv_kpc_to_arcsec/IDpixscale_banda)  # per smoothness and clumpiness
      print('smoothradius for segmap =',smoothradius_4segmap)

      IDnamefile=input_folder+'ID_'+str(IDgal)+'_'+banda_obj+'.fits' # write image filenames in the same format

      img=fits.open(IDnamefile)
      imagein=img[0].data # Remember that primary header is empty
      header=img[0].header 
      #WCS=astWCS.WCS(header, mode = "fits")
      print('Image shape =',imagein.shape)
      
      # *****************************************************************************
      # Define segmentation images and masks
      # *****************************************************************************
    
      segmap=np.ones(imagein.shape)
      # Here create maskgood nel caso in cui volessi applicare una maschera all'oggetto (per ora nessuna segmentation map !!!)
      positions = (imagein.shape[0]/2, imagein.shape[1]/2)
      # Estimate background first : 
      aperture_arcsec=0.3 ; aperture_bkg1=0.5 ; aperture_bkg2=0.8
  


      # PLOT IMAGES ( in case you decide to have a figure) :
      if ( (make_images==True) & (test_segmap_index==False) ) :
        #ax1=subpanels[count_image]
        imagestat1=calc_bg(imagein)
        back1=imagestat1[0]   #'background'
        sigmaback1=imagestat1[1]
        norm1=matplotlib.colors.Normalize(vmin=back1-2*sigmaback1, vmax=back1+5*sigmaback1, clip=False)
        axs[count_image].imshow(imagein, origin='lower', norm=norm1) #, cmap='cool')
        #ax1.set_title('Original')
        axs[count_image].text(0.05, 0.95, bande_all[jkk], horizontalalignment='left',verticalalignment='top', fontsize=14, color='w',transform=axs[count_image].transAxes, weight='bold')
        axs[count_image].set_xticks([])
        axs[count_image].set_yticks([])
      

      # DERIVE SEGMENTATION IMAGE 
      maskgood,binary_pawlik_4,binary_pawlik_4_complete,segmap_final=make_segmentation(segmap_detection_type,segmentation_type,imagein,IDpixscale_banda,input_folder,IDgal,banda_obj,galaxy_index,test_segmap_index,size_square_source,res,smoothradius_4segmap)
      

      # RIGA PER CREARE SORGENTI SECONDARIE DA RIMPIAZZARE CON SYNTHETIC SKY ANCHE DOVE E' TUTTO ZERO (immagini tagliate etc ...)
      problems_cut= imagein==0.
      problems_cut=problems_cut*1
      binary_pawlik_4_complete[problems_cut==1]=1
      #plt.imshow(binary_pawlik_4_complete,origin='lower')
      #plt.show()



      ###############################################################################
      # Calculate statistics image and derive bkg subtracted image                  #
      ###############################################################################

      background_all2 = maskgood < 0.5
      background_all2=1*background_all2
      image_bkg2=imagein*background_all2
      background_only2=image_bkg2[background_all2>0]
      
      #if ( (IDgal==76) &  (banda_obj=='f150w') ) :
      #  plt.imshow(imagein)
      #  plt.show()
      #  plt.imshow(maskgood)
      #  plt.show()
      # INTERNAL CHECK : 
      '''
      fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 6))
      ax1.imshow(imagein, origin='lower', interpolation='nearest',cmap='hot')
      ax1.set_title('Image ID '+str(IDgal)+' '+banda_obj)
      ax2.imshow(binary_pawlik_4, origin='lower', interpolation='nearest',cmap='hot')
      ax2.set_title('Segm object')
      ax3.imshow(binary_pawlik_4_complete, origin='lower', interpolation='nearest',cmap='hot')
      ax3.set_title('Segmentation complete')
      plt.show()
      '''
      #quit()

      # --------------------------------------------------------------
      imagestatBK=calc_bg(background_only2)
      backBK=imagestatBK[0]   #'background'
      sigmabackBK=imagestatBK[1]
      print('background mean and std =',backBK,sigmabackBK)
      
      # Subtract background
      image_new=imagein-backBK    # image_new is imagein background subtracted
      noisy_img = np.random.normal(0,sigmabackBK,imagein.shape)  
      #print('len(smoothRadius) =',len(smoothRadius))
      
      image_new_masked=maskgood*image_new
      image_new_conv=maskgood*0
      #singlegood=maskgood*1


      # -----------------------------------------------------------------------------------------
      
      if (calculate_area_SNpixel==True) :
        
        print('\n\nCalculate S/N per pixel (Lotz +2004) and standard way')
        area2=np.sum(maskgood)
        
        #continue_T=True 
        #if continue_T==True :
        try :
          only_source_4=imagein[ maskgood==1 ]
          only_back_4=imagein[ binary_pawlik_4_complete== 0]

          imagestat1=calc_bg(only_back_4)
          back1=imagestat1[0]   #'background'
          sigmaback1=imagestat1[1]
          
          bkgsub_source=only_source_4-back1

          sigma_average=np.median(bkgsub_source/sigmaback1)
          
          # bkg_pawlik_4,std_pawlik_4=background_pawlik(only_back_4)
          # denomin44=np.abs(only_source_4+std_pawlik_4**2)
          # #print('only_source_4 =',len(only_source_4),sum(only_source_4))
          # sommatoria_4=np.abs(only_source_4)/np.sqrt(denomin44)
          # SNpixel_Lotz04=sum(sommatoria_4)/len(only_source_4)
          # #print('bkg_pawlik_4,std_pawlik_4,sommatoria_4, ',bkg_pawlik_4,std_pawlik_4,sommatoria_4,len(only_source_4))
          # print('SNpixel_Lotz04 =',SNpixel_Lotz04*1)
          ##binary_pawlik,centroid=binary_detection_pawlik(imagein_smoothed_skysub*maskgood,std_pawlik)
          ##binary_pawlik[binary_pawlik>0]=1
          SNpixel_Lotz04=-9
          
          # SNpixel_img1=image_new*singlegood/sigmabackBK
          # #SNpixel_img2=imagepxm_all*singlegood/sigmabackBK
          # SNpixel_img1_X=SNpixel_img1[singlegood>0]
          # #SNpixel_img2_X=SNpixel_img2[imagepxm_all>0]
          # SNpixel2=np.median(SNpixel_img1_X)

          # SNpixel_img1=imagein*maskgood/sigmabackBK
          # #SNpixel_img2=imagepxm_all*singlegood/sigmabackBK
          # SNpixel_img1_X=SNpixel_img1[maskgood>0]
          # #SNpixel_img2_X=SNpixel_img2[imagepxm_all>0]
          # SNpixel2=np.median(SNpixel_img1_X)
          # SNpixel2=np.round(SNpixel2,1)
          # print('SNpixel standard =',SNpixel2)
          SNpixel2=np.round(sigma_average*1,3)
          area2=-9
          SNpixel_img1=-9
          
        except :
          SNpixel_Lotz04=-9
          SNpixel2=-9
          area2=-9
          SNpixel_img1=-9
  
      else :
        SNpixel2=-9 ; SNpixel_img1=-9
        SNpixel_Lotz04=-9
        area2=-9

      photoband_all.append(banda_obj_index)
      sn_per_pixel2.append(np.round(SNpixel2,1))
      sn_per_pixel3.append(np.round(SNpixel_Lotz04,1))
      area.append(0)

      # ------------------------------------------------------------------------------


    
      if replace_secondary_sources==True :
          print('\n\nReplace secondary sources')
        
        #try :
          check_result=False
          image_new_replaced=replace_secondary_sources_with_bkg(image_new,binary_pawlik_4,binary_pawlik_4_complete)
          imagein_replaced=replace_secondary_sources_with_bkg(imagein,binary_pawlik_4,binary_pawlik_4_complete)
          


          if check_result==True :
            f, ((ax1, ax2, ax3)) = plt.subplots(1, 3, sharey=True) #, sharex=True)
            imagestat1=calc_bg(image_new[maskgood==1])
            back1=imagestat1[0]   #'background'
            sigmaback1=imagestat1[1]
            norm1=matplotlib.colors.Normalize(vmin=back1-2*sigmaback1, vmax=back1+7*sigmaback1, clip=False)
            imagestat3=calc_bg(image_new[maskgood==0])
            #imagestat3=calc_bg(np.ravel(image_4_clumpiness))
            back3=imagestat3[0]   #'background'
            sigmaback3=imagestat3[1]
            norm3=matplotlib.colors.Normalize(vmin=back3-2*sigmaback3, vmax=back3+7*sigmaback3, clip=False)
            ax1.imshow(image_new, origin='lower', norm=norm1, cmap='hot')
            ax1.set_title('Image original')
            ax2.imshow(image_new_replaced, origin='lower', cmap='hot', norm=norm1)
            ax2.set_title('Secondary sources replaced')
            # Immagine originale con maschera
            #ax3.imshow(maskgood*image_4_clumpiness, origin='lower', norm=norm1, cmap='hot')
            #ax3.set_title(str(IDgal))
            #plt.title(str(IDgal))
            #plt.colorbar()
            ax3.imshow(image_new-image_new_replaced, origin='lower', norm=norm3, cmap='hot')
            ax3.set_title('Difference')
            # Segmentation map (Deblended) on clumps
            plt.subplots_adjust(wspace=0.2, hspace=0.01)
            #plt.tight_layout()
            #if (show_all==True) :
            namefile=output_folder+str(IDgal)+'_secondary_replaced_'+banda_obj+'.png'
            plt.savefig(namefile,dpi=80)
            plt.show()
        #except :
        #  print('Just keeping the original images. No replacement made.')
        #  image_new_replaced=image_new*1
        #  imagein_replaced=imagein*1
  



      # print('Now show images !!')
      #print('sum images =',np.sum(maskgood),np.sum(maskgood*imagein))
      #if ( (np.sum(maskgood)>0) & (np.sum(maskgood*imagein) > 0) & (test_segmap_index == True) ) :

      if ( (np.sum(maskgood)>0) & (np.sum(maskgood*imagein) > 0) ) :
        # Colormap normalization according to image statistics (background and rms)
        f = plt.figure(figsize=(15,3))

        #ax1, ax2, ax3, ax4 = f.add_subplot(1, 4, sharey=True, sharex=True)
        ax1 = plt.subplot(151)
        ax2 = plt.subplot(152, sharex = ax1, sharey = ax1)
        ax3 = plt.subplot(153, sharex = ax1, sharey = ax1)
        ax4 = plt.subplot(154, sharex = ax1, sharey = ax1)
        ax5 = plt.subplot(155, sharex = ax1, sharey = ax1)

        upper_ct=30

        imagestat1=calc_bg(imagein)
        back1=imagestat1[0]   #'background'
        sigmaback1=imagestat1[1]
        norm1=matplotlib.colors.Normalize(vmin=back1-2*sigmaback1, vmax=back1+upper_ct*sigmaback1, clip=False)
        ax1.imshow(imagein, origin='lower', norm=norm1, cmap='hot')
        ax1.set_title('Original')
     
        ax2.imshow(segmap_final, origin='lower', cmap='hot')
        ax2.set_title('Segmentation map')
  
        imagestat3=calc_bg(imagein)
        back3=imagestat3[0]   #'background'
        sigmaback3=imagestat3[1]
        norm3=matplotlib.colors.Normalize(vmin=back3-2*sigmaback3, vmax=back3+upper_ct*sigmaback3, clip=False)
        ax3.imshow(maskgood*imagein, origin='lower', norm=norm1, cmap='hot')
        ax3.set_title(str(IDgal)+' band '+banda_obj)
  
        imagestat4=calc_bg(image_new)
        back4=imagestat4[0]   #'background'
        sigmaback4=imagestat3[1]
        norm4=matplotlib.colors.Normalize(vmin=back4-3*sigmaback4, vmax=back4+upper_ct*sigmaback4, clip=False)
        ax4.imshow(image_new, origin='lower', norm=norm4, cmap='hot')
        ax4.set_title('bkg subtracted')

        ax5.imshow(image_new_replaced, origin='lower', norm=norm4, cmap='hot')
        ax5.set_title('no secondary')

        #plt.title(str(IDgal))
        #plt.colorbar()
        plt.subplots_adjust(wspace=0.01, hspace=0.01)
        plt.tight_layout()
        if (test_segmap_index == True) : plt.show()
        #if (show_all==True) : plt.show()
        figcutouts=output_folder+'ID_'+str(IDgal)+'_band_'+banda_obj+'_ok.png'
        plt.savefig(figcutouts,dpi=80)
        #sys.exit()
        plt.close(f)

      if test_segmap_index == True :
        print('This is just for checking segmentation region of the galaxy')
        continue

      if test_segmap_index==True :
        print('Saved segmentation map. Now continue.')
        continue

    










      # -------------- --------------- -------------- --------------- -----
      



      # CALCULATE ALL ASYMMETRIES 
      print('\n\nCalculate All Asymmetries (and asymmetry center)')
      #print('Unless the shape-asymmetry, which will be calculated afterwards !!!')
    
      if (calculate_asymmetry_center==True) : # Questo calcola anche l'asymmetry !!!!
      # The center of rotation is not defined a priori, but is measured through an iterative process whereby the value of the asymmetry is calculated at the initial central guess (usually the geometric center or light centroid) and then the asymmetry is calculated around this central guess using some fraction of a pixel difference. This is repeated until a global minimum is found (Conselice et al. 2000a).
        
        try : 
          imageX1=imagein_replaced*1
          dim9=round(imageX1.shape[0]/2)
      
          range9=4
          centri=np.arange(dim9-range9,dim9+range9,1)
          asymmetry_array=[] ; ik_all=[] ; ip_all=[] ; counterA=0
          support_asymm=imageX1>-1e20
          support2_asymm=support_asymm*1
          # centri=centri9*1  # centri8 for imageX2, otherwise centri9 for image_new
          for ik in centri :
            for ip in centri :
              # print('center =',ik,ip)
              # Questa specie di maschera la metto (support2_asymm), pero' non serve a niente !!!
              #_asymmetry1,_denominatore=asymmetry_simple(imageX2,support2_asymm,ik,ip)
              #_asymmetry1,_denominatore=asymmetry_simple(image_new,support2_asymm,ik,ip)
              # _asymmetry1,_denominatore=asymmetry_simple(imagein,support2_asymm,ik,ip)
              _asymmetry1,_denominatore=asymmetry_simple(imagein_replaced,maskgood,ik,ip)
              #asymmetry_last,area_best,xc_best,yc_best,denominatore_best=asymmetry_function2(ik,ip,imageX2,maskgood_smaller,show_result)
              asymmetry_array.append(_asymmetry1)
              ik_all.append(ik)
              ip_all.append(ip)
              counterA+=1
              #print('asymmetry last (numeratore) =',asymmetry_last)
              #print('denominatore best =',denominatore_best)
          
          #counter_vector=np.arange(1,counterA+1,1)
          asymmetry_array=np.array(asymmetry_array)
          #print(asymmetry_array)
          #quit()
          # asymmetry_array=asymmetry_array[(asymmetry_array>-1) & (asymmetry_array<10)]
          minimum_asymmetry=min(asymmetry_array)
          #print('minimum asymmetry =',minimum_asymmetry)
          indexmin=np.where(asymmetry_array<minimum_asymmetry+0.00001)[0]
          if len(indexmin)>=1 : 
            indexmin=indexmin[0]
          
          center_ax=ik_all[indexmin] ; center_ay=ip_all[indexmin]
          sizeKA=20 ; Rmax_guess=10
          which_calc='asymmetry'
          #background_asymmetry,background_smoothness_notgood=asymmetry_bkg_simple(sizeKA,maskgood,image_new,Rmax_guess,which_calc,Rmax_guess)  
          background_asymmetry,min_bkg_asymmetry,background_smoothness_notgood=asymmetry_bkg_simple(sizeKA,binary_pawlik_4_complete,imagein_replaced,Rmax_guess,which_calc,Rmax_guess)  
          print('background asymmetry and minimum =',background_asymmetry,min_bkg_asymmetry)
          asymmetry_final=minimum_asymmetry-background_asymmetry
          print('asymmetry center =',center_ax,center_ay)
          print('Asymmetry, asymm background, and asymmetry final =')
          print(minimum_asymmetry,min_bkg_asymmetry,asymmetry_final)
          # print('xc_asymmetry from statmorph =',morphsingle.xc_asymmetry)
          # print('yc_asymmetry =', morphsingle.yc_asymmetry)
          #print('asymmetry statmorph =', morphsingle.asymmetry)
          #print('asymmetry statmorph =', morph.asymmetry)
          # print('concentration statmorph =', morphsingle.concentration)
          #plt.plot(counter_vector,asymmetry_array,color='cyan',lw=1.2)
          asymmetry_R=np.round(asymmetry_final,3)
          #asymmetry_R=np.round(minimum_asymmetry,3)  # For now I do not save the background subtracted value !!!!

        except :
          asymmetry_R=0
          center_ax=imagein.shape[0]/2 ; center_ay=imagein.shape[1]/2
      else :
        asymmetry_R=0
        center_ax=imagein.shape[0]/2 ; center_ay=imagein.shape[1]/2
      
      #plt.imshow(maskgood,origin='lower')
      #plt.show()
      print('asymmetry center =',center_ax,center_ay)
      #sys.exit()
      asymmetry_all1.append(np.round(asymmetry_R,3))
      centerax_all.append(center_ax)
      centeray_all.append(center_ay)






      # CALCULATE RMAX 

      if calculate_Rmax==True :
        try :
          firstsegmap_snr=2 ; firstsegmap_npixels=5   
          threshold1 = photutils.detect_threshold(imagein, nsigma=firstsegmap_snr)
      
          npixelsV = firstsegmap_npixels*1  # minimum number of connected pixels
          segm1 = photutils.detect_sources(imagein, threshold1, npixelsV)
          # Although statmorph is designed to process all the sources labeled by the segmentation map, in this example we only focus on the main (largest) source found in the image.
          areas1=segm1.areas ; #print segm.areas
          # Keep only the largest segment
          label1 = np.argmax(segm1.areas) + 1  # Prende la regione con area piu' grande
          segmap = segm1.data == label1
      
          # CALCULATE STATISTICS ON BACKGROUND ALL
          object1 = segm1.data > 0.5
          print('Created segmap, now convert bool to int')
          #plt.imshow(object1)
          #plt.show()
          #sys.exit()
          object1=1*object1
          object1ok=imagein*object1
          radius_vv=np.arange(1.5,round(image_new.shape[0]/2-2),0.2)
          # radius_vv=radius_vv[::-1]
          array_inverso=np.arange(len(radius_vv)-1)
          array_inverso=array_inverso[::-1]
          for tt in array_inverso :
            hD,wD= image_new.shape
            #print(center_ax,center_ay)
            ap_tot_inner = create_circular_mask(hD, wD, center=(center_ax,center_ay),radius=radius_vv[tt])
            #ap_tot_inner = ap_tot_inner*1
            ap_tot_outer = create_circular_mask(hD, wD, center=[center_ax,center_ay],radius=radius_vv[tt+1])
            #ap_tot_outer = ap_tot_outer*1
            ap_annulus=np.logical_xor(ap_tot_inner,ap_tot_outer)
            #plt.imshow(ap_annulus*1)
            #plt.show()
            #print(ap_annulus)
            #quit()

            # aper_ann=CircularAnnulus((center_ax,center_ay),radius_vv[tt]-1,radius_vv[tt])
            # # If you want a mask from that aperture do :
            # mask_ann = aper_ann.to_mask(method='center')
            # data_weighted= mask_ann.to_image(imagein.shape) # Vedi se questo funziona !! 
            # # data_weighted = mask_ann.multiply(np.ones(imagein.shape))
            # ann_table = aperture_photometry(object1, aper_ann)
            # sum_ann=ann_table['aperture_sum'].value 
            image_annulus=imagein_replaced*binary_pawlik_4*ap_annulus
  
            if sum(np.ravel(image_annulus))>0 : 
              Rmax=(radius_vv[tt+1]+radius_vv[tt])/2.
              print('Found Rmax =',np.round(Rmax,1))
              # time.sleep(1)
              break
    
          #plt.imshow(object1)
          #plt.show()
          if (np.array(Rmax).size==0) : 
            print('Problem -- Rmax not found')
            sys.exit()
  
        except :
          Rmax=round(imagein.shape[0]/2.-2)
          print('why not entered')
          #quit()
      else :
        Rmax=round(imagein.shape[0]/2-1)

      Rmax_all.append(np.round(Rmax,1))
      
      # ----------- Go back to full source image -----------
    
    
    
    
      
    
      # ---------------------------------------------------------------------------------
      # Calculate concentration manually
      if (calculate_concentration==True):

        # OK I RISULTATI SONO GLI STESSI, MA QUINDI la concentration si calcola sulla bkg subtracted o su quella originale ??????
        try :
          image_4_conc=image_new_replaced*1
          #image_4_conc=image_new *1   # imagein or image_new
          # print('This needs the calculation of Rmax in the previous step !!!!')
          print('Calculate concentration (Use photutils aperture photometry')
          from photutils.aperture import CircularAperture,aperture_photometry
          max_radius_source=Rmax*1
          raggi_vector=np.arange(1,max_radius_source+0.2,0.2)
          #print('Rmax =',Rmax)
          
          maskedimage=image_4_conc*maskgood
          center_ax4=center_ax*1 ; center_ay4=center_ay*1
          
          hD,wD= image_4_conc.shape
          ap_tot = create_circular_mask(hD, wD, center=[center_ax4,center_ay4],radius=max_radius_source)
          ap_tot = ap_tot*1
          somma_totale=sum(np.ravel(image_4_conc*ap_tot)) 
          ##imagein_masked = imagein_smoothed_skysub.copy()
          ##imagein_masked[~maskD] = 0
          
          #  ap_tot=CircularAperture((center_ax,center_ay), r= max_radius_source)
          #  phot_table_tot = aperture_photometry(image_4_conc, ap_tot)
          #  somma_totale=phot_table_tot['aperture_sum'].value
          
          # print('somma totale source (way 1) =',somma_totale)
          #   somma_totale2=sum(np.ravel(maskedimage))
          #   print('somma totale source (way 2) =',somma_totale2)
          #   print("somma_totale2 e' piu affidabile !! NOOO")
          # Devi usare l'immagine background subtracted
          #print('concentration from stamorph =',morph.concentration,morphsingle.concentration)
          #print("concentration da morph evidentemente e' sbagliata !!")
          
          somme_parziali=[]
          for kk in np.arange(len(raggi_vector)) :
            apertureX = create_circular_mask(hD, wD, center=[center_ax4,center_ay4],radius=raggi_vector[kk])
            apertureX = apertureX*1
            somma_X=sum(np.ravel(image_4_conc*apertureX)) 
            # apertureX = CircularAperture((center_ax,center_ay), r= raggi_vector[kk])
            # phot_tableX = aperture_photometry(image_4_conc, apertureX)
            # somma_X=phot_tableX['aperture_sum'].value
            somme_parziali.append(somma_X)

          somme_parziali=np.array(somme_parziali)
          #somme_parziali_fraction=somme_parziali/somma_totale2
          somme_parziali_fraction=somme_parziali/somma_totale
          #print(np.ravel(somme_parziali_fraction))
          index20=np.where(somme_parziali_fraction>=0.2)[0]
          index80=np.where(somme_parziali_fraction>=0.8)[0]
          radius_20=raggi_vector[index20[0]]
          radius_80=raggi_vector[index80[0]]
          print('radius 20 and 80 in pixels =',radius_20,radius_80)
          concentration=5*np.log10(radius_80/radius_20)
          print('Handmade concentration =',concentration)
          #print('Perfetto !!! Uguale a statmorph !!!')
          concentration_all_singlevalue=np.round(concentration,3)
        except:
          concentration_all_singlevalue=-9.
      
      else :
        concentration_all_singlevalue=-9.
      
      if np.sum(np.ravel(maskgood)) <= 0 :
        concentration_all_singlevalue=-9.
      #plt.imshow(maskgood,origin='lower')
      #plt.show()
      #plt.imshow(binary_pawlik_4,origin='lower')
      #plt.show()
      concentration_all.append(np.round(concentration_all_singlevalue,3))
    













      # ---------------------------------------------------------------------------------
      if calculate_smoothness==True :
        print('\n\nCalculate smoothness')
        #quit()

        smooY=1  # np.sum(np.ravel(maskgood))

        if smooY==1 :
        #try :
          print('\nCalculate smoothness : --------------- --------------- ------------ -----------\n')
          print('image shape =',image_new_replaced.shape)
          #quit()
          #print('Pixel scale =',IDpixscale[iyy])
          #print('FWHM band =',IDfwhm[iyy])
          # Lo smoothradius e' lo stesso che uso dopo per la clumpiness !!!
          #smoothradiusX22=round(1*conv_kpc_to_arcsec/IDpixscale_banda) # 1 kpc is the size of the clumps we want to detect
          # Se per la smoothness usassi Rmax come smoothing_length anziche' Rmax/4 verrebbe un valore piu' alto
          #smoothradiusX22=Rmax/4.  # Ti ricordo che per la clumpiness handmade io uso Rmax !!! Mentre in Lotz 2004 e' scritto Rmax/4

          # Questo smoothradiusX22 verra' usato tale e quale come raggio per smussare l'immagine per il calcolo della clumpiness !!!!

          factor_smooth22=1

          #plt.imshow(image_new_replaced,origin='lower')
          #plt.show()
          
          print('smoothradius =',smoothradiusX22)
          print('smoothtype =',smoothtype)
          #quit()
          imagesmooth22=smooth(image_new_replaced,smoothtype,int(smoothradiusX22), factor_smooth22)  # Il 15 alla fine che vuol dire ??
      
          # Calculate the residual image 
          residuals_orig22=image_new_replaced-imagesmooth22  # image_4_clumpiness e' immagine background subtracted
          residuals22=residuals_orig22*1
          residuals22[residuals22 < 0] = 0
          
          # Questo perche' in Lotz+2004 c'e' scritto che per la computation della smoothness devi escludere le regioni centrali, entro 0.25*Rpetrosian
          
          #   apertureX22 = CircularAperture((center_ax,center_ay), r= smoothradiusX22)
          #   aper2222= apertureX22.to_mask(method='center')
          #   data_weighted22 = aper2222.to_image(image_new.shape)
          #   data_weighted22 = np.absolute(data_weighted22-1) # Prendo l'opposto della maschera
          #   # final_mask=binary_pawlik_4*data_weighted22 # Lotz+2004
          #   final_mask=binary_pawlik_4*1  # per momento facciamo cosi, mantieni tutto, anche il centro
          
          hD,wD=image_new_replaced.shape
          mask504=create_circular_mask(hD, wD, center=(center_ax,center_ay),radius=smoothradiusX22)
          mask505=~mask504
          #plt.imshow(binary_pawlik_4)
          #plt.show()
          #plt.imshow(final_mask)
      
          residuals22_1d=np.absolute(np.ravel(residuals22*maskgood))
          denomin_22=np.absolute(np.ravel(image_new_replaced*maskgood))

          if (  (sum(denomin_22)>0) & (sum(denomin_22)<99999) ) :
            smoothness22_source=sum(residuals22_1d)/sum(denomin_22)
            print('sum numeratore e denominatore =',sum(residuals22_1d),sum(denomin_22))
            print('Smoothness source =',smoothness22_source)
            #time.sleep(0.6)
          
            sizeKA33=skybox_arcsec/IDpixscale_banda
            which_calc='smoothness'
            background_asymmetry_notgood,min_as_back,background_smoothness_good=asymmetry_bkg_simple(sizeKA33,maskgood,imagein_replaced,Rmax,which_calc,smoothradiusX22)  
            print('smoothness background =',background_smoothness_good)
        
            smoothness22=smoothness22_source-background_smoothness_good   # For now it is not anymore background subtracted !!!
            print('Smoothness =',smoothness22)
            #time.sleep(0.1)
            #print('\n')
          else :
            smoothness22_source = -9
            smoothness22= -9.



          if show_images_masks_smoothness==True :
            try : 
              # VISUALIZE ORIGINAL IMAGE, CLUMPY PIXELS and deblended segmenttion map on detected clumps
              # Colormap normalization according to image statistics (background and rms)
              fi9, ((ax1, ax2, ax4), (ax7, ax5, ax6)) = plt.subplots(2, 3, sharey=True, sharex=True)
              
              upper_ct=10

              imagestat1=calc_bg(imagein[maskgood==1])
              back1=imagestat1[0]   #'background'
              sigmaback1=imagestat1[1]
              norm1=matplotlib.colors.Normalize(vmin=back1-2*sigmaback1, vmax=back1+upper_ct*sigmaback1, clip=False)
              
              imagestat3=calc_bg(image_new[maskgood==0])
              #imagestat3=calc_bg(np.ravel(image_4_clumpiness))
              back3=imagestat3[0]   #'background'
              sigmaback3=imagestat3[1]
              norm3=matplotlib.colors.Normalize(vmin=back3-2*sigmaback3, vmax=back3+upper_ct*sigmaback3, clip=False)
              
              #print(residuals22.shape)
              #imagestat4=calc_bg(residuals22[maskgood==0])
              #imagestat3=calc_bg(np.ravel(image_4_clumpiness))
              back4=np.median(residuals22_1d) # imagestat4[0]   #'background'
              sigmaback4=np.std(residuals22_1d) # imagestat4[1]
              norm4=matplotlib.colors.Normalize(vmin=back4-2*sigmaback4, vmax=back4+upper_ct*sigmaback4, clip=False)
              
    
              ax1.imshow(imagein, origin='lower', norm=norm1, cmap='hot')
              ax1.set_title('ID '+str(IDgal)+' (original)')
              
              ax2.imshow(image_new, origin='lower', cmap='hot', norm=norm3)
              ax2.set_title('bkg subtracted')
              
              # Immagine originale con maschera
              #ax3.imshow(maskgood*image_4_clumpiness, origin='lower', norm=norm1, cmap='hot')
              #ax3.set_title(str(IDgal))
              #plt.title(str(IDgal))
              #plt.colorbar()
              
              # Original - smoothed (= Residual)
              #imagestat4=calc_bg(residuals_orig)
              #back4=imagestat4[0]   #'background'
              #sigmaback4=imagestat4[1]
              #norm4=matplotlib.colors.Normalize(vmin=back4-2*sigmaback4, vmax=back4+5*sigmaback4, clip=False)
              ax4.imshow(image_new_replaced, origin='lower', norm=norm3, cmap='hot')
              ax4.set_title('replaced')
    
              # Segmentation map (Deblended) on clumps
              ax7.imshow(residuals22, origin='lower', cmap='hot', norm=norm4)
              ax7.set_title('residuals') 
              
              # Segmentation map (Deblended) on clumps
              #image_neg=image_4_clumpiness*1 
              #image_neg[pixelmap>0]=-10
              
              #CLUMPS_debl_original_visual=CLUMPS_debl_original*1
              #CLUMPS_debl_original_visual[CLUMPS_debl_original_visual==0]=-10
              ax5.imshow(image_new_replaced*maskgood, origin='lower', cmap='hot', norm=norm3) # Questo ti serve poi per rimuovere i nuclei
              # #ax6.imshow(segnuclei_regul*maskgood, origin='lower', cmap='hot')
              ax5.set_title('with maskgood')
              
              # Immagine originale senza clumps
              #ax5.imshow(pixelmap_neg*image_4_clumpiness, origin='lower', norm=norm3, cmap='hot')
              #image_neg2=pixelmap_neg*1 
              #image_neg2[pixelmap>0]=-10
              ax6.imshow(image_new_replaced*binary_pawlik_4, origin='lower', norm=norm3, cmap='hot')
              ax6.set_title('with binary_pawlik_4')
              
              plt.subplots_adjust(wspace=0.2, hspace=0.01)
              #plt.tight_layout()
              #if (show_all==True) :
              namefile=output_folder+str(IDlista[iyy])+'_smoothness_'+banda_obj+'.png'
              plt.savefig(namefile,dpi=80)

              plt.close(fi9)
            except :
              fig9, ((ax1, ax2)) = plt.subplots(1, 2)
              namefile=output_folder+str(IDlista[iyy])+'_smoothness_'+banda_obj+'.png'
              plt.savefig(namefile,dpi=80)
              plt.close(fig9)
  
        #except : smoothness22= 0
        else : 
          smoothness22= -9

      else :
        smoothness22= -9.
      
      if np.sum(np.ravel(maskgood)) <= 0 :
        smoothness22= -9.

      #print(np.round(smoothness22,3))
      #quit()
      smoothness_all.append(np.round(smoothness22,3))



















      # -------------- --------------- -------------- --------------- -----------------------
      
      if (calculate_gini==True) :
        
        gini1=1

        if gini1==1 :
        #try : 
          print('\n\nCalculate Gini parameter ----------- ')
          # FOLLOWING LOTZ ET AL. 2004
          
          # STEP 0) LOOK FOR THE PETROSIAN RADIUS
          # petrosian_radius=petrosian_radius(image_new,centroid_single)
  
          # -------------------------------------------------------------------------------
          #print('Calculate Gini manually :')
          #good=singlegood.ravel()
          
          segm_gal=image_new_masked[maskgood==1]
          image_new_1d=segm_gal.ravel()  # imagein o image_new ci devono essere (sono rispettivamente non e background subtracted)
          #print good[0]
          #image_new_1d=image_new_1d[good>0]
          image_new_1d=np.abs(image_new_1d)
          sorted_pixelvals1=np.sort(image_new_1d)
          #plt.imshow(segm_gal,origin='lower')
          #plt.show()
          # image_gini2=image_new[maskgood>0]
          # #image_gini=np.abs(image_gini)
          # #image_gini=np.abs(image_new2_reduced)*1
          # image_gini2=image_gini2.reshape(-1)
          # print('len image gini =',len(image_gini2))
          # sorted_pixelvals2 = np.sort(image_gini2)
          # -------------------------------------------------------------------------
          #sorted_pixelvals=sorted_pixelvals2*1  # The default was sorted_pixelvals1 , don't know why exactly !
          #sorted_pixelvals = np.sort(np.abs(image_new2_reduced))
          #print 'sorted_pixelvals =',sorted_pixelvals
          #print('len sorted pixvals =',sorted_pixelvals)
          n = len(sorted_pixelvals1)
          #print('n =',n)
          if n <= 1 or np.sum(sorted_pixelvals1) == 0 :
            print('[gini] Not enough data for Gini calculation.')
          indices = np.arange(1, n+1,1) # start at i=1
          #print('len indices =',len(indices),min(indices),max(indices))
          #print(indices[0:20])
          #print(sorted_pixelvals1[0:20])
          if (  (np.sum(sorted_pixelvals1)!=0) & (n>1) ) :
            gini = (np.sum( (2*indices-n-1) * sorted_pixelvals1)) / ( float(n)*float(n-1) * np.average(sorted_pixelvals1) )
          else : gini=-9 
          if show_figures_Gini==True :
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 6))
            ax1.imshow(image_new, origin='lower', interpolation='nearest',cmap='hot')
            ax1.set_title('Image original')
            ax2.imshow(image_new_conv, origin='lower', interpolation='nearest',cmap='hot')
            ax2.set_title('Image original conv')
            ax3.set_title('Segmentation galaxy\nfor GINI')
            ax3.imshow(image_new_masked, origin='lower', interpolation='nearest',cmap='hot')
            # plt.show()
            namefileF=output_folder+str(IDlista[iyy])+'_ginimap_distclumps_'+str(IDdist_clumps[iyy])+'_mag_'+str(IDmag[iyy])+'.png'
            plt.savefig(namefileF,dpi=80)
        #except :
        else :
          gini=-9.
      else :
        gini=-9.
      gini_all.append(np.round(gini,3))
      print('Gini (hand-made) =',np.round(gini,3)) # ,'\n') # ,morph.gini,morphsingle.gini)
  










      # -------------- --------------- -------------- --------------- -----------------------

      # Calculate M20 manually  (Calculate the M_20 coefficient as described in Lotz et al. (2004))
      '''
      if (calculate_m20_oldstyle==True):
        print('NOT USED THIS VERSION')
        # Commento se mi interessa solo la clumpiness, pero' funziona !!!!
        #m20_all[uu,rr]=0
        try :
          imageX1=image_new_replaced*maskgood
          #imageX1=image_new*maskgood     # or image_new2_reduced
          #imageX=image_new2_reduced*1
          dim9=int(imageX1.shape[0]/2)
          row_cc=imageX1[dim9,:]
          column_cc=imageX1[:,dim9]
          row_start=np.where(row_cc>0)[0]
          #print(row_start)
          column_start=np.where(column_cc>0)[0]
          #print('column_start =',column_start)
          if segmentation_type=='square' :
            imageX2=image_new[column_start[0]:column_start[-1]+1,row_start[0]:row_start[-1]+1]
          # fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 6))
          # ax1.imshow(image_new, origin='lower', interpolation='nearest',cmap='hot')
          # ax1.set_title('Image original')
          # ax2.imshow(imageX1, origin='lower', interpolation='nearest',cmap='hot')
          # ax2.set_title('Image masked')
          # ax3.imshow(imageX2, origin='lower', interpolation='nearest',cmap='hot')
          # ax3.set_title('Image cropped')
          # plt.show()
          #namefileF=output_folder+str(IDlista[iyy])+'_ginimap_distclumps_'+str(IDdist_clumps[iyy])+'_mag_'+str(IDmag[iyy])+'.png'
          #plt.savefig(namefileF,dpi=80)
          # if show_figures==True : plt.show()
          imageX=imageX2*1   # imageX1 or imageX2
          M = skimage.measure.moments(imageX, order=1) 
          yc = M[1, 0] / M[0, 0]
          xc = M[0, 1] / M[0, 0]
          # Note that we do not shift (yc, xc) by 0.5 pixels here, since # (yc, xc) is only used as input for other skimage functions.
          # Calculate second total central moment
          Mc = skimage.measure.moments_central(imageX, center=(yc, xc), order=2)
          second_moment_tot = Mc[0, 2] + Mc[2, 0]
          # Calculate threshold pixel value
          sorted_pixelvals = np.sort(imageX.flatten())
          flux_fraction = np.cumsum(sorted_pixelvals) / np.sum(sorted_pixelvals) 
          sorted_pixelvals_20 = sorted_pixelvals[flux_fraction >= 0.8]
          if len(sorted_pixelvals_20) == 0:
            # This can happen when there are very few pixels.
            print('Not enough data for M20 calculation')
          thresholdm20 = sorted_pixelvals_20[0]
          image_20 = np.where(imageX >= thresholdm20, imageX, 0.0)
          Mc_20 = skimage.measure.moments_central(image_20, center=(yc, xc), order=2)
          second_moment_20 = Mc_20[0, 2] + Mc_20[2, 0]
          if (second_moment_20 <= 0) | (second_moment_tot <= 0): 
            print('[m20] Negative second moment(s)')
          m20 = np.log10(second_moment_20 / second_moment_tot)
          print('m20 (by me and by statmorph) =',m20) #,morph.m20
          m20_all.append(np.round(m20,3))
        except :
          m20_all.append(-9)
      else :
        m20_all.append(-9)
        print('continue') # m20_all[uu,rr]=0
      '''   
      
  





      # ------------------------------------------------------------------------------
      if (calculate_m20==True):  # See Rodriguez-Gomez +2019
        print('\n\nCalculate M20')
        #m20_all[uu,rr]=0
        #'
        m20_1=np.sum(np.ravel(maskgood))
        #try :


        if m20_1>0 :
        #kkll=2
        #if kkll==2 :
          imageX1=imagein_replaced*maskgood # maskgood*imagein_replaced # image_new_replaced*maskgood
          #plt.imshow(imageX1,origin='lower')
          #plt.show()
          
          passo=10
          centersF_x=np.arange(round(center_ax)-passo,round(center_ax)+passo,0.5) 
          centersF_y=np.arange(round(center_ay)-passo,round(center_ay)+passo,0.5) 
          muTOT_final_value,x_center_muTOT,y_center_muTOT=minimize_mu(imageX1,centersF_x,centersF_y)
          mu20=M20_simple(imageX1,x_center_muTOT,y_center_muTOT)
          print('Mu_tot and mu20 and centers =',muTOT_final_value,mu20,x_center_muTOT,y_center_muTOT)
          m20a=mu20/muTOT_final_value ; 
          m20=np.log10(m20a)
          m20=np.round(m20,3)
          print('M20 final (handmade) =',m20)
          print('\n\n')
          #quit()
        #except :
        else :
          m20=-9.
      else :
        m20=-9.
      
      #if np.sum(np.ravel(maskgood))<= 0 : m20=-9.
      
      #quit()
      

      m20_all.append(m20)
      #print('M20 final =',m20_all)
      #time.sleep(0.1)
      
      #if counter989==20 :
      #  quit()
      #print(m20)
      #quit()
      print('Gini and M20 =',gini,m20)
      time.sleep(0.7)
         

      # *******************************************************************

      final_mask_galaxy_segm=maskgood*1
      mask_bkg_segmentation=final_mask_galaxy_segm==0
      mask_bkg_segmentation=mask_bkg_segmentation*1
      background_new=image_new[mask_bkg_segmentation==1]
      background_new_1d=np.ravel(background_new)
      #print('len background =',len(background_new_1d))
      #mask_galaxy_segmentation=mask_galaxy_segmentation*image_new
      #plt.imshow(image_new*mask_bkg_segmentation)
      #plt.show()
      medbackBK_new=np.median(background_new_1d)
      sigmabackBK_new=np.std(background_new_1d)
      image_background=image_new_masked*1 #*maskgood
      noisy_new=np.zeros(image_background.shape)
      noisy_new2=image_new*1

      for i1 in np.arange(image_new.shape[0]) :
        for j1 in np.arange(image_new.shape[1]) :
          if (image_background[i1,j1]==0) :
            #print('Go on')
            image_background[i1,j1]=random.gauss(medbackBK_new,sigmabackBK_new)
            #noisy_new2[i1,j1]=np.abs(noisy_new2[i1,j1])
          else :
            noisy_new2[i1,j1]=np.abs(random.gauss(medbackBK_new,sigmabackBK_new))
          noisy_new[i1,j1]=random.gauss(medbackBK_new,sigmabackBK_new)
      #maskones=np.ones(image_new.shape)
      #maskones[100:400,100:400]=1  # Puoi metterla o non metterla
      #plt.imshow(noisy_new2)
      #plt.show()

      # *******************************************************************


















      # --------------------------------------------------------------------     

      # STATMORPH ALL 

      if use_statmorph==True : 
        print('\n\nNOW GOING TO USE STATMORPH TO CALCULATE MORPHOLOGICAL PARAMETERS')
        
        print('iyy and jkk =',iyy,' ',jkk)
        print('        ')
        print('********')
        # Funziona, ma commento perche' non serve piu' : 
        size = size_psf*1   # on each side from the center
        sigma_psf = IDfwhm_banda[jkk]/IDpixscale_banda
        #  print('sigma_psf =',sigma_psf)
        y, x = np.mgrid[-size:size+1, -size:size+1]
        psf = np.exp(-(x**2 + y**2)/(2.0*sigma_psf**2))
        psf /= np.sum(psf)
        
        start = time.time()
  
        mask_ones=image_background>=-1e10
        mask_ones=mask_ones*1
  
        #maskgood_reverse=np.abs(maskgood-1)
        # print(maskgood)
        # print(mask_ones)
        mask_galaxy_segmentation=final_mask_galaxy_segm==1
        mask_galaxy_segmentation=mask_galaxy_segmentation*1
  
        #print(np.ravel(mask_galaxy_segmentation)[0:50])
        #plt.imshow(image_new)
        #plt.show()
        #plt.imshow(mask_galaxy_segmentation)
        #plt.show()
        #plt.imshow(noisy_new2)
        #plt.show()
        #plt.imshow(psf)
        #plt.show()
        
        #   print(min(np.ravel(image_new)),max(np.ravel(image_new)))
        #   print(min(np.ravel(mask_galaxy_segmentation)),max(np.ravel(mask_galaxy_segmentation)))
        #   print(min(np.ravel(noisy_new2)),max(np.ravel(noisy_new2)))
        #   #quit()
        #   print(min(np.ravel(imagein)),max(np.ravel(imagein)))
        #   print(min(np.ravel(noisy_img)),max(np.ravel(noisy_img)))
        #  weightmap=noisy_new2,
        #continueX=True
        #if continueX==True 
        #usa_x=1 
        #if usa_x==1 :
        #try :

        '''
        noisy_img2=noisy_img*1
        absnoise=np.abs(noisy_img2)
        validpixels= (absnoise>1e-10) # |  (imagein_replaced2<1)
        #print(nonvalid*1)
        validmask=validpixels*1
        nonvalidpixels=np.abs(validmask-1)
        #print(np.where(nonvalid==1)[0])
        noisy_img2[nonvalidpixels]=1.

        # -----------------------------------------------
        imagein_replaced2=imagein_replaced*1
        absimagein=np.abs(imagein_replaced2)
        validpixels= (absimagein>1e-8) # |  (imagein_replaced2<1)
        #print(nonvalid*1)
        validmask=validpixels*1
        nonvalidpixels=np.abs(validmask-1)
        #print(np.where(nonvalid==1)[0])
        imagein_replaced2[nonvalidpixels]=1e-6

        dove=noisy_img2[noisy_img2==0]
        dove3=imagein_replaced2[imagein_replaced2==0]
        fff= (imagein_replaced2>-1) & (imagein_replaced2<1)
        good_fff=fff*1
        #print('How many pixels good =',sum(np.ravel(good_fff)))

        #quit()
        print(dove3)

        imagein_replaced2=imagein_replaced2*1e5

        plt.imshow(imagein_replaced2,origin='lower')
        plt.show()
        source_morphs = statmorph.source_morphology(imagein_replaced2, mask_galaxy_segmentation, psf=psf,skybox_size=round(skybox_arcsec/IDpixscale_banda),gain=1 ) #,weightmap=noisy_img2,gain=1)
        morph = source_morphs[0]  # Scelgo il segmentation part numero 0
        #print('xc_asymmetry =', morph.xc_asymmetry)
        #quit()
        print('xc_asymmetry =', morph.xc_asymmetry)
        print('yc_asymmetry =', morph.yc_asymmetry)
        print('Gini =', morph.gini)
        print('M20 =', morph.m20)
        print('sn_per_pixel =', morph.sn_per_pixel)
        print('C =', morph.concentration)
        #print('concentration hand-made =', concentration)
        print('A =', morph.asymmetry)
        print('S =', morph.smoothness)
        print('flag [0=good, 1=bad] =', morph.flag)
        gini_mm=np.round(morph.gini,3)
        flag_mm=morph.flag
        m20_mm=np.round(morph.m20,3) 
        concentration_mm=np.round(morph.concentration,3)
        smoothness_mm=np.round(morph.smoothness,3)
        AAA=morph.asymmetry
        #asymmetry_mm=np.round(morph.asymmetry,3)
        gini_SM.append(gini_mm)
        flag_SM.append(flag_mm)
        concentration_SM.append(concentration_mm)
        smoothness_SM.append(smoothness_mm)
        ellipticity_SM.append(-9)
        elongation_SM.append(-9)
        #rpetro_SM.append(np.round(morph.rpetro_circ,3))
        #rhalf_SM.append(np.round(morph.rhalf_ellip,3))
        m20_SM.append(m20_mm)
        asymmetry_SM.append(np.round(AAA,3))
        rpetro_SM.append(np.round(morph.rpetro_circ,3))
        rhalf_SM.append(np.round(morph.rhalf_ellip,3))
        quit()
        '''
        

        ######
        do_SM=1
        if ( (iyy==2) & (jkk==1) ) :
          do_SM=0

        if do_SM==1 :
        #result = None
        #while result is None:


          try:
            noisy_img2=noisy_img*1
            absnoise=np.abs(noisy_img2)
            validpixels= (absnoise>1e-10) # |  (imagein_replaced2<1)
            #print(nonvalid*1)
            validmask=validpixels*1
            nonvalidpixels=np.abs(validmask-1)
            #print(np.where(nonvalid==1)[0])
            noisy_img2[nonvalidpixels]=1
            imagein_replaced2=imagein_replaced*1

            absimagein=np.abs(imagein_replaced2)
            validpixels= (absimagein>1e-10) # |  (imagein_replaced2<1)
            #print(nonvalid*1)
            validmask=validpixels*1
            nonvalidpixels=np.abs(validmask-1)
            #print(np.where(nonvalid==1)[0])
            imagein_replaced2[nonvalidpixels]=1e-6
            #plt.imshow(psf,origin='lower')
            #plt.show()
            print('\n\n')
            

            imagein_replaced2=imagein_replaced2*1e4
            noisy_img2=noisy_img2*1e6

            source_morphs = statmorph.source_morphology(imagein_replaced2, mask_galaxy_segmentation, psf=psf,skybox_size=round(skybox_arcsec/IDpixscale_banda),gain=1) #weightmap=noisy_img2)

            # source_morphs = statmorph.source_morphology(image_new, mask_galaxy_segmentation, psf=psf,skybox_size=skybox_pixel,weightmap=noisy_new2)
            # source_morphs = statmorph.source_morphology(imagein_replaced2, mask_galaxy_segmentation, psf=psf,skybox_size=round(skybox_arcsec/IDpixscale_banda),gain=1) #      ,weightmap=noisy_img2)
            # source_morphs = statmorph.source_morphology(image_new, mask_galaxy_segmentation, psf=psf, weightmap=noisy_new2,skybox_size=15)
            # source_morphs = statmorph.source_morphology(image_background, mask_galaxy_segmentation, psf=psf, weightmap=noisy_new,skybox_size=25)
            # source_morphs = statmorph.source_morphology(image_background, mask_ones, psf=psf, weightmap=noisy_new,skybox_size=25)
            # print('Time: %g s.' % (time.time() - start))
            #source_morphs = statmorph.source_morphology(image_background, maskones, weightmap=noisy_img, skybox_size=1)
            morph = source_morphs[0]  # Scelgo il segmentation part numero 0
            #print(source_morphs[1])
            #sys.exit()
            # from statmorph.utils.image_diagnostics import make_figure
            print('xc_asymmetry =', morph.xc_asymmetry)
            print('yc_asymmetry =', morph.yc_asymmetry)
            print('Gini =', morph.gini)
            print('M20 =', morph.m20)
            print('sn_per_pixel =', morph.sn_per_pixel)
            print('C =', morph.concentration)
            #print('concentration hand-made =', concentration)
            print('A =', morph.asymmetry)
            print('S =', morph.smoothness)
            print('flag [0=good, 1=bad] =', morph.flag)
            gini_mm=np.round(morph.gini,3)
            flag_mm=morph.flag
            m20_mm=np.round(morph.m20,3) 
            concentration_mm=np.round(morph.concentration,3)
            smoothness_mm=np.round(morph.smoothness,3)
            AAA=morph.asymmetry
            #asymmetry_mm=np.round(morph.asymmetry,3)
            gini_SM.append(gini_mm)
            flag_SM.append(flag_mm)
            concentration_SM.append(concentration_mm)
            smoothness_SM.append(smoothness_mm)
            ellipticity_SM.append(-9)
            elongation_SM.append(-9)
            sn_per_pixel_SM.append(np.round(morph.sn_per_pixel,3))
            #rpetro_SM.append(np.round(morph.rpetro_circ,3))
            #rhalf_SM.append(np.round(morph.rhalf_ellip,3))
            m20_SM.append(m20_mm)
            asymmetry_SM.append(np.round(AAA,3))
            rpetro_SM.append(np.round(morph.rpetro_circ,3))
            rhalf_SM.append(np.round(morph.rhalf_ellip,3))
            #gini_SM.append(np.round(morph.gini,3))
            #flag_SM.append(morph.flag)
            ##sn_per_pixel_SM.append(np.round(morph.sn_per_pixel,3))
            #asymmetry_SM.append(np.round(morph.asymmetry,3))
            #concentration_SM.append(np.round(morph.concentration,3))
            ##print(np.round(morph.smoothness,3))
            ##time.sleep(2)
            #smoothness_SM.append(np.round(morph.smoothness,3))
            #m20_SM.append(np.round(morph.m20,3))
            #ellipticity_SM.append(-9)
            #elongation_SM.append(-9)
            #rpetro_SM.append(np.round(morph.rpetro_circ,3))
            #rhalf_SM.append(np.round(morph.rhalf_ellip,3))
            # try :
            #   fig = make_figure(morph)
            #   fig.savefig('tutorial.png', dpi=150)
            #   plt.close(fig)
            '''
            print('xc_centroid =', morph.xc_centroid)
            print('yc_centroid =', morph.yc_centroid)
            #print('ellipticity =', morph.ellipticity_centroid)
            #print('elongation =', morph.elongation_centroid)
            #print('orientation =', morph.orientation_centroid)
            print('rpetro_circ =', morph.rpetro_circ)
            print('rpetro_ellip =', morph.rpetro_ellip)
            #print('rhalf_circ =', morph.rhalf_circ)
            #print('rhalf_ellip =', morph.rhalf_ellip)
            print('xc_asymmetry =', morph.xc_asymmetry)
            print('yc_asymmetry =', morph.yc_asymmetry)
            #print('r20 =', morph.r20)
            #print('r80 =', morph.r80)
            print('Gini =', morph.gini)
            print('M20 =', morph.m20)
            print('F(G, M20) =', morph.gini_m20_bulge)
            print('S(G, M20) =', morph.gini_m20_merger)
            print('sn_per_pixel =', morph.sn_per_pixel)
            print('C =', morph.concentration)
            print('A =', morph.asymmetry)
            print('S =', morph.smoothness)
            #print('sersic_amplitude =', morph.sersic_amplitude)
            #print('sersic_rhalf =', morph.sersic_rhalf)
            #print('sersic_n =', morph.sersic_n)
            #print('sersic_xc =', morph.sersic_xc)
            #print('sersic_yc =', morph.sersic_yc)
            #print('sersic_ellip =', morph.sersic_ellip)
            #print('sersic_theta =', morph.sersic_theta)
            print('sky_mean =', morph.sky_mean)
            print('sky_median =', morph.sky_median)
            print('sky_sigma =', morph.sky_sigma)
            print('flag [0=good, 1=bad] =', morph.flag)
            #sys.exit()
            '''
            # Use STATMORPH on single galaxy images :
            print('Use statmorph on single galaxy image (masked)')
            image_background2=image_new*maskgood # sicuro che e' giusto ?
            for i1 in np.arange(imagein.shape[0]) :
              for j1 in np.arange(imagein.shape[1]) :
                if (image_background2[i1,j1]==0) :
                  image_background2[i1,j1]=random.gauss(0,sigmabackBK)
            maskones=np.zeros(imagein.shape)
            maskones[maskgood==1]=1  # Puoi metterla o non metterla
            mask_ok=maskgood==1.
            mask_ok=mask_ok*1
            
            # connect

            result = AAA*1
            print('result =',result)
            time.sleep(0.7)
          except:

            maskones=np.zeros(imagein.shape)
            mask_ok=maskones*1
            gini_SM.append(-9)
            flag_SM.append(-9)
            sn_per_pixel_SM.append(-9)
            asymmetry_SM.append(-9)
            concentration_SM.append(-9)
            smoothness_SM.append(-9)
            m20_SM.append(-9)
            ellipticity_SM.append(-9)
            elongation_SM.append(-9)
            rpetro_SM.append(-9)
            rhalf_SM.append(-9)

            #pass


          #plt.imshow(image_background2)
          #plt.show()
          # PREVIOUS VERSION : (DAVA PROBLEMI SE METTO IMMAGINE BKG SUBTRACTED)
          # source_morphs_single = statmorph.source_morphology(image_background2, mask_ok, weightmap=noisy_img,gain=2, psf=psf,skybox_size=skybox_pixel)
          #mask_reverse=np.abs(singlegood-1)
          #source_morphs_single = statmorph.source_morphology(imagein, maskones,weightmap=noisy_img) # , skybox_size=1)
          '''
          print(min(np.ravel(image_new)),max(np.ravel(image_new)))
          print(min(np.ravel(mask_galaxy_segmentation)),max(np.ravel(mask_galaxy_segmentation)))
          print(min(np.ravel(noisy_new2)),max(np.ravel(noisy_new2)))
          #quit()
          print('\n')
          print(min(np.ravel(imagein)),max(np.ravel(imagein)))
          print(min(np.ravel(noisy_img)),max(np.ravel(noisy_img)))
          #  weightmap=noisy_new2,
          '''
          # source_morphs = statmorph.source_morphology(image_new, mask_galaxy_segmentation, psf=psf,skybox_size=skybox_pixel,weightmap=noisy_new2)
          
          '''
          source_morphs_single = statmorph.source_morphology(imagein, mask_ok, psf=psf,skybox_size=round(skybox_arcsec/IDpixscale_banda),weightmap=noisy_img)
          #source_morphs_single = statmorph.source_morphology(imagein, maskgood_debl_seg, gain=0.5)
          #source_morphs_single = statmorph.source_morphology(imagein, singlegood, gain=0.5)
          morphsingle = source_morphs_single[0]  # Scelgo il segmentation part numero 0
          #print 'Now we print some of the morphological properties just calculated\n'
          #print('ellipticity =', morphsingle.ellipticity_centroid)
          #print('elongation =', morphsingle.elongation_centroid)
          #print('rpetro_ellip =', morphsingle.rpetro_ellip)
          #print('rhalf_ellip =', morphsingle.rhalf_ellip)
          print('xc_centroid =', morphsingle.xc_centroid)
          print('yc_centroid =', morphsingle.yc_centroid)
          print('xc_asymmetry =', morphsingle.xc_asymmetry)
          print('yc_asymmetry =', morphsingle.yc_asymmetry)
          print('gini =', morphsingle.gini)
          print('M20 =', morphsingle.m20)
          print('sn_per_pixel =', morphsingle.sn_per_pixel)
          #print 'elongation =',morphsingle.elongation_centroid
          #print('r half =', morphsingle.rhalf_ellip)
          print('r half circ =', morphsingle.rhalf_circ)
          print('r petro =', morphsingle.rpetro_ellip)
          print('r petro circ =', morphsingle.rpetro_circ)
          print('C =', morphsingle.concentration)
          print('A =', morphsingle.asymmetry) # That was calculated on full source, not just on single galaxy !!!!!
          print('S =', morphsingle.smoothness)
          print('smoothness hand-made =',smoothness22)
          print('sky_mean =', morphsingle.sky_mean)
          print('sky_median =', morphsingle.sky_median)
          print('sky_sigma =', morphsingle.sky_sigma)
          print('flag [0=good, 1=bad]=', morphsingle.flag)
          '''
          #ellipticity_all[uu,rr]=np.round(morphsingle.ellipticity_centroid)
          #elongation_all[uu,rr]=np.round(morphsingle.elongation_centroid,2)
          #rpetro_all[uu,rr]=np.round(morphsingle.rpetro_ellip,2)
          #rhalf_all[uu,rr]=np.round(morphsingle.rhalf_ellip,2)
          #sn_per_pixel[uu,rr]=round(morphsingle.sn_per_pixel)
          #centroid_single=(morphsingle.xc_centroid,morphsingle.yc_centroid)
          #petrosian_radius_single=np.round(morphsingle.rpetro_circ,1)
          

          # gini_SM.append(np.round(morph.gini,3))
          # flag_SM.append(morph.flag)
          # #sn_per_pixel_SM.append(np.round(morph.sn_per_pixel,3))
          # asymmetry_SM.append(np.round(morph.asymmetry,3))
          # concentration_SM.append(np.round(morph.concentration,3))
          # #print(np.round(morph.smoothness,3))
          # #time.sleep(2)
          # smoothness_SM.append(np.round(morph.smoothness,3))
          # m20_SM.append(np.round(morph.m20,3))
          # ellipticity_SM.append(-9)
          # elongation_SM.append(-9)
          # sys.exit()
          # rpetro_SM.append(np.round(morph.rpetro_circ,3))
          # rhalf_SM.append(np.round(morph.rhalf_ellip,3))
          '''
          print('Va bene, poi li confronto con il mio metodo e vediamo quali sono giusti\n\n')
          #sys.exit()
          print('\ncompare ginis. Mine, morph and morphsingle =',)
          print(gini,np.round(morph.gini,3),np.round(morphsingle.gini,3))
          print('compare M20s. Mine, morph and morphsingle =',)
          print(m20,np.round(morph.m20,3),np.round(morphsingle.m20,3))
          #sys.exit()
          print('asymmetry center (handmade) =',center_ax,center_ay)
          print('Asymmetry , asymm background and asymmetry final (handmade) =')
          print(minimum_asymmetry,background_asymmetry,asymmetry_final)
          print('xc_asymmetry from statmorph =',morphsingle.xc_asymmetry)
          print('yc_asymmetry =', morphsingle.yc_asymmetry)
          print('asymmetry statmorph =', morphsingle.asymmetry)
          print('asymmetry statmorph previous =', morph.asymmetry)
          print('concentration hand-made =', concentration)
          print('concentration from stamorph =',morph.concentration,morphsingle.concentration)
          print("concentration da morph evidentemente e' sbagliata !!")
          #quit()
          print("Gini. Il secondo metodo di stima sembra migliore, ovvero piu' simile a quello che calcolo io a mano. Con statmorph su morphsingle.")
          time.sleep(sleep_time)
          #time.sleep(4)
          '''
          #print('smoothness SM =',smoothness_SM)
          #quit()
          # -----------------------------------------------------------------------------------
      #  #except :
        else :
          maskones=np.zeros(imagein.shape)
          mask_ok=maskones*1
          gini_SM.append(-9)
          flag_SM.append(-9)
          sn_per_pixel_SM.append(-9)
          asymmetry_SM.append(-9)
          concentration_SM.append(-9)
          smoothness_SM.append(-9)
          m20_SM.append(-9)
          ellipticity_SM.append(-9)
          elongation_SM.append(-9)
          rpetro_SM.append(-9)
          rhalf_SM.append(-9)
      else :
        maskones=np.zeros(imagein.shape)
        mask_ok=maskones*1
        gini_SM.append(-9)
        flag_SM.append(-9)
        sn_per_pixel_SM.append(-9)
        asymmetry_SM.append(-9)
        concentration_SM.append(-9)
        smoothness_SM.append(-9)
        m20_SM.append(-9)
        ellipticity_SM.append(-9)
        elongation_SM.append(-9)
        rpetro_SM.append(-9)
        rhalf_SM.append(-9)


# -----------------------------------------------------------------------------------------------------------

    

      if (calculate_shape_asymmetry_easier==True):
        print('\n\nCalculate shape asymmetry')
        shapeasy_1=1 


        if shapeasy_1==1 :
        #try :
          #print('\n Calculate the shape asymmetry...')
          # maskgood_smaller e' la versione reduced(tagliata) di binary_pawlik
          show_result=False
          shape_asymmetry_last,area_best,xc_best,yc_best=asymmetry_function2_shape(center_ax,center_ay,maskgood,maskgood,show_result)
          # print('shape asymmetry last (num) =',shape_asymmetry_last)
          shape_asymmetry=shape_asymmetry_last/(2.*np.sum(maskgood))
          # print('Shape asymmetry =',shape_asymmetry)
          shape_asymmetry=np.round(shape_asymmetry,3)
          #shape_asymmetry_R=np.round(shape_asymmetry,3)
          print('Shape asymmetry =',shape_asymmetry)
          if ( (shape_asymmetry > -5) & (shape_asymmetry < 5) ) :
            print('ok shape asymmetry')
          else : shape_asymmetry=-9
          #sys.exit()
          #print('Anche per la shape asymmetry devi togliere quella del background ???')
        #except :
        else :
          shape_asymmetry=-9.
      else : 
        shape_asymmetry=-9.

      if np.sum(np.ravel(maskgood)) <= 0 :
        shape_asymmetry=-9.
    
    
      # --------- -------- --------- --------- ---------- ---------
      
      #print(shape_asymmetry_all)
      #print(smoothness_all)
      #quit()
      shape_asymmetry_all.append(shape_asymmetry)
      #asymmetry_all2.append(-9)
      ID_all.append(IDgal)













      # CLUMPINESS

      # -----------------------------------------------------------------------------------------------------------

      if calculate_clumpiness==True :

        # INIZIO SECONDA PARTE CLUMPINESS CALCULATION 
        # Importo la tabella con info di dove sono i nuclei (fino a un max di 3), e per ciascun clump calcolo la distanza da quel nucleo
        #  deb = np.genfromtxt(file_nuclei_smaller,names=True,dtype=None)
        #  
        #  label1x=deb[deb['IDgal']==IDgal]['label1x']
        #  label1y=deb[deb['IDgal']==IDgal]['label1y']
        #  label2x=deb[deb['IDgal']==IDgal]['label2x']
        #  label2y=deb[deb['IDgal']==IDgal]['label2y']
        #  label3x=deb[deb['IDgal']==IDgal]['label3x']
        #  label3y=deb[deb['IDgal']==IDgal]['label3y']
        #  label4x=deb[deb['IDgal']==IDgal]['label4x']
        #  label4y=deb[deb['IDgal']==IDgal]['label4y']
        #  debbl=deb[deb['IDgal']==IDgal]['deblend']
        #  label_merger=deb[deb['IDgal']==IDgal]['merger_vis']
        #  merger_all[uu,rr]=label_merger*1
        #label_all=[label1x,label1y,label2x,label2y,label3x,label3y,label4x,label4y]
        # ADESSO DEVO SCEGLIERE LA REGIONE CHE è PIÙ VICINA AL CENTRO
        #print "Ho sostituito morph con morphsingle, spero che trovi lo stesso il centro, ovvero che il centroide della full source e' uguale al centroide della single galaxy"


        # ---------------------------------------------------------------------------------
      
        image_4_clumpiness=image_new*1  # image_4_clumpiness e' immagine background subtracted 
        #deblend_clumpmap=0  # Default was 1 !!!!
    
        # 3b) SMOOTHED IMAGE AND ORIGINAL-SMOOTHED
     
        ##############################
        # Create the smoothed image  #
        ##############################
    
        #print('Image should be background subtracted ????')
        # smoothradiusX=1*conv_kpc_to_arcsec/IDpixscale_banda[jkk] # 1 kpc is the size of the clumps we want to detect
        # smoothradiusX2=Rmax/5. # Come in pawlik credo, o in conselice ?
        # smoothradiusX3=Rmax*2
        #print('smoothradius 1 and 2 (pixels) and pixel-scale =',smoothradiusX,smoothradiusX2,smoothradiusX3,IDpixscale[iyy])
        # smoothradius_eff=smoothradiusX*1

        smoothradius_eff=smoothradiusX22*1 # Same smoothing radius used for the smoothness calculation
        factor_smooth=1
        # clumpiness con smoothing length= 1 kpc oppure Rmax/5 mi viene simile a smoothness da morphsingle. Interessante !!!
    
        #imagesmooth=smooth(image_4_clumpiness,smoothtype,smoothRadius[uu],15)  # Il 15 alla fine che vuol dire ??
        imagesmooth=smooth(image_4_clumpiness,smoothtype,int(smoothradius_eff), factor_smooth)  # Il 15 alla fine che vuol dire ??
    
        # Calculate the residual image
        residuals_orig=image_4_clumpiness-imagesmooth  # image_4_clumpiness e' immagine background subtracted
        residuals=residuals_orig*1
        residuals[residuals < 0]=0    # Per Conselice 2003 !!!!
        ###residuals=residuals*maskin  # QUESTO CI VA SE STAI APPLICANDO UNA MASCHERA, SE NON APPLICHI MASCHERE LASCIA PERDERE. 
        # print('End smoothed image and original-smoothed calculation')
        ############################
        # END     smoothed image  #
        ############################
        
        # 3c) Calculate the clumpiness (Conselice 2003/4)
        ############################ #########################################################
        # HERE CALCULATE THE CLUMPINESS (STANDARD METHOD, or CONSELICE 2004/3)
        ############################ #########################################################
        
        # METHOD 1) : (standard)
        if (clumps_method==0) :
          # Calculate good points for clumpiness estimation
          sigmaback=sigmabackBK*1
          # print('sigmaback =',sigmabackBK)
          #print('maschera esiste ? ',maschera_)
          #sigmaback=0.0024   # Scusa perche' il background e' fisso a 0.0024 ??? Perche' sono simulazioni. E' solo per le simulazioni.
          
          maschera_for_clumpiness=True  # Puoi scegliere maskgood che al momento e' un quadrato, oppure binary_pawlik !!
          #mask_4_clumps=maskgood*1
          mask_4_clumps=binary_pawlik_4*1
          #plt.imshow(mask_4_clumps,origin='lower')
          #plt.show()
          
          if (maschera_for_clumpiness==True):
            # SICURO CHE SERVONO TUTTE E DUE LE RIGHE pxm e goodPointsClump ?
            goodPointsClump= np.where( (residuals/sigmaback >= threshold_clumpiness) & (mask_4_clumps == 1)) 
            ctClump=len(goodPointsClump[0])
            pxm = ( (residuals/sigmaback >= threshold_clumpiness) & (mask_4_clumps ==1) )
            pxm_negative = ( (residuals/sigmaback < threshold_clumpiness) & (mask_4_clumps == 1) )
          else :
            goodPointsClump= np.where( (residuals/sigmaback >= threshold_clumpiness))  
            ctClump=len(goodPointsClump[0])
            pxm = (residuals/sigmaback >= threshold_clumpiness)
            pxm_negative = (residuals/sigmaback < threshold_clumpiness)
            print('Attenzione non stai mettendo una maschera. Non va bene')
            #sys.exit()
            #quit()
          
          print('ctClump (how many pixels in clumps) =',ctClump)
          # Derive the pixelmap
          pixelmap = pxm.astype(int)
          pixelmap_neg = pxm_negative.astype(int)
          thresholdnuclei= np.zeros(pixelmap.shape)
        
        else : print('Only method 0 implemented as of now') ; sys.exit()
        # print('Stoppa un momento qua -----')
          
    
        # ------- ------- -------- -------- -------- ------- ------- -------  -------
  
        if (clumps_method==0) :
          # 3d) Segmentation map on clumps
          if (ctClump > 0 ) :
            #kernel1b=np.array((makeGaussian(20, fwhm = IDfwhm[i]/IDpixscale[i] , center=None)))
            #kernel1b=kernel1b/np.sum(kernel1b)
            
            #   # PER TROVARE IL NUCLEO devo fare una segmentation map e poi rimuovere nucleo :
            #   # APPLICA UNA SEGMENTATION MAP ALLA PIXELMAP
            #   CLUMPS_segm = detect_sources(pixelmap, thresholdnuclei , npixels=npixels,connectivity=8) # connectivity puo essere 4 oppure 8 (vedi web x significato) # Beh per trovare il nucleo npixels puo' anche essere 1, no ???
            #   
            #   # FACCIO IL DEBLENDING IN BASE A QUELLO CHE HO SCRITTO IN UN FILE x ciascuna galassia:
            #   #CLUMPS_debl = deblend_sources(imagein, CLUMPS_segm, npixels=npixels, filter_kernel=kernel1b,relabel=True)
            #   CLUMPS_debl = deblend_sources(imagein*maskgood, CLUMPS_segm, npixels=npixels,relabel=True,contrast=0, nlevels=20,connectivity=8, mode='linear')
            #   
            #   CLUMPS_segm_original = detect_sources(pixelmap, thresholdnuclei , npixels=npixels,connectivity=8) # connectivity puo essere 4 oppure 8 (vedi web x significato)
            #   #CLUMPS_debl_original = deblend_sources(imagein, CLUMPS_segm_original, npixels=npixels, filter_kernel=kernel1b,relabel=True)
            #   CLUMPS_debl_original = deblend_sources(imagein*maskgood, CLUMPS_segm_original, npixels=npixels,relabel=True, contrast=0, nlevels=20,connectivity=8, mode='linear')   # 0.00001
            #   CLUMPS_segm_nucleus = detect_sources(pixelmap, thresholdnuclei , npixels=npixels,connectivity=8) # connectivity puo essere 4 oppure 8 (vedi web x significato)
            #   #CLUMPS_debl_nucleus = deblend_sources(imagein, CLUMPS_segm_nucleus, npixels=npixels, filter_kernel=kernel1b,relabel=True)
            #   CLUMPS_debl_nucleus = deblend_sources(imagein*maskgood, CLUMPS_segm_nucleus, npixels=npixels,relabel=True,contrast=0, nlevels=20,connectivity=8, mode='linear')
            
            # print("Per la threshold 1.5 ho deciso di non fare il deblending dei clumps, altrimenti mi da' un problema nel deblending. Da threshold 3 in poi invece sono sicuro che va tutto bene.")
            thresholdnuclei = photutils.detect_threshold(pixelmap, nsigma=0) 
            #npixels=npixels_original*1
            #print('npixels =',npixels)
            npixels=1 # Also 1 pixel is enough to make a clump
            #print('npixels =',npixels)
            # PER TROVARE IL NUCLEO devo fare una segmentation map e poi rimuovere nucleo :
            # APPLICA UNA SEGMENTATION MAP ALLA PIXELMAP
            CLUMPS_debl = detect_sources(pixelmap, thresholdnuclei, npixels=npixels) # connectivity puo essere 4 oppure 8 (vedi web x significato) # Beh per trovare il nucleo npixels puo' anche essere 1, no ???
            # FACCIO IL DEBLENDING IN BASE A QUELLO CHE HO SCRITTO IN UN FILE x ciascuna galassia:
            ##CLUMPS_debl = deblend_sources(imagein, CLUMPS_segm, npixels=npixels, filter_kernel=kernel1b,relabel=True)
            #CLUMPS_debl = deblend_sources(imagein*maskgood, CLUMPS_segm, npixels=npixels,relabel=True,contrast=0)
            CLUMPS_debl_original = detect_sources(pixelmap, thresholdnuclei, npixels=npixels) # connectivity puo essere 4 oppure 8 (vedi web x significato)
            ##CLUMPS_debl_original = deblend_sources(imagein, CLUMPS_segm_original, npixels=npixels, filter_kernel=kernel1b,relabel=True)
            #CLUMPS_debl_original = deblend_sources(imagein*maskgood, CLUMPS_segm_original, npixels=npixels,relabel=True, contrast=0)   # 0.00001
            # CLUMPS_debl_nucleus = detect_sources(pixelmap, thresholdnuclei , npixels=npixels,connectivity=8) # connectivity puo essere 4 oppure 8 (vedi web x significato)
            ##CLUMPS_debl_nucleus = deblend_sources(imagein, CLUMPS_segm_nucleus, npixels=npixels, filter_kernel=kernel1b,relabel=True)
            # CLUMPS_debl_nucleus = deblend_sources(imagein*maskgood, CLUMPS_segm_nucleus, npixels=npixels,relabel=True,contrast=0, nlevels=20,connectivity=8, mode='linear')
            if ( (deblend_clumpmap==1) ) :
              nlevelsX=4
              print('Also deblend segmentation map of clumps')
              # CLUMPS_debl = deblend_sources(imagein*maskgood, CLUMPS_debl, npixels=npixels,relabel=True,contrast=0, nlevels=nlevelsX,connectivity=8, mode='linear')
              # CLUMPS_debl_original = deblend_sources(imagein*maskgood, CLUMPS_debl_original, npixels=npixels,relabel=True, contrast=0, nlevels=nlevelsX,connectivity=8, mode='linear')
              CLUMPS_debl = deblend_sources(image_4_clumpiness, CLUMPS_debl, npixels=npixels,relabel=True,contrast=0, nlevels=nlevelsX,connectivity=8, mode='linear')
              CLUMPS_debl_original = deblend_sources(image_4_clumpiness, CLUMPS_debl_original, npixels=npixels,relabel=True, contrast=0, nlevels=nlevelsX,connectivity=8, mode='sinh') # exponential, sinh,linear
              # CLUMPS_debl_nucleus = deblend_sources(imagein*maskgood, CLUMPS_debl_nucleus, npixels=npixels,relabel=True,contrast=0, nlevels=nlevelsX,connectivity=8, mode='linear')
            mapclumps_all=CLUMPS_debl_original.data
            #if banda_obj=='f090w' :
            #  f55, ((ax1, ax2)) = plt.subplots(1, 2, sharey=True, sharex=True)
            #  ax1.imshow(image_4_clumpiness,origin='lower')
            #  ax1.set_title(str(IDgal),fontsize=20)
            #  ax2.imshow(mapclumps_all,origin='lower')
            #  plt.show()
            #  plt.close(f55)
            if remove_nucleus==True :
              print('\nRemove nuclei')
              #print(N1_x,N1_y,N2_x,N2_y)
              #quit()
              indice_mappa1=mapclumps_all[N1_y,N1_x]
              indice_mappa2=mapclumps_all[N2_y,N2_x]
              mapclumps_all[mapclumps_all==indice_mappa1]=0
              mapclumps_all[mapclumps_all==indice_mappa2]=0
          else :
            print('galaxy ID =',IDgal,' non ha clumps')
            CLUMPS_debl_original = np.zeros(image_4_clumpiness.shape)
            mapclumps_all=CLUMPS_debl_original.data


          
    
    
        # CONTINUATION OF CLUMPINESS CALCULATION 
  
        # 3e) visualization clumps and segmentation map on clumps
        
        if show_figures_clumps==True :
          try :
            # VISUALIZE ORIGINAL IMAGE, CLUMPY PIXELS and deblended segmenttion map on detected clumps
            # Colormap normalization according to image statistics (background and rms)
            f, ((ax1, ax2, ax4), (ax5, ax6, ax7)) = plt.subplots(2, 3, sharey=True, sharex=True)
            
            upper_ct=10

            imagestat1=calc_bg(image_4_clumpiness[maskgood==1])
            back1=imagestat1[0]   #'background'
            sigmaback1=imagestat1[1]
            norm1=matplotlib.colors.Normalize(vmin=back1-2*sigmaback1, vmax=back1+upper_ct*sigmaback1, clip=False)
            
            imagestat3=calc_bg(image_4_clumpiness[maskgood==0])
            #imagestat3=calc_bg(np.ravel(image_4_clumpiness))
            back3=imagestat3[0]   #'background'
            sigmaback3=imagestat3[1]
            norm3=matplotlib.colors.Normalize(vmin=back3-2*sigmaback3, vmax=back3+upper_ct*sigmaback3, clip=False)
      
            ax1.imshow(image_4_clumpiness, origin='lower', norm=norm1, cmap='hot')
            ax1.set_title('ID '+str(IDgal)+' (original)')
          
            ax2.imshow(maskgood, origin='lower', cmap='hot')
            ax2.set_title('Segmentation map')
            
            # Immagine originale con maschera
            #ax3.imshow(maskgood*image_4_clumpiness, origin='lower', norm=norm1, cmap='hot')
            #ax3.set_title(str(IDgal))
            #plt.title(str(IDgal))
            #plt.colorbar()
            
            # Original - smoothed (= Residual)
            imagestat4=calc_bg(residuals_orig)
            back4=imagestat4[0]   #'background'
            sigmaback4=imagestat4[1]
            norm4=matplotlib.colors.Normalize(vmin=back4-2*sigmaback4, vmax=back4+upper_ct*sigmaback4, clip=False)
            ax4.imshow(residuals, origin='lower', norm=norm4, cmap='hot')
            ax4.set_title('Orig - smoothed')
      
            # Segmentation map (Deblended) on clumps
            image_neg=image_4_clumpiness*1 
            image_neg[pixelmap>0]=-10
            
            #CLUMPS_debl_original_visual=CLUMPS_debl_original*1
            #CLUMPS_debl_original_visual[CLUMPS_debl_original_visual==0]=-10
            ax5.imshow(image_neg, origin='lower', cmap='hot', norm=norm3) # Questo ti serve poi per rimuovere i nuclei
            # #ax6.imshow(segnuclei_regul*maskgood, origin='lower', cmap='hot')
            ax5.set_title('Orig no clumps')
      
            # Immagine originale senza clumps
            #ax5.imshow(pixelmap_neg*image_4_clumpiness, origin='lower', norm=norm3, cmap='hot')
            image_neg2=pixelmap_neg*1 
            image_neg2[pixelmap>0]=-10
            ax6.imshow(image_neg2*image_4_clumpiness, origin='lower', norm=norm3, cmap='hot')
            ax6.set_title('masked')
            
            # Segmentation map (Deblended) on clumps
            ax7.imshow(mapclumps_all, origin='lower', cmap='hot', norm=norm1)
            ax7.set_title('Map clumps')      
            
            plt.subplots_adjust(wspace=0.2, hspace=0.01)
            #plt.tight_layout()
            #if (show_all==True) :
            namefile=output_folder+str(IDgal)+'_clumpmap_'+banda_obj+'.png'
            plt.savefig(namefile,dpi=80)
            #plt.show()
            plt.close(f)
          except :
            f, ((ax1, ax2)) = plt.subplots(1, 2)
            namefile=output_folder+str(IDgal)+'_clumpmap_'+banda_obj+'.png'
            plt.savefig(namefile,dpi=80)
            plt.close(f)
        

        #quit()
    
        # 3f) Analyze clumps properties :
        
        # Utilizzo STATMORPH su original image :
        #source_morphs = statmorph.source_morphology(imagein, maskgood, gain=1) #, mask=pixelmap_neg)  # maskgood e' la maschera che contiene la galassia scelta. Se e' un merger system ci saranno piu' parti, o anche se c'e' dell'emissione diffusa staccata.
        #print 'source_morphs =',len(source_morphs)
        #print('Analyze clumps properties !!!')    
        
        # OPTION 1 (VA INSIEME A N.3)
        
        #    import random
        #    image_new=imagein-backBK
        #    image_background=image_new*maskgood
        #    for i1 in np.arange(imagein.shape[0]) :
        #      for j1 in np.arange(imagein.shape[1]) :
        #        if (image_background[i1,j1]==0) :
        #          image_background[i1,j1]=random.gauss(0,sigmabackBK)
        #    # Calculate weightmap (NOISE map !!)
        #    noisy_img = np.random.normal(0, sigmabackBK, imagein.shape)
        #    ##noisy_img_clipped = np.clip(noisy_img, 0, 255)
        
        # OPTION 2 (GOES WITH N.5 later)
  
        #image_background[image_background==0]=random.gauss(backBK,sigmabackBK) #np.random.normal(backBK , sigmabackBK, imagein.shape)
        
        # N.1
        #source_morphs = statmorph.source_morphology(imagein*maskgood, maskgood, weightmap=noisy_img, psf=kernel1b)
        # N.2
        #source_morphs = statmorph.source_morphology(imagein, maskgood, gain=1)
        # N.3
        
        #    #limit1=int(len(imagein)/2-200)
        #    #limit2=int(len(imagein)/2+200)
        #    #maskones=maskgood*1
        #    ##maskones[222:282,212:272]=1
        #    source_morphs = statmorph.source_morphology(image_background, maskgood, weightmap=noisy_img, psf=kernel1b)
        #    # N.4  - non funziona
        #    #source_morphs = statmorph.source_morphology(image_background, np.ones(imagein.shape), weightmap=noisy_img, psf=kernel1b)
    

        # CONTINUATION OF CLUMPINESS CALCULATION 
  
        if (ctClump>0) :
          # Creo una tabella con la proprieta' dei clumps
          #print(CLUMPS_debl)
          
          props = SourceCatalog(image_4_clumpiness, CLUMPS_debl)
          #print('\n-----------------------')
          #print(props.centroid)
          #sys.exit()
          xcen=np.empty(len(CLUMPS_debl.labels)) ; ycen=np.empty(len(CLUMPS_debl.labels))
          dist=np.empty(len(CLUMPS_debl.labels)) ; dist2=np.empty(len(CLUMPS_debl.labels)) ; dist3=np.empty(len(CLUMPS_debl.labels)) ; dist4=np.empty(len(CLUMPS_debl.labels)) ; dist5=np.empty(len(CLUMPS_debl.labels)) ; dist6=np.empty(len(CLUMPS_debl.labels)) ;  _area=np.empty(len(CLUMPS_debl.labels)) ; _photo=np.empty(len(CLUMPS_debl.labels)) ; _peak=np.empty(len(CLUMPS_debl.labels))
          for l in CLUMPS_debl.labels :
            prop=props[l-1]
            position = (prop.xcentroid, prop.ycentroid)
            #print('position centroid =',position)
            xcen[l-1]=prop.xcentroid
            ycen[l-1]=prop.ycentroid
            #print position
            dist2[l-1] = math.hypot(position[0] - center_ax , position[1] - center_ay)
            dist[l-1] = math.hypot(position[0] - imagein.shape[0]/2, position[1] - imagein.shape[1]  /2  ) ; _area[l-1] = prop.area.value ; 
            test_data=prop.data_ma
            #print('area clump =',prop.area.value)
            #print(test_data)
            #print('Summed area =',test_data.sum())
            #sys.exit()
            _photo[l-1] = prop.data_ma.sum() # prop.source_sum ; 
            _peak[l-1] = prop.max_value
            #dist3[l-1]= math.hypot(position[0] - label1x, position[1] - label1y)
            #dist4[l-1]= math.hypot(position[0] - label2x, position[1] - label2y)
            #dist5[l-1]= math.hypot(position[0] - label3x, position[1] - label3y)
            #dist6[l-1]= math.hypot(position[0] - label4x, position[1] - label4y)
            
          #print 'distanze da centroide =',dist2
          #mlabels= np.column_stack((CLUMPS_debl.labels,dist,dist2,dist3,dist4,dist5,dist6,_area,_photo,_peak,xcen,ycen))
          mlabels= np.column_stack((CLUMPS_debl.labels,dist,dist2,_area,_photo,_peak,xcen,ycen))
  
          # print('mlabels =')
          # print('# label distance dist2 area photo peak xcen ycen')
          # for ikk in np.arange(len(mlabels)):  print(mlabels[ikk])
          
          # NO NEED TO REMOVE NUCLEUS NOW !!!!
          # 3g) Identify and remove nucleus :    
          # Devi far girare una volta il programma con tutti zeri nel file (nelle 6 colonne), e poi ti scrivi nella tabella tutte le coordinate di dove sono i nuclei. E poi lo fai girare una seconda volta con la tabella nuova
          '''
          def GetKey0(item):
            return item[7]   # TOTAL FLUX
          def GetKey1(item):
            return item[8]   # PEAK FLUX OF THE CLUMP
          def GetKey2(item): 
            # (criterio per ordinare le colonne)
            return item[2]   # 1 è la distanza dal centro dell'immagine # 2 e' la distanza dal   centroide calcolato con sersic profile fitting (in realta' con statmorph sulla   COMPONENTE 0 pero', se ci sono piu' componenti puo' essere un problema)
          def GetKey3(item):
            return item[3]
          def GetKey4(item):
            return item[4]
          def GetKey5(item):
            return item[5]
          def GetKey6(item):
            return item[6]
          
          print('label1x =',label1x)
          
          if ((label1x[0]==0)):  # Tanto se e' 0 il primo sono zero anche i successivi
            #mgood2=sorted(mlabels,key=GetKey2) # DISTANZA DAL centroide del clump
            #labelnucleus2=int(mgood2[0][0]) 
            mgood2=sorted(mlabels,key=GetKey0)
            labelnucleus2=int(mgood2[-1][0]) # Prendo l'ultimo elemento, cioe' il flusso totale maggiore
            # CASO 1) Rimuovo clump piu' vicino a centro immagine di default
            # Questo serve per ordinare le righe della tabella a seconda del valore di una colonna (come in topcat ..)
            xcenG=int(mgood2[-1][10])   # 10 e' il x-centroide
            ycenG=int(mgood2[-1][11])  # 11 e' il y centroide
            print('Spero x e y non siano scambiati anche qua !!!')
            #if (SB==1): # Per le Starbursts rimuovo anche il secondo nucleo piu luminoso (vediamo se funziona cosi' risparmio un po' di lavoro)
            #  mgood4=sorted(mlabels,key=GetKey0)
            #  labelnucleus4=int(mgood4[-2][0]) # Prendo l'ultimo elemento, cioe' il flusso totale maggiore
            #  # CASO 1) Rimuovo clump piu' vicino a centro immagine di default
            #  CLUMPS_debl.remove_labels(labelnucleus2,labelnucleus4)
            #  # Questo serve per ordinare le righe della tabella a seconda del valore di una colonna (come in topcat ..)
            #  CLUMPS_debl_nucleus.keep_labels(labelnucleus2,labelnucleus4)
            #else :
            CLUMPS_debl.remove_labels(labelnucleus2)
            CLUMPS_debl_nucleus.keep_labels(labelnucleus2)
    
          elif (label1x[0]==-1):
            print("Non c'e' alcun nucleo")
    
          else :
            # CASO 2) Uso i label dal file (pero' cosi' c'e' il rischio che se faccio la procedura in un modo diverso poi cambia questo)
            #if (label1 != 0) :
            #  CLUMPS_debl.remove_labels(label1)
            #if (label2 != 0) :
            #  CLUMPS_debl.remove_labels(label2)
            #if (label3 != 0) :
            #  CLUMPS_debl.remove_labels(label3) 
            allnuc=[]
            if (label1x[0]!=0) :   # RIMUOVO PRIMO NUCLEO
              mgood3=sorted(mlabels,key=GetKey3)
              labelnucleus3=int(mgood3[0][0]) 
              allnuc=np.append(allnuc,labelnucleus3)
              CLUMPS_debl.remove_labels(labelnucleus3)
            if (label2x[0]!=0) :   # RIMUOVO SECONDO NUCLEO
              mgood4=sorted(mlabels,key=GetKey4)
              labelnucleus4=int(mgood4[0][0]) 
              allnuc=np.append(allnuc,labelnucleus4)
              CLUMPS_debl.remove_labels(labelnucleus4)
            if (label3x[0]!=0) :   # RIMUOVO TERZO NUCLEO
              mgood5=sorted(mlabels,key=GetKey5)
              labelnucleus5=int(mgood5[0][0]) 
              allnuc=np.append(allnuc,labelnucleus5)
              CLUMPS_debl.remove_labels(labelnucleus5) 
            if (label4x[0]!=0) :   # RIMUOVO QUARTO NUCLEO
              mgood6=sorted(mlabels,key=GetKey6)
              labelnucleus6=int(mgood6[0][0]) 
              allnuc=np.append(allnuc,labelnucleus6)
              CLUMPS_debl.remove_labels(labelnucleus6) 
            CLUMPS_debl_nucleus.keep_labels(allnuc)
    
    
          # Costruisco una mappa estesa dei nuclei (che pero' e' approssimata eh, perche' le coordinate che do' io dei nuclei sono approssimative !! E tra l'altro sono numeri interi)
          map_nuclei=np.zeros(imagein.shape)
          RR=5   # 3 pixels
          if ((label1x[0]==0)):
            print('xcen, ycen =',xcenG,ycenG)
            print('masking the only nucleus present here')
            map_nuclei[ycenG-RR:ycenG+RR,xcenG-RR:xcenG+RR]=1
          else :  
            if (label1x[0]!=0):
              map_nuclei[label1x[0]-RR:label1x[0]+RR,label1y[0]-RR:label1y[0]+RR]=1
            if (label2x[0]!=0):
              map_nuclei[label2x[0]-RR:label2x[0]+RR,label2y[0]-RR:label2y[0]+RR]=1
            if (label3x[0]!=0):
              map_nuclei[label3x[0]-RR:label3x[0]+RR,label3y[0]-RR:label3y[0]+RR]=1
            if (label4x[0]!=0):
              map_nuclei[label4x[0]-RR:label4x[0]+RR,label4y[0]-RR:label4y[0]+RR]=1
          map_no_nuclei=np.absolute(map_nuclei-1)
          #map_no_nuclei.astype(int)
          #print map_no_nuclei
          
          #plt.imshow(imagein)
          #plt.show()
    
          #plt.imshow(imagein*map_no_nuclei)
          #plt.show()
    
          
          mapnucleus=CLUMPS_debl_nucleus.data
          mapnucleus[mapnucleus>0]=1   # Per far si che anche se e' un intero > 1 diventi 1 (cosi ho solo 1 e 0)
          imagenucleus=imagein*mapnucleus # Immagine moltiplicata per la maschera del nucleo (Rimane imagein solo dove c'e' il nucleo)
          #MASKf=CLUMPS_debl.data/labelnucleus
          # In questo modo dovrei avere tutti 1 per i labelgood e tutti 0 outside !
          #imageregion=imagein*MASKf
          #plt.imshow(imagenucleus,origin='lower')
          #plt.show()
          #sys.exit()
          
          mapclumps_all=CLUMPS_debl_original.data
          #plt.imshow(mapclumps_all,origin='lower')
          #plt.show()
          mapclumps_all[mapclumps_all>0]=1  # PER CREARE LA MASCHERA DI 1 dove ci sono i clumps
          imagepxm_all=imagein*mapclumps_all
  
          mapclumps=CLUMPS_debl.data
          #print(mapclumps)
          #plt.imshow(mapclumps,origin='lower')
          #plt.show()
          #sys.exit()
          mapclumps[mapclumps>0]=1  # PER CREARE LA MASCHERA DI 1 dove ci sono i clumps
          imagepxm=imagein*mapclumps   # IMMAGINE MOLTIPLICATA PER LA MASCHERA DEI CLUMPS (si vede imagein solo dove sono i clumps)
          nclumps=CLUMPS_debl.nlabels
          number_clumps[uu,rr]=nclumps*1
          # CALCOLO L'IMMAGINE DIFFERENZA, CIOÈ CHE RIMANE FUORI DAI CLUMPS (compreso il nucleo)
          difference = imagein*maskgood - imagepxm
          difference_original = imagein*maskgood- imagepxm_all
        else :
          difference= imagein*maskgood
          difference_original = difference*1
        '''
  
  
        # ------------------------------------------------------------------------------
    
    
        # 3h) Clumpiness derivation     # Devi sottrarre la clumpiness del background
        # CLUMPINESS derivation :::::        fluxtot=sum(imageregion)
  
        #print 'Clumpiness 1 e 2 sono calcolati su tutto il sistema interagente (non sulla galassia singola'
        if (ctClump > 0) : 
          #mapclumps_all=CLUMPS_debl_original.data
          #plt.imshow(mapclumps_all,origin='lower')
          #plt.show()
          mapclumps_all[mapclumps_all>0]=1  # PER CREARE LA MASCHERA DI 1 dove ci sono i clumps
          imagepxm_all=image_4_clumpiness*mapclumps_all # *maskgood # QUesta mappa contiene i flussi nei clumps 
          #plt.imshow(imagepxm_all,origin='lower')
          #plt.show()
  
          #plt.imshow(maskgood*imagein,origin='lower')
          #plt.show()
  
          #print('sum imagepxm =',np.sum(imagepxm_all))
          #print('clumpiness (nuclei inclusi) no rounded =',np.sum(imagepxm_all)/np.sum(maskgood*imagein))
          clumpiness_BB=np.round(np.sum(imagepxm_all)/np.sum(maskgood*image_4_clumpiness),3)
          print('\nClumpiness 1 =',clumpiness_BB)
          #print('clumpiness (nuclei inclusi) =',clumpiness[uu,rr])
          #print('smoothness (me and statmorph) =',smoothness22,morphsingle.smoothness,morph.smoothness)
          #print 'qui rar'
          #if (label1x[0]==-1) : # Non ci sono nuclei
          #  clumpiness2[uu,rr]=clumpiness[uu,rr]*1
          #else :
          #  clumpiness2[uu,rr]=np.round(np.sum(imagepxm_all-imagenucleus)/np.sum(maskgood*imagein-imagenucleus),4)
          #  #print '\nclumpiness =',clumpiness[uu,rr],clumpiness2[uu,rr]
          #print('clumpiness2 (senza nucleo) =',clumpiness2[uu,rr])
          clumpiness2_BB=0
          #sys.exit()
  
        else :
          clumpiness_BB=0
          clumpiness2_BB=0
          CLUMPS_debl_original=np.zeros(imagein.shape)
          CLUMPS_debl_nucleus=np.zeros(imagein.shape)
        #sys.exit()
  
        ############################ ############################ ############################
        # HERE ENDS THE CLUMPINESS (METODO STANDARD, o CONSELICE 2004/3)
        ############################ ############################ ############################
        '''
        # 3iB) CLUMPINESS derivation per single galaxies case :: fluxtot=sum(imageregion)
        if (ctClump > 0) : 
          print("clumpiness 3 fa rapporto flusso nei clumps e flusso totale (prendendo pero' solo singlegood)")
          #clumpiness3[uu,rr]=np.round(np.sum(maskgood_debl*imagepxm)/np.sum(maskgood_debl*imagein),4)
          clumpiness3[uu,rr]=np.round(np.sum(singlegood*imagepxm)/np.sum(singlegood*imagein),4)
          print('clumpiness3 =',clumpiness3[uu,rr])
          #print 'qui rar'
          if (label1x==-1) : # Non ci sono nuclei
            clumpiness4[uu,rr]=clumpiness[uu,rr]*1
          else :
            print('Clumpiness 4 fa rapporto flusso nei clumps e flusso totale senza nucleo (prendendo solo singlegood)')
            #clumpiness4[uu,rr]=np.round(np.sum(imagepxm*maskgood_debl)/np.sum(maskgood_debl*imagein-imagenucleus*maskgood_debl),4)
            clumpiness4[uu,rr]=np.round(np.sum(imagepxm*singlegood-imagepxm*imagenucleus)/np.sum(singlegood*imagein-singlegood*imagenucleus),4)
          print('clumpiness4 =',clumpiness4[uu,rr])
  
          if (label1x==-1) : # Non ci sono nuclei
            clumpiness5[uu,rr]=clumpiness[uu,rr]*1
          else :
            print('Clumpiness 5 fa rapporto flusso nei clumps e flusso totale (comprende tutto) - prende solo singlegood !!!')
            #clumpiness6[uu,rr]=np.round(np.sum(imagepxm*maskgood_debl)/np.sum(maskgood_debl*imagein-imagenucleus*maskgood_debl),4)
            clumpiness5[uu,rr]=np.round(np.sum(imagepxm*singlegood*map_no_nuclei)/np.sum(singlegood*imagein),4)
          print('clumpiness5 =',clumpiness5[uu,rr])
  
          #image_nuclei=imagein*map_nuclei
          if (label1x==-1) : # Non ci sono nuclei
            clumpiness6[uu,rr]=clumpiness[uu,rr]*1
          else :
            print('Clumpiness 6 fa rapporto flusso nei clumps e flusso totale senza nuclei (prende solo singlegood)')
            print('Quindi tra 4 e 6 quale essere la differenza ?')
            #clumpiness6[uu,rr]=np.round(np.sum(imagepxm*maskgood_debl)/np.sum(maskgood_debl*imagein-imagenucleus*maskgood_debl),4)
            clumpiness6[uu,rr]=np.round(np.sum(imagepxm*singlegood*map_no_nuclei)/np.sum(singlegood*imagein*map_no_nuclei),4)
          print('clumpiness6 =',clumpiness6[uu,rr])
  
        else :
          clumpiness3[uu,rr]=0 ; clumpiness4[uu,rr]=0 ; clumpiness5[uu,rr]=0 ; clumpiness6[uu,rr]=0
        '''
        #clumpiness3[uu,rr]=0
        #clumpiness4[uu,rr]=0
        #clumpiness5[uu,rr]=0
        #clumpiness6[uu,rr]=0
        #print('clumpiness 3, 4, 5, 6 =',clumpiness3[uu,rr],clumpiness4[uu,rr],clumpiness5[uu,rr],clumpiness6[uu,rr])
        #print '\nclumpiness 3 e 4=',clumpiness3[uu,rr],clumpiness4[uu,rr]
  
        #fin5.write(str(IDgal)+' '+str(np.round(clumpiness[uu,rr],3))+' '+str(np.round(clumpiness2[uu,rr],3))+' '+str(np.round(clumpiness3[uu,rr],3))+' '+str(np.round(clumpiness4[uu,rr],3))+' '+str(np.round(clumpiness5[uu,rr],3))+' '+str(np.round(clumpiness6[uu,rr],3))+'\n')
  
  
  
        # -------------------------------------------------------------
  
        
        # 3k) Visualization clumps :
        
        try :
          # For logarithmic scale :
          log_stretch = LogStretch(a=10000.0)
          def normalize(image):
            m, M = np.min(image), np.max(image)
            return (image-m) / (M-m)
          # plt.imshow(log_stretch(normalize(image)), origin='lower', cmap='gray')
          
          '''
          # f, ((ax1, ax2, ax4),(ax7, ax9, ax11)) = plt.subplots(2, 3, sharey=True, sharex=True,figsize=(12,10))
          f, ((ax1, ax2, ax4)) = plt.subplots(1, 3, figsize=(16,6))
          
          limit1=int(len(image_4_clumpiness)/2-3*morphsingle.rhalf_ellip)
          limit2=int(len(image_4_clumpiness)/2+3*morphsingle.rhalf_ellip)
          #limit1=int(len(image_4_clumpiness)/2-100)
          #limit2=int(len(image_4_clumpiness)/2+100)
    
          imagestat1=calc_bg(image_4_clumpiness[limit1:limit2,limit1:limit2])
          back1=imagestat1[0]   #'background'
          sigmaback1=imagestat1[1]
          norm1=matplotlib.colors.Normalize(vmin=back1-1*sigmaback1, vmax=back1+3*sigmaback1, clip=True)
          ax1.imshow(image_4_clumpiness, origin='lower', norm=norm1, cmap='hot') # norm=norm1
          ax1.set_title(str(IDgal))
    
          # Segmentation map on galaxy 
          #ax2.imshow(maskgood, origin='lower', cmap='hot')
          ax2.imshow(singlegood, origin='lower', cmap='hot')
          ax2.set_title('Segmap '+str(IDgal)) 
    
          # Original - smoothed (= Residual)
          imagestat4=calc_bg(residuals_orig)
          back4=imagestat4[0]   #'background'
          sigmaback4=imagestat4[1]
          norm4=matplotlib.colors.Normalize(vmin=back4-2*sigmaback4, vmax=back4+5*sigmaback4, clip=False)
          #ax4.imshow(residuals*maskgood, origin='lower', norm=norm4, cmap='hot')
          ax4.imshow(residuals*singlegood, origin='lower', norm=norm4, cmap='hot')
          ax4.set_title('Orig - smoothed')      
          
    
          # Segmentation map (Deblended) on clumps
          
          #   CLUMPS_debl_original_visual=CLUMPS_debl_original*1
          #   CLUMPS_debl_original_visual[CLUMPS_debl_original_visual==0]=-10
          #   ax7.imshow(CLUMPS_debl_original, origin='lower', cmap='hot')
          #   ax7.set_title('Segmap deb clumps '+str(IDgal),color='green',size=15)      
          #   
          #   ax9.imshow(singlegood*difference_original, norm=norm1, origin='lower', cmap='hot')
          #   ax9.set_title(str(i)+' (clumps black)',color='green',size=15) 
          #   # Segmentation map (Deblended) on clumps
          #   # Per momento non c'e' alcun nucleo ??
          #   ax11.imshow(CLUMPS_debl_nucleus, origin='lower', cmap='hot') # ,norm=norm1)
          #   ax11.set_title('nucleus',color='green',size=15) 
          
          # Part of image used by gini and M20 calculation (within 1.5 time rpetro)
          
          #zeros_temp=np.zeros(image_4_clumpiness.shape)
          #zeros_temp[xc_asym-rpetrocirc:xc_asym+rpetrocirc,yc_asym-rpetrocirc:yc_asym+rpetrocirc]=1
          #ax8.imshow(image_4_clumpiness*zeros_temp*maskgood, origin='lower', cmap='hot',norm=norm1)
          #ax8.set_title('image reduced',color='red') 
          
    
          #ax8.imshow(maskgood_debl_seg,origin='lower', cmap='hot',norm=norm1)
          #ax8.imshow(image_4_clumpiness*maskgood_debl,origin='lower', cmap='hot',norm=norm1)
          #ax8.set_title('single')
    
          # Immagine differenza
          #imagestat9=calc_bg(difference)
          #back9=imagestat9[0]   #'background'
          #sigmaback9=imagestat9[1]
          #norm9=matplotlib.colors.Normalize(vmin=back9-2*sigmaback9, vmax=back9+5*sigmaback9, clip=False)
    
          #ax10.imshow(difference, origin='lower', norm=norm1, cmap='hot')
          #ax10.set_title('difference',color='red')  
    
          # Immagine dei soli clumps
          #imagestat10=calc_bg(imagepxm)
          #back10=imagestat10[0]   #'background'
          #sigmaback10=imagestat10[1]
          #norm10=matplotlib.colors.Normalize(vmin=back10-2*sigmaback10, vmax=back10+5*sigmaback10, clip=False)
          #ax9.imshow(imagepxm*maskgood, origin='lower', norm=norm1, cmap='hot')
          #ax9.set_title('clumps',color='white') 
          #ax9.imshow(maskgood_debl*difference, norm=norm1, origin='lower', cmap='hot')
    
          #ax12.imshow(singlegood*difference, origin='lower', cmap='hot', norm=norm1)
          #ax12.set_title(str(i)+' (without nucleus)',color='green',size=15) 
          
          plt.xlim(limit1,limit2)
          plt.ylim(limit1,limit2)
          plt.subplots_adjust(wspace=0.02, hspace=0.02)
          #plt.tight_layout()
          #if ( (show_all==True)) :
          plt.show()
          '''
    
          # 3l) Per salvare dei fits files :
          #print('Stop before saving table')
          #sys.exit()
          if ((save_int==True) ):  # Se vuoi salvare per tutti
            # Save fits file 
            #outputFileName1=home_folder2+str(IDgal)+'_img.fits'
            #print('threshold =',threshold)
            outputFileName1=home_folder+str(IDgal)+'_clumps_'+str(deblend_clumpmap)+'_'+which_data+'.fits'    # '_thresh'+str(threshold[0])+'_.fits'
            outputFileName2=home_folder+str(IDgal)+'_clumpmap_'+str(deblend_clumpmap)+'_'+which_data+'.pdf'  # '_thresh'+str(threshold[0])+'_.pdf'
            astImages.saveFITS(outputFileName1,difference) #, imageWCS=WCS)  # Facendo vedere l'immagine differenza penso si vedano meglio i clumps
            f.savefig(outputFileName2)
            plt.close(f)
          time.sleep(sleep_time)
        except :
          time.sleep(sleep_time)
        clumpiness_all.append(clumpiness_BB)
        clumpiness2_all.append(clumpiness2_BB)
      else :
        clumpiness_all.append(-9)
        clumpiness2_all.append(-9)

          
  
      # ----------------------------------------------------------------------------------
    
    
      #print(clumpiness_all)
      #print('smoothness_SM =',smoothness_SM)
      #if counter989==8 :
      #  quit()
      #quit()
      gini_4_stamp=np.round(gini,2)
      m20_4_stamp=np.round(m20,2)
      concentration_4_stamp=np.round(concentration,2)
      shape_asymmetry_4_stamp=np.round(shape_asymmetry,2)
      clumpiness_4_stamp=np.round(clumpiness_BB,2)

      # NOW WRITE RESULTS IN STAMPS :
      if make_images==True :
        # WRITE RESULTS ON TOP OF CUTOUTS :
  
        axs[count_image].text(0.03, 0.2, 'G ='+str(gini_4_stamp), horizontalalignment='left',verticalalignment='top', fontsize=12, color='w',transform=axs[count_image].transAxes, weight='bold')
        axs[count_image].text(0.03, 0.1, r'M$_{20}$ ='+str(m20_4_stamp), horizontalalignment='left',verticalalignment='top', fontsize=12, color='w',transform=axs[count_image].transAxes, weight='bold')
        axs[count_image].text(0.97, 0.3, 'C ='+str(concentration_4_stamp), horizontalalignment='right',verticalalignment='top', fontsize=12, color='w',transform=axs[count_image].transAxes, weight='bold')
        axs[count_image].text(0.97, 0.2, r'A$_S$ ='+str(shape_asymmetry_4_stamp), horizontalalignment='right',verticalalignment='top', fontsize=12, color='w',transform=axs[count_image].transAxes, weight='bold')
        axs[count_image].text(0.97, 0.1, r'CL ='+str(clumpiness_4_stamp), horizontalalignment='right',verticalalignment='top', fontsize=12, color='w',transform=axs[count_image].transAxes, weight='bold')
        count_image+=1

    # Indentation of this part is correct  !!!!
    if ( (make_images==True) & (test_segmap_index==False) ) :
      plt.subplots_adjust(wspace=0.01) #, hspace=0.01)
      #plt.tight_layout()
      savefinalimage=output_folder+'image_cutouts_'+str(IDgal)+'_'+banda_obj+'.png'
      plt.savefig(savefinalimage, format='png', dpi=100)
      print('saved figure in '+savefinalimage)
      #if show_all==True : plt.show()
      plt.close(fig66)
       


  # 3m) Store clumpiness and clump properties in table (and then save)
  #print 'store clumpiness info in table'
  # print('\n\nSaving arrays into table (len vectors =)')

  #########################################################################
  # Store clumpiness and number of clumps in table/column to be saved :   #
  #########################################################################
  gini_vecSM=np.array(gini_SM)
  flag_vecSM=np.array(flag_SM)
  sn_per_pixel_vecSM=np.array(sn_per_pixel_SM)
  asymmetry_vecSM=np.array(asymmetry_SM)
  concentration_vecSM=np.array(concentration_SM)
  smoothness_vecSM=np.array(smoothness_SM)
  m20_vecSM=np.array(m20_SM)
  ellipticity_vecSM=np.array(ellipticity_SM)
  elongation_vecSM=np.array(elongation_SM)
  rpetro_vecSM=np.array(rpetro_SM)
  rhalf_vecSM=np.array(rhalf_SM)
  IDvector=np.array(ID_all)    
  clumpvector=np.array(clumpiness_all)
  clumpvector2=np.array(clumpiness2_all)
  #clumpvector3=np.array(clumpiness3)
  #clumpvector4=np.array(clumpiness4)
  #clumpvector5=np.array(clumpiness5)
  #clumpvector6=np.array(clumpiness6)
  #print(clumpiness)
  #print(len(clumpvector3))
  photoband_vector=np.array(photoband_all)
  Rmax_vector=np.array(Rmax_all)
  # nclumpvector=np.array(number_clumps)
  gini_vector=np.array(gini_all)
  m20_vector=np.array(m20_all)
  #print(len(gini_vector),len(m20_vector))
  # distlotz_vector=np.array(distlotz_all)
  #sn_per_pixel_vector1=np.array(sn_per_pixel)
  sn_per_pixel_vector2=np.array(sn_per_pixel2)
  sn_per_pixel_vector3=np.array(sn_per_pixel3)
  area_vector=np.array(area)
  #print(len(sn_per_pixel_vector),len(sn_per_pixel_vector2),len(sn_per_pixel_vector3),len(area_vector))
  # ginim20a_vector=np.array(GiniM20_1)
  #ginim20b_vector=np.array(GiniM20_2)
  asymmetry_vector1=np.array(asymmetry_all1)
  asymmetry_vector2=np.array(asymmetry_all2)
  concentration_vector=np.array(concentration_all)
  #print(len(asymmetry_vector1),len(asymmetry_vector2),len(concentration_vector))
  smoothness_vector=np.array(smoothness_all)
  #print(len(smoothness_vector))
  #merger_vector=np.array(merger_all)
  ellipticity_vector=np.array(ellipticity_all)
  elongation_vector=np.array(elongation_all)
  rpetro_vector=np.array(rpetro_all)
  rhalf_vector=np.array(rhalf_all)
  shapeasy_vector=np.array(shape_asymmetry_all)
  centerax_all_vector=np.array(centerax_all)
  centeray_all_vector=np.array(centeray_all)
  #print(len(shapeasy_vector))
  mylist=np.array([int(IDgal)]) # Per la colonna degli ID mi serve
  # listarea=np.array(areagood)
  # listdist=np.array(distancegood)
  # listphoto=np.array(photometry)
  # listpeak=np.array(peakgood)
  #print len(clumpiness),len(listarea),len(listdist),len(listphoto),len(listpeak)
  # if (silent != True):
  #  print 'clumpiness 1 e 2, area, dist, photo, peak =',clumpvector,clumpvector2,listarea, listdist, listphoto, listpeak
  # galaxy_table=np.hstack((mylist,clumpvector,clumpvector2,listarea,listdist,listphoto,listpeak)) #,nc_nodeblvec,nclumpvector))
  if (silent != True): print('')

  #print(len(mylist))
  #print(len(IDvector))
  #print(len(clumpvector))
  #quit()
  print('len several vectors =')
  print(len(IDvector))
  print(len(photoband_vector))
  print(len(clumpvector))
  print(len(smoothness_SM))
  print(len(sn_per_pixel_vector2))
  print(len(rpetro_vecSM))
  print(len(sn_per_pixel_vecSM))
  print(len(centeray_all_vector))

  galaxy_table=np.hstack((IDvector,photoband_vector,clumpvector,clumpvector2,Rmax_vector,gini_vector,m20_vector,concentration_vector,asymmetry_vector1,centerax_all_vector,centeray_all_vector,shapeasy_vector,smoothness_vector,sn_per_pixel_vector2,sn_per_pixel_vector3,area_vector,flag_vecSM,gini_vecSM,m20_vecSM,concentration_vecSM,asymmetry_vecSM,smoothness_vecSM,sn_per_pixel_vecSM,rpetro_vecSM)) # ,ginim20a_vector,distlotz_vector)) #,nc_nodeblvec,nclumpvector))
  #galaxy_table=np.ravel(galaxy_table)
  print('shape galaxy table =',galaxy_table.shape)
  
  ##print('galaxy table =',galaxy_table)
  ##assert(len(galaxy_table)==22)
  #final_table=np.append(final_table,galaxy_table)
  ##print('len final table =',len(final_table))
  ##time.sleep(15)

  #########################################################################
  # FINE storing IDgal, clumpiness and clumps number into a table         #
  #########################################################################
  
  #print('OK JOOO')
  #quit()
  print('OK GOOOD\n')
  #quit()
    
  ############################################
  # Save table with clumpiness values :      #
  ############################################
  if (test_segmap_index==True):
    print('Stop here - No table to save')
    sys.exit()

  print('Save final table')
  print('counter989 =',counter989)
  inputlist=np.array(['# IDgal','band','c1','c2','Rmax','gini','m20','C','A1','Axc','Ayc','shapeasy','S','SNpixel1','SNpixel2','Area','SM_flag','SM_gini','SM_m20','SM_C','SM_A','SM_S','SM_SNpixel','SM_Rpetro\n'])   #   # VEDI A INIZIO PROGRAMMA
  #final_table=galaxy_table.reshape((counter989,-1))
  final_table=galaxy_table.reshape((-1,counter989))
  final_table=final_table.transpose()
  if (saveoutput==True) :
    with open(filenameoutput, "w") as myfile :
     for itN in inputlist :
      myfile.write("%s " % itN)
     np.savetxt(myfile, final_table ,delimiter='   ', fmt='%s')
  print('Stored information in table =',filenameoutput)
  print('--- END ---\n\n')
  # END OF THIS LOOP
  time.sleep(5)
  #quit()
  #print('asse x completo per visualizzazione =')
  #print(asse_x_completo)
