#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################
__author__ = "Antonello Calabro"
__affiliation__= "INAF - Osservatorio Astronomico di Roma"
#__copyright__ = "Copyright 2007, The Cogent Project"
#__credits__ = []
__license__ = "GPL"
#__version__ = "1.0.1"
__maintainer__ = "Antonello Calabro"
__email__ = "antonello.calabro@inaf.it"
__status__ = "Production"
##########################################

# OPEN TABLE FILE WITH RESULTS

# GINI PARAMETER
# M20 PARAMETER 
# CONCENTRATION
# SHAPE ASYMMETRY
# SMOOTHNESS
# CLUMPINESS

import os,sys,time
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator
# Import necessary modules

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy,math,glob,random
import pandas as pd

from gauss_fit_routines import gaussian,moments,fitgaussian

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
redshift=6
_conv_kpc_to_arcsec=cosmo.arcsec_per_kpc_proper(redshift)
_conv_arcsec_to_kpc=1./_conv_kpc_to_arcsec
conv_kpc_to_arcsec=_conv_kpc_to_arcsec.value
conv_arcsec_to_kpc=_conv_arcsec_to_kpc.value
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,FixedLocator
from scipy import signal,ndimage,stats

from auxiliary import import_table,region,create_circular_mask,fit_gauss,replace_secondary_sources_with_bkg,calc_bg,rmsdiff,cartesian_product,makeGaussian,petrosian_radius,binary_detection_pawlik,background_pawlik,RMAX_calc,smoothed_mask_skysub,asymmetry_function,asymmetry_function2,asymmetry_simple,Cutout2D_mine,minimize_mu,asymmetry_verysimple_4bkg,M20_simple,asymmetry_function2_shape,asymmetry_function_bkg,bkg_search,asymmetry_bkg_simple
from auxiliary import make_segmentation

sleep_time=0.01
cwd=os.getcwd()+'/'

# __*****************************************************************__

calculate_parameters=True # Default is True
make_images=False 
segmap_detection_type='automatic'  # automatic or indices (in the second case you have to write all the indices of the segmentation region corresponding to the galaxy itself)
test_segmap_index=False  # Put True if you have to select the segmentation region of the galaxies


# IMPORTA TABELLA CON I FILE DA PROCESSARE
# START PROGRAM
# LOOP ON SIMULATIONS
# OPEN PHOTOMETRIC CATALOG
# START LOOP ON IMAGES 
# PLOT IMAGES
# Define segmentation images and masks
# DERIVE SEGMENTATION IMAGE
# NOW GOING TO USE STATMORPH TO CALCULATE MORPHOLOGICAL PARAMETERS
# CONTINUATION OF CLUMPINESS CALCULATION
# Save table with clumpiness values

which_subset='subset_EOR'    # it's just a name that you give to distinguish among different subsets
segmentation_type='Pawlik'  #  'square', 'photutils', 'Pawlik'      # 'petrosian_based' or 'photutils' # This is the segmentation for the gini calculation
size_square_box_small=0.4 ; size_square_box_large=0.7  # arcsec
skybox_pixel=15 ; skybox_arcsec=0.5 ; size_psf=5
pixel_scale_SW=0.031 # arcsec/pixel
pixel_scale_LW=0.063 # arcsec/pixel

# 0) INITIAL PARAMETERS CLUMPINESS
# CLUMPINESS measurement
# DEFINISCO LISTA immagini su cui calcolare la clumpiness
# 3b) SMOOTHED IMAGE AND ORIGINAL-SMOOTHED
# 3c) Calcolo la clumpiness (Conselice 2003/4)
# 3c) Calcolo la clumpiness (My method)  # (PSF gauss fit visualization)
# 3d) Segmentation map on clumps
# 3e) visualization clumps e segmentation map on clumps
# 3f) Analyze clumps properties  
# Use STATMORPH on original image 
# Use STATMORPH on single galaxy image
# INIZIO SECONDA PARTE CLUMPINESS CALCULATION
# CONTINUATION OF CLUMPINESS CALCULATION
# CALCULATE ALL ASYMMETRIES
# Calculate M20 manually
# Calculate concentration manually


# OPEN PHOTOMETRIC CATALOG : 

input_folder=cwd+"data/"
file_catalog=input_folder+"assembled_catalog.cat"
table_catalog=np.genfromtxt(file_catalog,names=True,dtype=None)
IDobj_all=table_catalog['ID'] ; ra_all=table_catalog['RA'] ; dec_all=table_catalog['DEC']
if which_subset=='subset_EOR' :
    ID_selected=[65,76,79,80,91,104,113,125,129] # 65,76,79,80,91,104,113,125,129]
    IDlista=np.array(ID_selected)

# OPEN GALAXY INDICES SEGMENTATION MAPS : 
IN_catalog=input_folder+"galaxy_segmap_indices.txt"
indices_catalog=np.genfromtxt(IN_catalog,names=True,dtype=None)

output_folder=cwd+'results/'
output_folder_segmap=cwd+'segmaps/'
filenameoutput=output_folder+"results_subset_"+which_subset+".txt"

bande_all=np.array(['f090','f115','f150','f200','f277','f356','f444'])
IDpixscale=[pixel_scale_SW,pixel_scale_SW,pixel_scale_SW,pixel_scale_SW,pixel_scale_LW,pixel_scale_LW,pixel_scale_LW]
IDfwhm_banda=[0.034,0.040,0.050,0.066,0.091,0.115,0.145] # Taken from https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-predicted-performance/nircam-point-spread-functions






# OPEN TABLE FILE WITH RESULTS : 

column_names=np.array(['# IDgal','ID2','band','c1','c2','Rmax','gini','m20','C','A1','shapeasy','A2','S','ellipticity','elongation','rpetro','rhalf','SNpixel','SNpixel2','SNpixel3','Area\n']) # VEDI A INIZIO PROGRAMMA
results_table=np.genfromtxt(filenameoutput,names=True,dtype=None)

IDgal_all=results_table['IDgal']
bands_all=results_table['band']
clumpiness=results_table['c1']
gini=results_table['gini']
m20=results_table['m20']
C=results_table['C']
S=results_table['S']
shapeasy=results_table['shapeasy']

# This should be copied from simulations results :
mag=np.array([26,27,28])
uncertainty_gini=np.array([0.05,0.05,0.05])
uncertainty_m20= np.array([0.08,0.11,0.16])
uncertainty_C=   np.array([0.15,0.3,0.4])
uncertainty_shA= np.array([0.07,0.09,0.11])
uncertainty_S=   np.array([0.04,0.06,0.07])
uncertainty_cl=  np.array([0.02,0.04,0.07])





for IDgal in ID_selected :
  print('\n\nStarting analyzing galaxy '+str(IDgal))
  
  selection=np.where(IDgal_all==IDgal)[0]

  sel_clumpiness=clumpiness[selection]
  sel_gini=gini[selection]
  sel_m20=m20[selection]
  sel_C=C[selection]
  sel_S=S[selection]
  sel_shapeasy=shapeasy[selection]
  
  # Styling
  plt.style.use("seaborn-darkgrid")
  plt.rcParams["font.family"] = "Avenir"
  plt.rcParams["font.size"] = 16

  # Create Blank Figure
  fig = plt.figure(figsize=(10, 12))
  # Create 4x4 Grid
  gs = fig.add_gridspec(nrows=3, ncols=2) # , height_ratios=[2, 1], width_ratios=[2, 1])

  # Create Three Axes Objects
  ax1 = fig.add_subplot(gs[0, 0])
  ax2 = fig.add_subplot(gs[0, 1])
  ax3 = fig.add_subplot(gs[1, 0])
  ax4 = fig.add_subplot(gs[1, 1])
  ax5 = fig.add_subplot(gs[2, 0])
  ax6 = fig.add_subplot(gs[2, 1])
  
  # Plot Data

  # GINI PARAMETER 
  condition= (sel_gini>=0) & (sel_gini<=1)
  ax1.scatter(bande_all[condition],  sel_gini[condition],marker='s',facecolor='blue',edgecolor='k',lw=1,s=200)
  ax1.errorbar(bande_all[condition], sel_gini[condition],yerr=uncertainty_gini[1],markersize=0,marker='s',color='blue',capsize=5,lw=1,ls='None')
  
  for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(18) 
  for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)
  majY=0.2   ;  minY=0.05
  majorLocatorY = MultipleLocator(majY)
  minorLocatorY = MultipleLocator(minY)    
  ax1.yaxis.set_major_locator(majorLocatorY)
  ax1.yaxis.set_minor_locator(minorLocatorY)
  ax1.tick_params(axis='both',bottom=True,top=True,right=True,which='both')
  ax1.set_xlabel(' ',size=0) 
  ax1.set_ylabel('gini',size=21.5)
  #ax1.legend(loc=legplace,prop={'size':legsize})



  # M20 PARAMETER  
  condition= (sel_m20<=0) & (sel_m20>=-3)
  ax2.scatter(bande_all[condition],  sel_m20[condition],marker='s',facecolor='blue',edgecolor='k',lw=1,s=200)
  ax2.errorbar(bande_all[condition], sel_m20[condition],yerr=uncertainty_m20[1],markersize=0,marker='s',color='blue',capsize=5,lw=1,ls='None')
  
  for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(18) 
  for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)
  majY=0.2   ;  minY=0.05
  majorLocatorY = MultipleLocator(majY)
  minorLocatorY = MultipleLocator(minY)    
  ax2.yaxis.set_major_locator(majorLocatorY)
  ax2.yaxis.set_minor_locator(minorLocatorY)
  ax2.tick_params(axis='both',bottom=True,top=True,right=True,which='both')
  ax2.set_xlabel(' ',size=0) 
  ax2.set_ylabel('m20',size=21.5)


  # CONCENTRATION
  condition= (sel_C>=0) & (sel_C<=5)
  ax3.scatter(bande_all[condition],  sel_C[condition],marker='s',facecolor='blue',edgecolor='k',lw=1,s=200)
  ax3.errorbar(bande_all[condition], sel_C[condition],yerr=uncertainty_C[1],markersize=0,marker='s',color='blue',capsize=5,lw=1,ls='None')
  
  for tick in ax3.xaxis.get_major_ticks():
    tick.label.set_fontsize(18) 
  for tick in ax3.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)
  majY=0.5   ;  minY=0.1
  majorLocatorY = MultipleLocator(majY)
  minorLocatorY = MultipleLocator(minY)    
  ax3.yaxis.set_major_locator(majorLocatorY)
  ax3.yaxis.set_minor_locator(minorLocatorY)
  ax3.tick_params(axis='both',bottom=True,top=True,right=True,which='both')
  ax3.set_xlabel(' ',size=0) 
  ax3.set_ylabel('concentration',size=21.5)
  


  # SHAPE ASYMMETRY
  condition= (sel_shapeasy>=-4) & (sel_shapeasy<=5)
  ax4.scatter(bande_all[condition],  sel_shapeasy[condition],marker='s',facecolor='blue',edgecolor='k',lw=1,s=200)
  ax4.errorbar(bande_all[condition], sel_shapeasy[condition],yerr=uncertainty_shA[1],markersize=0,marker='s',color='blue',capsize=5,lw=1,ls='None')
  
  for tick in ax4.xaxis.get_major_ticks():
    tick.label.set_fontsize(18) 
  for tick in ax4.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)
  majY=0.2   ;  minY=0.05
  majorLocatorY = MultipleLocator(majY)
  minorLocatorY = MultipleLocator(minY)    
  ax4.yaxis.set_major_locator(majorLocatorY)
  ax4.yaxis.set_minor_locator(minorLocatorY)
  ax4.tick_params(axis='both',bottom=True,top=True,right=True,which='both')
  ax4.set_xlabel(' ',size=0) 
  ax4.set_ylabel('shape asymmetry',size=21.5)



  # SMOOTHNESS
  condition= (sel_S>=0) & (sel_S<=1)
  ax5.scatter(bande_all[condition],  sel_S[condition],marker='s',facecolor='blue',edgecolor='k',lw=1,s=200)
  ax5.errorbar(bande_all[condition], sel_S[condition],yerr=uncertainty_S[1],markersize=0,marker='s',color='blue',capsize=5,lw=1,ls='None')
  
  for tick in ax5.xaxis.get_major_ticks():
    tick.label.set_fontsize(18) 
  for tick in ax5.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)
  majY=0.2   ;  minY=0.05
  majorLocatorY = MultipleLocator(majY)
  minorLocatorY = MultipleLocator(minY)    
  ax5.yaxis.set_major_locator(majorLocatorY)
  ax5.yaxis.set_minor_locator(minorLocatorY)
  ax5.tick_params(axis='both',bottom=True,top=True,right=True,which='both')
  ax5.set_xlabel(' ',size=0) 
  ax5.set_ylabel('smoothness',size=21.5)



  # CLUMPINESS
  condition= (sel_clumpiness>=0) & (sel_clumpiness<=1)
  ax6.scatter(bande_all[condition],  sel_clumpiness[condition],marker='s',facecolor='blue',edgecolor='k',lw=1,s=200)
  ax6.errorbar(bande_all[condition], sel_clumpiness[condition],yerr=uncertainty_cl[1],markersize=0,marker='s',color='blue',capsize=5,lw=1,ls='None')
  
  for tick in ax6.xaxis.get_major_ticks():
    tick.label.set_fontsize(18) 
  for tick in ax6.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)
  majY=0.2   ;  minY=0.05
  majorLocatorY = MultipleLocator(majY)
  minorLocatorY = MultipleLocator(minY)    
  ax6.yaxis.set_major_locator(majorLocatorY)
  ax6.yaxis.set_minor_locator(minorLocatorY)
  ax6.tick_params(axis='both',bottom=True,top=True,right=True,which='both')
  ax6.set_xlabel(' ',size=0) 
  ax6.set_ylabel('clumpiness',size=21.5)



  plt.suptitle('ID galaxy = '+str(IDgal),fontsize=25,weight='bold')
  plt.tight_layout()
  savefigname=output_folder+'IDgal_'+str(IDgal)+'_summary.png'
  plt.savefig(savefigname,dpi=80,fomrat='png')
  plt.close(fig)
































print('Stop here')
sys.exit()
# ---------------------------------------------------------------------

stop_here=False
clumps_method = 0   # 0=standard (Conselice 2003/4)
npixels = 5  # minimum number of connected pixels
npixels_original = 5

# PROPERTIES :
calculate_area_SNpixel=True # E' importante per capire se i risultati possono essere affidabili o no !!
replace_secondary_sources=True
calculate_Rmax=True
calculate_smoothness=True
calculate_asymmetry_center=True # contains also asymmetry calculation
calculate_shape_asymmetry_easier=True
calculate_concentration=True
calculate_gini=True
calculate_m20=True
# def asymmetry_function2
# Calculate Shape Asymmetry
# Calculate Asymmetry
calculate_shape_asymmetry_old=False
calculate_asymmetry_oldstyle=False # Di default perche' l'asimmetria e' calcolata con calculate_asymmetry_center
calculate_m20_oldstyle=False

show_figures=False
show_figures_Gini=False
show_figures_clumps=True             # True
show_image_shape_asymmetry=False
show_images_masks_smoothness=False   # True
plot_asymmetry=False


# 3g) Identify and remove nucleus
# 3h) Clumpiness derivation
# 3i) Derivation of clumpiness, ellipticity and elongation for single (central) galaxies in merging system / pair
# Calculate asymmetry manually
# Calculate smoothness manually
# 3iB) CLUMPINESS derivation per single galaxies case
# 3k) Visualization clumps
# 3l) Per salvare dei fits files
# 3m) Store clumpiness and clump properties in table (and then save)

# ---------------- ---------------- ---------------- -------------- -------------
show_all=True ; cutstamp_new=True ; #  maschera=True # Non metto una maschera quando calcolo la clumpiness
firstsegmap_snr=3   ;   firstsegmap_npixels=1    
firstsegmap_snr_single=2  ;  firstsegmap_npixels_single=50

####--------------------------------------------------------------------------------------------