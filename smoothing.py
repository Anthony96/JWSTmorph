from scipy import signal
from scipy import ndimage
from skimage.morphology import disk
import numpy as np  
from astropy.convolution import Box2DKernel,Gaussian2DKernel,convolve

def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    import numpy as np

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    
    return np.round(np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2),5)

#SIZE=30
#FWHM=15
#tableG=makeGaussian(SIZE,fwhm=FWHM)

def gkern(kernlen=21, nsig=3):
    """Returns a 2D Gaussian kernel array."""

    interval = (2*nsig+1.)/(kernlen)
    x = np.linspace(-nsig-interval/2., nsig+interval/2., kernlen+1)
    kern1d = np.diff(st.norm.cdf(x))
    kernel_raw = np.sqrt(np.outer(kern1d, kern1d))
    kernel = kernel_raw/kernel_raw.sum()
    return kernel

def gkern2(kernlen=21, nsig=3):
    """Returns a 2D Gaussian kernel array."""

    # create nxn zeros
    inp = np.zeros((kernlen, kernlen))
    # set element at the middle to one, a dirac delta
    inp[kernlen//2, kernlen//2] = 1
    # gaussian-smooth the dirac, resulting in a gaussian filter mask
    return fi.gaussian_filter(inp, nsig)

def makeGaussian0(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    
    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


def gauss_kern1(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = scipy.mgrid[-size:size+1, -sizey:sizey+1]
    g = scipy.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.max()

def gauss_kern2(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

# ALTRO MODO PER CREARE UN KERNEL GAUSSIANO :::
#from astropy.convolution import Gaussian2DKernel
#sigma = 3.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))   # FWHM = 3
#kernel = Gaussian2DKernel(sigma, x_size=5, y_size=5)
#kernel.normalize()


def smooth(imagein,smoothtype,smoothRadius,RR) :
  selem = disk(RR)
  # CHOOSE THE SMOOTHED IMAGE
  if (smoothtype=='gauss'):
    # CON ndimage :
    img_gauss = ndimage.filters.gaussian_filter(imagein, smoothRadius, mode='reflect')  
    # mode : {reflect, constant, nearest, mirror, wrap}, #prima avevo messo reflect
    imagesmooth=img_gauss
  elif (smoothtype=='uniform'): 
    img_uniform = ndimage.filters.uniform_filter(imagein, smoothRadius)   #size=11 
    imagesmooth=img_uniform
  elif (smoothtype=='median'):
    # MEDIAN FILTER (should better preserve the edges)
    img_median = ndimage.median_filter(imagein, smoothRadius)
    imagesmooth=img_median
  elif (smoothtype=='percentile'):
    # SMOOTHING WITH SCIKIT-IMAGE
    percentile_result = rank.mean_percentile(imagein, selem=selem, p0=.05, p1=.95)
    percentile_result= percentile_result*np.sum(imagein)/np.sum(percentile_result)
    imagesmooth=percentile_result
  elif (smoothtype=='bilateral'):
    bilateral_result = rank.mean_bilateral(imagein, selem=selem, s0=500, s1=500)
    bilateral_result= bilateral_result*np.sum(imagein)/np.sum(bilateral_result)
    imagesmooth=bilateral_result
  elif (smoothtype=='normal'):
    normal_result = rank.mean(imagein, selem=selem)
    normal_result= normal_result*np.sum(imagein)/np.sum(normal_result)
    imagesmooth=normal_result
  elif (smoothtype=='convolve'):
    # Other methods for smoothing : # KERNEL AND CONVOLVE 2D
    # make some kind of kernel, there are many ways to do this...
    dimk=10
    t = 1 - np.abs(np.linspace(-1, 1, dimk))
    kernel = t.reshape(dimk, 1) * t.reshape(1, dimk)
    kernel /= kernel.sum()   # kernel should sum to 1!  :)   #print kernel
    # convolve 2d the kernel 
    img_gaus2 = signal.convolve2d(imagein, kernel, mode='same')
    imagesmooth=img_gaus2
  elif (smoothtype=='boxcar'):
    box_2D_kernel = Box2DKernel(boxsize)
    smoothed_data_boxcar = convolve(imagein, box_2D_kernel)
    imagesmooth = smoothed_data_boxcar
  return imagesmooth

# Gaussian filter 
'''
input : array_like
Input array to filter.
sigma : scalar or sequence of scalars
Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.
order : {0, 1, 2, 3} or sequence from same set, optional
The order of the filter along each axis is given as a sequence of integers, or as a single number. An order of 0 corresponds to convolution with a Gaussian kernel. An order of 1, 2, or 3 corresponds to convolution with the first, second or third derivatives of a Gaussian. Higher order derivatives are not implemented
output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {reflect, constant, nearest, mirror, wrap}, optional
The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to constant. Default is reflect
cval : scalar, optional
Value to fill past edges of input if mode is constant. Default is 0.0
truncate : float
Truncate the filter at this many standard deviations. Default is 4.0.
Returns:  
gaussian_filter : ndarray
Returned array of same shape as input.
'''  