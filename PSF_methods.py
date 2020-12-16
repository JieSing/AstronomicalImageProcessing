# -*- coding: utf-8 -*-
from scipy import fftpack
import numpy as np
from astropy.convolution import Gaussian2DKernel, AiryDisk2DKernel, Tophat2DKernel
import skimage.restoration as rst

def deconvolve_airydisk(Data,Nstd=3,iterations=10):
    """
    Deconvolves the data with an airy disk of a given standard deviation using 
    Richardson-Lucy deconvolution.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        Nstd: Number of standard deviations of the airy disk
        iterations: Number of iterations of the Richardson–Lucy algorithm
    
    Returns:
        The deconvolved image (array) with the same shape 
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
    if isinstance(Data[0][0],(int,float)) is not True:
                raise Exception('Data values must be integers or floats')
    if isinstance(Nstd,(int,float)) is not True:
                raise Exception('Nstd must be an integer or a float')
    if isinstance(iterations,(int,float)) is not True:
                raise Exception('iterations must be an integer or a float')
                
    psf = AiryDisk2DKernel(Nstd)
    deconvolved_image = rst.richardson_lucy(Data, psf.array, iterations=iterations,clip=False)
    return deconvolved_image

def deconvolve_gaussian(Data,Nstd=3,iterations=10):
    """
    Deconvolves the data with a Gaussian of a given standard deviation using 
    Richardson-Lucy deconvolution.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        Nstd: Number of standard deviations of the airy disk
        iterations: Number of iterations of the Richardson–Lucy algorithm
    
    Returns:
        The deconvolved image (array) with the same shape 
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
    if isinstance(Data[0][0],(int,float)) is not True:
                raise Exception('Data values must be integers or floats')
    if isinstance(Nstd,(int,float)) is not True:
                raise Exception('Nstd must be an integer or a float')
    if isinstance(iterations,(int,float)) is not True:
                raise Exception('iterations must be an integer or a float')
                
    psf = Gaussian2DKernel(Nstd)
    deconvolved_image = rst.richardson_lucy(Data, psf.array, iterations=iterations,clip=False)
    return deconvolved_image

def deconvolve_tophat(Data,Nstd=3,iterations=10):
    """
    Deconvolves the data with a top hat kernel of a given standard deviation using 
    Richardson-Lucy deconvolution.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        Nstd: Number of standard deviations of the airy disk
        iterations: Number of iterations of the Richardson–Lucy algorithm
    
    Returns:
        The deconvolved image (array) with the same shape 
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
    if isinstance(Data[0][0],(int,float)) is not True:
                raise Exception('Data values must be integers or floats')
    if isinstance(Nstd,(int,float)) is not True:
                raise Exception('Nstd must be an integer or a float')
    if isinstance(iterations,(int,float)) is not True:
                raise Exception('iterations must be an integer or a float')
                
    psf = Tophat2DKernel(Nstd)
    deconvolved_image = rst.richardson_lucy(Data, psf.array, iterations=iterations,clip=False)
    return deconvolved_image
