# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:09:17 2020

@author: Harry Anthony
"""

from astropy.io import fits
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import random
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import scipy.ndimage as snd
from photutils import CircularAnnulus, fit_2dgaussian, CircularAperture, aperture_photometry, EllipticalAperture, EllipticalAnnulus, ApertureMask, data_properties
from matplotlib.colors import LogNorm
from astropy.wcs import WCS
from scipy import stats  
from PSF_methods import *

from photutils import Background2D, MedianBackground
from photutils import detect_threshold, detect_sources
from photutils import deblend_sources
from photutils import detect_sources
from photutils import centroid_sources, centroid_com
from photutils import source_properties

def Show_image(Data):
    """
    The method displays the image data in a 2D plot.
    
    Parameters:
        Data: 2D array of pixel data from the CCD        
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
                
    Data_log = np.log(np.log(np.log(np.log(Data)))) #Logarithm put on data to make stars more visible
    norm = ImageNormalize(stretch=SqrtStretch()) #Normalisation for the image
    plt.imshow(Data_log[::-1],'Greys_r',norm=norm)


def Locate_galaxies(Data,N_std,mask=None):
    global Image_no_noise, centres, Image_no_noise_watershed
    """
    The method locates galaxies in data above the threshold (mean + Nstd standard
    deviations) and returns their centres.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        N_std: Number of standard deviations above the mean the threshold value is   
        mask: 2D boolean array to mask 'bad pixels'
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
    if isinstance(N_std,(int,float)) is not True:
                raise Exception('N_std must be a positive integer or a float')
    if N_std <= 0:
                raise Exception('N_std must be a positive integer or a float')

    #Threshold from which below is considered noise
    bkg_estimator = MedianBackground()
    bkg = Background2D(Data[::-1], (50, 50), filter_size=(3, 3),bkg_estimator=bkg_estimator)
    #The threshold depends on the number of standard deviations
    threshold = bkg.background + (N_std * bkg.background_rms)

    #Isolates the shapes of the galaxies using a segmentation algorithm
    Image_no_noise = detect_sources(Data[::-1], threshold, npixels=6,mask=mask[::-1])

    #Close by galaxies are resolved using a watershed algorithm
    Image_no_noise_watershed = deblend_sources(Data[::-1], Image_no_noise, npixels=6,
                                nlevels=12,contrast=0.001)
    
    #The centres of the galaxies are then estimated
    source_prop = source_properties(Data[::-1], Image_no_noise_watershed)

    #First an initial guess of the position is made
    x_init = [obt.xcentroid.value for obt in source_prop]
    y_init = [obt.ycentroid.value for obt in source_prop]

    #The centre of mass is then calculated usingn moments
    x, y = centroid_sources(Data[::-1], x_init, y_init, box_size=21,
                        centroid_func=centroid_com)  
    
    #Append centres of the galaxies to array centres
    centres = [[y[i],x[i]] for i in range(0,len(x))]
    
    return centres
    
    
def Image_no_noise_array(Data,N_std,mask=None):
    """
    The method shows the shapes of the isolated galaxies calculated using the algorithm.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        N_std: Number of standard deviations above the mean the threshold value is   
        mask: 2D boolean array to mask 'bad pixels'
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
    if isinstance(N_std,(int,float)) is not True:
                raise Exception('N_std must be a positive integer or a float')
    if N_std <= 0:
                raise Exception('N_std must be a positive integer or a float')
                
    Locate_galaxies(Data,N_std,mask)
    cmap = Image_no_noise.make_cmap(seed=123)
    plt.imshow(Image_no_noise.data,cmap=cmap)
    plt.show()
    return Image_no_noise


def Image_no_noise_watershed_array(Data,N_std,mask=None):
    """
    The method shows the shapes of the isolated galaxies calculated using the 
    watershed algorithm.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        N_std: Number of standard deviations above the mean the threshold value is   
        mask: 2D boolean array to mask 'bad pixels'
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
    if isinstance(N_std,(int,float)) is not True:
                raise Exception('N_std must be a positive integer or a float')
    if N_std <= 0:
                raise Exception('N_std must be a positive integer or a float')
                
    Locate_galaxies(Data,N_std,mask)
    cmap = Image_no_noise_watershed.make_cmap(seed=123)
    plt.imshow(Image_no_noise_watershed.data,cmap=cmap)
    plt.show()
    return Image_no_noise_watershed
    
def Show_centres(Data,N_std,mask=None):
    global Image_no_noise, centres
    """
    This method plots the measured centres of the galaxies against the image
    of the galaxy.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        N_std: Number of standard deviations above the mean the threshold value is 
        mask: 2D boolean array to mask 'bad pixels'
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
    if isinstance(N_std,(int,float)) is not True:
                raise Exception('N_std must be a positive integer or a float')
    if N_std <= 0:
                raise Exception('N_std must be a positive integer or a float')
    
    Locate_galaxies(Data,N_std,mask)
    Show_image(Data)
    for x in range(0,len(centres)):
        plt.plot(centres[x][1],centres[x][0],marker='+',color='cyan')
        
        
def Centre_on_galaxy(Data,c,Data_error):
    global Galaxy_centred_Data, Masking, galaxy_size, centres, Error, Image_no_noise
    """
    Crops the image to a 300x300 image with galaxy with centre index c at its
    centre (150,150).
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        c: Centre value index
        Data_error: 2D array of the error on pixel
        
    Returns:
        Galaxy_centred_data: A 2D array of pixel values with galaxy with
        centre index c at its centre.
        Masking: A 2D boolean array of the masked pixels.
        Error: A 2D array of the errors associated with each pixel.
    """
    if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
    if isinstance(Data_error,(list,np.ndarray)) is not True:
            raise Exception("Data_error must be an array")
    if len(Data) != len(Data_error):
            raise Exception("Data and Data_error must be the same shape")
    if isinstance(c,(int)) is not True:
                raise Exception('c must be a positive integer')
    if c < 0:
                raise Exception('c must be a positive integer')
    
    #Crops the image to 300x300 pixels
    Galaxy_centred_Data = Data[::-1][int(centres[c][0]-150):int(centres[c][0
                              ]+150)][:,int(centres[c][1]-150):int(centres[c][1]+150)]

    #This code adds to the mask array where another galaxy is located which is not
    #being studied.
    Masking = np.zeros(Galaxy_centred_Data.shape,dtype=bool) #Create a 2D boolean array
    #All stars with a difference centre index are masked
    Masking = Image_no_noise_watershed.data[int(centres[c][0]-150):int(centres[c][0]+150)][:,int(
            centres[c][1]-150):int(centres[c][1]+150)] != c+1
    #Galaxy size is calculated by taking the masked pixels away from the total image
    galaxy_size = Galaxy_centred_Data.shape[0]*Galaxy_centred_Data.shape[1] - Masking.sum()
    
    Error = np.zeros(Galaxy_centred_Data.shape,dtype=bool) #Create a 2D error array
    Error = Data_error[int(centres[c][0]-150):int(centres[c][0]+150)][:,int(
            centres[c][1]-150):int(centres[c][1]+150)] != c+1
    
    return Galaxy_centred_Data, Masking, galaxy_size, Error
    

def Calculate_2Dstd(Data,Masking):
        """
        Calculates the radial standard deviation of a galaxy by fitting a 
        2D Gaussian kernel to the data.
        
        Parameters:
            Data: 2D array of pixel data from the CCD
            Masking: 2D boolean array of the masked pixels
            
        Returns:
            The 2D standard deviation of the galaxy
        """
        if isinstance(Data,(list,np.ndarray)) is not True:
            raise Exception("Data must be an array")
        if isinstance(Masking,(list,np.ndarray)) is not True:
            raise Exception("Masking must be an array")
        if len(Data) != len(Masking):
            raise Exception("Data and Masking must be the same shape")
            
        #Applys a 2D Gaussian fit to the data
        Gauss_fit = fit_2dgaussian(Data,mask=Masking)
        #Calculate the radial standard deviation
        std = np.sqrt(Gauss_fit.y_stddev.value**2 + Gauss_fit.x_stddev.value**2)
        if std > 30: #Prevents stds being too large
            std = 29
        return std
    
    
def Pixel_rejection_annulus(N_std):    
    """
    Creates a mask for pixels in the annulus which are above the threshold of
    the mean of pixels inside the annulus plus a given number of standard deviations
    (N_std).
        
    Parameters:
        N_std: Number of standard deviations above the mean the threshold value is 
        
    Returns:
        A Boolean mask array for the annulus 
    """
    if isinstance(N_std,(int,float)) is not True:
        raise Exception('N_std must be a positive integer or a float')
    if N_std <= 0:
        raise Exception('N_std must be a positive integer or a float')
                
    #The shape of the annulus is derived
    Annulus_mask = aperture_and_annulus[1].to_mask(method='center')
    Annulus_shape = Annulus_mask.to_image(shape=((300, 300)))
    
    #The pixels in the annulus are then isolated
    Pixels_in_annulus = Galaxy_centred_Data*Annulus_shape
    
    #If the pixel values are above the threshold, they are masked
    Annulus_mask = Pixels_in_annulus > np.mean(Pixels_in_annulus[
                Pixels_in_annulus != 0]) + N_std*np.std(Pixels_in_annulus[
                        Pixels_in_annulus != 0])
    
    return Annulus_mask


def Counts_single_galaxy_circular_aperture(Data,c,N_std,Data_error,visualise=False):
    global count_list, aperture_and_annulus, aperture_measurement, centres
    """
    Calculates the number of counts of galaxy with centre index c using a circular
    aperture and annulus.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        c: Centre value index
        N_std: Number of standard deviations above the mean the threshold value is   
        Data_error: 2D array of the error on pixel
        Visualise: Boolean value which will make the function plot the aperture on
        top of the image of the galaxy if true.
        
    Returns:
        Counts of the galaxy
        Local background of the galaxy
        Error associated with counts of the galaxy
    """
    if isinstance(visualise,(bool)) is not True:
        raise Exception('visualise must be a boolean')
    
    Galaxy_centred_Data, Masking, galaxy_size, Error = Centre_on_galaxy(Data,c,
                                                                Data_error)
    
    if galaxy_size > 6:  #Galaxies must be greater than 6 pixels to ensure it is not noise
        std = Calculate_2Dstd(Galaxy_centred_Data,Masking)
        ##Create the aperature and annulus for the star
        aperture_and_annulus = [CircularAperture((150,150), r=2*std),
                            CircularAnnulus((150,150), r_in=3*std, r_out=5*std)]        
        Annulus_mask = Pixel_rejection_annulus(N_std)        
        #Measure the counts in the aperature and annulus
        aperture_measurement = aperture_photometry(Galaxy_centred_Data, aperture_and_annulus,
                                       subpixels=10,mask=Annulus_mask,error=Error)    
        
        #Counts and background are calculated
        counts_with_background = aperture_measurement[0][3]
        background = aperture_measurement[0][5]/aperture_and_annulus[1].area
        counts_error = aperture_measurement[0][4]
        background_error = aperture_measurement[0][6]

        if visualise == True:
                #Plot the aperature and annulus
                fig, axes = plt.subplots()
                Data_log = np.log(Galaxy_centred_Data) #Logarithm put on data to make stars more visible
                norm = ImageNormalize(stretch=SqrtStretch()) #Normalisation for the image
                axes.imshow(Data_log,'Greys_r',norm=norm)
                aperture_and_annulus[0].plot(color='Chartreuse')
                aperture_and_annulus[1].plot(color='Crimson')
            
        #Galaxy_counts = counts_with_background - background * aperture_and_annulus[0].area
        return (counts_with_background - background * aperture_and_annulus[0
                    ].area)/(aperture_and_annulus[0
                    ].area),counts_error + background_error * aperture_and_annulus[
                            0].area, background
    else:
        return 'Galaxy too small','Galaxy too small', 'Galaxy too small'
    
    
def Counts_single_galaxy_elliptical_aperture(Data,c,N_std,Data_error,visualise=False):
    global count_list, aperture_and_annulus, aperture_measurement, centres
    """
    Calculates the number of counts of galaxy with centre index c using an elliptical
    aperture and annulus.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        c: Centre value index
        N_std: Number of standard deviations above the mean the threshold value is   
        Data_error: 2D array of the error on pixel
        Visualise: Boolean value which will make the function plot the aperture on
        top of the image of the galaxy if true.
        
    Returns:
        Counts of the galaxy
        Local background of the galaxy
        Error associated with counts of the galaxy
    """
    if isinstance(visualise,(bool)) is not True:
        raise Exception('visualise must be a boolean')
                
    Galaxy_centred_Data, Masking, galaxy_size, Error = Centre_on_galaxy(Data,c,
                                                                    Data_error)
    
    if galaxy_size > 6:  #Galaxies must be greater than 6 pixels to ensure it is not noise
        #Measure the features of the galaxy
        Galaxy_measurement = data_properties(Galaxy_centred_Data,mask=Masking)
        #Create the aperature and annulus for the star
        aperture_and_annulus = [EllipticalAperture((150,150),
                                Galaxy_measurement.semimajor_axis_sigma.value * 2 * np.sqrt(2*np.log(2))
                              ,Galaxy_measurement.semiminor_axis_sigma.value * 2 * np.sqrt(2*np.log(2))
                              ,theta=Galaxy_measurement.orientation.value), 
            EllipticalAnnulus((150,150), Galaxy_measurement.semimajor_axis_sigma.value * 3 * np.sqrt(2*np.log(2)),
                              Galaxy_measurement.semimajor_axis_sigma.value * 5 * np.sqrt(2*np.log(2)),
                              Galaxy_measurement.semiminor_axis_sigma.value * 5 * np.sqrt(2*np.log(2)),
                              theta=Galaxy_measurement.orientation.value)]    
            
        Annulus_mask = Pixel_rejection_annulus(N_std)  
        #Measure the counts in the aperature and annulus
        aperture_measurement = aperture_photometry(Galaxy_centred_Data, aperture_and_annulus,
                                       subpixels=10,mask=Annulus_mask,error=Error)    
        
        #Counts and background are calculated
        counts_with_background = aperture_measurement[0][3]
        background = aperture_measurement[0][5]/aperture_and_annulus[1].area
        counts_error = aperture_measurement[0][4]
        background_error = aperture_measurement[0][6]

        if visualise == True:
                #Plot the aperature and annulus
                fig, axes = plt.subplots()
                Data_log = np.log(Galaxy_centred_Data) #Logarithm put on data to make stars more visible
                norm = ImageNormalize(stretch=SqrtStretch()) #Normalisation for the image
                axes.imshow(Data_log,'Greys_r',norm=norm)
                aperture_and_annulus[0].plot(color='Chartreuse')
                aperture_and_annulus[1].plot(color='Crimson')

        return counts_with_background - background * aperture_and_annulus[0
                    ].area,counts_error + background_error * aperture_and_annulus[
                            0].area, background
    
    else:
        return 'Galaxy too small','Galaxy too small','Galaxy too small'
  
def Deconvolve_with_psf(Data,psf=None):
    """
    Deconvolves data with a psf using the functions from PSF_methods.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        psf: Psf to deconvolve the data with; it can have values 'None', 'Airy',
        'Gaussian' or 'Top_hat'.
        
    Returns:
        A 2D array containing the deconvolved data
    """
    print('yo')
    print(psf)
    print(str(psf))
    if str(psf) != 'None' and str(psf) != 'Airy' and str(psf) != 'Gaussian' and str(psf) !='Top_hat':
        raise Exception("psf can have values 'None', 'Airy','Gaussian','Top_hat'") 
        
    if psf == 'Airy':
        Deconvolved_data = deconvolve_airydisk(Data,3,10)
    if psf == 'Gaussian':
        Deconvolved_data = deconvolve_gaussian(Data,3,10)
    if psf == 'Top_hat':
        Deconvolved_data = deconvolve_tophat(Data,3,10)
    
    return Deconvolved_data
          
def Counts_all_galaxies_circular_aperture(Data,N_std,Data_error,mask=None,Deconvolve=False,psf=None):
    global centres, counts_list, centre_list, error_list, background_list, Image_no_noise
    """
    Calculates the number of counts of all the galaxies in the image with 
    circular apertures using aperture photometry.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        N_std: Number of standard deviations above the mean the threshold value is   
        Data_error: 2D array of the error on pixel
        mask: 2D boolean array to mask 'bad pixels'
        Deconvolve: Boolean value which will make the function deconvolve the
        image with a psf before locating the galaxy centres if true.
        psf: Psf to deconvolve the data with; it can have values 'None', 'Airy',
        'Gaussian' or 'Top_hat'.
        
    Returns:
        Counts_list: List of number of counts of each galaxy.
        Centre_list: List of centres of each galaxy.
        error_list: List of errors associated with the number of counts of each 
        galaxy.
        background_list: List of local background values for each galaxy.
    """
    if isinstance(Deconvolve,(bool)) is not True:
        raise Exception('Deconvolve must be a boolean')
    
    #Locate the centres of the galaxies
    if Deconvolve == True:
        centres = Locate_galaxies(Deconvolve_with_psf(Data,psf),N_std,mask)
    else:
        centres = Locate_galaxies(Data,N_std,mask)
    
    #Initialise arrays to input our results into
    counts_list = []
    centre_list = []
    error_list = []
    background_list = []
    
    #Loop through all centres
    for c in range(0,len(centres)):
        print(str(c/len(centres))+' % complete')
        galaxy_count, count_error, background = Counts_single_galaxy_circular_aperture(
                Data,c,N_std,Data_error,visualise=False)
        if galaxy_count != 'Galaxy too small':
            counts_list.append(galaxy_count)
            error_list.append(count_error)
            centre_list.append([centres[c][1],centres[c][0]])
            background_list.append(background)
    
    return counts_list,centre_list,error_list, background_list


def Counts_all_galaxies_elliptical_aperture(Data,N_std,Data_error,mask=None,Deconvolve=False,psf=None):
    global centres, counts_list, centre_list, error_list, background_list, Image_no_noise
    """
    Calculates the number of counts of all the galaxies in the image with 
    circular apertures using aperture photometry.
    
    Parameters:
        Data: 2D array of pixel data from the CCD
        N_std: Number of standard deviations above the mean the threshold value is   
        Data_error: 2D array of the error on pixel
        mask: 2D boolean array to mask 'bad pixels'
        Deconvolve: Boolean value which will make the function deconvolve the
        image with a psf before locating the galaxy centres if true.
        psf: Psf to deconvolve the data with; it can have values 'None', 'Airy',
        'Gaussian' or 'Top_hat'.
        
    Returns:
        Counts_list: List of number of counts of each galaxy.
        Centre_list: List of centres of each galaxy.
        error_list: List of errors associated with the number of counts of each 
        galaxy.
        background_list: List of local background values for each galaxy.
    """
    if isinstance(Deconvolve,(bool)) is not True:
        raise Exception('Deconvolve must be a boolean')
    
    #Calculate the centres of the galaxies
    if Deconvolve == True:
        centres = Locate_galaxies(Deconvolve_with_psf(Data,psf),N_std,mask)
    else:
        centres = Locate_galaxies(Data,N_std,mask)
    
    #Initialise arrays to input our results into
    counts_list = []
    centre_list = []
    error_list = []
    background_list = []
    for c in range(0,len(centres)):
        print(str(c/len(centres))+' % complete')
        galaxy_count, count_error, background = Counts_single_galaxy_elliptical_aperture(
                Data,c,N_std,Data_error,visualise=False)
        if galaxy_count != 'Galaxy too small':
            counts_list.append(galaxy_count)
            error_list.append(count_error)
            centre_list.append([centres[c][1],centres[c][0]])
            background_list.append(background)
    
    return counts_list,centre_list,error_list, background_list

