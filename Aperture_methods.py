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




def Show_image(Data):
    """Plot the Image"""
    Data_log = np.log(np.log(np.log(np.log(Data))))                #Logarithm function on each data to make objects more visible
    norm = ImageNormalize(stretch=SqrtStretch())                   #Normalisation for the image
    plt.imshow(Data_log[::-1],'Greys_r',norm=norm)




def Locate_galaxies(Data,N_std):
    """Identify the centre location of the star"""
    global Image_no_noise, centres 
    Data_log = np.log(np.log(np.log(np.log(Data))))                #Logarithm function makes the bright sources more distinct
    threshold = np.mean(Data_log) + N_std*np.std(Data_log)         #Threshold value from which below is considered noise
    Image_no_noise, length = snd.label(Data_log[::-1] > threshold, np.ones((3,3)))
    centres = snd.center_of_mass(Data_log[::-1], Image_no_noise,   #Calculates the centres of the stars
                np.arange(1,length+1,1)) 
    
    return centres
    

  
 
def Show_centres(Data,N_std):
    """Plot the centre of the bright object""" 
    Locate_galaxies(Data,N_std)
    Show_image(Data)
    for x in range(0,len(centres)):
        plt.plot(centres[x][1],centres[x][0],marker='+',color='cyan')
 
        
        
def Centre_on_galaxy(Data,c,Data_error):
    """Crop the image to get a 300x300 galaxy image """
    global Galaxy_centred_Data, Masking, galaxy_size, centres, Error
    # Crop image
    Galaxy_centred_Data = Data[::-1][int(centres[c][0]-150):int(centres[c][0
                              ]+150)][:,int(centres[c][1]-150):int(centres[c][1]+150)] 

    """This code adds the location of another star which is not
    being studied to the mask array."""
    Masking = np.zeros(Galaxy_centred_Data.shape,dtype=bool)    #Create a 2D mask array
    Masking = Image_no_noise[int(centres[c][0]-150):int(centres[c][0]+150)][:,int(
            centres[c][1]-150):int(centres[c][1]+150)] != c+1
    galaxy_size = Galaxy_centred_Data.shape[0]*Galaxy_centred_Data.shape[1] - Masking.sum()
    
    Error = np.zeros(Galaxy_centred_Data.shape,dtype=bool) #Create a 2D error array
    Error = Data_error[int(centres[c][0]-150):int(centres[c][0]+150)][:,int(
            centres[c][1]-150):int(centres[c][1]+150)] != c+1
    
    return Galaxy_centred_Data, Masking, galaxy_size, Error
    


def Calculate_2Dstd(Data,Masking):
        """Calculate the radial standard deviation"""
        Gauss_fit = fit_2dgaussian(Data,mask=Masking)                       #Apply a 2D Gaussian fit to the data
        std = np.sqrt(Gauss_fit.y_stddev.value**2 + Gauss_fit.x_stddev.value**2)
        if std > 30:                                                        #Prevents stds being too large
            std = 29
            
        return std
    
  
def Pixel_rejection_annulus(N_std):    
    """ Put a mask on the pixel tha exceeds the given std dev """
    Annulus_mask = aperture_and_annulus[1].to_mask(method='center')
    Annulus_shape = Annulus_mask.to_image(shape=((300, 300)))
    Pixels_in_annulus = Galaxy_centred_Data*Annulus_shape
    Annulus_mask = Pixels_in_annulus > np.mean(Pixels_in_annulus[
                Pixels_in_annulus != 0]) + N_std*np.std(Pixels_in_annulus[
                        Pixels_in_annulus != 0])
    
    return Annulus_mask


def Counts_single_galaxy_circular_aperture(Data,c,N_std,Data_error,visualise=False):
    """ Calculate number of counts of a galaxy """
    global count_list, aperture_and_annulus, aperture_measurement, centres
    
    centres = Locate_galaxies(Data,N_std)
    Galaxy_centred_Data, Masking, galaxy_size, Error = Centre_on_galaxy(Data,c,
                                                                Data_error)
    
    #Reject the galaxy that is smaller than 6 pixel as they can't be differ from noise
    if galaxy_size > 6:  
        std = Calculate_2Dstd(Data,Masking)
        #Create the aperture and annulus for the galaxy
        aperture_and_annulus = [CircularAperture((150,150), r=2*std),
                            CircularAnnulus((150,150), r_in=3*std, r_out=5*std)]        
        Annulus_mask = Pixel_rejection_annulus(N_std)        
        #Measure the counts in the aperature and annulus
        aperture_measurement = aperture_photometry(Galaxy_centred_Data, aperture_and_annulus,
                                       subpixels=10,mask=Annulus_mask,error=Error)    
        
        counts_with_background = aperture_measurement[0][3]
        background = aperture_measurement[0][5]/aperture_and_annulus[1].area
        counts_error = aperture_measurement[0][4]
        background_error = aperture_measurement[0][6]

        if visualise == True:
                #Plot the aperture and annulus
                fig, axes = plt.subplots()
                Data_log = np.log(Galaxy_centred_Data) #Logarithm put on data to make stars more visible
                norm = ImageNormalize(stretch=SqrtStretch()) #Normalisation for the image
                axes.imshow(Data_log,'Greys_r',norm=norm)
                aperture_and_annulus[0].plot(color='Chartreuse')
                aperture_and_annulus[1].plot(color='Crimson')
            
        #Galaxy_counts = counts_with_background - background * aperture_and_annulus[0].area
        return counts_with_background - background * aperture_and_annulus[0
                    ].area,counts_error + background_error * aperture_and_annulus[
                            0].area, background
    else:
        return 'Galaxy too small','Galaxy too small', 'Galaxy too small'
    
    
def Counts_single_galaxy_elliptical_aperture(Data,c,N_std,Data_error,visualise=False):
    global count_list, aperture_and_annulus, aperture_measurement, centres
    
    centres = Locate_galaxies(Data,N_std)
    Galaxy_centred_Data, Masking, galaxy_size, Error = Centre_on_galaxy(Data,c,
                                                                    Data_error)

    if galaxy_size > 6:  
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
        return counts_with_background - background * aperture_and_annulus[0
                    ].area,counts_error + background_error * aperture_and_annulus[
                            0].area, background
    
    else:
        return 'Galaxy too small','Galaxy too small','Galaxy too small'
      
          
def Counts_all_galaxies_circular_aperture(Data,N_std,Data_error):
    global centres, counts_list
    
    counts_list = []
    centre_list = []
    error_list = []
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


def Counts_all_galaxies_elliptical_aperture(Data,N_std,Data_error):
    global centres, counts_list, centre_list, error_list, background_list
    
    print(len(centres))
    
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

