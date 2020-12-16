# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 12:40:04 2020

@author: Harry Anthony
"""

from Aperture_methods import *

FITSdata = fits.open('C:\\Users\\Harry\\Documents\\My documents\\Third year\
\\Laboratory work\\Astronomical image processing\\Padded_data.fits',
                    mode='readonly')

Header = FITSdata[0].header
Data = FITSdata[0].data

Mean = np.mean(Data) # Calculate the mean
std = np.std(Data)
#Data = np.pad(Data, 150, mode='constant',constant_values=Mean)

#hdu = fits.PrimaryHDU(Data)
#hdul = fits.HDUList([hdu])
#hdul.writeto('Padded_data.fits')
#%%

#Data_clipped = Data[len(Data)-2770:len(Data)-2940+100][:,1900+100:2050+100] 
 
#Data_clipped = Data[len(Data)-3200:len(Data)-2950]
Data_clipped = Data[2950:3050][:,2100:2200]

norm = ImageNormalize(stretch=SqrtStretch())
#plt.imshow(np.log(Data_clipped[::-1]))
plt.imshow(np.log(np.log(np.log(np.log(Data_clipped[::-1])))),'Greys_r',norm=norm)

Galaxy_centred_data = Data_clipped

#%%

Data_clipped = Data[1630:1690][:,960:1020]

norm = ImageNormalize(stretch=SqrtStretch())
#plt.imshow(np.log(Data_clipped[::-1]))
plt.imshow(np.log(np.log(np.log(np.log(Data_clipped[::-1])))),'Greys_r',norm=norm)

Galaxy_centred_data = Data_clipped.copy()
#%%

Data_log = np.log(Data) #Logarithm makes the bright sources more distinct
#Threshold from which below is considered noise
threshold = np.mean(Data_log) + np.std(Data_log) 
Image_no_noise, length = snd.label(Data_log[::-1] > threshold, np.ones((3,3)))
centres = snd.center_of_mass(Data_log[::-1], Image_no_noise, 
                np.arange(1,length+1,1)) #Calculates the centres of the stars
true_centres = [[len(Data)-centres[x][0],centres[x][1]] for x in range(
            0,len(centres))]
Galaxy_centred_Data = Data[::-1][int(centres[165][0]-150):int(centres[165][0]+150
                          )][:,int(centres[165][1]-150):int(centres[165][1]+150)]

#This code adds to the mask array where another star is located which is not
#being studied.
Masking = np.zeros(Galaxy_centred_Data.shape,dtype=bool) #Create a 2D mask array
Masking = Image_no_noise[int(centres[165][0]-150):int(centres[165][0]+150)][:,int(
            centres[165][1]-150):int(centres[165][1]+150)] != 166
galaxy_size = Galaxy_centred_Data.shape[0]*Galaxy_centred_Data.shape[1] - Masking.sum()

Galaxy_centred_Data = Galaxy_centred_Data[120:180][:,120:180]
    
#%%
from photutils import DAOStarFinder
from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils import find_peaks

psf_model = IntegratedGaussianPRF(sigma=3)
print(psf_model)

#a = DAOStarFinder(Mean+std,fwhm=6)
#sources = a(Galaxy_centred_data)
#positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
#apertures = CircularAperture(positions, r=4.)
#apertures.plot(color='blue', lw=1.5, alpha=0.5)

plt.imshow(Galaxy_centred_Data)

#tbl = find_peaks(Galaxy_centred_data, Mean+std, box_size=11)
#tbl['peak_value'].info.format = '%.8g'  # for consistent table output
#print(tbl)  # print only the first 10 peaks

iraffind = IRAFStarFinder(threshold=1*std,
                          fwhm=1*gaussian_sigma_to_fwhm)

from photutils.psf import BasicPSFPhotometry
my_photometry = IterativelySubtractedPSFPhotometry(
     finder= iraffind, group_maker=DAOGroup(5),
     bkg_estimator=MMMBackground(),psf_model = psf_model,niters=20,
     fitter=LevMarLSQFitter(),fitshape = (11,11))

result_tab = my_photometry(Galaxy_centred_Data)
residual_image = my_photometry.get_residual_image()

print(result_tab)
x,y = result_tab['x_fit'],result_tab['y_fit']

plt.plot(x,y,'kx')

plt.imshow(residual_image)
plt.colorbar()

#%%

from ..photometry import DAOPhotPSFPhotometry

# create a fittable 2d model
my_psf_model = IntegratedGaussianPRF(sigma=3.0)

# create a finder object
my_finder = DAOStarFinder(threshold=5., fwhm=3.0*gaussian_sigma_to_fwhm)

# create a grouping object
my_grouping = DAOGroup(crit_separation=1.5*3.0*gaussian_sigma_to_fwhm)

# create a background object
my_bkg = MMMBackground()

# create a fitter object
my_fitter = LevMarLSQFitter()

# create psf photometry object
my_photometry = DAOPhotPSFPhotometry(finder=my_finder, group_maker=my_grouping,
                                   bkg_estimator=my_bkg, psf_model=my_psf_model,
                                   fitter=my_fitter, fitshape=(7,7))

# get photometry results
photometry_results = my_photometry(Galaxy_centred_Data) 

#%%
plt.figure(1)
plt.imshow(Galaxy_centred_Data)
plt.colorbar()
plt.figure(2)
plt.imshow(residual_image,'Greys_r')
plt.colorbar()
#plt.imshow(Galaxy_centred_Data,'Greys_r',norm=LogNorm())

#%%

#Data_clipped = Data[2950:3050][:,2100:2200]

#norm = ImageNormalize(stretch=SqrtStretch())
#plt.imshow(np.log(Data_clipped[::-1]))
#plt.imshow(np.log(np.log(np.log(np.log(Data_clipped[::-1])))),'Greys_r',norm=norm)

#Galaxy_centred_data = Data_clipped.copy()

from Aperture_methods import *
Show_centres(Data_clipped,0.7)
centres = Locate_galaxies(Data_clipped,0.7)

#%%

from astropy.convolution import Gaussian2DKernel, AiryDisk2DKernel, convolve_fft

Galaxy_centred_Data, Masking, galaxy_size = Centre_on_galaxy(Data_clipped,0)
std = Calculate_2Dstd(Galaxy_centred_Data,Masking)
psf = AiryDisk2DKernel(std)
plt.imshow(psf)
convolved_image = convolve_fft(Galaxy_centred_Data, psf, boundary='wrap')
plt.figure(2)
plt.imshow(convolved_image)
plt.figure(3)
plt.imshow(Galaxy_centred_Data)

#%%

from astropy.convolution import Gaussian2DKernel, AiryDisk2DKernel, convolve_fft

Masking = np.zeros(Data_clipped.shape,dtype=bool) #Create a 2D mask array
Masking = Image_no_noise != 0+1
galaxy_size = Data_clipped.shape[0]*Data_clipped.shape[1] - Masking.sum()
    
#%%
plt.figure(1)
plt.imshow(Data_clipped[::-1])
plt.figure(2)
plt.imshow(Masking)

#%%    
#from scipy.signal import deconvolve
from scipy import fftpack
from astropy.convolution import Gaussian2DKernel, AiryDisk2DKernel, convolve_fft, Tophat2DKernel
#from pysap.plugins.astro.deconvolution.deconvolve import sparse_deconv_condatvu
#from scikitimage import restoration
from scipy.signal import convolve2d as conv2

from skimage import color, data, restoration

Masking = np.zeros(Data_clipped.shape,dtype=bool)
std = Calculate_2Dstd(Data_clipped,Masking)
psf = AiryDisk2DKernel(2)
#psf = Tophat2DKernel(2)
plt.imshow(psf.array)
convolved_image = convolve_fft(Data_clipped, psf, boundary='wrap')
plt.figure(2)
plt.imshow(convolved_image)
plt.figure(3)
plt.imshow(Data_clipped)
plt.figure(4)
#plt.imshow(Data_clipped - convolved_image)
deconvolved_RL = restoration.richardson_lucy(Data_clipped, psf.array, iterations=3000)
#plt.imshow(deconvolved_RL,vmin=Data_clipped.min(), vmax=Data_clipped.max(),norm=norm)
plt.imshow(deconvolved_RL.real)
#delta_x_psf=5 # number of pixels from the edges
#xmin, xmax = -psf.shape[1]-delta_x_psf, -delta_x_psf
#ymin, ymax = delta_x_psf, delta_x_psf+psf.shape[0]
#convolved_image[xmin:xmax, ymin:ymax] = psf.array/psf.array.max()*10

#fig = plt.figure(figsize=(8,12))
#i_plot = fig.add_subplot(111)
#plt.imshow(np.log10(convolved_image+1e-3), vmin=-1, vmax=1.0, origin='lower')#, cmap=plt.cm.viridis)
#plt.colorbar()

#def deconvolve(star, psf):
#    star_fft = fftpack.fftshift(fftpack.fftn(star))
#    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
#    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft/psf_fft)))

#deconv_data = sparse_deconv_condatvu(Data, psf.data, n_iter=3000)
#deconvolved = deconvolve(Data_clipped,psf)
#plt.figure(4)
#plt.imshow(deconvolved)
#%%

from Aperture_methods import *
Show_centres(convolved_image,0.7)
centres = Locate_galaxies(convolved_image,0.7)

plt.figure(2)
Show_centres(Data_clipped,0.7)

#%%

#import k2
import phot
import pylab as py

fn = 'kplr060018142-2014044044430_lpd-targ.fits'
f = pyfits.open('kplr060018142-2014044044430_lpd-targ.fits')
image = f[1].data['flux'][0]
prf = k2.loadPRF(fn)
out = phot.psffit(prf, image, (33, 23), scale=50, dframe=7, verbose=True)

py.figure()
py.subplot(132)
py.imshow(out[0])
py.title('Best-fit Model PRF')
py.colorbar()
py.subplot(131)
py.imshow(out[1])
py.title('Observed Data')
py.clim([out[0].min(), out[0].max()])
py.colorbar()
py.subplot(133)
py.imshow(out[1] - out[0])
py.title('Data - Model')
py.colorbar()