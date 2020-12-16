# -*- coding: utf-8 -*-

"""
Created on Wed Dec 16 11:24:46 2020

@author: Harry Anthony
"""

from astropy.io import fits
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from pylab import figure, cm
from matplotlib.colors import LogNorm

""" Opens the file, reads the header and data, and close the file"""

FITSdata = fits.open('C:\\Users\\Harry\\Documents\\My documents\\Third year\
\\Laboratory work\\Astronomical image processing\\A1_Mosaic.fits',
                    mode='readonly')

Header = FITSdata[0].header
Data = FITSdata[0].data
FITSdata.close()            

#%%
#Visualise the data as an image
norm = ImageNormalize(stretch=SqrtStretch())
plt.imshow(Data[::-1],'Greys_r',norm=norm)

#%%
#Plot the histogram of the image
plt.hist(Data.flatten(),bins=1000)
plt.xlim('Value')
plt.ylim('Relative frequency')

#%%
#Finds the brightest pixel in the image

maximum = 0
pos_x = 0
pos_y = 0
for y in range(0,len(Data)):
    for x in range(0,len(Data[y])):
        if float(Data[y][x]) > float(maximum):
            minimum = Data[y][x]                    #Finds the maximum value of data
            pos_x = x
            pos_y = y
            

#%%
#Calculates image statistics     

print(np.mean(Data))
print(np.std(Data))
print(Data.shape)

#%%

#Crops the edges of the image (Noise)

cropped_data = []        #The remaining pixels

y = 0
while y < 4406-424:
    y_array = Data[424+y]
    cropped_data.append(y_array[242:2378])
    y = y + 1
    
#%%
#The percentage of the image that are croppped
cropped_data = np.array(cropped_data)
print((cropped_data.shape[0]*cropped_data.shape[1])/(Data.shape[0]*Data.shape[1]))

#%%
#Calculate the centres of galaxy
import scipy.ndimage as snd
Image_no_noise, length = snd.label(cropped_data[::-1] > 25000, np.ones((3,3)))                  # Isolate the stellar object from noise
centres = snd.center_of_mass(cropped_data[::-1], Image_no_noise, np.arange(1,length+1,1))       # Calculate the centre of the stellar object

#%%
#Mask the blooming stars
mean = np.median(cropped_data)
new_data = cropped_data.copy()

#Mask the central star with the mean value of image
for y in range(0,3982):
    new_data[y][1180:1210] = np.full(30,mean)
for y in range(0,39):
    new_data[y][860+6*y:1460-6*y] = np.full(1460-860-12*y,mean) 
for y in range(2670-105,2900+135):
    new_data[y][1050-100:1350+100] = np.full(1350-1050+200,mean) 
for y in range(2670,2900):
    new_data[y][1120:1290] = np.full(1290-1120,mean) 
for y in range(0,27):
    new_data[y][782:803] = np.full(803-782,mean)
    

#Mask the next largest star
for y in range(2775,3000):
    new_data[y][500-19:570+19] = np.full(570-500+38,mean)
    
#Mask the next star
for y in range(2250,2420):
    new_data[y][700-19:760+19] = np.full(760-700+38,mean)
    
#Mask the next star
for y in range(1780,1940):
    new_data[y][640-19:690+19] = np.full(690-640+38,mean)
    
#Mask the next star
for y in range(3280,3380):
    new_data[y][1870-19:1920+19] = np.full(1920-1870+38,mean)
    

#%%
#Visualise the image
    
new_data_brighter = np.log(np.log(np.log(np.log(new_data))))
plt.imshow(new_data_brighter[::-1],'Greys_r')

#%%
#Masks the unusually bright stars 

mean = np.median(Data)
masked_data = new_data.copy()

Image_no_noise, length = snd.label(new_data[::-1] > 25000, np.ones((3,3)))
centres = snd.center_of_mass(new_data[::-1], Image_no_noise, np.arange(1,length+1,1))
true_centres = [[len(new_data)-centres[x][0],centres[x][1]] for x in range(0,len(centres))]


#Finds the bright stars and replace them with the mean value of the image
for p in range(2,len(centres)):
    y = int(np.around(true_centres[p][0],0))
    x = int(np.around(true_centres[p][1],0))
    z = 0
    
    while z < 60:
        w = 0
        while w < 60:
            masked_data[y+z-30][x+w-30] = mean
            w = w + 1
        z = z + 1

#%%
#Visualise the final image     
norm = ImageNormalize(stretch=SqrtStretch())
new_image = np.log(np.log(np.log(np.log(masked_data))))
plt.imshow(new_image[::-1], 'Greys_r',norm=norm)

#%%

hdu = fits.PrimaryHDU(masked_data)
hdul = fits.HDUList([hdu])
hdul.writeto('Mosaic_no_blooming20.fits')
