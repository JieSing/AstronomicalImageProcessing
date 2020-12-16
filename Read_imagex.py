# -*- coding: utf-8 -*-

from astropy.io import fits
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

""" Opens the file, reads the header and data, and close the file"""

FITSdata = fits.open('H:\Imperial\EXPT LAB\Y3\A1\mosaic.fits',
                    mode='readonly')

Header = FITSdata[0].header
Data = FITSdata[0].data
FITSdata.close()            

#%%
#Visualise the data as an image
plt.imshow(Data[::-1])

#%%
#Plot the histogram of the image
plt.hist(Data.flatten(),bins=1000)

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

New = []        #The remaining pixels

y = 0
while y < 4406-424:
    y_array = Data[424+y]
    New.append(y_array[242:2378])
    y = y + 1
    
#%%
#The percentage of the image that are croppped
New = np.array(New)
print((New.shape[0]*New.shape[1])/(Data.shape[0]*Data.shape[1]))

#%%
#Write the cropped data as Fits file 

plt.imshow(New)
hdu = fits.PrimaryHDU(New)
hdul = fits.HDUList([hdu])
hdul.writeto('Mosaic_no_blooming9.fits')

#%%

#Calculate the centres of galaxy
import scipy.ndimage as snd
Image_no_noise, length = snd.label(Data[::-1] > 25000, np.ones((3,3)))                  # Isolate the stellar object from noise
centres = snd.center_of_mass(Data[::-1], Image_no_noise, np.arange(1,length+1,1))       # Calculate the centre of the stellar object


#%%
#Mask the blooming stars
mean = np.median(Data)
new_data = Data.copy()
meanz = np.full(30,mean)

#Mask the central star with the mean value of image
for y in range(0,3982):
    new_data[y][1180:1210] = meanz
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
    
norm = ImageNormalize(stretch=SqrtStretch())
new_data_brighter = np.log(new_data)
plt.imshow(new_data_brighter,'Greys_r',norm=norm)

#%%
#Masks the unusually bright stars 

mean = np.median(Data)
new_y = new_data.copy()

Image_no_noise, length = snd.label(new_data[::-1] > 25000, np.ones((3,3)))
centres = snd.center_of_mass(new_data[::-1], Image_no_noise, np.arange(1,length+1,1))
true_centres = [[3982-centres[x][0],centres[x][1]] for x in range(0,len(centres))]


#Finds the bright stars and replace them with the mean value of the image
for p in range(2,len(centres)):
    y = int(np.around(true_centres[p][0],0))
    x = int(np.around(true_centres[p][1],0))
    print(str(y)+' '+str(x))
    z = 0
    print(p)
    
    while z < 60:
        w = 0
        while w < 60:
#            print(y+z-30)
#            print(x+w-30)
            new_y[y+z-30][x+w-30] = mean
            w = w + 1
        z = z + 1

#%%
#Visualise the final image
        
from pylab import figure, cm
from matplotlib.colors import LogNorm

new_image = np.log(new_y)
plt.imshow(new_image, 'Greys_r',norm=norm)
        
