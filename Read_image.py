# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:25:56 2020

@author: Harry Anthony
"""

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
#Calculate image statistics     

print(np.mean(Data))
print(np.std(Data))
print(Data.shape)



#%%

#Crops the edges of the image (Noise)

New = []        #Crops the image

y = 0
while y < 4406-424:
    y_array = Data[424+y]
    New.append(y_array[242:2378])
    y = y + 1
    
#%%
#    what percentage of the image are croppped
New = np.array(New)

print((New.shape[0]*New.shape[1])/(Data.shape[0]*Data.shape[1]))

#%%
#Write as Fits file 

plt.imshow(newy)
hdu = fits.PrimaryHDU(newy)
hdul = fits.HDUList([hdu])
hdul.writeto('Mosaic_no_blooming9.fits')

#%%

#calculate the centres of galaxy
import scipy.ndimage as snd
Image_no_noise, length = snd.label(Data[::-1] > 25000, np.ones((3,3)))                  # Isolate the stellar object from noise
centres = snd.center_of_mass(Data[::-1], Image_no_noise, np.arange(1,length+1,1))       # Calculate the centre of the stellar object




#%%


# masking the blooming stars
mean = np.median(Data)
newd = Data.copy()
meanz = np.full(30,mean)

#shape of the central star is replaced by the mean of image masking
for y in range(0,3982):
    newd[y][1180:1210] = meanz
for y in range(0,39):
    newd[y][860+6*y:1460-6*y] = np.full(1460-860-12*y,mean) 
for y in range(2670-105,2900+135):
    newd[y][1050-100:1350+100] = np.full(1350-1050+200,mean) 
for y in range(2670,2900):
    newd[y][1120:1290] = np.full(1290-1120,mean) 
for y in range(0,27):
    newd[y][782:803] = np.full(803-782,mean)
    
# next few largest star
#Next largest star
for y in range(2775,3000):
    newd[y][500-19:570+19] = np.full(570-500+38,mean)
    
#Next star
for y in range(2250,2420):
    newd[y][700-19:760+19] = np.full(760-700+38,mean)
    
#Next star
for y in range(1780,1940):
    newd[y][640-19:690+19] = np.full(690-640+38,mean)
    
#Next star
for y in range(3280,3380):
    newd[y][1870-19:1920+19] = np.full(1920-1870+38,mean)
    

#%%
    
#visualise the final image
#    change image variab;e newdy
    
norm = ImageNormalize(stretch=SqrtStretch())
newdy = np.log(newd)
plt.imshow(newdy,'Greys_r',norm=norm)

#%%
#mask the stars considered brightest pixels
#really bright ones
#lots of blooming
#finds the bright stars
#replace them with a mean

mean = np.median(Data)
newy = newd.copy()

Image_no_noise, length = snd.label(newd[::-1] > 25000, np.ones((3,3)))
centres = snd.center_of_mass(newd[::-1], Image_no_noise, np.arange(1,length+1,1))
true_centres = [[3982-centres[x][0],centres[x][1]] for x in range(0,len(centres))]


for p in range(2,len(centres)):
    y = int(np.around(true_centres[p][0],0))
    x = int(np.around(true_centres[p][1],0))
    print(str(y)+' '+str(x))
    z = 0
    print(p)
    
    while z < 60:
        w = 0
        while w < 60:
            #print(y+z-30)
            #print(x+w-30)
            newy[y+z-30][x+w-30] = mean
            w = w + 1
        z = z + 1

#%%
#        visualise image
from pylab import figure, cm
from matplotlib.colors import LogNorm

newyy = np.log(newy)
plt.imshow(newyy, 'Greys_r',norm=norm)
        

    