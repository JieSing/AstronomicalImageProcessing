# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 10:06:28 2020

@author: jsy18
"""

#stellar_creator

from astropy.io import fits
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import scipy as sp
from scipy import signal
from scipy import optimize

FITSdata = fits.open("H:\Imperial\EXPT LAB\Y3\A1\Mosaic_no_blooming9.fits", mode='readonly')
"""blooming9_fits pixels : 2136x3982 """

Header = FITSdata[0].header
Data = FITSdata[0].data

# Visualise the data as an image
plt.figure(1)
plt.imshow(Data[::-1])      
plt.xlabel('x')
plt.ylabel('y')  
plt.show()


#%%

x = np.linspace(-2,2,10)                        
y = np.linspace(-2,2,10)    
xx, yy = np.meshgrid(x,y)                                   #Creates a 2D square grid    
d = np.sqrt(xx**2 + yy**2)                                  #Values to input into 2D gaussian function
print('d=',d)

# Intialising sigma and mu
sigma = 0.6
mu = 0.0
  

gauss = 20000*np.exp(-( (d-mu)**2 / ( 2 * sigma**2 ) ) ) + np.mean(Data)    # Calculate Gaussian array 
                                                                            # intensity need to be larger than the mean of diagram to be obvious

plt.figure(2)                                                               # Visualise the cluster object
h = plt.contourf(xx,yy,gauss)
plt.show()

# add cluster objects onto the picture randomly
mark=0                                      
while mark < 9:                                                             # no. of cluster objects to add
    a = np.random.randint(0,len(Data)-gauss.shape[0])                       # row
    b = np.random.randint(0,len(Data[0])-gauss.shape[1])                    # column
    print(a,b)
    for i in range(0, gauss.shape[0]):
        for j in range(0, gauss.shape[1]):            
            if gauss[i,j] > Data[a+i, b+j]:                                 # To avoid adding dark pixels on top of bright objects
                Data[a+i, b+j] = gauss[i,j]                                 # Subst the gaussian array into the Data array (based on matrix notation: x row y column)
    mark+=1
    print(mark)

            
# Plot the final image
print('data=',Data)
plt.figure(3)
plt.imshow(Data[::-1])      
plt.xlabel('x')
plt.ylabel('y')  


# Write the final data as a FITS file
hdu = fits.PrimaryHDU(Data)
hdul = fits.HDUList([hdu])
hdul.writeto('Mosaic_no_blooming9_addobject.fits')
