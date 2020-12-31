# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 10:40:12 2020

@author: Harry Anthony
"""
from Aperture_methods import *
from PSF_methods import *


#Open the data and the error file
FITSdata = fits.open('C:\\Users\\Harry\\Documents\\My documents\\Third year\
\\Laboratory work\\Astronomical image processing\\Mosaic_no_blooming9.fits',
                    mode='readonly')

Header = FITSdata[0].header
Data = FITSdata[0].data
FITSdata.close()

FITSdata = fits.open('C:\\Users\\Harry\\Documents\\My documents\\Third year\
\\Laboratory work\\Astronomical image processing\\poisson_error (1).fits',
                    mode='readonly')

Header_error = FITSdata[0].header
Data_error = FITSdata[0].data
FITSdata.close()


Mean = np.mean(Data) # Calculate the mean
Data = np.pad(Data, 300, mode='constant',constant_values=Mean) #Pad the data to stop edge effects
Mean_error = np.mean(Data_error)
Data_error = np.pad(Data_error, 300, mode='constant',constant_values=Mean_error) #Pad the data to stop edge effects

#%%
#This cell creates a mask array for the data

masked = np.full((Data.shape[0],Data.shape[1]),True,dtype=bool) #Create a 2D boolean array

#Mask the edges of the image
y = 0
while y < len(Data)-600:
    masked[300+y][300:len(Data[0])-300] = np.full(len(Data[0])-600,False)
    y = y + 1

#Mask the central star
for y in range(0,39):
    masked[y+424][860+6*y+242:1460-6*y+242] = np.full(1460-860-12*y,True) 
for y in range(2470-105,2900+135):
    masked[y+424][1050-100+242:1500+100+242] = np.full(1500-1050+200,True) 
for y in range(2670,2900):
    masked[y+424][1120+242:1290+242] = np.full(1290-1120,True) 
for y in range(2670,2900):
    masked[y+424][1120+242:1290+242] = np.full(1290-1120,True) 
for y in range(0,27):
    masked[y+424][782+242:803+242] = np.full(803-782,True)
    
#Next largest star
for y in range(2775,3000):
    masked[y+424][500-19+242:570+19+242] = np.full(570-500+38,True)
    
#Next star
for y in range(2250,2420):
    masked[y+424][700-19+242:760+19+242] = np.full(760-700+38,True)
    
#Next star
for y in range(1780,1940):
    masked[y+424][640-19+242:690+19+242] = np.full(690-640+38,True)
    
#Next star
for y in range(3280,3380):
    masked[y+424][1870-19+242:1920+19+242] = np.full(1920-1870+38,True)

#%%
#Display the image
plt.figure(1)
Show_image(Data)
plt.show()

#Plot the stars located by the method
plt.figure(2)
Show_centres(Data,2,masked)

#%%
#This cell plots the shapes of the galaxies calculated with and without the watershed algorithm

plt.figure(1)
Image_no_noise_array(Data,2,masked)
plt.figure(2)
Image_no_noise_watershed_array(Data,2,masked)


#%%
#Calculate the counts of each of the galaxies detected

counts_list,centre_list,error_list,background_list = Counts_all_galaxies_elliptical_aperture(
        Data,2,Data_error,mask)

#Convert the counts to magnitude
MAGZPT  = 2.530E+01
MAGZRR  = 2.000E-02 
calibrated_mag = MAGZPT - 2.5* np.log10(counts_list)
calibrated_mag_err = MAGZRR + 2.5*(np.array(error_list)/np.array(counts_list))

#Write the result into a textfile
myfile=open('Counts_measurement_circular_survey.txt','a')
myfile.write('x'+','+'y'+','+'Magnitude'+','+'Magnitude error'+','+'background'+'\n')
for x in range(0,len(counts_list)):
    myfile.write(str(centre_list[x][0]-300+242)+','+str(len(Data)-centre_list[x][1]-300+424
                                                       )+','+str(calibrated_mag[x])+','+str(calibrated_mag_err[x]
                                                                               )+','+str(background_list[x])+'\n')

myfile.close()

#%%
#Read the results from a textfile
x,y,ra,dec,calibrated_mag,calibrated_mag_err,background = np.loadtxt('Catalogue_elliptical_aperture_and_psf.txt', dtype=str, skiprows=1,delimiter=',', unpack=True)
calibrated_mag_err = calibrated_mag_err.astype(np.float)
calibrated_mag = calibrated_mag.astype(np.float)


#%%
"""
This section of code plots the magnitude plot (ln(m) vs m) for the measured data
"""

#Sort the magnitudes in order
calibrated_mag_with_err = [[calibrated_mag[x],calibrated_mag_err[x]] for x in range(0,len(calibrated_mag))]
calibrated_mag_with_err.sort(key=lambda x: x[0]) #Sort by the zeroth element

#Initialise variables to append our results in to
Nm = []
M = []
M_err = []

#Iterate through all the magnitudes
x = 0 
while x < len(calibrated_mag_with_err):
    cumulative = 0 #Counts the number of galaxies with magnitude below m
    for j in range(0,len(calibrated_mag_with_err)):
        global Next_val
        if calibrated_mag_with_err[j][0] <= calibrated_mag_with_err[x][0]:
            cumulative = cumulative + 1 
        else:
            Next_val = j 
            break
    #Append the results to the array
    Nm.append(cumulative)
    M.append(calibrated_mag_with_err[x][0])
    M_err.append(calibrated_mag_with_err[x][1])
    #Ends the loop if it has iterated through all magnitudes (increases efficiency)
    if x == len(calibrated_mag_with_err)-1:
        break
    #While loop to skip through galaxies with the same magnitude (increases efficiency)
    while calibrated_mag_with_err[x][0] == calibrated_mag_with_err[x+1][0]:
        x = x+1
    x = x + 1
    
#Error associated with Nm can be modelled as Poisson
Nm_err = np.sqrt(Nm)

#Plot the total N(m) against m graph
plt.figure(1)
plt.plot(M,np.log10(Nm),'kx')
plt.fill_between(M, np.log10(Nm)-Nm_err/Nm, np.log10(Nm)+Nm_err/Nm,
                 capstyle='round',lw=0.5,color = (255/256,127/256,80/256), alpha = 1)
plt.xlabel('Magnitude M')
plt.ylabel('log$_{10}$(N(M))')
plt.grid()
plt.xlim(6)
plt.ylim(0)
plt.show()

#Plot only the linear region of the N(m) against m graph
plt.figure(2)

#Initialise arrays to hold the data
M_clipped = []
M_err_clipped = []
Nm_clipped = []
Nm_err_clipped = []

#Append values only in the linear range
for t in range(0,len(M)):
    if M[t] > 11 and M[t] < 14:
        M_clipped.append(M[t])
        M_err_clipped.append(M_err[t])
        Nm_clipped.append(Nm[t])
        Nm_err_clipped.append(Nm_err[t])

#Plot the data in the linear region
plt.plot(M,np.log10(Nm),'kx')
plt.fill_between(M_clipped, np.log10(Nm_clipped)-np.array(Nm_err_clipped
                 )/np.array(Nm_clipped), np.log10(Nm_clipped)+np.array(Nm_err_clipped
                           )/np.array(Nm_clipped), color = (255/256,127/256,80/256
                                     ),capstyle='round',lw=5, alpha=1)

#Calculate the linear line of best fit for the data
fit,cov = np.polyfit(M_clipped, np.log10(Nm_clipped),1, w=np.array(Nm_clipped
                       )/np.array(Nm_err_clipped), cov=True )
linebestfit= np.poly1d(fit)
h = np.arange(0,20,1)
plt.plot(h,linebestfit(h), color='turquoise',label='Experimental value trend line')

#Calculate the gradient and the error
gradient = str(("%." + str(3) + "e") % fit[0]) #Round gradient to 3sf
gradient_err = str(("%." + str(3) + "e") % np.sqrt(cov[0][0]))
plt.text(13,2.6,'Gradient = '+gradient + ' +/- '+ gradient_err)

plt.xlabel('Magnitude M')
plt.ylabel('log$_{10}$(N(M))')
plt.grid()
plt.xlim(11,14)
plt.legend()
plt.show()

#%%
#This cell places the galaxies into 0.5 magnitude bins

lognm = np.log10(Nm)
#Combine the magnitude and logarithm
combo = [[M[t],lognm[t]] for t in range(0,len(M))]
x = 7.5 #Start at magnitude 7.5
bins = []
#Appends the average in the bin into the array bins
while x < 22:
    av = [] 
    for y in range(0,len(lognm)):
        if x - 0.5 <= combo[y][0] and x >= combo[y][0]:
            av.append(combo[y][1])
    bins.append(np.mean(av))
    x = x + 0.5

#%%
#This cells plots the ln(N(m)) vs m plot with 0.5 magnitude bins

#Graph embellishments
params = {
   'axes.labelsize': 20,
   'axes.titlesize': 20,
   'font.size': 20,
   'font.family': 'serif',
   'legend.fontsize': 20,
   'xtick.labelsize': 20,
   'ytick.labelsize': 20,
   'figure.figsize': [10, 7]
   }
plt.rcParams.update(params)
fam = {'fontname':'Times New Roman'}

#Errors are of the form
errors = (10**(np.array(bins))*np.log(10)*(1/(2*np.array(bins))))

#Plot the graph using log scales
x = np.arange(7.5,22,0.5)
plt.plot(x,10**(np.array(bins)),'kx',label='Experimental values')
plt.errorbar(x,10**(np.array(bins)),yerr=errors,fmt='kx',ms=10,capsize=5)
plt.yscale('log') #Place a log scale on the y-axis
plt.xlabel('Magnitude m')
plt.ylabel('Cumulative number of galaxies N(m)')
plt.ylim(0)
plt.grid()
plt.xlim(8,22)
plt.ylim(1,10**5)


