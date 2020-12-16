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
#Display the image
plt.figure(1)
Show_image(Data)
plt.show()

#Plot the stars located by the method
plt.figure(2)
Show_centres(Data,0.7)

#%%
#Calculate the counts of each of the galaxies detected

counts_list,centre_list,error_list,background_list = Counts_all_galaxies_circular_aperture(
        Data,0.7,Data_error)

#Convert the counts to magnitude
MAGZPT  = 2.530E+01
MAGZRR  = 2.000E-02 
calibrated_mag = MAGZPT - 2.5* np.log10(counts_list)
calibrated_mag_err = MAGZRR + 2.5*(np.array(error_list)/np.array(counts_list))

#Write the result into a textfile
myfile=open('Counts_measurement_circular_survey5.txt','a')
myfile.write('x'+','+'y'+','+'Magnitude'+','+'Magnitude error'+','+'background'+'\n')
for x in range(0,len(counts_list)):
    myfile.write(str(centre_list[x][0]-300+242)+','+str(len(Data)-centre_list[x][1]-300+424)+','+str(calibrated_mag[x])+','+str(calibrated_mag_err[x])+','+str(background_list[x])+'\n')

myfile.close()

#%%
#Read the results from a textfile
x,y,calibrated_mag,calibrated_mag_err,background = np.loadtxt('Counts_measurement_circular_survey5.txt', dtype=str, skiprows=1,delimiter=',', unpack=True)
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
#plt.fill_betweenx(np.log10(Nm), np.array(M)-np.array(M_err), np.array(M
#                  )+np.array(M_err),capstyle='round',lw=0.5, color = (255/256,127/256,80/256), alpha = 1)
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
#plt.fill_betweenx(np.log10(Nm_clipped), np.array(M_clipped)-np.array(M_err_clipped
#                  ), np.array(M_clipped)+np.array(M_err_clipped), color = (
#                          255/256,127/256,80/256), capstyle='round',lw=0.5,alpha = 1)


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
plt.xlim(12,14)
plt.ylim(2.25,3.25)
plt.legend()
plt.show()
