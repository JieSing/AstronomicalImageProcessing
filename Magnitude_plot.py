# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 10:40:12 2020

@author: Harry Anthony
"""

from Aperture_methods import *

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


Mean = np.mean(Data)                                                            # Calculate the mean
Data = np.pad(Data, 300, mode='constant',constant_values=Mean)                  #Pad the data to eliminate edge effects
Mean_error = np.mean(Data_error)
Data_error = np.pad(Data_error, 300, mode='constant',constant_values=Mean_error) #Pad the data to eliminate edge effects

#%%

plt.figure(1)                                                                   #Display the image
Show_image(Data)
plt.show()


plt.figure(2)                                                                   #Plot the stars located by the code
Show_centres(Data,0.7)

#%%
counts_list,centre_list,error_list, background_list = Counts_all_galaxies_elliptical_aperture(
        Data,0.7,Data_error)


#Use the data to estimate the calibrated magnitude of the star
"""
EXPTIME = 720.000
MAGZPT  = 2.530E+01
MAGZRR  = 2.000E-02 
print(calibrated_mag)
"""

MAGZPT  = 2.530E+01
MAGZRR  = 2.000E-02 
calibrated_mag = MAGZPT - 2.5* np.log10(counts_list)                            #The counts are transformed into relative amplitudes 
calibrated_mag_err = 2.5*(np.array(error_list)/np.array(counts_list))

myfile=open('Counts_measurement_elliptical_survey.txt','a')                     # Write the galactic survey into a text file
myfile.write('x'+','+'y'+','+'Magnitude'+','+'Magnitude error'+','+'background'+'\n')
for x in range(0,len(counts_list)):
    myfile.write(str(centre_list[x][0]-300+242)+','+str(len(Data)-centre_list[x][1]-300+424)+','+str(calibrated_mag[x])+','+str(calibrated_mag_err[x])+','+str(background_list[x])+'\n')

myfile.close()

#%%

x,y,calibrated_mag,calibrated_mag_err,background = np.loadtxt('Counts_measurement_elliptical_survey.txt', dtype=str, skiprows=1,delimiter=',', unpack=True)
calibrated_mag = calibrated_mag.astype(np.float)
calibrated_mag_err = calibrated_mag_err.astype(np.float)


#%%

calibrated_mag_with_err = [[calibrated_mag[x],calibrated_mag_err[x]] for x in range(0,len(calibrated_mag))]
calibrated_mag_with_err.sort(key=lambda x: x[0])

Nm = []
M = []
M_err = []
#Cumulative_err = []
x = 0
while x < len(calibrated_mag_with_err):
    #print(x)
    cumulative = 0
    for j in range(0,len(calibrated_mag_with_err)):
        global Next_val
        if calibrated_mag_with_err[j][0] <= calibrated_mag_with_err[x][0]:
            cumulative = cumulative + 1
            
        else:
            Next_val = j
            break
    Nm.append(cumulative)
    M.append(calibrated_mag_with_err[x][0])
    M_err.append(calibrated_mag_with_err[x][1])
    if x == len(calibrated_mag_with_err)-1:
        break
    while calibrated_mag_with_err[x][0] == calibrated_mag_with_err[x+1][0]:
        x = x+1
    x = x + 1
    

plt.plot(M,np.log10(Nm),'kx')
plt.fill_betweenx(np.log10(Nm), np.array(M)-np.array(M_err), np.array(M)+np.array(M_err), color = (255/256,127/256,80/256), zorder = 2, alpha = 0.9)
plt.xlabel('Magnitude M')
plt.ylabel('log$_{10}$(N(M))')
plt.grid()

plt.figure(2)
plt.plot(M,np.log10(Nm),'kx')
plt.xlabel('Magnitude M')
plt.ylabel('log$_{10}$(N(M))')
plt.grid()
plt.xlim(12,14)

M_clipped = []
M_err_clipped = []
Nm_clipped = []

for t in range(0,len(M)):
    if M[t] > 11 and M[t] < 14:
        M_clipped.append(M[t])
        M_err_clipped.append(M_err[t])
        Nm_clipped.append(Nm[t])

fit2,cov2 = np.polyfit(M_clipped, np.log10(Nm_clipped), 1, cov=True )
linebestfit= np.poly1d(fit2)
h = np.arange(0,20,1)
plt.plot(h,linebestfit(h), color='turquoise',label='Experimental value trend line')
plt.fill_betweenx(np.log10(Nm_clipped), np.array(M_clipped)-np.array(M_err_clipped), np.array(M_clipped)+np.array(M_err_clipped), color = (255/256,127/256,80/256), zorder = 2, alpha = 0.9)
plt.legend()
plt.ylim(1,3)
print(fit2)
print(cov2)
plt.text(13,2,'Gradient = '+str(fit2[0]))
