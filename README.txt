Computational Image Processing: A Deep Galaxy Survey
built with Python by Harry Anthony and Jie Sing Yoo

Aim: To determine the number counts of galaxies in the local universe.



The main branch should include the following files:

1) Read_image.py
This code loads the original image into a data array, with each element corresponds the number of counts of a pixel. 
It plots the histogram of the number of counts, and calculate the mean and standard deviation of the distribution.
It removes the noise, bright stars and blooming effect from the original image, and save the result as a FITS file.

2) Cluster_generator.py
This code generates 2D gaussian plot as stars on the image for investigation on the effect of stars.

3) Aperture_methods.py
Method to calculate the number of counts of galaxy with the circular and elliptical aperture is defined in this python script.

4) PSF_methods.py
Method for Richardson-Lucy deconvolution of data is defined in this script.

5) Magnitude_plot.py
This code is used to run and store the results of galaxy survey




How to use the code
The code can only be executed in Spyder IDE, as it includes certain features(e.g. ‘#%%’) that could only be run using the Spyder IDE.
Run 'pip install photutils' in C drive in the command console to install photutils library.
Install scikit-image (https://scikit-image.org/docs/stable/install.html)

Firstly, run the 'mosaic.fits' file in read_image.py. The final data with the noise and bright stars removed would be stored as 'Mosaic_no_blooming20.fits'.

To add artificial bright objects onto the image, run Mosaic_no_blooming20.fits in Cluster_generator.py 

To conduct galaxy survey and store the results, run Mosaic_no_blooming20.fits in magnitudeplot.py 
To use a circular aperture instead of an ellitical one, change 'Counts_all_galaxies_elliptical_aperture' to 'Counts_all_galaxies_circular_aperture' in magnitudeplot.py.

"""""""""NEED TO ADD PSF METHODS DECONVOLUTION
