Computational Image Processing: A Deep Galaxy Survey
built with Python by Harry Anthony and Jie Sing Yoo

Aim: To determine the number counts of galaxies in the local universe.



The main branch should include the following files:

1) Read_image.py
This file loads the original image into a data array, with each element corresponds the number of counts of a pixel. 
It plots the histogram of the number of counts, and calculate the mean and standard deviation of the distribution.
It removes the noise, bright stars and blooming effect from the original image, and save the result as a FITS file.

2) Cluster_generator.py
This file adds 2D Gaussian stars to a image, in order to investigate the impact of non-uniformly distributed stars on the galaxy catalogue statistics.

3) Aperture_methods.py
This file contains several methods which are be used to detect galaxies and measure their number of counts using either circular or elliptical apertures. 

4) PSF_methods.py
This file contains several methods for Richardson-Lucy deconvolution of the data against three different point source functions: an airy disk, a 2D Gaussian and a 2D top hat function.

5) Magnitude_plot.py
This code utilises the methods from Aperture_methods.py and PSF_methods.py to create galaxy surveys and write them into a text file (examples are shown in branch 'Catalogues'). The code also reads these galaxy catalgoue text files and uses the data to create a plot of magnitude m against the number of sources below the given magnitude N(m).

The code also contains two branches:
- fits-file (Sontaining the astronomical data used by the code)
- Catalogue (Containing the galaxy surveys the code has created)



How to use the code
The code can only be executed in Spyder IDE, as it includes certain features(e.g. ‘#%%’) that could only be run using the Spyder IDE.

The following libraries are required:
  - photutils
    Run 'pip install photutils' in C drive in the command console to install photutils library. (https://photutils.readthedocs.io/en/stable/install.html)
  - scikit-image
    Install scikit-image (https://scikit-image.org/docs/stable/install.html)

Firstly, run the 'mosaic.fits' file in read_image.py. The final data with the noise and bright stars removed would be stored as 'Mosaic_no_blooming20.fits'.

To add artificial bright objects onto the image, run Mosaic_no_blooming20.fits in Cluster_generator.py 

To conduct galaxy survey and store the results, run Mosaic_no_blooming20.fits in magnitudeplot.py 
To use a circular aperture instead of an ellitical one, change 'Counts_all_galaxies_elliptical_aperture' to 'Counts_all_galaxies_circular_aperture' in magnitudeplot.py.

To recover the original image that was blurred by a point source function (psf), the 'psf_methods.py' file contains functions that deconvolves the image with a several different point source functions through a Richardson-Lucy algorithm. These methods can be applied to the galaxy catalogue by changing the parameters deconvolve to 'true' and psf to 'Airy', 'Gaussian' or 'Top_hat' in the functions 'Counts_all_galaxies_elliptical_aperture' and 'Counts_all_galaxies_circular_aperture'.
