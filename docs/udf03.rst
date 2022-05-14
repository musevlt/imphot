.. _UDF03:

The photometry of MUSE field UDF03
==================================

Muse field UDF03 is a 1x1 arc-minute field centered at Right Ascension
03:32:45, and declination -27:47:12, within the Hubble UDF. The HST
image of this region, seen through the F606W filter and resampled onto
the pixel grid of the MUSE images of UDF03, is shown below.

.. image:: _static/imphot/hst_udf03.jpeg

This field contains no stars or other bright point sources, so its
photometric parameters can only be fit using the global image fitting
method. When a fit is performed using the F606W HST image, the results
are as follows::

  % fit_photometry hst_F606W_for_UDF03.fits wfc_F606W_UDF03.fits --fix_beta=2.8 --hardcopy=jpeg
  # MUSE observation ID              Method    Flux    FWHM    beta      Flux  x-offset  y-offset
  #                                           scale     (")            offset       (")       (")
  #--------------------------------- ------  ------  ------  ------  --------  --------  --------
                     wfc_F606W_UDF03  image  0.9798  0.6685  2.8000   0.05663  -0.01618   0.01768

This recorded the following plot of the fitted images and their residuals:

.. image:: _static/imphot/udf03_image_fit.jpeg

The residual image is virtually empty, implying that a good fit was
achieved.

FWHM versus wavelength
----------------------

When the above fit was performed on images with the response curves of
the HST F606W, F775W, F814P, and F850LP filters, the fitted FWHMs of
the PSF had the values shown in the following plot.

.. image:: _static/imphot/udf03_fwhms_vs_lambda.png

The FWHMs closely follow a straight line. The plotted line is the best
fit line through them.
