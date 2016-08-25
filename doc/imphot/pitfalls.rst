.. _pitfalls:

Avoiding common pitfalls
========================

Always exclude stars
--------------------

Traditionally, photometric parameters have been derived by performing
PSF fits to bright stars. This is a good method for finding the PSF
dimensions of a single exposure, or for aligning multiple exposures
that are taken over a few hours. However bright stars usually have
enough proper motion that they should not be used to align multiple
images that are taken more than a few weeks apart. Similarly, stars
should not be used to align one exposure with an observation from
another telescope at a different epoch. Using stars to align images
generally leads to other sources in the field not being well aligned.

If many observations of a field are taken over a long period, their
relative positions must be aligned before the images can be summed to
generate a higher SNR image. The alignment procedure is generally
automated, using either cross-correlation, or an image fitting
procedure like the one implemented by the `imphot` module. It is
critically important that any stars in the images be masked out before
this procedure is performed. Stars are usually the brightest and
narrowest objects in a field, so their positions bias any automated
measurement of the image misalignment. The result of summing images
that have been mis-aligned in this way, is that all sources except the
star are smeared slightly along the proper motion direction of the
star.

When the photometric fitting procedures provided by the `imphot`
module are applied to an image that has been affected by the above
problem, a couple of common symptoms are seen:

1. Measurements of the PSF of the bright star in the image yield a
   slightly narrower PSF than measurements of the PSF of other sources
   in the image, because the other sources have been smeared along the
   proper motion direction of the star.
2. The flux scale needed to scale the pixel values of the MUSE image
   to match those of the equivalent HST image, is higher for the star
   than for other sources. Again this is because the fluxes of the
   other sources are spread over a wider area by the smearing, whereas
   the exposures of the star have been optimally added because of
   their excellent alignment.

To avoid these problems, stars should always be excluded from
photometric fits, both when constructing summed images of multiple
exposures, and when analyzing the summed images.  When using the
:ref:`fit_photometry<fit_photometry>` script, or the equivalent python
functions in the `imphot` module, this can be done by passing a ds9
region file to the procedure. This ds9 region file should contain
circular regions centered on each star, each marked as *Exclude*
regions by ds9. Usually a radius of 2 to 3 arcsec is sufficient for
these regions. For details on how to create compatible region files,
and how to pass them to the :ref:`fit_photometry<fit_photometry>`
script, see the :ref:`regions` section.

Exclude bright QSOs
-------------------

Bright QSOs (Quasi-Stellar-Objects) are point-sources that look
similar to stars. However unlike most stars, they don't have
measurable proper motions. Also unlike stars, their fluxes are often
significantly variable on both short and long timescales.  Any image
that contains a variable source will yield a poor fit when the
photometric image-fitting method is applied to it. This is because
variable sources break the assumption that any differences between a
MUSE image and a corresponding HST image, are solely caused by
instrumental and atmospheric effects.

It is tempting to use a QSO as a stable position reference for
aligning multiple exposures. However any jitter in the fitted centroid
of this source, due to noise or atmospheric effects, then becomes
applied to all sources in the image. The end result will be that the
effective PSFs of all other sources in the image will be slightly
wider than they could be.  Similarly, any systematic position shift
due to a localized plate distortion at the position of the QSO, ends
up being applied to the whole image.  The result of this is that other
sources in the image will be slightly misaligned, when compared to
images from other telescopes, such as the HST.

For the reasons discussed above it is better to base the alignment of
an image on an average misalignment, measured from all of the sources
in the image, than on the centroid of a single QSO. To do this, the
QSO should be excluded from any automated alignment procedure. For the
procedures implemented by the `imphot` module, this can be done using
ds9 region files, as described in the previous section on stars.

QSOs should also be excluded when fitting for the photometric
parameters of images that have already been aligned, since their flux
variations will otherwise degrade the quality of the fit.

Exclude variable sources
------------------------

For the reasons already discussed in the above sections, any sources
with variable fluxes or positions should be excluded using ds9 region
files, both when aligning images before summing, and when fitting for
the PSF of an image. Weak variable sources may not have much affect on
the fits, but it is easy to mask them out, so this is worth
trying. The quickest way to find out if there are any problematic
sources of this form, is to perform a photometry fit without
specifying a region file, request that the residual image be displayed
(eg. using the :ref:`--display<plotting_options>` option of
:ref:`fit_photometry<fit_photometry>`, or by passing the
``display=True`` to the `~muse_analysis.imphot.fit_image_photometry()`
function), then see if anything significant remains in the residual
image of the fit.
