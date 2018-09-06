===================================================
The *imphot* python module and command-line scripts
===================================================

:author: Martin Shepherd

The ``imphot`` (Imaging Photometry) python module provides a set of
functions and scripts that can be used to determine the photometric
parameters of individual MUSE images. These parameters include the
width of the point-spread function, the scaling and zero offset of the
flux in the image, and the pointing error of the observation. The
photometric parameters of a MUSE image are estimated by comparing it
to an HST image of the same region of sky. Two algorithms are
provided:

1. An overall image fitting algorithm is provided that uses the HST
   image as an estimate of the flux distribution on the sky. It uses
   the HST image to reproduce the MUSE image by convolving the HST
   image with a PSF that best reproduces the resolution of the MUSE
   image, scales and offsets the HST image with numbers that best
   reproduce the fluxes of the MUSE image, and shifts the HST image by
   an amount that best lines up features in the two images. The PSF
   parameters, calibration factors and pointing offsets that best
   reproduce the MUSE image, are the outputs of this algorithm.

2. A star fitting algorithm is also provided. This fits a Moffat PSF
   profile to a given star in both the MUSE image and the HST
   image. The FWHM and beta parameters that are fitted to the star in
   the MUSE observation, are reported as the fitted PSF parameters of
   the image, whereas the calibration errors and pointing errors of
   the MUSE image are determined by comparing the total fluxes, flux
   offsets and positions of the star fits in the two images.

The imphot module provides scripts that perform these algorithms,
either separately or simultaneously. It also provides python functions
and objects that can be used to perform these fits as part of external
scripts. In both cases, multi-process iterators are available for
efficiently fitting multiple MUSE exposures of a single field in
parallel. The fitting procedures also provide the option of plotting
the results of each fit, as it completes, and of generating either
verbose or summarized textual reports about the fitted parameters.

General Usage Instructions
==========================

The following are the main steps involved in fitting photometry to
one or more exposures of a single field.

- Preparatory steps:

  1. For each of the following HST WFC imaging filters, F606W, F775W,
     F814W and F850LP, obtain an HST image of the target field.

  2. Obtain one or more MUSE exposures of the target field. These
     should be MUSE cubes.

  3. For each of the HST images, use the filter curve of that image to
     extract a MUSE image with an equivalent spectral response from
     each MUSE cube, as described in section:

     :ref:`make_wideband_image`

  4. Resample the HST image(s) onto the pixel grid of the MUSE images,
     as described in section:

     :ref:`regrid_hst_to_muse`.

  5. Identify any bright stars and QSOs in the target field and use
     ds9 to create a region file that excludes all of them. If you are
     interested in an individual star or QSO in the field, also create
     a separate region-file for it that selects a circular area around
     it of about 3 arcsec radius.

- Fitting steps:

  1. Run the image and/or star fitting procedures on each of the
     extracted MUSE images of the field. If the region contains any
     stars, it is usually best to exclude them, as described in the
     :ref:`pitfalls<pitfalls>` section.

  2. Repeat the fitting procedure as many times as desired, such as
     experimenting with star fits, or changing which photometric
     parameters are to be held fixed at specified values, or
     restricting the fit to operate on different parts of the field,
     to check that the results are consistent from one part of the
     field to another.

Tutorials
=========

* :ref:`tutorial`

* :ref:`pitfalls`

* :ref:`udf_photometry`

Resources
=========

* :ref:`hst_filter_curves`

* :ref:`udf_region_files`

* :ref:`docstrings`

* :ref:`whitepapers`

Provided command-line scripts
=============================

The following scripts are provided for fitting photometry parameters
from the command line:

  :ref:`regrid_hst_to_muse<regrid_hst_to_muse>`
     Resample an HST image onto the same grid of pixels as a MUSE image.
  :ref:`make_wideband_image<make_wideband_image>`
     Assemble a MUSE image with a given spectral response curve
     from a MUSE cube.
  :ref:`fit_photometry<fit_photometry>`
     Fit for the photometry parameters of one or more MUSE images of a field.

A categorized list of python functions and objects
==================================================

In addition to the command-line scripts listed above, python functions
and objects are provided for invoking any of the preparatory
procedures and fitting procedures from custom python scripts.  The
source-code documentation for each of these functions and classes can
be found in section, :ref:`docstrings`. The following is a categorized
index of links into that section.

Image Fitting
-------------
  `~muse_analysis.imphot.fit_image_photometry()`
     Fit for the photometric parameters of a single MUSE image, by
     simulating the effects of different seeing conditions,
     calibration errors on pointing errors on a template HST image.
  `~muse_analysis.imphot.FitImagePhotometryMP()`
     An iterator which uses multiple processes to run
     `~muse_analysis.imphot.fit_image_photometry()` on a sequence of MUSE exposures of a
     field.
  `~muse_analysis.imphot.FittedImagePhotometry`
     The object in which photometry results are returned by
     `~muse_analysis.imphot.fit_image_photometry()` and `~muse_analysis.imphot.FitImagePhotometryMP()`.

Star Fitting
------------
  `~muse_analysis.imphot.fit_star_photometry()`
    Fit for the photometric parameters of a single MUSE image, by
    fitting a Moffat PSF profile to a star in both the MUSE image and
    a corresponding HST image.
  `~muse_analysis.imphot.FitStarPhotometryMP()`
     An iterator which uses multiple processes to run
     `~muse_analysis.imphot.fit_star_photometry()` on multiple MUSE exposures of a
     field.
  `~muse_analysis.imphot.FittedStarPhotometry`
     The object in which photometry results are returned by
     `~muse_analysis.imphot.fit_star_photometry()` and `~muse_analysis.imphot.FitStarPhotometryMP()`.

Combined image and star fitting
-------------------------------
  `~muse_analysis.imphot.fit_image_and_star_photometry()`
    Fit for the photometric parameters of a single MUSE image,
    returning a tuple of results from both
    `~muse_analysis.imphot.fit_image_photometry()` and
    `~muse_analysis.imphot.fit_star_photometry()`.
  `~muse_analysis.imphot.FitImageAndStarPhotometryMP()`
     An iterator which uses multiple processes to run
     `~muse_analysis.imphot.fit_image_and_star_photometry()` on a sequence of
     MUSE exposures of a field.

Region handling
---------------
  `~muse_analysis.imphot.ds9regions`
     A module for reading simple ds9 region files.

Miscellaneous
-------------
  `~muse_analysis.imphot.rescale_hst_like_muse()`
     Given an HST image, convert its flux units from electrons s\
     :sup:`-1` to the flux units of a MUSE observation (usually
     1e-20 erg cm\ :sup:`-2` s\ :sup:`-1` Angstrom\ :sup:`-1`).
  `~muse_analysis.imphot.regrid_hst_like_muse()`
     Resample an HST image onto the pixel grid of a MUSE image.
  `~muse_analysis.imphot.image_grids_aligned()`
     Return ``True`` if the pixel coordinates of two MPDAF
     images are the same.
  `~muse_analysis.imphot.FittedPhotometry`
     The base class of `muse_analysis.imphot.FittedImagePhotometry` and
     `muse_analysis.imphot.FittedStarPhotometry`. This contains the fitted photometric
     parameters and functions to report them, without the extra details
     that are specific to the two fitting methods.
  `~muse_analysis.imphot.FittedValue`
     A single fitted photometric value, which includes the
     value, its uncertainty and a bool value that indicates whether this
     parameter was held fixed during the least-squares fit, or allowed
     to vary.
  `~muse_analysis.imphot.HstFilterInfo`
     An object that contains the filter characteristics of an HST image.
  `~muse_analysis.imphot.UserError`
     An exception which is raised for errors that are due to incorrect
     user-input. In the provided scripts, errors of this type are caught
     and reported without showing a traceback of the call stack.
  `~muse_analysis.imphot.WorkerError`
     An exception forwarded from a secondary worker process.
  `~muse_analysis.imphot.ImphotArgumentParser`
     A command-line argument parser for receiving user arguments for the
     `~muse_analysis.imphot.fit_image_photometry()`,
     `~muse_analysis.imphot.fit_star_photometry()`, `~muse_analysis.imphot.FitImagePhotometryMP()`,
     `~muse_analysis.imphot.FitStarPhotometryMP()` and
     `~muse_analysis.imphot.FitImageAndStarPhotometryMP()` functions.
  `~muse_analysis.imphot.extract_function_args()`
     This takes either a dictionary of keyword/value pairs, or an
     ``argparse.Namespace`` object, returned by
     ``argparse.ArgumentParser.parse_args()``, and creates a new
     dictionary that just contains the subset of the keyword/value pairs
     that are valid arguments of a specified function.

Page index
----------
.. toctree::
   :maxdepth: 2

   regrid_hst_to_muse.rst
   make_wideband_image.rst
   fit_photometry.rst
   tutorial.rst
   pitfalls.rst
   docstrings.rst
   whitepapers.rst
   udf_photometry.rst
